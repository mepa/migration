#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <utility>

#include "function.hh"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf.h>

using namespace std;

inline double sq(double x) { return x * x; }
inline double cub(double x) { return x * x * x; }

const double pc = 3.0e18;
const double kpc = 1.0e3 * pc;
const double M_sun = 2.0e33;
const double G = 6.67e-8;
const double yr = 3.14e7;
const double two_pi = 2.0 * M_PI;
const double four_pi = 4.0 * M_PI;

const double R_e = 500.0 * pc;
const double n = 1.5;
const double mass = 1.0e+09 * M_sun;

const double r_min = 1.0 * pc;
const double r_max = 50.0 * kpc;

const double mass_halo = 1.0e+10 * M_sun;
const double concentration = 3.5;
const double r_vir = 10.0 * kpc;

const double M_min = 100.0 * M_sun;
const double M_max = 1.0e+05 * M_sun;

size_t N = 1000;
size_t num_pts = 100;
size_t step = N / num_pts; 

struct NFW
{
  double mass, r_vir, c, r_s, rho_s;

  NFW(double _mass, double _r_vir, double _c) : 
    mass(_mass), r_vir(_r_vir), c(_c), r_s(r_vir / c), 
    rho_s(mass * cub(c) * (1.0 + c) / (4.0 * M_PI * cub(r_vir) * ((1.0 + c) * log(1.0 + c) - c))) { }

  double rho(double r) const { return rho_s / ((r / r_s) * sq(1.0 + r / r_s)); }

  double M(double r) const 
  {
    return four_pi * rho_s * cub(r_s) * (- log(r_s) * r - log(r_s) * r_s - r + log(r + r_s) * r_s + log(r + r_s) * r) / (r + r_s);
  }  
};


struct Model
{
  Model(double _R_e, double _n, double _mass) : 
    R_e(_R_e), n(_n), inverse_n(1.0 / n), mass(_mass) // Take p(1.5) ~ 2/3
  { 
    p = 2.0 / 3.0; // for n = 1.5: p = 1.0 - 0.6097 * inverse_n + 0.05463 * sq(inverse_n) ~ 0.617813;
    half_p = 1.0 / 3.0;
    a = (3.0 - p) * n;
    gamma_fcn = gsl_sf_gamma(a);
    normalization = mass / ( four_pi * cub(R_e) * n * gamma_fcn );

    if (a <= 0) 
      {
	cout << "gamma function undefined" << endl;
	abort();
      } 
  }
  
  double R_e; // R_e actually means R_0 here
  double n;
  double inverse_n;
  double mass;
  double normalization;
  double p;
  double half_p;
  double a;
  double gamma_fcn;

  double rho(double r) const
  {
    double x = r / R_e;

    return normalization * pow(x, -p) * exp(- pow(x, inverse_n) );
  }

  double M(double r) const
  {
    double x = r / R_e;
    //const double gamma_inc_lower = gamma_fcn - gsl_sf_gamma_inc(a, pow(x, inverse_n));
    const double approx_gamma_inc_lower 
      = (-20 * x - 8 * x * pow(x, p) - 30 * pow(x, half_p)) / (8.0 * exp(pow(x, p))) 
      + 15 * sqrt(M_PI) * erf(pow(x, half_p)) / 8.0;

    return four_pi * normalization * cub(R_e) * n * approx_gamma_inc_lower;
  }
};


vector<double> Project(const vector<double> & r, const vector<double> & rho)
{
  vector<double> Sigma(r.size(), 0.0);

  for(size_t i = 0; i < r.size(); i ++)
    {
      Sigma[i] += 2.0 * sqrt(r[i + 1] * r[i + 1] - r[i] * r[i]) * 0.5 * (rho[i + 1] + rho[i]);

      for(size_t j = i + 1; j < r.size() - 1; j ++)
	{
	  const double r1 = r[j];
	  const double r2 = r[j + 1];
	  const double int_1 = rho[j] * 2.0 * r1 / sqrt(r1 * r1 - r[i] * r[i]);
	  const double int_2 = rho[j + 1] * 2.0 * r2 / sqrt(r2 * r2 - r[i] * r[i]);

	  Sigma[i] += 0.5 * (r2 - r1) * (int_1 + int_2);
	}
    }

  return Sigma;
}


struct Profile 
{
  Model model;
  NFW nfw;
  
  vector<double> r, Sigma_baryon, rho_baryon, mass_baryon, Sigma_NFW, rho_NFW, rho_total, mass_total;

  double mass_nucleus;
  vector<double> evolving_mass_profile;
  
  Profile(const Model & _model, const NFW & _nfw,
	  double r_min, double r_max, size_t c) : 
    model(_model), nfw(_nfw), r(c), Sigma_baryon(c), rho_baryon(c), mass_baryon(c), Sigma_NFW(c), rho_NFW(c), rho_total(c), mass_total(c),
    mass_nucleus(0.0), evolving_mass_profile(c, 0.0) 
  {
    for(size_t i = 0; i < c; i ++)
      {
	r[i] = r_min * exp(double(i) * log(r_max / r_min) / double(c - 1));
	rho_baryon[i] = model.rho(r[i]);
	rho_NFW[i] = nfw.rho(r[i]); 
	rho_total[i] = rho_baryon[i] + rho_NFW[i];
      }
    Sigma_baryon = Project(r, rho_baryon);
    Sigma_NFW = Project(r, rho_NFW);

    const double r0 = r.front();
 
    const double central_mass_baryon = model.M(r0);
    const double central_mass_nfw = nfw.M(r0);

    cout << "# central_mass_baryon, central_mass_nfw: " << central_mass_baryon / M_sun << " "
	 << central_mass_nfw / M_sun << endl; 
  
    for(size_t i = 0; i < c; i ++)
      {
	mass_baryon[i] = model.M(r[i]);
	mass_total[i] = mass_baryon[i] + nfw.M(r[i]);
      }
  }

  size_t size() const { return r.size(); } 

  double randomRadius() const
  {  
    const double partial_mass = drand48() * (mass_baryon.back() - mass_baryon.front()) + mass_baryon.front();
  
    double r_min = r.front();
    double r_max = r.back();

    const double tiny_delta_r = 1.0e-6 * r_min;

    while(r_max - r_min > tiny_delta_r)
      {
	const double r_middle = 0.5 * (r_min + r_max);
	if(partial_mass < model.M(r_middle)) r_max = r_middle;
	else r_min = r_middle;
      }
    return 0.5 * (r_min + r_max);
  }
  
  double getMassBaryon() const { return mass_baryon.back(); }
  double getMassTotal() const { return mass_total.back(); }

};


const Model model(R_e, n, mass);

const NFW nfw(mass_halo, r_vir, concentration);

const Profile profile(model, nfw, r_min, r_max, N);


ostream & operator<< (ostream & out, const Profile & profile)
{
  for(size_t i = 0; i < profile.size(); i ++)
    {
      out << profile.r[i] << " " << profile.Sigma_baryon[i] << " " << profile.Sigma_NFW[i] << " " << profile.rho_baryon[i] << " " << profile.rho_NFW[i] << " " << profile.mass_baryon[i] << " " << profile.mass_total[i] << endl;
    }
  return out;
}

double timeMigration(double cluster_mass, double r, const Profile & profile)
{
  const double R = max<double> (r, profile.r.front());

  const double ln_Lambda = 5.0;

  const size_t i = size_t(floor(double(profile.r.size()) * log(R / profile.r.front()) / log(profile.r.back() / profile.r.front())));

  const double evolving_mass = profile.evolving_mass_profile[i] + 
    (profile.evolving_mass_profile[i + 1] - profile.evolving_mass_profile[i]) * (R - profile.r[i]) / (profile.r[i + 1] - profile.r[i]);

  const double evolving_density = ((profile.evolving_mass_profile[i + 1] - profile.evolving_mass_profile[i]) 
				   / (profile.r[i + 1] - profile.r[i])) / (four_pi * sq(R));

  const double mass = profile.nfw.M(r) + profile.mass_nucleus + max<double> (profile.model.M(r), evolving_mass);

  const double rho = profile.nfw.rho(r) + max<double> (profile.model.rho(R), evolving_density);
  
  return 0.5 * sqrt(mass) * pow(R, 1.5) / (cluster_mass * ln_Lambda * sqrt(G)) + 0.125 * pow(mass, 1.5) / (cluster_mass * ln_Lambda * pow(r, 1.5) * sqrt(G) * M_PI * rho);
}

double timeDisruption(double cluster_mass, double r, const Profile & profile)
{
  const double radius_factor = 1.0; // circular orbit?
  const double f_dis = 0.2;
  const double gamma = 0.62;
  const double radius_pericenter = max<double> (radius_factor * r, profile.r.front());

 if(radius_pericenter < profile.r.front() || radius_pericenter >= profile.r.back()) 
    {
      cout << "out of range: " << radius_pericenter << endl;
      abort();
    }

  const size_t i = size_t(floor(double(profile.r.size()) * log(radius_pericenter / profile.r.front()) / log(profile.r.back() / profile.r.front())));

  const double evolving_mass = profile.evolving_mass_profile[i] + 
    (profile.evolving_mass_profile[i + 1] - profile.evolving_mass_profile[i]) * (radius_pericenter - profile.r[i]) / (profile.r[i + 1] - profile.r[i]);

  const double evolving_density = ((profile.evolving_mass_profile[i + 1] - profile.evolving_mass_profile[i]) 
				   / (profile.r[i + 1] - profile.r[i])) / (four_pi * sq(radius_pericenter));

  const double mass = max<double> (evolving_mass, profile.model.M(radius_pericenter)) + profile.nfw.M(radius_pericenter) + profile.mass_nucleus;

  const double rho = max<double> (profile.model.rho(radius_pericenter), evolving_density) + profile.nfw.rho(radius_pericenter);

  const double Omega = sqrt(G * mass / cub(radius_pericenter));

//   // OLD
//   const double t_0 = f_dis / Omega;
//   const double t_4 = 1.355 * yr * pow(10.0, 4.0 * gamma) * (1.0 / gamma) * pow(t_0 / yr, 0.967);
//   const double t_dis = t_4 * pow(cluster_mass / (1.0e4 * M_sun), gamma);

  // NEW
  const double dOmega = - 1.5 * Omega / r + 0.5 * Omega * four_pi * sq(r) * rho / mass;
  const double t_0 = f_dis / abs(r * dOmega);
  const double t_dis = t_0 * pow(cluster_mass / M_sun, gamma);

  return t_dis;
}

int func (double t, const double y[], double f[],
	  void *params)
{
  const pair<double, Profile const *> & data = 
    * static_cast<pair<double, Profile const *> *> (params);

  const double r = y[0];
  const double mass = y[1];

  if (isnan(r) || r < 0.0) 
    {
      return GSL_FAILURE;
    }

  const double t_mig = timeMigration(mass, r, * data.second);
  const double t_dis = timeDisruption(mass, r, * data.second);

  f[0] = - y[0] / t_mig;

  f[1] = - y[1] / t_dis;

  return GSL_SUCCESS;
}

struct Element
{
  double mass, r, dM_dr;
  
  Element(double _mass, double _r, double _dM_dr) : mass(_mass), r(_r), dM_dr(_dM_dr) { }
};

struct DepositedMass : public vector<Element> 
{
  double disrupt_radius, disrupt_mass, center_mass, time;

  DepositedMass() : vector<Element> (), disrupt_radius(0.0), disrupt_mass(0.0), center_mass(0.0), time(0.0) { }
};

void decayOrbit(double cluster_mass_init, double r_init, const Profile & profile, DepositedMass & dM_dr)
{
  dM_dr.clear();

  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk4;

  gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
     
  const pair<double, Profile const *> data(cluster_mass_init, &profile);

  gsl_odeiv_system sys = {func, 0, 2, const_cast<pair<double, const Profile *> *> (&data)};
     
  double t = 0.0;
    
  double h = yr;
  
  double y[2] = { r_init, cluster_mass_init }, y_err[1];

  double dydt_in[2], dydt_out[2];

  GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);

  dM_dr.push_back(Element(cluster_mass_init, r_init, dydt_in[1] / dydt_in[0]));

  while (true)
    {
      const double y_save[2] = { y[0], y[1] };
     
      int status;

      do 
	{
	  status = gsl_odeiv_step_apply (s, t, h, 
					     y, y_err, 
					     dydt_in, 
					     dydt_out, 
					     &sys);

	  if (status != GSL_SUCCESS)
	    {
	      h *= 0.5;
	      y[0] = y_save[0]; 
	      y[1] = y_save[1];

	    }
	}
      while (status != GSL_SUCCESS);

      dydt_in[0] = dydt_out[0];
      dydt_in[1] = dydt_out[1];
      
      t += h;

      const double epsilon_step = 0.01;

      h = min<double> (- epsilon_step * y[0] / dydt_in[0], - epsilon_step * y[1] / dydt_in[1]);

      const double r = y[0];
      const double cluster_mass = y[1];

      dM_dr.push_back(Element(cluster_mass, r, dydt_out[1] / dydt_out[0]));
  
      if(r < profile.r.front()) 
	{
	  dM_dr.disrupt_radius = r;
	  dM_dr.center_mass = cluster_mass;
	  dM_dr.time = t;
	  break;
	}

      if(cluster_mass < 100.0 * M_sun) 
	{
	  dM_dr.disrupt_radius = r;
	  dM_dr.disrupt_mass = cluster_mass;
	  dM_dr.time = t;
	  break;
	}
    }

  gsl_odeiv_step_free (s);
}

struct Cluster 
{
  double mass, r;
  Cluster(double _mass, double _r) : mass(_mass), r(_r) { }
};

ostream & operator<< (ostream & out, const Cluster & cluster) 
{
  out << cluster.mass << " " << cluster.r;
  return out;
}

struct Clusters : public vector<Cluster> 
{ 
  Clusters(double M_min, double M_max, double slope, const Profile & profile)
  {
    double sum_mass = 0.0;

    while(sum_mass < profile.getMassBaryon()) 
      {
	const Cluster cluster(genMass(M_min, M_max, slope), profile.randomRadius());

	sum_mass += cluster.mass;

	push_back(cluster);
      }

    const double factor = profile.getMassBaryon() / sum_mass;

    for(vector<Cluster>::iterator i = begin(); i != end(); i ++)
      {
	i->mass *= factor;
      }
  }
  
  double genMass(double M_min, double M_max, double slope)
  {
    const double x = drand48();
    if(slope != 2.0) return pow(pow(M_min, 2.0 - slope) 
				- pow(M_min, 2.0 - slope) * x 
				+ pow(M_max, 2.0 - slope) * x, 1.0 / (2.0 - slope));
    else return M_min * exp(x * log(M_max / M_min));
  }

  double mass() const 
  {
    double sum = 0.0;
    for(vector<Cluster>::const_iterator i = begin(); i != end(); i ++) sum += i->mass;
    return sum;
  }
};

bool operator< (const Cluster & cluster_1, const Cluster & cluster_2)
{
  return timeMigration(cluster_1.mass, cluster_1.r, profile) < timeMigration(cluster_2.mass, cluster_2.r, profile);
}

struct Center
{
  vector<double> mass, mass_excess, mass_excess_actual, time, time_actual;
};

int main()
{
  cout << "# R_e, n, mass, r_min, r_max, mass_halo, concentration, r_vir, M_min, M_max, N, num_pts: " 
       << R_e / pc << ", " 
       << n << ", " 
       << mass / M_sun << ", " 
       << r_min / pc << ", " 
       << r_max / pc << ", " 
       << mass_halo / M_sun << ", " 
       << concentration << ", " 
       << r_vir / pc << " "
       << M_min / M_sun << ", " 
       << M_max / M_sun <<  ", "
       << N << ", " << num_pts << endl;
 
  ofstream profile_out("profile.dat");
  profile_out << profile;

  //const double M_cl_min = 1.0e1 * M_sun;
  //const double M_cl_max = 1.0e4 * M_sun;
  //vector<double> M_cl(100);
  //for (size_t i = 0; i < M_cl.size(); i ++)
  //  {
  //    M_cl[i] = M_cl_min * exp(double(i) * log(M_cl_max / M_cl_min) / double(100 - 1));
  //    double time_dissociation = timeDisruption(M_cl[i], 5.0 * pc, profile);
  //    cout << M_cl[i] / M_sun << " " << time_dissociation / yr << endl;
  //  }
  //abort();
  //cout << time_dissociation / yr << endl;

  Clusters clusters(M_min, M_max, 2.0, profile);

  //sort(clusters.begin(), clusters.end());

  vector<double> time_mig(clusters.size(), 0.0);
  for(size_t i = 0; i < clusters.size(); i ++) 
    {
      time_mig[i] = timeMigration(clusters[i].mass, clusters[i].r, profile);
    }

  cout << "# generated " << clusters.size() << " clusters" << " "
       << "with total mass " << clusters.mass() / M_sun << endl;

ofstream cluster_out("data/script/cluster_sph_nonordered_1.0e+05_high_z.out");
  copy(clusters.begin(), clusters.end(), ostream_iterator<Cluster> (cluster_out, "\n"));

  vector<double> initial_mass(profile.r.size(), 0.0), final_mass(profile.r.size(), 0.0), final_mass_disrupt(profile.r.size(), 0.0);
  vector<double> mass_center(clusters.size(), 0.0), mass_center_excess(clusters.size(), 0.0), mass_center_excess_actual(clusters.size(), 0.0);
 
  double mass_disrupted = 0.0, mass_deposited = 0.0, mass_nucleus = 0.0;
  double time_total = 0.0;

  const double radius_center = 10.0 * pc;
  double time_center_total = 0.0;

  Center center;

  for(size_t i = 0; i < clusters.size(); i ++)
    {
      const double r_cluster = clusters[i].r;
      for(size_t j = 0; j < profile.r.size(); j ++)
	if (profile.r[j] > r_cluster) 
	  {
	    initial_mass[j] += clusters[i].mass;
	  }

      DepositedMass dM_dr;

      decayOrbit(clusters[i].mass, clusters[i].r, profile, dM_dr);
  
      reverse(dM_dr.begin(), dM_dr.end());

      for(size_t j = 1; j < dM_dr.size(); j ++)
	{
	  const double r = 0.5 * (dM_dr[j - 1].r + dM_dr[j].r);

	  const double delta_mass = (dM_dr[j].r - dM_dr[j - 1].r) * 0.5 * (dM_dr[j].dM_dr + dM_dr[j - 1].dM_dr);
 
	  for(size_t k = 0; k < profile.r.size(); k ++)
	    if (profile.r[k] > r) 
	      {
	  	final_mass[k] += delta_mass;
	      }

	  mass_deposited += delta_mass;
	}
	  
      for(size_t j = 0; j < profile.r.size(); j ++)
	if (profile.r[j] > dM_dr.disrupt_radius)
	  {
	    final_mass_disrupt[j] += dM_dr.disrupt_mass; //doesn't include mass if radii equal
	  }
      
      mass_nucleus += dM_dr.center_mass;
      mass_disrupted += dM_dr.disrupt_mass;

      const_cast<Profile &> (profile).mass_nucleus = mass_nucleus;
      
      transform(final_mass.begin(), final_mass.end(), final_mass_disrupt.begin(), const_cast<Profile &> (profile).evolving_mass_profile.begin(), plus<double> ());

      if (dM_dr.disrupt_radius < radius_center)  
	{
	  for(size_t j = 0; j < profile.r.size(); j ++)
	    if (profile.r[j] < radius_center)
	      {
		mass_center[i] = profile.mass_nucleus + profile.evolving_mass_profile[j];
		mass_center_excess[i] = mass_center[i] - profile.model.M(profile.r[j]);
		mass_center_excess_actual[i] = mass_center[i] - initial_mass[j];
	      }

	  if (dM_dr.time > time_center_total) time_center_total = dM_dr.time;

	  center.mass.push_back(mass_center[i]);
	  center.mass_excess.push_back(mass_center_excess[i]);
	  center.mass_excess_actual.push_back(mass_center_excess_actual[i]);
	  center.time.push_back(time_mig[i]);
	  center.time_actual.push_back(dM_dr.time);
	}

      if (dM_dr.time > time_total) time_total = dM_dr.time;
     }

ofstream center_out("data/script/center_sph_nonordered_1.0e+05_high_z.out");
  for (size_t i = 0; i < center.mass.size(); i ++) 
    center_out << center.time[i] << " " 
	       << center.time_actual[i] << " " 
	       << center.mass[i] << " " 
	       << center.mass_excess[i] << " " 
	       << center.mass_excess_actual[i] << endl;

  cout << "# mass_nucleus " << mass_nucleus / M_sun << " " << mass_nucleus / clusters.mass() << endl;
   
  cout << "# mass_disrupted " << mass_disrupted / M_sun << " " << mass_disrupted / clusters.mass() << endl;

  cout << "# mass_deposited " << mass_deposited / M_sun << " " << mass_deposited / clusters.mass() << endl;

  cout << "# before after difference relative "
       << clusters.mass() / M_sun << " " 
       << (mass_nucleus + mass_disrupted + mass_deposited) / M_sun << " "
       << fabs(mass_nucleus + mass_disrupted + mass_deposited - clusters.mass()) / M_sun << " " 
       << fabs(mass_nucleus + mass_disrupted + mass_deposited - clusters.mass()) / clusters.mass() << " " 
       << endl;

  vector<double> final_mass_total(profile.r.size(), 0.0); 
  for(size_t i = 0; i < profile.r.size(); i ++)
    {
      final_mass_total[i] = mass_nucleus + final_mass_disrupt[i] + final_mass[i]; 
    }

  vector<double> initial_rho(profile.r.size()), final_rho(profile.r.size()), final_rho_disrupt(profile.r.size());

  for(size_t i = 1; i < final_mass_total.size() - 1; i += 1)
    {
      initial_rho[i] = ((initial_mass[i + 1] - initial_mass[i - 1]) / (profile.r[i + 1] - profile.r[i - 1])) / (four_pi * sq(0.5 * (profile.r[i + 1] + profile.r[i - 1])));
      
      final_rho[i] = ((final_mass_total[i + 1] - final_mass_total[i - 1]) / (profile.r[i + 1] - profile.r[i - 1])) / (four_pi * sq(0.5 * (profile.r[i + 1] + profile.r[i - 1])));

      final_rho_disrupt[i] = ((final_mass_disrupt[i + 1] - final_mass_disrupt[i - 1]) / (profile.r[i + 1] - profile.r[i - 1])) / (four_pi * sq(0.5 * (profile.r[i + 1] + profile.r[i - 1])));
    }
  initial_rho[0] = initial_rho[1];
  initial_rho[profile.r.size() - 1] = initial_rho[profile.r.size() - (2 * 1)];
  final_rho[0] = final_rho[1];
  final_rho[profile.r.size() - 1] = final_rho[profile.r.size() - (2 * 1)];
  final_rho_disrupt[0] = final_rho_disrupt[1];
  final_rho_disrupt[profile.r.size() - 1] = final_rho_disrupt[profile.r.size() - (2 * 1)];

  const vector<double> initial_Sigma = Project(profile.r, initial_rho);
  const vector<double> final_Sigma = Project(profile.r, final_rho);
  const vector<double> final_Sigma_disrupt = Project(profile.r, final_rho_disrupt);

  vector<double> final_Sigma_nuc(final_mass_total.size(), 0.0);
  const double r_nucleus = 2.0 * pc;  
  for(size_t i = 0; i < final_mass_total.size(); i ++)
    {
      final_Sigma_nuc[i] += 0.5 * mass_nucleus * r_nucleus / (M_PI * pow(sq(r_nucleus) + sq(profile.r[i]), 1.5));
    }

  vector<double> final_Sigma_total(final_mass_total.size(), 0.0);
  for(size_t i = 0; i < final_mass_total.size(); i ++)
    {
      final_Sigma_total[i] = final_Sigma[i] + final_Sigma_nuc[i];
    }

  double actual_time = 0.0;
  for(size_t i = 0; i < center.time_actual.size(); i ++)
      if (center.time_actual[i] > actual_time) actual_time = center.time_actual[i];

  cout << "# (in central " << radius_center / pc << " pc) mass, mass_excess, time, time_actual: " << center.mass.back() / M_sun << " " << center.mass_excess.back() / M_sun << " " << center.time.back() / yr << " " << actual_time / yr << endl;

  double radius_upturn = 0.0;
  double mass_upturn = 0.0;
  double mass_upturn_excess = 0.0;
  double mass_upturn_excess_actual = 0.0;
  size_t index = 0;
  while(final_Sigma_total[index] - profile.Sigma_baryon[index] > 0.0)
    {
      radius_upturn = profile.r[index];
      mass_upturn = final_mass_total[index];
      mass_upturn_excess = final_mass_total[index] - profile.model.M(profile.r[index]);
      mass_upturn_excess_actual = final_mass_total[index] - initial_mass[index];

      index ++;
    }

  cout << "# (at upturn) radius, mass, mass_excess, mass_excess_actual: " << radius_upturn / pc << " " 
       << mass_upturn / M_sun << " " <<  mass_upturn_excess / M_sun << " " << mass_upturn_excess_actual / M_sun << endl;
  
  cout << "# time_total, time_center_total: " << time_total / yr << " " << time_center_total / yr << endl;

  
  for(size_t i = 0; i < final_mass_total.size(); i += step)
    {
      cout << profile.r[i] << " " 
	   << profile.Sigma_NFW[i] << " " 
	   << profile.Sigma_baryon[i] << " " 
	   << initial_Sigma[i] << " " 
	   << final_Sigma_total[i] << " " 
	   << final_Sigma_nuc[i] << " "
	   << final_Sigma_disrupt[i] << " " 
	   << final_Sigma[i] << endl; 
     }
  
  return 0;
}
