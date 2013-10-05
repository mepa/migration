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

using namespace std;

inline double sq(double x) { return x * x; }
inline double cub(double x) { return x * x * x; }

const double pc = 3.0e18;
const double kpc = 1.0e3 * pc;
const double M_sun = 2.0e33;
const double G = 6.67e-8;

const double yr = 3.14e7;

const double four_pi = 4.0 * M_PI;
const double two_pi = 2.0 * M_PI;

struct NFW
{
  double mass, r_vir, c, r_s, rho_s;

  NFW(double _mass, double _r_vir, double _c) : 
    mass(_mass), r_vir(_r_vir), c(_c), r_s(r_vir / c),
    rho_s(mass * c * c * c * (1.0 + c) / (4.0 * M_PI * r_vir * r_vir * r_vir * ((1.0 + c) * log(1.0 + c) - c))) { }

  double rho(double r) const { return rho_s / ((r / r_s) * sq(1.0 + r / r_s)); }

  double M(double r) const 
  {
    return 4.0 * M_PI * rho_s * cub(r_s) * (- log(r_s) * r - log(r_s) * r_s - r + log(r + r_s) * r_s + log(r + r_s) * r) / (r + r_s);
  }
};

struct Model
{
  Model(double _R_e, double _n, double _mass) : R_e(_R_e), n(_n), inverse_n(1.0 / n), mass(_mass), 
    normalization(calcNormalization()) { }
  
  double R_e;
  double n;
  double inverse_n;
  double mass;
  double normalization;

  double Sigma(double R) const
  {
    return normalization * exp(- pow(R / R_e, inverse_n));
  }

  double Sigma_deriv(double R) const 
  {
    return - Sigma(R) * inverse_n * pow(R / R_e, inverse_n) / R;
  }

  double rho(double r) const
  {
    const double R_max = 1.0e3 * R_e;
    
    if(r >= R_max) return 0.0;

    const size_t c = 1000;

    const double R_min = r * exp(log(R_max / r) / double(c - 1));
 
    vector<double> R(c);

    for(size_t i = 0; i < c; i ++) 
      R[i] = R_min * exp(double(i) * log(R_max / R_min) / double(c - 1));

    double s = 0.0;
    
    for(size_t i = 0; i < c - 1; i ++)
      {
	const double R_1 = R[i];
	const double R_2 = R[i + 1];
	const double I_1 = Sigma_deriv(R_1) / sqrt(R_1 * R_1 - r * r);
	const double I_2 = Sigma_deriv(R_2) / sqrt(R_2 * R_2 - r * r);

	s += 0.5 * (R_2 - R_1) * (I_1 + I_2);
      }

    s += (log(R_min + sqrt(R_min * R_min - r * r)) - log(r)) * Sigma_deriv(r);

    return - s / M_PI;
  }

 private:

  double Sigma_unnormalized(double R) const
  {
    return exp(- pow(R / R_e, 1.0 / n));
  }
  double calcNormalization() const
  {
    const double R_min = 1.0e-3 * R_e;
    const double R_max = 1.0e3 * R_e;
    const size_t c = 1000;

    vector<double> R(c);
    
    for(size_t i = 0; i < c; i ++) 
      R[i] = R_min * exp(double(i) * log(R_max / R_min) / double(c - 1));

    double s = 0.0;

    for(size_t i = 0; i < c - 1; i ++)
      {
	const double R_1 = R[i];
	const double R_2 = R[i + 1];
	const double M_1 = 2.0 * M_PI * R_1 * Sigma_unnormalized(R_1);
	const double M_2 = 2.0 * M_PI * R_2 * Sigma_unnormalized(R_2);

	s += 0.5 * (R_2 - R_1) * (M_1 + M_2);
      }

    return mass / s;
  }
};



struct Profile 
{
  Model model;
  NFW nfw;
  
  vector<double> r, rho_baryon, rho_total, Sigma, mass_baryon, mass_total, Omega, dOmega;

  Profile(const Model & _model, const NFW & _nfw,
	  double r_min, double r_max, size_t c) : 
    model(_model), nfw(_nfw), r(c), rho_baryon(c), rho_total(c), Sigma(c), mass_baryon(c), mass_total(c), Omega(c), dOmega(c) 
  {
    for(size_t i = 0; i < c; i ++)
      {
	r[i] = r_min * exp(double(i) * log(r_max / r_min) / double(c - 1));
	rho_baryon[i] = model.rho(r[i]);
	rho_total[i] = rho_baryon[i] + nfw.rho(r[i]);
	Sigma[i] = model.Sigma(r[i]);
      }

    const double rho0 = rho_baryon[0];
    const double rho1 = rho_baryon[1];
    const double r0 = r[0];
    const double r1 = r[1];
 
    const double central_mass_baryon = 4.0*M_PI*rho0*(-log(r1)+log(r0))*r0*r0*r0/(-3.0*log(r1)+3.0*log(r0)-log(rho1)+log(rho0));
    const double central_mass_nfw = nfw.M(r0);

    cout << "# central_mass_baryon = " << central_mass_baryon / M_sun << endl;
    cout << "# central_mass_nfw = " << central_mass_nfw / M_sun << endl; 
  
    double mass_baryon_sum = central_mass_baryon;
  
    for(size_t i = 0; i < c; i ++)
      {
	if(i > 0) mass_baryon_sum += 0.5 * (r[i] - r[i - 1]) * 4.0 * M_PI * 
		    (rho_baryon[i] * r[i] * r[i] + rho_baryon[i - 1] * r[i - 1] * r[i - 1]);

	mass_baryon[i] = mass_baryon_sum;
	mass_total[i] = mass_baryon[i] + nfw.M(r[i]);
	Omega[i] = sqrt(G * mass_total[i] / (r[i] * r[i] * r[i]));
	if(i > 0) dOmega[i] = - (Omega[i] - Omega[i - 1]) / (r[i] - r[i - 1]); // MINUS SIGN
	if(i == 1) dOmega[0] = dOmega[1];
      }
  }

  size_t size() const { return r.size(); } 

  double randomRadius() const
  {  
    const double partial_mass = drand48() * (mass_baryon.back() - mass_baryon.front()) + mass_baryon.front();

    size_t i = 0;

    while(mass_baryon[i] < partial_mass) i ++;

    return r[i - 1] + (partial_mass - mass_baryon[i - 1]) * (r[i] - r[i - 1]) / (mass_baryon[i] - mass_baryon[i - 1]);
  }
  
  double getMassBaryon() const { return mass_baryon.back(); }
  double getMassTotal() const { return mass_total.back(); }

  double interp(const vector<double> & variable, double radius) const
  {
    if(radius < r.front()) return variable.front();
    if(radius >= r.back()) return variable.back();
    size_t i = 0;
    while(r[i] < radius) i ++;
    return variable[i - 1] + (radius - r[i - 1]) * 
      (variable[i] - variable[i - 1]) / (r[i] - r[i -1]);
  }
};

ostream & operator<< (ostream & out, const Profile & profile)
{
  for(size_t i = 0; i < profile.size(); i ++)
    {
      out << profile.r[i] << " " << profile.rho_baryon[i] << " " << profile.rho_total[i] << " " << profile.Sigma[i] << " " << profile.mass_baryon[i] << " "<< profile.mass_total[i] << " " << profile.Omega[i] << " " << profile.dOmega[i] << endl;
    }
  return out;
}

struct SplinedProfile : public Profile
{
  PowerLaw rho_baryon_spline, rho_total_spline, Sigma_spline, mass_baryon_spline, mass_total_spline, Omega_spline, dOmega_spline;

  SplinedProfile(const Profile & profile) : 
    Profile(profile), rho_baryon_spline(r, rho_baryon), rho_total_spline(r, rho_total), Sigma_spline(r, Sigma), 
    mass_baryon_spline(r, mass_baryon), mass_total_spline(r, mass_total), 
    Omega_spline(r, Omega), dOmega_spline(r, dOmega) { }
					   
};

double timeMigration(double cluster_mass, double r, const SplinedProfile & profile)
{
  const double ln_Lambda = 5.0;
  const double mass = profile.mass_total_spline(r);
  const double rho = profile.rho_total_spline(r);

  return 0.5 * sqrt(mass) * pow(r, 1.5) / (cluster_mass * ln_Lambda * sqrt(G)) + 0.125 * pow(mass, 1.5) / (cluster_mass * ln_Lambda * pow(r, 1.5) * sqrt(G) * M_PI * rho);
}

double timeDisruption(double cluster_mass, double r, const SplinedProfile & profile)
{
  const double radius_factor = 0.333;

  const double gamma = 0.62;

  const double t_0 = 0.1 / profile.Omega_spline(radius_factor * r); 

  const double t_4 = 1.355 * yr * pow(10.0, 4.0 * gamma) * (1.0 / gamma) * pow(t_0 / yr, 0.967);

  const double t_dis = t_4 * pow(cluster_mass / (1.0e4 * M_sun), gamma);
      
  return t_dis;
}

int func (double t, const double y[], double f[],
	  void *params)
{
  const pair<double, SplinedProfile const *> & data = 
    * static_cast<pair<double, SplinedProfile const *> *> (params);

  const double r = y[0];
  const double mass = y[1];

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

vector<Element> decayOrbit(double cluster_mass_init, double r_init, const SplinedProfile & profile)
{
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk4;

  gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
     
  const pair<double, SplinedProfile const *> data(cluster_mass_init, &profile);

  gsl_odeiv_system sys = {func, 0, 2, const_cast<pair<double, const SplinedProfile *> *> (&data)};
     
  double t = 0.0;
    
  double h = yr;
  
  double y[2] = { r_init, cluster_mass_init }, y_err[1];

  double dydt_in[2], dydt_out[2];
     
  GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);
     
  vector<Element> dM_dr;

  while (true)
    {
      int status = gsl_odeiv_step_apply (s, t, h, 
					 y, y_err, 
					 dydt_in, 
					 dydt_out, 
					 &sys);
     
      if (status != GSL_SUCCESS) abort();
     
      dydt_in[0] = dydt_out[0];
      dydt_in[1] = dydt_out[1];
                
      t += h;
     
      h = min<double> (- 0.1 * y[0] / dydt_in[0], - 0.1 * y[1] / dydt_in[1]);

      const double r = y[0];
      const double cluster_mass = y[1];

      dM_dr.push_back(Element(cluster_mass, r, dydt_out[1] / dydt_out[0]));
      
      if(r < 2.0 * profile.r.front()) break;

      //if(cluster_mass < 100.0 * M_sun) break;
      if(cluster_mass < 10000.0 * M_sun) break;
    }
     
  gsl_odeiv_step_free (s);
 
  return dM_dr;
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
    //    const double c_s = 1.0e6;
    double sum_mass = 0.0;

    while(sum_mass < profile.getMassBaryon()) 
      {
	const Cluster cluster(genMass(M_min, M_max, slope), profile.randomRadius());

	sum_mass += cluster.mass;

// 	if(cluster.mass < profile.getMassBaryon() - profile.interp(profile.mass, cluster.r) &&
//  	   c_s * profile.interp(profile.Omega, cluster.r) / (M_PI * G * profile.interp(profile.Sigma, cluster.r)) < 1.0 &&
//  	   cluster.mass < pow(c_s, 4.0) / (G * G * profile.interp(profile.Sigma, cluster.r)) &&
//  	   (M_PI / 6.0) * c_s * c_s * c_s / (sqrt(G * G * G) * sqrt(profile.interp(profile.rho, cluster.r))))
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

int main()
{
  const double R_e = 500.0 * pc;
  const double n = 1.5;
  const double mass = 1.0e10 * M_sun;

  const Model model(R_e, n, mass);

  const double r_min = pc;
  const double r_max = 50.0 * kpc;
  
  const double mass_halo = 1.0e11 * M_sun;
  const double concentration = 10;
  const double r_vir = 50 * kpc;

  const double M_min = 10000.0 * M_sun; //100.0 * M_sun;
  const double M_max = 1.0e6 * M_sun; //2.0e5 * M_sun;

  cout << "# R_e, n, mass, r_min, r_max, mass_halo, concentration, r_vir: " << R_e / pc << ", " << n << ", " << mass / M_sun << ", " << r_min / pc << ", " << r_max / pc << ", " << mass_halo / M_sun << ", " << concentration << ", " << r_vir / pc << endl; 

  cout << "# M_min, M_max: " << M_min / M_sun << ", " << M_max / M_sun << endl;

  const NFW nfw(mass_halo, r_vir, concentration);

  const Profile profile(model, nfw, r_min, r_max, 1000);
  
  // for(size_t i = 0; i < 100; i ++) cout << profile.randomRadius() / pc << endl;

  ofstream profile_out("profile.dat");
  profile_out << profile;


//   double total_mass_baryon_Sigma = 0.0;
//   double total_mass_baryon_rho = 0.0;
//   double total_mass_NFW_rho = 0.0;
//   for (size_t i = 0; i < profile.r.size() - 1; i ++)
//     {
//       const double r1 = profile.r[i];
//       const double r2 = profile.r[i+1];

//       const double int_1_sigma = two_pi * r1 * profile.model.Sigma(r1);
//       const double int_2_sigma = two_pi * r2 * profile.model.Sigma(r2);

//       const double int_1_rho = four_pi * sq(r1) * profile.model.rho(r1);
//       const double int_2_rho = four_pi * sq(r2) * profile.model.rho(r2);

//       const double int_1_rho_NFW = four_pi * sq(r1) * profile.nfw.rho(r1);
//       const double int_2_rho_NFW = four_pi * sq(r2) * profile.nfw.rho(r2);

//       //cout << profile.model.Sigma(r1) << " " << profile.model.Sigma(r2) << " " << profile.model.rho(r1) << " " << profile.model.rho(r2) << " " << profile.nfw.rho(r1) << " " << profile.nfw.rho(r2) << endl;

//       total_mass_baryon_Sigma += 0.5 * (r2 - r1) * (int_1_sigma + int_2_sigma);
//       total_mass_baryon_rho   += 0.5 * (r2 - r1) * (int_1_rho + int_2_rho);
//       total_mass_NFW_rho      += 0.5 * (r2 - r1) * (int_1_rho_NFW + int_2_rho_NFW);
//       //cout << total_mass_baryon_rho / M_sun << " " << total_mass_NFW_rho / M_sun << endl;
      
//     }
//   cout << total_mass_baryon_Sigma / M_sun << " " << total_mass_baryon_rho / M_sun << " " << total_mass_NFW_rho / M_sun << endl;
//   abort();

//   const double factor = 2.0;
//   for (size_t i = 0; i < profile.r.size(); i ++)
//     {
//       const double radius = profile.r[i];
//       const double mass_enclosed = profile.mass_baryon[i] + profile.nfw.M(radius);
//       const double velocity_dispersion = sqrt(G * mass_enclosed / (factor * radius));
//       cout << radius / pc << " " << velocity_dispersion / 1.0e5 << endl;
//     }
//   abort();
  
  const SplinedProfile splined_profile(profile);

  const Clusters clusters(M_min, M_max, 2.0, profile);

  cout << "# generated " << clusters.size() << " clusters" << endl;
  cout << "# with total mass " << clusters.mass() / M_sun << endl;

  double mass_nucleus = 0.0;
  vector<double> final_mass(profile.r.size(), 0.0);
  
  for(size_t i = 0; i < clusters.size(); i ++)
    {

      vector<Element> dM_dr = decayOrbit(clusters[i].mass, clusters[i].r, splined_profile);

      reverse(dM_dr.begin(), dM_dr.end());

      if(dM_dr.front().r < 2.0 * profile.r.front()) mass_nucleus += dM_dr.front().mass;

      //      cout << dM_dr.back().r / pc << " " << dM_dr.back().mass / M_sun << " " << clusters[i].mass / M_sun << endl;

    //   const size_t initial_index = size_t(floor(double(profile.r.size()) * log(dM_dr.front().r / profile.r.front()) / log(profile.r.back() / profile.r.front())));

//       for(size_t k = initial_index; k < profile.r.size(); k ++)
// 	{
// 	  final_mass[k] += dM_dr.front().mass;
// 	}

      for(size_t j = 1; j < dM_dr.size(); j ++)
	{
	  const double r = 0.5 * (dM_dr[j - 1].r + dM_dr[j].r);

	  const double delta_mass = (dM_dr[j].r - dM_dr[j - 1].r) * 0.5 * (dM_dr[j].dM_dr + dM_dr[j - 1].dM_dr);

	  const size_t index = size_t(floor(double(profile.r.size()) * log(r / profile.r.front()) / log(profile.r.back() / profile.r.front())));
	  
	  for(size_t k = index; k < profile.r.size(); k ++)
	    {
	      final_mass[k] += delta_mass;
	    }
	}
    }

  vector<double> final_rho(profile.r.size());
  for(size_t i = 1; i < final_mass.size(); i ++)
    {
      final_rho[i] = ((final_mass[i] - final_mass[i - 1]) / (profile.r[i] - profile.r[i - 1])) / (4.0 * M_PI * sq(0.5 * (profile.r[i] + profile.r[i - 1])));
    }
  final_rho[0] = final_rho[1];

  const double r_nucleus = 0.5 * pc;
  for(size_t i = 0; i < final_mass.size(); i ++)
    {
      final_rho[i] += (mass_nucleus / (240.0 * r_nucleus * r_nucleus * r_nucleus)) * exp(- sqrt(profile.r[i] / r_nucleus));

	// exp(-0.5 * sq(profile.r[i] / r_nucleus)) / (cub(sqrt(2.0 * M_PI) * r_nucleus));
    }

  const vector<double> final_sigma_baryon = Project(profile.r, final_rho);

  cout << "# mass_nucleus " << mass_nucleus / M_sun << " " << mass_nucleus / clusters.mass() << endl;
  
  //cout << "# " << mass_nucleus / clusters.mass() << endl;

  vector<double> NFW_rho(final_mass.size());
  for(size_t i = 0; i < final_mass.size(); i ++)
    {
      NFW_rho[i] = profile.nfw.rho(profile.r[i]);
    }  
  const vector<double> NFW_Sigma = Project(profile.r, NFW_rho);

  for(size_t i = 0; i < final_mass.size(); i ++)
    {
      //      cout << profile.r[i] << " " << profile.mass[i] << " " << final_mass[i] << endl;

      cout << profile.r[i] << " " << profile.rho_baryon[i] << " " << final_rho[i] << " " << profile.nfw.rho(profile.r[i]) << " " << profile.Sigma[i] << " " << final_sigma_baryon[i] << " " << NFW_Sigma[i] << endl;
    }

  return 0;
}
