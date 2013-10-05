#!/usr/local/bin/perl -w

# NOTE: Replace 'tables_RunX_NAME'
#       Replace ../Data/RunX_NAME/

# Specify model
$galaxy = $ARGV[0];
$redshift = $ARGV[1];

$M_sun = 2.0e33;
$yr = 3.14e7;

if ($galaxy eq "disk") 
{
    $R_e = sprintf("%.1f", 1000.0);
    $n = sprintf("%.1f", 1.0);
    $mass = sprintf("%.1e", 5.0e9);
    $mass_halo = sprintf("%.1e", 5.0e10);
} 
else 
{
    $R_e = sprintf("%.1f", 500.0);
    $n = sprintf("%.1f", 1.5);
    $mass = sprintf("%.1e", 1.0e9);
    $mass_halo = sprintf("%.1e", 1.0e10);
}

if (($redshift eq "high") && ($mass_halo == 5.0e10))
{
    $concentration = sprintf("%.1f", 3.5);
    $r_vir = sprintf("%.1f", 20.0);
}
if (($redshift eq "low") && ($mass_halo == 5.0e10))
{
    $concentration = sprintf("%.1f", 15.0);
    $r_vir = sprintf("%.1f", 100.0);
}
if (($redshift eq "high") && ($mass_halo == 1.0e10))
{
    $concentration = sprintf("%.1f", 3.5);
    $r_vir = sprintf("%.1f", 10.0);
}
if (($redshift eq "low") && ($mass_halo == 1.0e10))
{
    $concentration = sprintf("%.1f", 15.0);
    $r_vir = sprintf("%.1f", 60.0);
}

$outfile = 'tables_RunX_NAME';
open(OFILE, ">> $outfile");

print OFILE "\\multirow\{4\}\{\*\}\{$galaxy $redshift-z\}\n";

for ($Ms = 1.0e4; $Ms < 1.0e8; $Ms=$Ms*10) 
{
    $M_max = sprintf("%.1e", $Ms);

    $infile = 'center_'.$galaxy.'_'.$M_max.'_'.$redshift.'_z.out';
    open(IFILE, "< ../Data/RunX_NAME/$infile");
    @ilines = <IFILE>;
    close(IFILE);
 
    $length = @ilines;
    @last = split(/ /, $ilines[$length - 1]);
    $last_mass = $last[2];
    $last_mass_M_sun = $last_mass / $M_sun;
    $half_mass = 0.5 * $last_mass;

    foreach $iline (@ilines) 
    {
	@this_line = split(/ /, $iline);
	$this_mass = $this_line[2];

	if ($this_mass < $half_mass)
	{
	 $half_time = $this_line[0];  
	}
	else
	{
	    last;
	}

    }
    $log_M_max = sprintf("%d", log($M_max));
    $log_half_time_yr = sprintf("%.2f", log($half_time / $yr) / log(10));
    $log_last_mass_M_sun = sprintf("%.2f", log($last_mass_M_sun) / log(10));
    $log_last_mass_over_mass = sprintf("%.2f", log($last_mass_M_sun / $mass) / log(10));
    $log_last_mass_over_M_max = sprintf("%.2f", log($last_mass_M_sun / $M_max) / log(10));
 
    print OFILE "\& \$$log_M_max\$ \& \$$log_half_time_yr\$ \& \$$log_last_mass_M_sun\$ \& \$$log_last_mass_over_mass\$ \& \$$log_last_mass_over_M_max\\\\\n";
}

close(OFILE);

