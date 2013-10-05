#!/usr/local/bin/perl -w

$N = 1000;
$num_pts = 100;

$galaxy = $ARGV[0];
$redshift = $ARGV[1];

$code =  'main_'.$galaxy.'_unordered.cc';
open(CFILE, "<$code");
@clines = <CFILE>;
close(CFILE);

open(MFILE, "<Makefile");
@mlines = <MFILE>;
close(MFILE);

if ($galaxy eq "disk") {
    $R_e = sprintf("%.1f", 1000.0);
    $n = sprintf("%.1f", 1.0);
    $mass = sprintf("%.1e", 5.0e9);
    $mass_halo = sprintf("%.1e", 5.0e10);
} else {
    $R_e = sprintf("%.1f", 500.0);
    $n = sprintf("%.1f", 1.5);
    $mass = sprintf("%.1e", 1.0e9);
    $mass_halo = sprintf("%.1e", 1.0e10);
}

if (($redshift eq "high") && ($mass_halo == 5.0e10)){
    $concentration = sprintf("%.1f", 3.5);
    $r_vir = sprintf("%.1f", 20.0);
}
if (($redshift eq "low") && ($mass_halo == 5.0e10)){
    $concentration = sprintf("%.1f", 15.0);
    $r_vir = sprintf("%.1f", 100.0);
}
if (($redshift eq "high") && ($mass_halo == 1.0e10)){
    $concentration = sprintf("%.1f", 3.5);
    $r_vir = sprintf("%.1f", 10.0);
}
if (($redshift eq "low") && ($mass_halo == 1.0e10)){
    $concentration = sprintf("%.1f", 15.0);
    $r_vir = sprintf("%.1f", 60.0);
}

$M_min = sprintf("%.1f", 100.0);
$r_min = sprintf("%.1f", 1.0);
$r_max = sprintf("%.1f", 50.0);

$g = substr($galaxy, 0, 1);
$r = substr($redshift, 0, 1);

for ($Ms = 1.0e4; $Ms < 1.0e8; $Ms=$Ms*10) 
{
    $M_max = sprintf("%.1e", $Ms);

    $e = substr($M_max, 6, 1);
    $exe = $g.$r.$e.'_4';

    $outfile = 'main_'.$galaxy.'_unordered_'.$redshift.'_'.$M_max.'.cc';
    open(OFILE, ">$outfile");  

    $center_out = 'center_'.$galaxy.'_'.$M_max.'_'.$redshift.'_z.out';
    $center_yes_out = 'center_yes_'.$galaxy.'_'.$M_max.'_'.$redshift.'_z.out';
    $center_no_out = 'center_no_'.$galaxy.'_'.$M_max.'_'.$redshift.'_z.out';
    
    $line_number = 0;
    foreach $line (@clines) 
    {
	if ($line_number < 50) 
	{
	    if ($line =~ /const double R_e =/) {
		$line = "const double R_e = $R_e * pc;\n";
	    }
	    if ($line =~ /const double n =/) {
		$line = "const double n = $n;\n";
	    }
	    if ($line =~ /const double mass =/) {
		$line = "const double mass = $mass * M_sun;\n";
	    }
	    if ($line =~ /const double r_min =/) {
		$line = "const double r_min = $r_min * pc;\n";
	    }
	    if ($line =~ /const double r_max =/) {
		$line = "const double r_max = $r_max * kpc;\n";
	    }
	    if ($line =~ /const double mass_halo =/) {
		$line = "const double mass_halo = $mass_halo * M_sun;\n";
	    }
	    if ($line =~ /const double concentration =/) {
		$line = "const double concentration = $concentration;\n";
	    }
	    if ($line =~ /const double r_vir =/) {
		$line = "const double r_vir = $r_vir * kpc;\n";
	    }
	    if ($line =~ /const double M_min =/) {
		$line = "const double M_min = $M_min * M_sun;\n";
	    }
	    if ($line =~ /const double M_max =/) {
		$line = "const double M_max = $M_max * M_sun;\n";
	    }
	    if ($line =~ /size_t N =/) {
		$line = "size_t N = $N;\n";
	    }
	    if ($line =~ /size_t num_pts =/) {
		$line = "size_t num_pts = $num_pts;\n";
	    }
	}

	if ($line =~ /ofstream center_out\(/) {
	    $line = "ofstream center_out(\"../Data/Run4_Unordered/$center_out\");\n";
	}
	if ($line =~ /ofstream center_yes_out\(/) {
	    $line = "ofstream center_yes_out(\"../Data/Run4_Unordered/$center_yes_out\");\n"; 
	} 
	if ($line =~ /ofstream center_no_out\(/) {
	    $line = "ofstream center_no_out(\"../Data/Run4_Unordered/$center_no_out\");\n";  
	}

	$line =~ s/\r\n/\n/g;
	#$line =~ s/e\+0/e/g;
	#$line =~ s/e\+/e/g;
	print OFILE $line;
	$line_number++;
    }
    close(OFILE);

    open(OOFILE, ">Makefile");
    foreach $line (@mlines) {
	if ($line =~ /main_/) {
	    $line = " $outfile\n";
	}
	if ($line =~ /RESULT =/) {
	    $line = "RESULT = $exe\n";
	}
	$line =~ s/\r\n/\n/g;
	print OOFILE $line;
    }
    close (OOFILE);

    $datafile = 'profile_'.$galaxy.'_unordered_'.$mass.'_'.$mass_halo.'_'.$M_max.'_'.$redshift.'_z_low_res.dat';
    print `make\n`;
    print `$exe > ../Data/Run4_Unordered/$datafile &\n`;
    print `cp -p $outfile ../Data/Run4_Unordered/ &\n`;  
}



