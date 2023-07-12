#!/usr/bin/env perl

# get CHARMM energy from PDB file 
#
# http://mmtsb.scripps.edu/doc/enerCHARMM.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   enerCHARMM.pl [options] [PDBfile]\n";
  printf STDERR "options: [-out total,bonds,angles,ureyb,dihedrals,impropers,\n";
  printf STDERR "               vdwaals,elec,gb,sasa,asp]\n";
  printf STDERR "         [-all]\n";
  printf STDERR "         [-oneline]\n";
  printf STDERR "         [-noout]\n";
  printf STDERR "         [-charge] [-weight] [-volume]\n";
  printf STDERR "         [-dipole] [-mass] [-oxyz]\n";
  printf STDERR "         [-dist Selection Selection]\n";
  printf STDERR "         [-angle Selection Selection Selection]\n";
  printf STDERR "         [-dihedral Selection Selection Selection Selection]\n";
  printf STDERR "         [-psf PSFfile CRDfile]\n";
  printf STDERR "         [-mol2 MOL2file]\n";
  printf STDERR "         [-par CHARMMparams]\n";
  printf STDERR "         [-l min:max[=...]] [-self]\n";
  printf STDERR "         [-log logFile] [-cmd logFile]\n";
  printf STDERR "         [-custom file]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use GenUtil;
use Molecule;
use CHARMM;

my %par;

my $logFile;
my $cmdlog;

my $oneline=0;

my $inpfile="-";
my $needsasa=0;

my @olist;
push(@olist,"total");

my $sellist;
my $selfe=0;

my $psffile;
my $crdfile;
my $mol2file;
my $all=0;

my $charge=0;
my $weight=0;
my $dipole=0;
my $volume=0;
my $dist=0;
my $angle=0;
my $dihedral=0;
my $mass=0;
my $oxyz=0;

my $selection1;
my $selection2;
my $selection3;
my $selection4;

my $customfile;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=(shift @ARGV);
  } elsif ($ARGV[0] eq "-oneline") {
    shift @ARGV;
    $oneline=1;
  } elsif ($ARGV[0] eq "-all") {
    shift @ARGV;
    $all=1;
  } elsif ($ARGV[0] eq "-charge") {
    shift @ARGV;
    $charge=1;
  } elsif ($ARGV[0] eq "-volume") {
    shift @ARGV;
    $volume=1;
  } elsif ($ARGV[0] eq "-weight") {
    shift @ARGV;
    $weight=1;
  } elsif ($ARGV[0] eq "-dipole") {
    shift @ARGV;
    $dipole=1;
  } elsif ($ARGV[0] eq "-dist") {
    shift @ARGV;
    $selection1=shift @ARGV;
    $selection2=shift @ARGV;
    $dist=1;
  } elsif ($ARGV[0] eq "-angle") {
    shift @ARGV;
    $selection1=shift @ARGV;
    $selection2=shift @ARGV;
    $selection3=shift @ARGV;
    $angle=1;
  } elsif ($ARGV[0] eq "-dihedral") {
    shift @ARGV;
    $selection1=shift @ARGV;
    $selection2=shift @ARGV;
    $selection3=shift @ARGV;
    $selection4=shift @ARGV;
    $dihedral=1;
  } elsif ($ARGV[0] eq "-mass") {
    shift @ARGV;
    $mass=1;
  } elsif ($ARGV[0] eq "-oxyz") {
    shift @ARGV;
    $oxyz=1;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $sellist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-self") {
    shift @ARGV;
    $selfe=1;
  } elsif ($ARGV[0] eq "-out") {
    shift @ARGV;
    my $line=shift @ARGV;
    $needsasa=1 if ($line=~/sasa/);
    @olist=split(/,/,$line);
  } elsif ($ARGV[0] eq "-noout") {
    shift @ARGV;
    @olist=();
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-psf") {
    shift @ARGV;
    $psffile=shift @ARGV;
    $crdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-mol2") {
    shift @ARGV;
    $mol2file=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    $inpfile=(shift @ARGV);
    $done=1;
  }
}

my $charmm=&CHARMM::new($logFile,$cmdlog);

$charmm->loadParameters(%par);

if (defined $psffile) {
  $charmm->setupFromPSF($psffile,$crdfile);
} elsif (defined $mol2file) {
  $charmm->setupFromMol2($mol2file);
} else {
  $charmm->setupFromPDB($inpfile);
}

$charmm->setupEnergy();
$charmm->noeRestraints();

if (defined $customfile) {
  foreach my $c ( split(/:/,$customfile))  {
    if (&GenUtil::checkFile($c)) {
      my $custom=&GenUtil::readData(&GenUtil::getInputFile($c));
      $charmm->stream($custom);
    }
  }
}

if ($charge) {
  my $chg=$charmm->getTotalCharge();

  $charmm->finish();

  printf "%10.5f\n",$chg;
} elsif ($weight) {
  my $wgt=$charmm->getTotalWeight();
  $charmm->finish();
  printf "%12.5f amu = %10.5e kg\n",$wgt,$wgt*1.66053886E-27;
} elsif ($volume) {
  my ($vol1,$vol2)=$charmm->analyzeMolecularVolume();
  $charmm->finish();
  printf "%lf %lf\n",$vol1,$vol2;
} elsif ($dipole) {
  my ($dipole,$dipx,$dipy,$dipz)=$charmm->getDipole($mass,$oxyz);
  $charmm->finish(); 
  printf "dipole: %10.5f D, x: %10.5f D, y: %10.5f D, z: %10.5f D\n",
    $dipole,$dipx,$dipy,$dipz;
} elsif ($dist) {
  my $sel1="SL1";
  $charmm->defineSelection($selection1,$sel1);
  my $sel2="SL2";
  $charmm->defineSelection($selection2,$sel2);
  my $result=$charmm->analyzeDistance(selection1=>$sel1,selection2=>$sel2,mass=>$mass);
  printf "distance: %lf\n",$result;
} elsif ($angle) {
  my $sel1="SL1";
  $charmm->defineSelection($selection1,$sel1);
  my $sel2="SL2";
  $charmm->defineSelection($selection2,$sel2);
  my $sel3="SL3";
  $charmm->defineSelection($selection3,$sel3);
  my $result=$charmm->analyzeAngle(selection1=>$sel1,selection2=>$sel2,selection3=>$sel3,mass=>$mass);
  printf "angle: %lf\n",$result;
} elsif ($dihedral) {
  my $sel1="SL1";
  $charmm->defineSelection($selection1,$sel1);
  my $sel2="SL2";
  $charmm->defineSelection($selection2,$sel2);
  my $sel3="SL3";
  $charmm->defineSelection($selection3,$sel3);
  my $sel4="SL4";
  $charmm->defineSelection($selection4,$sel4);
  my $result=$charmm->analyzeDihedral(selection1=>$sel1,selection2=>$sel2,selection3=>$sel3,selection4=>$sel4,mass=>$mass);
  printf "dihedral: %lf\n",$result;
} elsif ($#olist>=0) {
  my $ener=$charmm->getEnergy($sellist,$selfe);

  if ($needsasa) {
   $charmm->solvAccessSurf();
   my $sasa=$charmm->getSASAOutput();
   $ener->{sasa}=$sasa->{energy};
  }

  $charmm->finish();
  if ($all) {
    foreach my $k ( sort keys %{$ener} ) {
      printf "%-20s %15.4f\n",$k,$ener->{$k};
    }
  } else {  
    foreach my $o ( @olist ) {
    my $sum=0.0;
    my @addlist=($o=~/([+-]*)([^+-]+)/g);
    while (@addlist) {
      my $sgn=shift @addlist;
      my $a=shift @addlist;
      my $mult=1.0;
      if (defined $ener->{$a}) {
	if ($sgn eq "-") {
	  $mult=-1;
	}
	$sum+=$ener->{$a}*$mult;
      }
    }
    if ($#olist>0) {
      if ($oneline) {
	printf "%15.4f ",$sum;
      } else {
	printf "%-20s %15.4f\n",$o,$sum;
      }
    } else {
      printf "%15.4f\n",$sum;
    }
  }

  printf "\n" if ($oneline && $#olist>0);
  }
} else {
  $charmm->finish();
}

exit 0;

