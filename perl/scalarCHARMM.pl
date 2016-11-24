#!/usr/bin/env perl

# get scalar quantity from CHARMM
#
# Michael Feig, 2004, MSU

sub usage {
  printf STDERR "usage:   scalarCHARMM.pl [options] name [PDBfile]\n";
  printf STDERR "options: [-par CHARMMoptions]\n";
  printf STDERR "         [-psf PSFfile CRDfile]\n";
  printf STDERR "         [-mol2 MOL2file]\n";
  printf STDERR "         name: charge, radius, mass, type\n";
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

my $inpfile="-";
my $base="";
my $name="";
my $psffile;
my $crdfile;
my $mol2file;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
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
    $name=(shift @ARGV);
    $inpfile=(shift @ARGV);
    $done=1;
  }
}

my $charmm=&CHARMM::new();
$charmm->loadParameters(%par);

if (defined $psffile) {
  $charmm->setupFromPSF($psffile,$crdfile);
} elsif (defined $mol2file) {
  $charmm->setupFromMol2($mol2file);
} else {
  $charmm->setupFromPDB($inpfile);
}

my $res=$charmm->getScalar($name);
$charmm->finish();

foreach my $n (@{$res}) {
  printf "%f\n",$n;
}

exit 0;

