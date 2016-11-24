#!/usr/bin/env perl

# calculates the qscore for a protein structure with respect
# to a reference
#
# http://mmtsb.scripps.edu/doc/qscore.pl.html
# 2000, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   qscore.pl [options] refPDB [cmpPDB]\n";
  printf STDERR "options: [-all]\n";
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
use Analyze;

my ($pdb1,$pdb2);

my $all=0;
my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-all") {
    shift @ARGV;
    $all=1;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $pdb1 = shift @ARGV;
    $pdb2 = shift @ARGV;
    $done=1;
  }
}

my $refmol=Molecule::new($pdb1);
my $analyze=Analyze::new($refmol);

my $cmpmol=Molecule::new();
$cmpmol->readPDB($pdb2);
my $qsc=$analyze->qscore($cmpmol);

if ($all) {
  printf STDOUT "all:     %f\n",$qsc->{all};
  printf STDOUT " short:  %f\n",$qsc->{short};
  printf STDOUT " medium: %f\n",$qsc->{medium};
  printf STDOUT " long:   %f\n",$qsc->{long};
} else {
  printf STDOUT "%f\n",$qsc->{all};
}

