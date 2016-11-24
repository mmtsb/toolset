#!/usr/bin/env perl

# calculates the minimum distance between two residues in
# a given protein structure
#
# http://mmtsb.scripps.edu/doc/mindist.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage: mindist.pl residue1 residue2 pdbFile\n";
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

my $fname;
my $res1;
my $res2;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $res1=shift @ARGV;
    $res2=shift @ARGV;
    $fname=shift @ARGV;
    $done=1;
  }
}

my $mol=Molecule::new($fname);

my ($c1,$r1)=($res1=~/([A-Za-z]*)([0-9]+)/);
my ($c2,$r2)=($res2=~/([A-Za-z]*)([0-9]+)/);

my $ir=$mol->getResidue($r1,$c1);
my $jr=$mol->getResidue($r2,$c2);

if (!defined $ir || !defined $jr) {
  printf STDOUT "N/A\n";
} else {
  my $mind=$mol->minDistance($ir,$jr);
  printf STDOUT "%f\n",$mind;
}

