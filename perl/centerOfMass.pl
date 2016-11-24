#!/usr/bin/env perl

# calculates the center of mass for a protein structure
#
# 2005, Michael Feig, MSU

sub usage {
  printf STDERR "usage: centerOfMass.pl [pdbFile]\n";
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

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $fname = shift @ARGV;
    $done=1;
  }
}

my $mol=Molecule::new();
$mol->readPDB($fname);

my ($cx,$cy,$cz)=$mol->centerOfMass();

printf STDOUT "%f %f %f\n",$cx,$cy,$cz;

