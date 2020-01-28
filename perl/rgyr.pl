#!/usr/bin/env perl


# calculates the radius of gyration for a protein structure
#
# http://mmtsb.scripps.edu/doc/rgyr.pl.html
# 2000, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage: rgyr.pl -cofm x y z -caonly [pdbFile]\n";
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

my $fname;

my $cofmx=undef;
my $cofmy=undef;
my $cofmz=undef;

my $caonly=0;
my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-caonly") {
    shift @ARGV;
    $caonly=1;
  } elsif ($ARGV[0] eq "-cofm") {
    shift @ARGV;
    $cofmx=shift @ARGV;
    $cofmy=shift @ARGV;
    $cofmz=shift @ARGV;
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

printf STDOUT "%f\n",&Analyze::radiusOfGyration($mol,$caonly,$cofmx,$cofmy,$cofmz);

