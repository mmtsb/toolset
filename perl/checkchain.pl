#!/usr/bin/env perl

# checks a lattice chain to be valid MONSSTER input
# 
# http://mmtsb.scripps.edu/doc/checkchain.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage: checkchain.pl [-v] [chainFile]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use SICHO;

&usage() if (@ARGV && ($ARGV[0] eq "-help" || $ARGV[0] eq "-h"));

my $verbose=0;
if (@ARGV && $ARGV[0] eq "-v") {
  shift @ARGV;
  $verbose=1;
}

my $chain=SICHO::new();
$chain->readChain(shift @ARGV);
my $violate=$chain->checkMONSSTER($verbose);

if ($violate>0) {
  printf "%d violation(s) found in chain\n",$violate;
} elsif ($violate<0) {
  printf "No MONSSTER chain available\n";
} else {
  printf "Valid MONSSTER chain\n";
}


