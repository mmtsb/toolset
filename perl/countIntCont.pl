#!/usr/bin/env perl

# count the number of heavy-atom contacts at the interface between chains
#
# http://mmtsb.scripps.edu/doc/countIntCont.html
# 2003, John Karanicolas, Baker group, UW

sub usage {
  printf STDERR "usage:   countIntCont.pl [options] [PDBfile]\n";
  printf STDERR "options: [-distcut num] [-chainA id] [-chainB id]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;
use Molecule;

my $fname="-";
my $distcut=5;
my $chainA;
my $chainB;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-distcut") {
    shift @ARGV;
    $distcut=shift @ARGV;
  } elsif ($ARGV[0] eq "-chainA") {
    shift @ARGV;
    $chainA=shift @ARGV;
  } elsif ($ARGV[0] eq "-chainB") {
    shift @ARGV;
    $chainB=shift @ARGV;
  } else {
    $fname = shift @ARGV;
  }
}

my $mol=Molecule::new($fname);
my $contacts=$mol->countInterfaceContacts($distcut,$chainA,$chainB);

if ($#{$contacts} >= 0) {
    foreach my $nc ( @{$contacts} ) {
	printf STDOUT "%d atomic contacts between chains %s and %s\n",
   	  $nc->{AAcontacts}, $nc->{cidA}, $nc->{cidB};
    }
}

