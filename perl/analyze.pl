#!/usr/bin/env perl

# analyzes a protein structure with an external function
#

sub usage {
  printf STDERR "usage: analyze.pl -function file [pdbFile]\n";
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
use Math::Trig;

my $fname;

my $function;
my $done=0;

while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-function") {
    shift @ARGV;
    $function=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $fname = shift @ARGV;
    $done=1;
  }
}

die "please provide external function" if (!defined $function);
die "cannot read external function" if (!-r $function);

require "$function";

my $mol=Molecule::new();
$mol->readPDB($fname);

&start if (defined &start);
my @data=&analyze($mol);

printf STDOUT "%s\n",join(" ",@data) if (defined $data[0]);

&end if (defined &end);

