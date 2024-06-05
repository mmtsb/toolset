#!/usr/bin/env perl

# analyzes a protein structure with an external function
#

sub usage {
  printf STDERR "usage: analyze.pl -function file [pdbFile]\n";
  printf STDERR "options: [-readseg]\n";
  printf STDERR "         [-splitseg]\n";
  printf STDERR "         [-par a[:b...]]\n";
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
my $ignoreseg=1;
my $splitseg=0;
my @parameters=();

while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-function") {
    shift @ARGV;
    $function=shift @ARGV;
  } elsif ($ARGV[0] eq "-readseg") {
    shift @ARGV;
    $ignoreseg=0;
  } elsif ($ARGV[0] eq "-splitseg") {
    shift @ARGV;
    $ignoreseg=0;
    $splitseg=1;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    @parameters=split(/:/,shift @ARGV);
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
$mol->readPDB($fname,ignoreseg=>$ignoreseg,splitseg=>$splitseg);

&start if (defined &start);
my @data=&analyze($mol,@parameters);

printf STDOUT "%s\n",join(" ",@data) if (defined $data[0]);

&end if (defined &end);

