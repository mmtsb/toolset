#!/usr/bin/env perl

# runs MONSSTER program to evaluate lattice model energies
#
# http://mmtsb.scripps.edu/doc/latticeener.pl.html
# 2002, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   latticeener.pl [options] [seqFile chainFile | pdbFile]\n";
  printf STDERR "options: [-par stiff=val,short=val,central=val,kdcore=val,\n";
  printf STDERR "               softcore=val,pair=val,hbond=val,short=val,\n";
  printf STDERR "               burial=val,multibody=val,threebody=val,temp=val]\n";
  exit 1;
}

require 5.004;

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use GenUtil;
use Sequence;
use Molecule;
use SICHO;
use MONSSTER;

my %simpar=();
my $run=1;
my $file1;
my $file2;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-const") {
    shift @ARGV;
    $simpar{central}=0.5;
    $simpar{tsteps}=1;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      $simpar{$key}=$val;

    }
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    $file1=shift @ARGV if (!defined $file1);
    $file2=shift @ARGV if (!defined $file2);
  }
}

my $seq;
my $inpchain=SICHO::new();  

if (defined $file1 && defined $file2) {
  $seq=Sequence::new();
  $seq->readMONSSTER($file1);
  $inpchain->readChain($file2);
} else {
  my $mol=Molecule::new($file1);
  $mol->selectChain("");
  $seq=Sequence::new($mol);
  $seq->secFromDSSP($mol);
  $inpchain->genMONSSTERFromAllAtom($mol);
}  

my $nres=$#{$seq->{sequence}}+1;
my $nchain=$#{$inpchain->{sidechain}}+1;

my $monsster=MONSSTER::new($seq);
$monsster->setParameters(%simpar);
my ($energy,$eshort,$epair,$eburial,$ecorrect,$ecent)=$monsster->energy($inpchain);

printf "%f %f %f %f %f %f\n",$energy,$eshort,$epair,$eburial,$ecorrect,$ecent;

