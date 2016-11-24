#!/usr/bin/env perl

# runs MONSSTER program to perform lattice-based 
# protein simulations
#
# http://mmtsb.scripps.edu/doc/latticesim.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   latticesim.pl [options] seqFile\n";
  printf STDERR "options: [-chain file | -rnd | -pdb file] [-g gridsize]\n";
  printf STDERR "         [-par tsteps=val,ncycle=val,icycle=val,\n";
  printf STDERR "               stiff=val,short=val,central=val,kdcore=val]\n";
  printf STDERR "         [-l force min:max[=min:max ...]]\n";
  printf STDERR "         [-d force res1:res2[=res1:res2 ...]]\n";
  printf STDERR "         [-sa temp] [-const temp]\n";
  printf STDERR "         [-norun]\n";
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

my ($chainfile, $pdbfile);
my $rndchain=1;
my $gridsize;
my ($fraglist,$limforce);
my ($drestforce, $drestlist);
my %simpar;
my $seqfile;
my $run=1;

PARLOOP:
while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-chain") {
    shift @ARGV;
    $chainfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-rnd") {
    shift @ARGV;
    $rndchain=1;
  } elsif ($ARGV[0] eq "-pdb") {
    shift @ARGV;
    $pdbfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-g") {
    shift @ARGV;
    $gridsize=shift @ARGV;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $limforce=shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-d") {
    shift @ARGV;
    $drestforce=shift @ARGV;
    $drestlist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-sa") {
    shift @ARGV;
    $simpar{temp}=(shift @ARGV).":1.0"
  } elsif ($ARGV[0] eq "-const") {
    shift @ARGV;
    $simpar{temp}=(shift @ARGV);
    $simpar{temp}.=":$simpar{temp}";
    $simpar{central}=0.5;
    $simpar{tsteps}=1;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      $simpar{$key}=$val;
    }
  } elsif ($ARGV[0] eq "-norun") {
    shift @ARGV;
    $run=0;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    $seqfile=shift @ARGV;
    last PARLOOP;
  }
}

my $seq=Sequence::new();
$seq->readMONSSTER($seqfile);
my $nres=$#{$seq->{sequence}}+1;

my $inpchain=SICHO::new(gridsize=>$gridsize);
if (defined $chainfile) {
  $inpchain->readChain($chainfile);
} elsif (defined $pdbfile) {
  my $mol=Molecule::new($pdbfile);
  $mol->selectChain("");
  $inpchain->genMONSSTERFromAllAtom($mol,fraglist=>$fraglist);
} else {
  $inpchain->genRandomMONSSTER($#{$seq->{sequence}}+1);
}
my $nchain=$#{$inpchain->{sidechain}}+1;

die "sequence and chain lengths do not match (seq: $nres, chain: $nchain)"
  if ($nres != $nchain-2);

my $monsster=MONSSTER::new($seq);
$monsster->setParameters(%simpar);

if (defined $fraglist) {
  $seq->setValidResidues($fraglist,1);
  $monsster->
    setPositionalRestraints(&GenUtil::gradForceList($seq->listFromValid(),$limforce));
}

$monsster->setDistanceRestraints($drestforce, $drestlist) 
  if (defined $drestforce && defined $drestlist);

if ($run) {
  $monsster->run($inpchain);
} else {
  $monsster->setup($inpchain);
}

