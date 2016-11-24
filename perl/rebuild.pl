#!/usr/bin/env perl

# rebuilds all-atom PDB from SICHO chain
#
# http://mmtsb.scripps.edu/doc/rebuild.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   rebuild.pl [options] [ -pdb file | seq chain ]\n";
  printf STDERR "options: [-l refpdb min:max[=min:max...]] [-fixca]\n";
  printf STDERR "         [-r resolution] [-o offsetx offsety offsetz]\n";
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
use Sequence;
use SICHO;

my $resolution=1.45;
my $offsetx=50.0;
my $offsety=50.0;
my $offsetz=50.0;
my ($fraglist, $refpdb);
my ($inpchain, $inpseq);
my $chainpdb;
my $havepdb=0;
my $fixca;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-r") {
    shift @ARGV;
    $resolution=shift @ARGV;
  } elsif ($ARGV[0] eq "-o") {
    shift @ARGV;
    $offsetx=shift @ARGV;
    $offsety=shift @ARGV;
    $offsetz=shift @ARGV;
  } elsif ($ARGV[0] eq "-fixca") {
    shift @ARGV;
    $fixca=1;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $refpdb=shift @ARGV;
    $fraglist=shift @ARGV;
  } elsif ($ARGV[0] eq "-pdb") {
    shift @ARGV;
    $havepdb=1;
    $chainpdb=shift @ARGV;
    $resolution=0;
    $offsetx=$offsety=$offsetz=0;
  } elsif ($ARGV[0] =~ /^-.+/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    $inpseq=shift @ARGV;
    $inpchain=shift @ARGV;
  }
}

unless ((defined $inpseq && &GenUtil::checkFile($inpseq) && 
	 defined $inpchain && &GenUtil::checkFile($inpchain)) || $havepdb) {
  print STDERR "need sequence/chain or chain pdb file\n";
  &usage();
}

my $chain=SICHO::new(resolution=>$resolution, 
		     offsetx=>$offsetx, offsety=>$offsety, offsetz=>$offsetz);
my $sequence;

if ($havepdb) {
  my $chainmol=Molecule::new();
  $chainmol->readPDB($chainpdb);
  $chainmol->translate("GENERIC");
  $chain->fromMolecule($chainmol);
  $sequence=Sequence::new($chainmol);
} else {
  $chain->readChain($inpchain);
  $sequence=Sequence::new();
  $sequence->readMONSSTER($inpseq);
}

my $rebmol=Molecule::new();
$rebmol->rebuildFromSICHO($sequence,$chain,$fraglist,$refpdb,$fixca);
$rebmol->writePDB(\*STDOUT,ssbond=>0);

