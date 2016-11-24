#!/usr/bin/env perl

# extracts lattice chain from MONSSTER trajectory
#
# http://mmtsb.scripps.edu/doc/xtracttra.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   xtracttra.pl tstep:cycle [trajectoryFile]\n";
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
use SICHO;
use Molecule;
use Sequence;

&usage() if ((@ARGV && ($ARGV[0] eq "-help" || $ARGV[0] eq "-h")) || !(@ARGV));

my $pdb=0;
my $seqfile;
if ((@ARGV && ($ARGV[0] eq "-pdb"))) {
  $pdb=1;
  shift @ARGV;
  $seqfile=shift @ARGV;
}

my @f=split(/:/,shift @ARGV);
my $tstp=$f[0];
my $cyc=$f[1];

my $trafile=shift @ARGV;

my $sicho=SICHO::new();
$sicho->readMONSSTERTraj($trafile,$tstp,$cyc);

die "cannot find conformation" 
  if ($#{$sicho->{sidechain}}<0);

if ($pdb) {
  my $outmol=Molecule::new();
  my $seq;
  $seq=Sequence::new();
  $seq->readMONSSTER($seqfile);
  $sicho->{resolution}=1.0;
  $sicho->{offset}->{xcoor}=0.0;
  $sicho->{offset}->{ycoor}=0.0;
  $sicho->{offset}->{zcoor}=0.0;
  $outmol->fromSICHO($seq,$sicho);
  $outmol->writePDB(\*STDOUT,ssbond=>0);
} else {
  $sicho->writeChain(\*STDOUT);
}

