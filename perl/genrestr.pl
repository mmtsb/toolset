#!/usr/bin/env perl

# generates MONSSTER restraint file
#
# http://mmtsb.scripps.edu/doc/genrestr.pl.html
# 2000, Michael Feig, Brooks group, The Scripps Research Institute
#

sub usage {
  printf STDERR "usage:   genrestr.pl [option]\n";
  printf STDERR "options: [-r min:max_force[=...]]\n";
  printf STDERR "         [-l force min:max[=...]]\n";
  printf STDERR "         [-seq file]\n";   
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
use MONSSTER;
use Sequence;

my $restlist;
my $seqfile;
my $exclmode;
my $force;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-r") {
    shift @ARGV;
    $restlist=&GenUtil::fragListFromOption(shift @ARGV);
    $exclmode=0;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    my @f=split(/:/,shift @ARGV);
    $force=shift @ARGV;
    $restlist=&GenUtil::fragListFromOption(shift @ARGV);
    $exclmode=1;
  } elsif ($ARGV[0] eq "-seq") {
    shift @ARGV;
    $seqfile=shift @ARGV;
  } else {
    &usage();
  }
}

die "no list specified"
  unless (defined $restlist && $#{$restlist}>=0);

die "need sequence with -l option"
  if ($exclmode && !defined $seqfile);

my $monsster;

if (defined $seqfile) {
  my $seq=Sequence::new();
  $seq->readMONSSTER($seqfile);
  if ($exclmode) {
    $seq->setValidResidues($restlist,$exclmode);
    $restlist=&GenUtil::gradForceList($seq->listFromValid(),$force);
  }
  $monsster=MONSSTER::new($seq);
} else {
  $monsster=MONSSTER::new();
}

$monsster->setPositionalRestraints($restlist);
$monsster->writePosRestraints(\*STDOUT);
