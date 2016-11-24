#!/usr/bin/env perl

# performs least squares fit between two protein structures
#
# http://mmtsb.scripps.edu/doc/lsqfit.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   lsqfit.pl [options] refPDB cmpPDB\n";
  printf STDERR "options: [-l min:max[=...]]\n";
  printf STDERR "         [-x min:max[=...]]\n";
  printf STDERR "         [-sel cab|ca|cb|p|heavy]\n";
  printf STDERR "         [-nowarn]\n";
  printf STDERR "         [-resnumonly] [-useseg]\n";
  printf STDERR "         [-s min:max[=...] min:max[=...]]\n";
  printf STDERR "         [-align fasta]\n";
  printf STDERR "         [-mem]\n";
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

my $fraglist;
my $exclmode;
my $selmode="cab";

my ($refpdb,$cmppdb);

my $warn=1;
my $resnumonly=undef;
my $useseg=0;

my $alignfile;

my $mem;

my $slist=();

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-sel") {
    shift @ARGV;
    $selmode=shift @ARGV;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
    $exclmode=0;
  } elsif ($ARGV[0] eq "-x") {
    shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
    $exclmode=1;
  } elsif ($ARGV[0] eq "-s") {
    shift @ARGV;
    my $trec={};
    $trec->{fraglist}=&GenUtil::fragListFromOption(shift @ARGV);
    $trec->{sellist}=&GenUtil::fragListFromOption(shift @ARGV);
    push(@{$slist},$trec);
    $exclmode=0;
  } elsif ($ARGV[0] eq "-nowarn") {
    shift @ARGV;
    $warn=0;
  } elsif ($ARGV[0] eq "-useseg") {
    shift @ARGV;
    $useseg=1;
  } elsif ($ARGV[0] eq "-resnumonly") {
    shift @ARGV;
    $resnumonly=1;
  } elsif ($ARGV[0] eq "-mem") {
    shift @ARGV;
    $mem=1;
  } elsif ($ARGV[0] eq "-align") {
    shift @ARGV;
    $resnumonly=1;
    $warn=0;
    $alignfile=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-.+/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $refpdb = shift @ARGV;
    $cmppdb = shift @ARGV;
    $done=1;
  }
}

my $refmol=Molecule::new($refpdb);
my $cmpmol=Molecule::new();
$cmpmol->readPDB($cmppdb);

my $analyze=Analyze::new($refmol);

if ($#{$slist}<0) {
  $cmpmol->setValidResidues($fraglist,$exclmode)
    if (defined $fraglist);
  
  $analyze->lsqfit($cmpmol,$selmode,$warn,$resnumonly,$alignfile,$useseg,0,$mem);
  $cmpmol->writePDB("-",translate=>"CHARMM22");
} else {
  my $newmol=Molecule::new();
  foreach my $s ( @{$slist} ) {
    $cmpmol->setValidResidues($s->{fraglist},$exclmode);
    $analyze->lsqfit($cmpmol,$selmode,$warn,$resnumonly,$alignfile,$useseg,0,$mem);
    $cmpmol->setValidResidues($s->{sellist},$exclmode);
    $newmol->merge($cmpmol);
  }
  $newmol->writePDB("-",translate=>"CHARMM22");
}

