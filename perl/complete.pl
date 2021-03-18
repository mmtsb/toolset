#!/usr/bin/env perl
#
# adds missing atoms to a PDB structure 
#
# http://mmtsb.scripps.edu/doc/complete.pl.html
# 2000, Michael Feig, Brooks group, TSRI
# 2006, Michael Feig, MSU

sub usage {
  printf STDERR "usage:   complete.pl [options] [PDBfile]\n";
  printf STDERR "options: [-hsd list] [-hse list] [-param 19|22|27] [-blocked]\n";
  printf STDERR "         [-log file] [-cmd file]\n";
  printf STDERR "         [-nosicho] [-nofixca]\n";
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
use CHARMM;

my %par = (
  hsd => "",
  hse => "",
  param => 36,
  blocked => 0
);

my $fname;

my $logfile;
my $cmdfile;

my $fraglist;
my $sicho=1;
my $fixca=1;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-hsd") {
    shift @ARGV;
    $par{hsd}=shift @ARGV;
  } elsif ($ARGV[0] eq "-hse") {
    shift @ARGV;
    $par{hse}=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-param") {
    shift @ARGV;
    $par{param}=shift @ARGV;
  } elsif ($ARGV[0] eq "-nosicho") {
    shift @ARGV;
    $sicho=0;
  } elsif ($ARGV[0] eq "-nofixca") {
    shift @ARGV;
    $fixca=0;
  } elsif ($ARGV[0] eq "-blocked") {
    shift @ARGV;
    $par{blocked}=1;
  } else {
    $fname=shift @ARGV;
  }
}

my $mol=Molecule::new();
$mol->readPDB($fname);
$mol->translate("CHARMM22");

$mol->markClashes();
foreach my $c ( @{$mol->{chain}}) {
  foreach my $a ( @{$c->{atom}}) {
    if ($a->{atomname} eq "CA" && !$a->{valid}) {
      $a->{valid}=1;
    }
  }
}

my $tmol=$mol->clone(1);

my $ss=$tmol->getSSBonds();

my $nmol=$tmol->completeWater();
$nmol->completeResidue($sicho,$fixca);
$nmol->fixCOO() if (!$par{blocked});
$nmol->translate("CHARMM22");
$nmol->fixHistidine($par{hsd},$par{hse});
$nmol->generateSegNames();

my $charmm=&CHARMM::new($logfile,$cmdfile);
$charmm->loadParameters(%par);
$charmm->setupFromMolecule($nmol);

my $chmoutpdb="t$$.out.pdb";
$charmm->writePDB($chmoutpdb);
my $outmol=Molecule::new();
$outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($par{param}),chainfromseg=>1);

#$outmol->markClashes();
#my $omol=$outmol->clone(1);

$outmol->setSSBonds($ss);
$outmol->writePDB("-",translate=>"CHARMM22");
&GenUtil::remove($chmoutpdb);
  
$charmm->finish();
