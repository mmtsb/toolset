#!/usr/bin/env perl
#
# run MD simulations with CHARMM
# http://mmtsb.scripps.edu/doc/mdCHARMM.pl.html
# 2001, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:    mdCHARMM.pl [options] PDBfile\n";
  printf STDERR "options:  [-restart filename]\n";
  printf STDERR "          [-enerout file] [-trajout file]\n";
  printf STDERR "          [-restout file] [-final pdb]\n";
  printf STDERR "          [-crdout file]\n";
  printf STDERR "          [-nochain]\n";
  printf STDERR "          [-elog file] [-log file]\n";
  printf STDERR "          [-cmd file]\n";
  printf STDERR "          [-stdout tag]\n";
  printf STDERR "          [-par CHARMMparams]\n";
  printf STDERR "          [-psf PSFfile CRDfile]\n";
  printf STDERR "          [-mol2 MOL2file]\n";
  printf STDERR "          [-boxsizefromrestart file]\n";
  printf STDERR "          [-l [ca|cb|cab|cabp|heavy] force refpdb|self min:max[=...]]\n";
  printf STDERR "          [-cons [ca|cb|cab|cabp|heavy] refpdb|self min:max[_force][=...]]\n";
  printf STDERR "          [-rmsd [ca|cb|cab|cabp|heavy] refpdb refval min:max[_force][=...]]\n";
  printf STDERR "          [-hmcm chainFile min:max[_force][=...]]\n";
  printf STDERR "          [-custom file]\n";
  printf STDERR "          [-comp PDBfile]\n";
  printf STDERR "          [-charmmexec charmmexec]\n";
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
# gb         =>  1,
 shake      =>  1,
 dynsteps   =>  1000,
 dyntemp    =>  298,
 dynens     =>  "NVT",
 dyneqfrq   =>  100,
 dynoutfrq  =>  50
); 

my $pdbfile;
my $inpmode="start";
my $restfile;

my $elog;
my $charmmlog;
my $cmdlog;

my $charmmexec;

my $trajout;
my $enerout;
my $restout;
my $finalpdb;
my $crdout;

my $cons=();

my $psffile;
my $crdfile;
my $mol2file;
my $pdbcomp;

my $customfile;

my $nochainoutput=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-restart") {
    shift @ARGV;
    $restfile=shift @ARGV;
    $inpmode="restart";
  } elsif ($ARGV[0] eq "-boxsizefromrestart") {
    shift @ARGV;
    $restfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $charmmlog=shift @ARGV;
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=shift @ARGV;
  } elsif ($ARGV[0] eq "-elog") {
    shift @ARGV;
    $elog=shift @ARGV;
  } elsif ($ARGV[0] eq "-enerout") {
    shift @ARGV;
    $enerout=shift @ARGV;
  } elsif ($ARGV[0] eq "-trajout") {
    shift @ARGV;
    $trajout=shift @ARGV;
  } elsif ($ARGV[0] eq "-restout") {
    shift @ARGV;
    $restout=shift @ARGV;
  } elsif ($ARGV[0] eq "-final") {
    shift @ARGV;
    $finalpdb=shift @ARGV;
  } elsif ($ARGV[0] eq "-crdout") {
    shift @ARGV;
    $crdout=shift @ARGV;
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-charmmexec") {
    shift @ARGV;
    $charmmexec=shift @ARGV;
  } elsif ($ARGV[0] eq "-nochain") {
    shift @ARGV;
    $nochainoutput=1;
  } elsif ($ARGV[0] eq "-stdout") {
    shift @ARGV;
    my $tag=shift @ARGV;
    $trajout=$tag.".traj.crd";
    $restout=$tag.".restart";
    $finalpdb=$tag.".final.pdb";
    $elog=$tag.".elog";
    $charmmlog=$tag.".log";
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-cons") {
    shift @ARGV;
    my $c={};
    if ($ARGV[0] =~ /^(ca|cb|cab|cabp|heavy)$/) {
      $c->{sel}=shift @ARGV;
    } else {
      $c->{sel}="heavy";
    }
    if ($ARGV[0] eq "self") {
      shift @ARGV;
      $c->{type}="self";
    } else {
      $c->{type}="ref";
      $c->{reffile}=shift @ARGV;
    }
    $c->{list}=&GenUtil::fragListFromOption(shift @ARGV);
    $c->{exclmode}=0;
    push(@{$cons},$c);
  } elsif ($ARGV[0] eq "-rmsd") {
    shift @ARGV;
    my $c={};
    if ($ARGV[0] =~ /^(ca|cb|cab|cabp|heavy)$/) {
      $c->{sel}=shift @ARGV;
    } else {
      $c->{sel}="heavy";
    }
    $c->{type}="rmsd";
    $c->{reffile}=shift @ARGV;
    $c->{refe}=shift @ARGV;
    $c->{list}=&GenUtil::fragListFromOption(shift @ARGV);
    $c->{orient}=1;
    $c->{exclmode}=0;
    push(@{$cons},$c);
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    my $c={};
    if ($ARGV[0] =~ /^(ca|cb|cab|cabp|heavy)$/) {
      $c->{sel}=shift @ARGV;
    } else {
      $c->{sel}="heavy";
    }

    $c->{force}=shift @ARGV;

    if ($ARGV[0] eq "self") {
      shift @ARGV;
      $c->{type}="self";
    } else {
      $c->{type}="ref";
      $c->{reffile}=shift @ARGV;
    }

    $c->{list}=&GenUtil::fragListFromOption(shift @ARGV);
    $c->{exclmode}=1;
    push(@{$cons},$c);
  } elsif ($ARGV[0] eq "-hmcm") {
    shift @ARGV;
    my $c={};
    $c->{type}="hmcm";
    $c->{reffile}=shift @ARGV;
    $c->{list}=&GenUtil::fragListFromOption(shift @ARGV);
    push(@{$cons},$c);
  } elsif ($ARGV[0] eq "-psf") {
    shift @ARGV;
    $psffile=shift @ARGV;
    $crdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-mol2") {
    shift @ARGV;
    $mol2file=shift @ARGV;
  } elsif ($ARGV[0] eq "-comp") {
    shift @ARGV;
    $pdbcomp=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "unknown option %s\n",shift @ARGV;
    &usage();
  } else {
    $pdbfile=shift @ARGV;
  }
}

$CHARMM::exec=$charmmexec if (defined $charmmexec);

my $charmm=&CHARMM::new($charmmlog,$cmdlog);
$charmm->setEnergyLogFile($elog) 
  if (defined $elog);

$charmm->loadParameters(%par);

if (defined $psffile) {
  $charmm->setupFromPSF($psffile,$crdfile);
} elsif (defined $mol2file) {
  $charmm->setupFromMol2($mol2file);
} else {
  $charmm->setupFromPDB($pdbfile);
}

if (!defined $charmm->{par}->{boxsize} &&
    !defined $charmm->{par}->{boxx} &&
    defined $restfile) {
  my @box=&CHARMM::boxsizeFromRestart($restfile);
  $charmm->setParameter(boxx=>$box[0], boxy=>$box[1], boxz=>$box[2], boxshape=>$box[3]);
}

if (defined $pdbcomp) {
  $charmm->loadReference($pdbcomp,9001.0);
}

$charmm->setupEnergy();

$charmm->setupRestraints(1.0,$cons) if ($#{$cons}>=0);
$charmm->noeRestraints();
$charmm->shake();

if (defined $customfile) {
  foreach my $c ( split(/:/,$customfile))  {
    if (&GenUtil::checkFile($c)) {
      my $custom=&GenUtil::readData(&GenUtil::getInputFile($c));
      $charmm->stream($custom);
    }
  }
}

$charmm->runDynamics(($inpmode eq "restart")?$restfile:undef,$restout,$trajout,$enerout,
		     dyntwin=>0.02*$charmm->{par}->{dyntemp});

$charmm->logMDEnergy("MD") if (defined $elog);

if (defined $crdout) {
  my $chmoutcrd=lc "crd$$-out";
  $charmm->writeCRD($chmoutcrd);
  my $ofile=&GenUtil::getOutputFile($crdout);
  my $ifile=&GenUtil::getInputFile($chmoutcrd);
  while (<$ifile>) {
    print $ofile $_;
  }
  close $ifile;
  close $ofile;
  &GenUtil::remove($chmoutcrd);
}

if (defined $finalpdb) {
  my $chmoutpdb="chmout.pdb";
  $charmm->writePDB($chmoutpdb);
  my $outmol=Molecule::new();
  $outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),chainfromseg=>(defined $psffile || $nochainoutput)?0:1);
  $outmol->setSSBonds($charmm->{molecule}->getSSBonds());
  $outmol->writePDB($finalpdb,translate=>"CHARMM22");
  &GenUtil::remove($chmoutpdb);
}

$charmm->finish();
