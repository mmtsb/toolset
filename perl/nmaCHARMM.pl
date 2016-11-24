#!/usr/bin/env perl

# carry out normal mode analysis from PDB file 
#
# 2005, Michael Feig, MSU
#

sub usage {
  printf STDERR "usage:   nmaCHARMM.pl [options] [PDBfile]\n";
  printf STDERR "options: [-modes value]\n";
  printf STDERR "         [-xtract prefix]\n";
  printf STDERR "         [-minimize steps]\n";
  printf STDERR "         [-minout file]\n";
  printf STDERR "         [-block]\n";
  printf STDERR "         [-quasi dcdfile] [-temp value]\n";
  printf STDERR "         [-quasel selection]\n";
  printf STDERR "         [-ic] [-icscale value]\n";
  printf STDERR "         [-ics angle=value,bond=value,dihedral=value modefile]\n";
  printf STDERR "         [-dihedsample modefile prefix from to delta]\n";
  printf STDERR "         [-entropy]\n";
  printf STDERR "         [-psf PSFfile CRDfile]\n";
  printf STDERR "         [-par CHARMMparams]\n";
  printf STDERR "         [-log logFile] [-cmd logFile]\n";
  printf STDERR "         [-custom file]\n";
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

my %par;

my $logFile;
my $cmdlog;

my $inpfile="-";

my $psffile;
my $crdfile;

my $customfile;

my $nmodes;
my $xtract=undef;
my $minimize=0;
my $minout;
my $block=0;
my $negmodes=0;
my $ic=0;
my $icscale=1.0;
my %icspar={ "angle"    => 0.0,
             "bond"     => 0.0,
             "dihedral" => 0.0 };
my $icsfile;
my $icsprefix;
my $icsfrom;
my $icsto;
my $icsdelta;

my $quasi;
my $selq;
my $quasel;
my $qnatom;

my $entropy=0;
my $temperature=300.0;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=(shift @ARGV);
  } elsif ($ARGV[0] eq "-modes") {
    shift @ARGV;
    $nmodes=shift @ARGV;
  } elsif ($ARGV[0] eq "-negmodes") {
    shift @ARGV;
    $negmodes=1;
  } elsif ($ARGV[0] eq "-xtract") {
    shift @ARGV;
    $xtract=shift @ARGV;
  } elsif ($ARGV[0] eq "-minimize") {
    shift @ARGV;
    $minimize=shift @ARGV;
  } elsif ($ARGV[0] eq "-minout") {
    shift @ARGV;
    $minout=shift @ARGV;
  } elsif ($ARGV[0] eq "-block") {
    shift @ARGV;
    $block=1;
  } elsif ($ARGV[0] eq "-ic") {
    shift @ARGV;
    $ic=1;
  } elsif ($ARGV[0] eq "-icscale") {
    shift @ARGV;
    $icscale=shift @ARGV;
  } elsif ($ARGV[0] eq "-ics") {
    shift @ARGV;
    foreach my $s ( split(/,/,shift @ARGV) )  {
      my @sp=split(/=/,$s);
      $icspar{$sp[0]}=$sp[1];
    }
    $icsfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-dihedsample") {
    shift @ARGV;
    $icsfile=shift @ARGV;
    $icsprefix=shift @ARGV;
    $icsfrom=shift @ARGV;
    $icsto=shift @ARGV;
    $icsdelta=shift @ARGV;
  } elsif ($ARGV[0] eq "-entropy") {
    shift @ARGV;
    $entropy=1;
  } elsif ($ARGV[0] eq "-quasi") {
    shift @ARGV;
    $quasi=shift @ARGV;
  } elsif ($ARGV[0] eq "-quasel"){
    shift @ARGV;
    $quasel=shift @ARGV;
  } elsif ($ARGV[0] eq "-temp") {
    shift @ARGV;
    $temperature=shift @ARGV;
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-psf") {
    shift @ARGV;
    $psffile=shift @ARGV;
    $crdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    $inpfile=(shift @ARGV);
    $done=1;
  }
}

my $charmm=&CHARMM::new($logFile,$cmdlog);

$charmm->loadParameters(%par);

if (defined $psffile) {
  $charmm->setupFromPSF($psffile,$crdfile);
} else {
  $charmm->setupFromPDB($inpfile);
}

$charmm->setupEnergy();

if (defined $customfile && &GenUtil::checkFile($customfile)) {
  my $custom=&GenUtil::readData(&GenUtil::getInputFile($customfile));
  $charmm->stream($custom);
}

if (defined $icsfile && -r $icsfile) {
  my $icsamplitude=$icsfrom;
  my $icsinx=1;
  do {
    if (defined $icsprefix) {
      $icspar{dihedral}=$icsamplitude;
    }

    $charmm->_sendCommand("coor copy comp");
    $charmm->nmaICSample($icsfile,$icscale,%icspar);

    my $chmoutpdb=lc "pdb$$-out";
    $charmm->writePDB($chmoutpdb);
    
    my $outmol=Molecule::new();
    $outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),
		     chainfromseg=>(defined $psffile)?0:1);
  
    $outmol->setSSBonds($charmm->{molecule}->getSSBonds());

    if (defined $icsprefix) {
      my $outname=sprintf("%s.%d.pdb",$icsprefix,$icsinx++);
      $outmol->writePDB($outname,translate=>"CHARMM22");
      $icsamplitude+=$icsdelta;
      $charmm->_sendCommand("coor copy");
    } else {
      $outmol->writePDB("-",translate=>"CHARMM22");
    }
    &GenUtil::remove($chmoutpdb);
  } while (defined $icsprefix && $icsamplitude<=$icsto+0.00001);
    
  $charmm->finish();
  exit 0;
} 

if ($minimize>0 && !defined $quasi) {
  if ($charmm->{par}->{sdsteps}>0) {
    $charmm->minimizeSD();
  }
  
  $charmm->minimize(minsteps=>$minimize);

  if (defined $minout) {
    my $chmoutpdb=lc "pdb$$-out";
    $charmm->writePDB($chmoutpdb);

    my $outmol=Molecule::new();
    $outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),
                     chainfromseg=>(defined $psffile)?0:1);

    $outmol->setSSBonds($charmm->{molecule}->getSSBonds());
    $outmol->writePDB($minout,translate=>"CHARMM22");
    &GenUtil::remove($chmoutpdb);
  }
}

if ($entropy && !defined $quasi) {
  my $cmol=$charmm->{molecule};
  my $natoms=0;
  foreach my $c ( @{$cmol->{chain}})  {
    $natoms+=$#{$c->{atom}}+1;
  } 
  $nmodes=$natoms*6 if (!defined $nmodes);
  my ($vib,$rot,$trans)=$charmm->getEntropy($nmodes);
  printf STDOUT "%f %f %f\n",$vib,$rot,$trans;
  $charmm->finish();
} else {
  $nmodes=10 if (!defined $nmodes);

  my $nma;
  if ($quasi) { 
      if (defined $quasel){
	  $selq="SLQ";
	  $charmm->defineSelection($quasel,$selq);
    $charmm->_getCHARMMOutput("PERLDONE\n  \n")
    if ($charmm->{_lastOutput} eq "");
    foreach my $l (split(/\n/,$charmm->{_lastOutput}) ) {
      if ($l=~/(\d+) atoms have been selected out of/){
        $qnatom=$1;
      }
    }
      }
      else{
	  $selq="SLQ";
	  $charmm->defineSelection("all",$selq);
    $qnatom=undef;
      }
    $nma=$charmm->quasiHarmonicAnalysis($quasi,$nmodes,$temperature,selq=>$selq,qnatom=>$qnatom);
  } else {
    $nma=$charmm->normalModeAnalysis($nmodes,$block,$ic,$icscale,$xtract,$negmodes);
  }

  $charmm->finish();

  for (my $in=0; $in<=$#{$nma->{frequencies}}; $in++) {
    printf STDOUT "%d %f %f\n",$in+1,$nma->{frequencies}->[$in],$nma->{transrot}->[$in];
  }

  if (defined $xtract && (!defined $ic || !$ic)) {
   my $imode=1;
   for (my $in=0; $in<=$#{$nma->{modes}}; $in++) {
    if ($nma->{transrot}->[$in]<20 && ($nma->{frequencies}->[$in]>0.01 || $negmodes)) {
      my $mode=$nma->{modes}->[$in];
      open OUT,">$xtract.$imode";
      foreach my $m ( @{$mode->{vec}} ) {
        printf OUT "%f %f %f\n",$m->{x},$m->{y},$m->{z};
      }
      close OUT;
      $imode++;
    }
   }
  }
}

exit 0;

