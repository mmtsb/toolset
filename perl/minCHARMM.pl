#!/usr/bin/env perl

# minimize structure from PDB file using CHARMM
#
# http://mmtsb.scripps.edu/doc/minCHARMM.pl.html
# 2000, Michael Feig, Brooks group, The Scripps Research Institute

sub usage {
  printf STDERR "usage:   minCHARMM.pl [options] [PDBfile]\n";
  printf STDERR "options: [-par CHARMMparams]\n";
  printf STDERR "         [-psf PSFfile CRDfile]\n";
  printf STDERR "         [-mol2 MOL2file]\n";
  printf STDERR "         [-crdout]\n";
  printf STDERR "         [-nochain] [-splitseg]\n";
  printf STDERR "         [-l [ca|cb|cab|heavy] force self|refpdb min:max[=...]]\n";
  printf STDERR "         [-cons [ca|cb|cab|heavy] self|refpdb min:max[_force][=...]]\n";
  printf STDERR "         [-hmcm chainFile min:max[_force][=...]]\n";
  printf STDERR "         [-rmsd [ca|cb|cab|heavy|cap|cabp] refpdb refval min:max[_force][=...]]\n";
  printf STDERR "         [-custom file[:file]]\n";
  printf STDERR "         [-comp PDBfile]\n";
  printf STDERR "         [-log logFile] [-elog energyLogFile]\n";
  printf STDERR "         [-cmd logFile]\n";
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

my $cons=();

my $logFile;
my $cmdlog;
my $energyFile;

my $inpfile="-";
my $base="";

my %par;

my $psffile;
my $crdfile;
my $mol2file;
my $pdbcomp;

my $crdout=0;
my $nochainoutput=0;
my $splitseg=0;

my $customfile;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-cons") {
    shift @ARGV;
    my $c={};
    if ($ARGV[0] =~ /^(ca|cb|cab|cabp|cap|heavy)$/) {
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
    if ($ARGV[0] =~ /^(ca|cb|cab|cabp|cap|heavy)$/) {
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
    if ($ARGV[0] =~ /^(ca|cb|cab|cabp|cap|heavy)$/) {
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
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=(shift @ARGV);
  } elsif ($ARGV[0] eq "-elog") {
    shift @ARGV;
    $energyFile=(shift @ARGV);
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
  } elsif ($ARGV[0] eq "-crdout") {
    shift @ARGV;
    $crdout=1;
  } elsif ($ARGV[0] eq "-nochain") {
    shift @ARGV;
    $nochainoutput=1;
  } elsif ($ARGV[0] eq "-splitseg") {
    shift @ARGV;
    $splitseg=1;
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    $inpfile=(shift @ARGV);
    $done=1;
  }
}

my $charmm=&CHARMM::new($logFile,$cmdlog);
$charmm->setEnergyLogFile($energyFile) 
  if (defined $energyFile);

$charmm->loadParameters(%par);

if (defined $psffile) {
  $charmm->setupFromPSF($psffile,$crdfile);
} elsif (defined $mol2file) {
  $charmm->setupFromMol2($mol2file);
} else {
  $charmm->setupFromPDB($inpfile);
}

if (defined $pdbcomp) {
  $charmm->loadReference($pdbcomp,9001.0);
}

$charmm->setupEnergy();
$charmm->shake();

$charmm->noeRestraints();

$charmm->setupRestraints(1.0,$cons)
  if ($#{$cons}>=0);

if (defined $customfile) {
  foreach my $c ( split(/:/,$customfile))  {
    if (&GenUtil::checkFile($c)) {
      my $custom=&GenUtil::readData(&GenUtil::getInputFile($c));
      $charmm->stream($custom);
    }
  }
}

if ($charmm->{par}->{sdsteps}>0) {
  $charmm->minimizeSD();
  $charmm->logEnergy("SD");
}

if ($charmm->{par}->{minsteps}>0) {
  $charmm->minimize();
  $charmm->logEnergy("Min final");
}

if ($crdout) {
  my $chmoutcrd=lc "crd$$-out";
  $charmm->writeCRD($chmoutcrd);
  my $ofile=&GenUtil::getOutputFile("-");
  my $ifile=&GenUtil::getInputFile($chmoutcrd);
  while (<$ifile>) {
    print $ofile $_;
  }
  close $ifile;
  close $ofile;
  &GenUtil::remove($chmoutcrd);
} else {
  my $chmoutpdb=lc "pdb$$-out";
  $charmm->writePDB($chmoutpdb);
  if (!-r $chmoutpdb) {
     sleep(20);
  } else {
     sleep(5);
  }
  my $outmol=Molecule::new();
  $outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),chainfromseg=>(defined $psffile || $nochainoutput)?0:1,splitseg=>($splitseg)?1:0);

  $outmol->setSSBonds($charmm->{molecule}->getSSBonds());
  $outmol->writePDB("-",translate=>"CHARMM22");
  &GenUtil::remove($chmoutpdb);
}

$charmm->finish();

exit 0;

