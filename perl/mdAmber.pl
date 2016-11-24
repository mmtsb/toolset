#!/usr/bin/env perl

# minimize structure from PDB file using Amber
#
# http://mmtsb.scripps.edu/doc/mdAmber.pl.html
# 2000, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   mdAmber.pl [options] PDBfile\n";
  printf STDERR "options: [-par AmberParams]\n";
  printf STDERR "         [-restart file] [-restout file]\n";
  printf STDERR "         [-trajout file] [-final file]\n";
  printf STDERR "         [-partop file]\n";
  printf STDERR "         [-log file] [-elog file]\n";
  printf STDERR "         [-l [ca|cb|cab|heavy] force refpdb min:max[=...]]\n";
  printf STDERR "         [-cons [ca|cb|cab|heavy] refpdb min:max[_force][=...]]\n";
  printf STDERR "         [-[no]translate]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use IO::Handle;
use IO::File;

use GenUtil;
use Molecule;
use Amber;

my %par;

my $translate=1;

my $cons=undef;

my $logFile;
my $energyFile;

my $inpfile="-";
my $base="";

my $inpmode="start";
my $restfile;

my $trajout;
my $enerout;
my $restout;
my $finalpdb;
my $partop;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-restart") {
    shift @ARGV;
    $restfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-trajout") {
    shift @ARGV;
    $trajout=shift @ARGV;
  } elsif ($ARGV[0] eq "-restout") {
    shift @ARGV;
    $restout=shift @ARGV;
  } elsif ($ARGV[0] eq "-partop") {
    shift @ARGV;
    $partop=shift @ARGV;
  } elsif ($ARGV[0] eq "-final") {
    shift @ARGV;
    $finalpdb=shift @ARGV;
  } elsif ($ARGV[0] eq "-cons") {
    shift @ARGV;
    my $c={};
    if ($ARGV[0] =~ /^(ca|cb|cab|heavy)$/) {
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
    $c->{exclmode}=0;
    $c->{list}=&GenUtil::fragListFromOption(shift @ARGV);
    $cons=$c;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    my $c={};
    if ($ARGV[0] =~ /^(ca|cb|cab|heavy)$/) {
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
    $cons=$c;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-translate") {
    shift @ARGV;
    $translate=1;
  } elsif ($ARGV[0] eq "-notranslate") {
    shift @ARGV;
    $translate=0;
  } elsif ($ARGV[0] eq "-elog") {
    shift @ARGV;
    $energyFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    $inpfile=(shift @ARGV);
    $done=1;
  }
}

my $ambpartop=lc "$$.top";
my $ambinpcoor=lc "$$.inp.x";
my $ambinpcmd=lc "$$.inp.cmd";
my $amboutcoor=lc "$$.out.x";
my $amboutlog=lc "$$.out.log";
my $ambrefcoor=lc "$$.ref.x";
my $ambtrjout=lc "$$.trj.x";

my $mol=Molecule::new();
$mol->readPDB($inpfile);
$mol->fixHistidine($par{hsd},$par{hse});

my @map;
foreach my $c ( @{$mol->{chain}} ) {
  push(@map,$mol->renumber(1,$c->{id}));
}

if ($translate) {
  $mol->translate("AMBER");
}

my $amber=Amber::new();

$par{restrain}=(defined $cons)?1:0;

if (!defined $partop || !-r $partop || !defined $restfile || !-r $restfile) {
  $amber->generateTopCoor($mol,$ambpartop,$ambinpcoor,%par);
  system "cp $ambpartop $partop" if (defined $partop && !-r $partop);
}

$partop=$ambpartop if (!defined $partop);

my $amberinp=new IO::File;
$amberinp->open(">$ambinpcmd");

$amber->genInputDynamics($amberinp,(defined $restfile),(defined $trajout),%par);

if (defined $cons) {
  my $refmol=Molecule::new();
  $refmol->readAmber($partop,$ambinpcoor);
  for (my $ic=0; $ic<=$#{$mol->{chain}}; $ic++) {
    $refmol->numberReset($map[$ic],$mol->{chain}->[$ic]->{id});
  }

  if ($cons->{type} eq "ref" && defined $cons->{reffile}) {
    $refmol->zapCoordinates();
    $refmol->fillCoorFromPDB($cons->{reffile});
  }
  $refmol->writeAmber($ambrefcoor);

  if (defined $cons->{exclmode} && $cons->{exclmode}) {
    $refmol->setValidResidues($cons->{list},1);
    $cons->{list}=&GenUtil::gradForceList($refmol,$cons->{force});
  }	
  $amber->genInputRestraints($amberinp,$cons);
}

undef $amberinp;

$amber->runSander(input=>$ambinpcmd, partop=>$partop,
                  inpcoor=>(defined $restfile && -r $restfile)?$restfile:$ambinpcoor, 
		  outcoor=>$amboutcoor,log=>$amboutlog,restraint=>(defined $cons)?$ambrefcoor:undef,
		  trajout=>(defined $trajout)?$ambtrjout:undef);

system "cp $amboutlog $logFile" 
  if (defined $logFile);

system "cp $amboutcoor $restout" 
  if (defined $restout);

system "cp $ambtrjout $trajout" 
  if (defined $trajout);

$amber->logEnergy($amboutlog, $energyFile) 
  if (defined $energyFile);

if (defined $finalpdb) {
  my $outmol=Molecule::new();
  $outmol->readAmber($partop,$amboutcoor);

  for (my $ic=0; $ic<=$#{$mol->{chain}}; $ic++) {
    $outmol->numberReset($map[$ic],$mol->{chain}->[$ic]->{id});
  }

  $outmol->writePDB($finalpdb,translate=>"CHARMM22");
}

&GenUtil::remove($ambinpcmd);
&GenUtil::remove($ambpartop);
&GenUtil::remove($ambinpcoor);
&GenUtil::remove($amboutlog);
&GenUtil::remove($amboutcoor);
&GenUtil::remove("leap.out");
&GenUtil::remove($ambrefcoor);
&GenUtil::remove($ambtrjout);

exit 0;



