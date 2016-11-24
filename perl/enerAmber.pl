#!/usr/bin/env perl

# minimize structure from PDB file using Amber
#
# http://mmtsb.scripps.edu/doc/enerAmber.pl.html
# 2000, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   enerAmber.pl [options] PDBfile\n";
  printf STDERR "options: [-out total,bonds,angles,dihedral,impropers,vdwaals,elec,gb,sasa]\n";
  printf STDERR "         [-oneline]\n";
  printf STDERR "         [-[no]translate]\n";
  printf STDERR "         [-par AmberParams]\n";
  printf STDERR "         [-log logFile]\n";
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

my $oneline=0;
my @olist;
push(@olist,"total");

my $translate=1;

my $logFile;

my $inpfile="-";
my $base="";

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-oneline") {
    shift @ARGV;
    $oneline=1;
  } elsif ($ARGV[0] eq "-out") {
    shift @ARGV;
    my $line=shift @ARGV;
    @olist=split(/,/,$line);
  } elsif ($ARGV[0] eq "-translate") {
    shift @ARGV;
    $translate=1;
  } elsif ($ARGV[0] eq "-notranslate") {
    shift @ARGV;
    $translate=0;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
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

my $mol=Molecule::new();
$mol->readPDB($inpfile);
$mol->fixHistidine($par{hsd},$par{hse});
#$mol->selectChain("");

if ($translate) {
  $mol->translate("AMBER");
}

my $amber=Amber::new();

$par{minsteps}=0;
$par{sdsteps}=0;
$amber->generateTopCoor($mol,$ambpartop,$ambinpcoor,%par);

my $amberinp=new IO::File;
$amberinp->open(">$ambinpcmd");

$amber->genInputMinimize($amberinp);

undef $amberinp;

$amber->runSander(input=>$ambinpcmd,partop=>$ambpartop,
                  inpcoor=>$ambinpcoor,outcoor=>$amboutcoor,
                  log=>$amboutlog);

system "cp $amboutlog $logFile" 
  if (defined $logFile);

my $ener=$amber->getEnergy($amboutlog);

foreach my $o ( @olist ) {
  my $sum=0.0;
  my @addlist=($o=~/([+-]*)([^+-]+)/g);
  while (@addlist) {
    my $sgn=shift @addlist;
    my $a=shift @addlist;
    my $mult=1.0;
    if (defined $ener->{$a}) {
      if ($sgn eq "-") {
	$mult=-1;
      }
      $sum+=$ener->{$a}*$mult;
    }
  }
  if ($#olist>0) {
    if ($oneline) {
      printf "%15.4f ",$sum;
    } else {
      printf "%-20s %15.4f\n",$o,$sum;
    }
  } else {
    printf "%15.4f\n",$sum;
  }
}

printf "\n" if ($oneline && $#olist>0);

&GenUtil::remove($ambinpcmd);
&GenUtil::remove($ambpartop);
&GenUtil::remove($ambinpcoor);
&GenUtil::remove($amboutlog);
&GenUtil::remove($amboutcoor);
&GenUtil::remove("leap.out");

exit 0;



