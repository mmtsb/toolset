#!/usr/bin/env perl
#
# run MD simulations with NAMD
# 2007, Michael Feig, Michigan State University
#

sub usage {
  printf STDERR "usage:    mdNAMD.pl [options] PSFfile PDBfile\n";
  printf STDERR "options:  [-restout tag] [-restart tag] [-final pdb] [-trajout DCDname]\n";
  printf STDERR "          [-coor file] [-ext file]\n";
  printf STDERR "          [-elog file] [-log file] [-cmd file]\n";
  printf STDERR "          [-par NAMDparams] [-first step]\n";
  printf STDERR "          [-custom file]\n";
  printf STDERR "          [-consref file]\n";
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
use NAMD;

my %par = ( 
 shake      =>  1,
 dynsteps   =>  1000,
 dyntemp    =>  298,
 dynens     =>  "NPT",
 dynoutfrq  =>  500
); 

my $pdbfile;

my $elog;
my $namdlog;
my $cmdlog;

my $restname;
my $outname;
my $trajout;
my $finalpdb;
my $coorfile;
my $extfile;
my $tmdfile;
my $consfile;

my $psffile;

my $customfile;

my $firststep=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $namdlog=shift @ARGV;
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=shift @ARGV;
  } elsif ($ARGV[0] eq "-elog") {
    shift @ARGV;
    $elog=shift @ARGV;
  } elsif ($ARGV[0] eq "-trajout") {
    shift @ARGV;
    $trajout=shift @ARGV;
  } elsif ($ARGV[0] eq "-restout") {
    shift @ARGV;
    $outname=shift @ARGV;
  } elsif ($ARGV[0] eq "-restart") {
    shift @ARGV;
    $restname=shift @ARGV;
  } elsif ($ARGV[0] eq "-final") {
    shift @ARGV;
    $finalpdb=shift @ARGV;
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-coor") {
    shift @ARGV;
    $coorfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-ext") {
    shift @ARGV;
    $extfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-first") {
    shift @ARGV;
    $firststep=shift @ARGV;
  } elsif ($ARGV[0] eq "-consref") {
    shift @ARGV;
    $consfile=shift @ARGV;
    $par{cons}=1;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "unknown option %s\n",shift @ARGV;
    &usage();
  } else {
    $psffile=shift @ARGV;
    $pdbfile=shift @ARGV;
  }
}

$outname="run" if (!defined $outname);

my $namd=&NAMD::new();

$namd->setParameter(%par);
$namd->fixParameters();
$namd->runDynamics($namdlog,$cmdlog,$psffile,$pdbfile,$coorfile,$extfile,$trajout,$outname,$restname,$elog,$firststep,$customfile,$consfile);

if (defined $finalpdb && -r "$outname.coor") {
  my $mol=Molecule::new();
  $mol->readPDB($pdbfile);
  
  my $coorfile=&GenUtil::getInputFile("$outname.coor");
  binmode $coorfile;
  
  my $buffer="";
  read($coorfile,$buffer,4);
  my $natoms=unpack("L",$buffer);

  my $nn=0;
  foreach my $c ( @{$mol->{chain}}) {
    foreach my $a ( @{$c->{atom}}) {
      read($coorfile,$buffer,24);
      my ($xval,$yval,$zval)=unpack("d*",$buffer);
      $a->{xcoor}=$xval;
      $a->{ycoor}=$yval;
      $a->{zcoor}=$zval;
      $nn++;
    }
  }
  if ($nn != $natoms) {
    printf STDERR "%d atoms read from coordinate file, %d atoms expected from header\n",$nn,$natoms;
  }
  close $coorfile;
  undef $coorfile;
  $mol->writePDB($finalpdb,translate=>"CHARMM22");
}

