#!/usr/bin/env perl

# check NOE restraints for a given structure with CHARMM
#
# http://mmtsb.scripps.edu/doc/noecheck.pl.html
# 2002, Michael Feig, Brooks group, The Scripps Research Institute

sub usage {
  printf STDERR "usage:   noecheck.pl [options] PDBfile\n";
  printf STDERR "options: [-par noerest=file,xnoerest=file]\n";
  printf STDERR "         [-log logFile] [-cmd logFile]\n";
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

my $logFile;
my $cmdlog;

my $inpfile="-";
my $base="";

my %par = ( 
 param   =>  22
); 

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      $par{$key}=(defined $val)?$val:1;
    }
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=(shift @ARGV);
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
$charmm->setupFromPDB($inpfile);
$charmm->noeRestraints();
my $noev=$charmm->getNOEAnalysis();

foreach my $n ( @{$noev} ) {
  my $a1=sprintf("%s:%d:%s",$n->{seg1},$n->{res1},$n->{at1});
  my $a2=sprintf("%s:%d:%s",$n->{seg2},$n->{res2},$n->{at2});
  printf "%-5d %-20s %-20s %5.2f %5.2f %5.2f %s\n",
  $n->{inx},$a1,$a2,$n->{rmin},$n->{rmax},$n->{ractual},
  ($n->{violation}>0.001 || $n->{violation}<-0.001)?sprintf("%5.2f",$n->{violation}):"ok";
}

$charmm->finish();

exit 0;

