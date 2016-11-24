#!/usr/bin/env perl

# generate PSF with CHARMM 
#
# 2006, Michael Feig, Michigan State University
#

sub usage {
  printf STDERR "usage:   genPSF.pl [PDBfile]\n";
  printf STDERR "         [-par CHARMMparams] [-crdout file]\n";
  printf STDERR "         [-xplor]\n";
  printf STDERR "         [-log file] [-cmd file]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use POSIX ":sys_wait_h";

use Fcntl;

use GenUtil;
use Molecule;
use CHARMM;

my %par;

my $logFile;
my $cmdlog;

my $pdbfile;
my $crdout;

my $xplor=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=(shift @ARGV);
  } elsif ($ARGV[0] eq "-crdout") {
    shift @ARGV;
    $crdout=shift @ARGV;
  } elsif ($ARGV[0] eq "-xplor") {
    shift @ARGV;
    $xplor=1;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    $pdbfile=shift @ARGV;
  }
}

my $charmm=CHARMM::new($logFile,$cmdlog);
$charmm->loadParameters(%par);
$charmm->setupFromPDB($pdbfile);
$charmm->writePSF("$$.psf",$xplor);
$charmm->writeCRD($crdout) if (defined $crdout);
$charmm->finish();

open INP,"$$.psf";
while(<INP>) {
  print;
}
close INP;

system "rm $$.psf";

exit 0;

  
