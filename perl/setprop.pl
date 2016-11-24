#!/usr/bin/env perl
#
# calculate missing properties in ensemble 
#
# http://mmtsb.scripps.edu/doc/setprop.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   setprop.pl [options] tag proptag [value]\n";
  printf STDERR "options: [-dir datadir]\n";
  printf STDERR "         [-remove]\n";
  printf STDERR "         [-at index]\n";
  printf STDERR "         [-f file]\n";
  printf STDERR "         [-inx index]\n";
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
use Ensemble;

my $tag;
my $proptag;
my $dir=".";
my $inpfile;
my $start=1;
my $value;
my $inx=0;
my $remove=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-at") {
    shift @ARGV;
    $start=shift @ARGV;
  } elsif ($ARGV[0] eq "-f") {
    shift @ARGV;
    $inpfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-remove") {
    shift @ARGV;
    $remove=1;
  } elsif ($ARGV[0] eq "-inx") {
    shift @ARGV;
    $inx=(shift @ARGV)-1;
  } elsif ($ARGV[0] =~/^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    $tag=shift @ARGV
      if (!defined $tag);
    $proptag=shift @ARGV
      if (!defined $proptag);
    $value=shift @ARGV
      if (!defined $value);
  }    
}

&usage() 
  if (!defined $tag || !defined $proptag);

my $ens=Ensemble->new($tag,$dir);

my $n=$start;

if ($remove) {
  $ens->removeProp($proptag);
} elsif (defined $value) {
  $ens->setProp($proptag,$n,$value);
} else {
  $inpfile="-" if (!defined $inpfile);
  my $inp=&GenUtil::getInputFile($inpfile);
  while (<$inp>) {
    chomp;
    s/^[ \t]+//;
    my @f=split(/[ \t\n]+/);
    $ens->setProp($proptag,$n++,$f[$inx]);
  }
  close $inp;
}

$ens->save();
