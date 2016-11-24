#!/usr/bin/env perl

# generates MONSSTER input file
#
# http://mmtsb.scripps.edu/doc/geninput.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   geninput.pl [options]\n";
  printf STDERR "options: [-d force res1:res2[=res1:res2 ...]]\n";
  printf STDERR "         [-par random=val,ncycle=val,icycle=val,\n";
  printf STDERR "               tsteps=val,resmin=val,resmax=val,\n";
  printf STDERR "               temp=first[:last],softcore=val,\n";
  printf STDERR "               central=val,stiff=val,pair=val,\n";
  printf STDERR "               kdcore=val,hbond=val,short=val,\n";
  printf STDERR "               burial=val,multibody=val,threebody=val]\n";
  printf STDERR "         [-seq file]\n";
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
use MONSSTER;
use Sequence;

my ($force, $fraglist);
my %par=();

my $seqfile;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-d") {
    shift @ARGV;
    $force=shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      $par{$key}=$val;
    }
  } elsif ($ARGV[0] eq "-seq") {
    shift @ARGV;
    $seqfile=shift @ARGV;
  } else {
    &usage();
  }
}

my $monsster;

if (defined $seqfile) {
  my $seq=Sequence::new();
  $seq->readMONSSTER($seqfile);
  $monsster=MONSSTER::new($seq);
} else {
  $monsster=MONSSTER::new();
}

$monsster->setDistanceRestraints($force,$fraglist) 
  if (defined $force && defined $fraglist);

foreach my $n (keys %par) {
  $monsster->setParameters($n => $par{$n});
}

$monsster->writeInput(\*STDOUT);
