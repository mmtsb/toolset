#!/usr/bin/env perl
#
# extract a list of files from ensemble structures
#
# http://mmtsb.scripps.edu/doc/showcluster.pl.html
# 2001, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   showcluster.pl [options] tag\n";
  printf STDERR "options: [-dir workdir]\n";
  printf STDERR "         [-ctag tag]\n";
  exit 1;
}

require 5.004;

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use Cluster;
use Ensemble;
use GenUtil;

my $workdir=".";
my $tag;

my $cluster;

my $ctag;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $workdir=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-ctag") {
    shift @ARGV;
    $ctag=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    die "unknown option $ARGV[0]";
  } else {
    $tag=shift @ARGV
      if (!defined $tag);
  }
}

&usage() if (!defined $tag);

my $ens=Ensemble->new($tag,$workdir);

$ctag=$ens->{tag} if (!defined $ctag);

my $cluster=Cluster::new();

die "no clusters available"
  if (!-r "$ens->{dir}/$ctag.cluster");

$cluster->readFile("$ens->{dir}/$ctag.cluster");

my $list=$cluster->allClusters();

die "no clusters found"
  if (!defined $list || $#{$list}<0);

foreach my $l ( @{$list} ) {
  for (my $i=1; $i<=$l->{level}; $i++) {
    print '  ';
  }
  printf "%-12s %5d %5d\n",$l->{tag},$#{$l->{element}}+1,$#{$l->{subcl}}+1;
}

exit 0;

