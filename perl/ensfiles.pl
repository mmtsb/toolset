#!/usr/bin/env perl
#
# extract a list of files from ensemble structures
#
# http://mmtsb.scripps.edu/doc/ensfiles.pl.html
# 2001, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   ensfiles.pl [options] tag\n";
  printf STDERR "options: [-cluster level]\n";
  printf STDERR "         [-ctag tag]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-list file]\n";
  printf STDERR "         [-sort proptag]\n";
  printf STDERR "         [-prop proptag]\n";
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

my $listfile;

my $proptag;
my @xtags;
my $sort=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-cluster") {
    shift @ARGV;
    $cluster=shift @ARGV;
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $workdir=shift @ARGV;
  } elsif ($ARGV[0] eq "-list") {
    shift @ARGV;
    $listfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-sort") {
    shift @ARGV;
    @xtags=split(/,/,shift @ARGV);
    $proptag=shift @xtags;
    $sort=1;
  } elsif ($ARGV[0] eq "-prop") {
    shift @ARGV;
    @xtags=split(/,/,shift @ARGV);
    $proptag=shift @xtags;
    $sort=0;
  } elsif ($ARGV[0] eq "-ctag") {
    shift @ARGV;
    $ctag=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] =~ /^-/) {
    die "unknown option $ARGV[0]";
  } else {
    $tag=shift @ARGV
      if (!defined $tag);
  }
}

&usage() if (!defined $tag);

my $ens=Ensemble->new($tag,$workdir);

my $filelist=();
if (!defined $cluster) {
  my $list;
  if (defined $listfile) {
    open INP,$listfile;
    while (<INP>) {
      chomp;
      s/^[ \t]+//;
      my @f=split(/[ \t\n]+/);
      push(@{$list},$f[0]);
    }
    close INP;
  }

  $filelist=$ens->fileList($list);
} else {
  my $cl=Cluster::new();

  $ctag=$ens->{tag} if (!defined $ctag);

  die "no clusters available"
    if (!-r "$ens->{dir}/$ctag.cluster");

  $cl->readFile("$ens->{dir}/$ctag.cluster");

  $filelist=$cl->fileList($cluster);
}

if (defined $filelist) {
  if (defined $proptag) {
    if ($sort) {
     foreach my $f ( @{$ens->getSortedList($proptag,$filelist,@xtags)} ) {
       printf "%s %f ",$f->{name},$f->{val};      
       foreach my $x ( @{$f->{xtag}} ) {
         printf "%f ",$x;
       }
       printf "\n";
     } 
    } else {
     foreach my $f ( @{$ens->getList($proptag,$filelist,@xtags)} ) {
       printf "%s %f ",$f->{name},$f->{val};      
       foreach my $x ( @{$f->{xtag}} ) {
         printf "%f ",$x;
       }
       printf "\n";
     } 
    }
  } else {
    foreach my $f ( @{$filelist} ) { #sort { $a->{inx}<=>$b->{inx} } @{$filelist} ) {
      printf "%s %f\n",$f->{name},$f->{dist};
    }
  }
} else {
  die "no files found";
}

exit 0;

