#!/usr/bin/env perl
#
# get best clusters
#
# http://mmtsb.scripps.edu/doc/bestcluster.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   bestcluster.pl [options] tag\n";
  printf STDERR "options: [-dir datadir]\n";
  printf STDERR "         [-level num]\n";
  printf STDERR "         [-ctag tag]\n";
  printf STDERR "         [-prop tag[+tag...]]\n";
  printf STDERR "         [-size]\n";
  printf STDERR "         [-crit avg|avglow|avgcent|best|best#|median]\n";
  printf STDERR "         [-limstd multiple]\n";
  printf STDERR "         [-lowest]\n";
  printf STDERR "         [-xlowest tags]\n";
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
use Cluster;

my $tag;
my $dir=".";

my $ctag;

my $level=-1;

my $proptag="etot";

my $size=0;

my $crit="avg";
my $limstd=2.0;

my $lowest=0;
my $xtags;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-level") {
    shift @ARGV;
    $level=shift @ARGV;
  } elsif ($ARGV[0] eq "-prop") {
    shift @ARGV;
    $proptag=shift @ARGV;
  } elsif ($ARGV[0] eq "-size") {
    shift @ARGV;
    $size=1;
  } elsif ($ARGV[0] eq "-lowest") {
    shift @ARGV;
    $lowest=1;
  } elsif ($ARGV[0] eq "-xlowest") {
    shift @ARGV;
    $lowest=1;
    $xtags=shift @ARGV;
  } elsif ($ARGV[0] eq "-ctag") {
    shift @ARGV;
    $ctag=shift @ARGV;
  } elsif ($ARGV[0] eq "-crit") {
    shift @ARGV;
    $crit=shift @ARGV;
  } elsif ($ARGV[0] eq "-limstd") {
    shift @ARGV;
    $limstd=shift @ARGV;
  } elsif ($ARGV[0] =~/^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    $tag=shift @ARGV
      if (!defined $tag);
  }    
}

&usage() 
  if (!defined $tag);

my $ens=Ensemble->new($tag,$dir);

$ctag=$ens->{tag} if (!defined $ctag);

my $cluster=Cluster::new();

die "no clusters available"
  if (!-r "$ens->{dir}/$ctag.cluster");

$cluster->readFile("$ens->{dir}/$ctag.cluster");  

my $clist=$cluster->clusterList($level);

die "no clusters found"
  if (!defined $clist || $#{$clist}<0);

my $res=();
foreach my $c ( @{$clist} ) {
  my $plist=$ens->getPropList($proptag,$c->{element});
  my $rec={};
  ($rec->{score},$rec->{sdev},$rec->{nscore})=&_getScore($plist,$crit,$limstd);
  $rec->{cluster}=$c;
  $rec->{size}=$#{$c->{element}}+1;
  push (@{$res},$rec);
}

my @slist;

if ($size) {
  @slist=sort { $b->{size} <=> $a->{size} || $a->{score} <=> $b->{score} } @{$res};
} else {
  @slist=sort { $a->{score} <=> $b->{score} } @{$res};
}

foreach my $c ( @slist ) {
  if ($lowest) {
    my $flist=$cluster->fileList($c->{cluster}->{tag});
    my $plist= (defined $xtags) ? $ens->getSortedList($proptag,$flist,split(/,/,$xtags)) : $ens->getSortedList($proptag,$flist);
    printf "%-12s %5d %5d %12.4f %9.4f %9.4f %10.5f",
      $c->{cluster}->{tag},$c->{size},
	$c->{nscore},$c->{score},$c->{sdev},$c->{sdev}/sqrt($c->{nscore}),
	  $plist->[0]->{val};
    if (defined $xtags) {
      foreach my $xval ( @{$plist->[0]->{xtag}} ) {
	printf " %10.5f",$xval;
      }
    }
    printf " %s\n",$plist->[0]->{name};
  } else {
    printf "%-12s %5d %5d %12.4f %9.4f %9.4f\n",
      $c->{cluster}->{tag},$c->{size},
      $c->{nscore},$c->{score},$c->{sdev},($c->{nscore}!=0.0)?$c->{sdev}/sqrt($c->{nscore}):0.0;
  }
}

exit 0;

sub _getScore {
  my $list=shift;
  my $crit=shift;
  my $limstd=shift;

  my $score;
  my $sdev;
  my $nscore;

  my $olist=();
  my $lastlimit;

  if ($crit eq "avg") {
    ($score,$sdev,$nscore)=&GenUtil::average($list);
  } elsif ($crit eq "avglow") {
    ($olist,$lastlimit,$score,$sdev,$nscore)=&GenUtil::limCore($list,upper=>1,lower=>0,mult=>$limstd);
  } elsif ($crit eq "avgcent") {
    ($olist,$lastlimit,$score,$sdev,$nscore)=&GenUtil::limCore($list,upper=>1,lower=>1,mult=>$limstd);
  } elsif ($crit=~/^best/) {
    my $nbest=$crit;
    $nbest=~s/best//;
    $nbest=1 if ($nbest eq "");

    my @slist=sort { $a->{val} <=> $b->{val} } @{$list};

    if ($nbest==1) {
      $score=$slist[0]->{val};
      $sdev=0.0;
      $nscore=1;
    } else {
      $olist=();    
      for (my $i=0; $i<=$#slist && $i<$nbest; $i++) {
	push(@{$olist},$slist[$i]);
      }
      ($score,$sdev,$nscore)=&GenUtil::average($olist);
    }
  } elsif ($crit eq "median") {
    my @slist=sort { $a->{val} <=> $b->{val} } @{$list};
    $score=$slist[$#slist/2]->{val};
    $sdev=0.0;
    $nscore=1;
  }

  return ($score,$sdev,$nscore);
}

