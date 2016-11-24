#!/usr/bin/env perl
#
# calculate missing properties in ensemble 
#
# http://mmtsb.scripps.edu/doc/getprop.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   getprop.pl [options] tag\n";
  printf STDERR "options: [-dir datadir]\n";
  printf STDERR "         [-show]\n";
  printf STDERR "         [-prop prop[+prop][,prop]]\n";
  printf STDERR "         [-cluster tag] [-ctag tag]\n";
  printf STDERR "         [-score avg|avglow|avgcent|median]\n";
  printf STDERR "         [-limit low|high|cent]\n";
  printf STDERR "         [-liminx index]\n";
  printf STDERR "         [-limstd multiple]\n";  
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

my $show=0;
my @proplist;

my $score;
my $limstd=3.0;

my $cluster;
my $ctag;

my $limit;
my $liminx=1;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-show") {
    shift @ARGV;
    $show=1;
  } elsif ($ARGV[0] eq "-prop") {
    shift @ARGV;
    @proplist=split(/,/,shift @ARGV);
  } elsif ($ARGV[0] eq "-cluster") {
    shift @ARGV;
    $cluster=shift @ARGV;
  } elsif ($ARGV[0] eq "-ctag") {
    shift @ARGV;
    $ctag=shift @ARGV;
  } elsif ($ARGV[0] eq "-score") {
    shift @ARGV;
    $score=shift @ARGV;
  } elsif ($ARGV[0] eq "-limit") {
    shift @ARGV;
    $limit=shift @ARGV;
  } elsif ($ARGV[0] eq "-liminx") {
    shift @ARGV;
    $liminx=shift @ARGV;
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
  if (!defined $tag || (!$show && $#proplist<0));

my $ens=Ensemble->new($tag,$dir);

if ($show) {
  my $ptags=$ens->getPropTags();
  print "available properties are:\n";
  my $n=0;
  for my $p ( @{$ptags} ) {
    print "\n" if ((($n++)%4==0) && $n>1);
    print " $p";
  } 
  print "\n";
} else {
  my $cllist;
  if (defined $cluster) {
    $ctag=$ens->{tag} unless (defined $ctag);
    
    my $cl=Cluster::new();

    die "no clusters available"
      if (!-r "$ens->{dir}/$ctag.cluster");
    $cl->readFile("$ens->{dir}/$ctag.cluster");  

    $cllist=$cl->fileList($cluster);
    die "no cluster elements found"
      if (!defined $cllist);
  }

  if (!defined $score) {
    my $pdat=();
    for (my $i=1; $i<=$ens->{par}->{runs}; $i++) {
      $pdat->[$i]=();
    }

    my $p0list;
    for (my $i=0; $i<=$#proplist; $i++) {
      my $plist=$ens->getPropList($proplist[$i]);
      $p0list=$plist if ($i==($liminx-1));
      if (defined $plist) {
	foreach my $p (@{$plist}) {
	  $pdat->[$p->{inx}]->[$i]=$p->{val};
	}
      }
    }

    my %limitshow;
    if ($limit eq "low" && defined $p0list) {
      my ($olist,$lastlimit,$pavg,$psdev,$pn)=&GenUtil::limCore($p0list,upper=>1,lower=>0,mult=>$limstd);
      foreach my $o ( @{$olist} ) {
	$limitshow{$o->{inx}}=1;
      }
    } elsif ($limit eq "cent" && defined $p0list) {
      my ($olist,$lastlimit,$score,$sdev,$nscore)=&GenUtil::limCore($p0list,upper=>1,lower=>1,mult=>$limstd);
      foreach my $o ( @{$olist} ) {
	$limitshow{$o->{inx}}=1;
      }
    } elsif ($limit eq "high" && defined $p0list) {
      my ($olist,$lastlimit,$score,$sdev,$nscore)=&GenUtil::limCore($p0list,upper=>0,lower=>1,mult=>$limstd);
      foreach my $o ( @{$olist} ) {
	$limitshow{$o->{inx}}=1;
      }
    }

    if (defined $cllist) {
      foreach my $p ( sort { $a->{inx}<=>$b->{inx} } @{$cllist} ) {
	_showData($p->{inx},$pdat,$#proplist) if (!defined $limit || $limitshow{$p->{inx}});
      }
    } else {
      for (my $i=1; $i<=$ens->{par}->{runs}; $i++) {
	_showData($i,$pdat,$#proplist) if (!defined $limit || $limitshow{$i});
      }
    }
  } else {
    for (my $i=0; $i<=$#proplist; $i++) {
      my $plist=$ens->getPropList($proplist[$i],$cllist);
      die "no property values found for $proplist[$i]" if (!defined $plist);
      my ($score,$sdev,$nscore)=&_getScore($plist,$score,$limstd);
      printf " %f",$score;
    }
    printf "\n";
  }
}

exit 0;

sub _showData {
  my $inx=shift;
  my $pdat=shift;
  my $nplist=shift;
  
  printf "%d ",$inx;
  for (my $j=0; $j<=$nplist; $j++) {
    printf "%f ",$pdat->[$inx]->[$j];
  }
  printf "\n";
}

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

