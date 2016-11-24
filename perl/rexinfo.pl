#!/usr/bin/env perl
#
# Extract replica exchange data
# (show clients, show conditions, report data by client, report data by condition)
#
# http://mmtsb.scripps.edu/doc/rexinfo.pl.html
# 2001, Michael Feig, Brooks group, TSRI
# 2001, John Karanicolas, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   rexinfo.pl [options]\n";
  printf STDERR "options: [-dir workdir]\n";
  printf STDERR "         [-inx from:to] [-step value]\n";
  printf STDERR "         [-condfile file]\n";
  printf STDERR "         [-clients]\n";
  printf STDERR "         [-conds]\n";
  printf STDERR "         [-byclient clientid]\n";
  printf STDERR "         [-bytemp temp]\n";
  printf STDERR "         [-bycond condindex]\n";
  printf STDERR "         [-clientsel clientid items]\n";
  printf STDERR "         [-condsel condindex items]\n";
  printf STDERR "         [-best]\n";
  printf STDERR "         [-files]\n";
  printf STDERR "         [-rankall] [-rank from:to] [-ranklast num]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use ReXServer;

my $showclients;
my $showconds;
my $clientid;
my $condinx;
my $dir=".";
my $condfile;
my $temp;
my $sellist;
my $files=0;
my $rank=0;
my $rankfrom;
my $rankto;
my $ranklast=100;

my $from=1;
my $to=9999999999;
my $step=1;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-condfile") {
    shift @ARGV;
    $condfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-clients") {
    shift @ARGV;
    $showclients=1;
  } elsif ($ARGV[0] eq "-conds") {
    shift @ARGV;
    $showconds=1;
  } elsif ($ARGV[0] eq "-byclient") {
    shift @ARGV;
    $clientid=shift @ARGV;
  } elsif ($ARGV[0] eq "-bytemp") {
    shift @ARGV;
    $temp=shift @ARGV;
  } elsif ($ARGV[0] eq "-bycond") {
    shift @ARGV;
    $condinx=shift @ARGV;
  } elsif ($ARGV[0] eq "-clientsel") {
    shift @ARGV;
    $clientid=shift @ARGV;
    $sellist=shift @ARGV;
  } elsif ($ARGV[0] eq "-condsel") {
    shift @ARGV;
    $condinx=shift @ARGV;
    $sellist=shift @ARGV;
  } elsif ($ARGV[0] eq "-best") {
    shift @ARGV;
    $condinx=0;
  } elsif ($ARGV[0] eq "-rankall") {
    shift @ARGV;
    $rank=1;
  } elsif ($ARGV[0] eq "-rank" ) {
    shift @ARGV;
    ($rankfrom,$rankto)=split(/:/,shift @ARGV);
    $rank=1;
  } elsif ($ARGV[0] eq "-ranklast") {
    shift @ARGV;
    $ranklast=shift @ARGV;
    $rankfrom=-1;
    $rankto=9999999;
    $rank=1;
  } elsif ($ARGV[0] eq "-files") {
    shift @ARGV;
    $files=1;
  } elsif ($ARGV[0] eq "-step") {
    shift @ARGV;
    $step=shift @ARGV;
  } elsif ($ARGV[0] eq "-inx") {
    shift @ARGV;
    my @f=split(/:/,shift @ARGV);
    $from=$f[0];
    $to=($#f>0)?$f[1]:$f[0];

  } elsif ($ARGV[0]=~/^-/) {
    printf STDERR "unknown option %s\n",shift @ARGV;
    &usage();
  }
}

my $rex=ReXServer->new(0,$dir);
$condfile=$dir."/rexserver.cond" if (!defined $condfile);
$rex->setup($condfile);
$rex->readData();


if (defined $temp) {
  my $diff=1E99;
  for (my $i=0; $i <= $#{$rex->{cond}}; $i++) {
    if ((! defined $diff) || (abs($rex->{cond}->[$i]->{temp} - $temp) < $diff)) {
      $condinx=$i;
      $diff=abs($rex->{cond}->[$i]->{temp} - $temp);
    }
  }
}

if (defined $showclients && $showclients) {
  die "client list not available"
    if ($#{$rex->{clientid}}<0);
  print join(" ",@{$rex->{clientid}})."\n";
} elsif (defined $showconds && $showconds) {
  my $cdat=$rex->getClientData();
  my $onecond=$rex->{cond}->[0];

  my $maxlen=0;
  for (my $i=0; $i <= $#{$rex->{cond}}; $i++) {
    my $thiscond=$rex->{cond}->[$i];
    for (my $j=0; $j < $rex->nUmbrellas(); $j++) {
      my $b=$thiscond->{bias}->[$j];
      my $t="";
      foreach my $kpar ( keys %{$b}) {
	$t.="$kpar: ".$b->{$kpar}."  ";
      }
      $maxlen=length($t) if (length($t)>$maxlen);
    }
  }

  printf "%-5s %-9s","INDEX","TEMP";
  for (my $j=0; $j < $rex->nUmbrellas(); $j++) {
    printf " %-*s",$maxlen,"BIAS ".($j+1);
  }
  printf "\n";

  for (my $i=0; $i <= $#{$rex->{cond}}; $i++) {
    my $thiscond=$rex->{cond}->[$i];
    printf "%-5d %-9.4f", $i, $thiscond->{temp};
    for (my $j=0; $j < $rex->nUmbrellas(); $j++) {
      my $b=$thiscond->{bias}->[$j];
      my $t="";
      foreach my $kpar ( keys %{$b}) {
	$t.="$kpar: ".$b->{$kpar}."  ";
      }
      printf " %-*s",$maxlen,$t;
    }
    printf "\n";
  }
} elsif (defined $clientid) {
  if ($clientid eq "all") {
    $clientid=join(":",sort @{$rex->{clientid}});
  }
  if (!defined $sellist) {
    if (!$files) {
      printf "Client %s\n",$clientid;
      printf "%-4s %-4s %-4s %-4s %-9s %-12s %-8s %-12s",
      "TRUN","MODE","RUN","COND","TEMP","ENERGY","RMSD","VOLUME";
      for (my $j=0; $j<$rex->nUmbrellas(); $j++) {
	printf " %-19s","BIAS ".($j+1);
      }
      printf "\n";
    }
    $sellist="trun,mode,run,cond,temp,ener,rmsd,volume,bias";
  }

  my @selarr=split(/,/,$sellist);

  for (my $i=$from; $i<=$rex->{trun} && $i<=$to; $i+=$step) {
    my @cdat=();
    foreach my $cidv ( split(/:/,$clientid) ) {
      my $tcdat=$rex->getClientData($cidv,$i);
      die "cannot get client data record for $cidv"
        if (!defined $tcdat);

      push(@cdat,$tcdat);
    }
    &showdata($rex,$i,$clientid,\@cdat,\@selarr,$files,$dir);
  }
} elsif (defined $condinx) {
  if (!defined $sellist) {
    if (!$files) {
      printf "Condition $condinx\n";
      printf "%-4s %-4s %-4s %-6s %-9s %-12s %-8s %-12s",
      "TRUN","MODE","RUN","CLIENT","TEMP","ENERGY","RMSD","VOLUME";
      for (my $j=0; $j<$rex->nUmbrellas(); $j++) {
	printf " %-19s","BIAS ".($j+1);
      }
      printf "\n";
    }
    $sellist="trun,mode,run,clientid,temp,ener,rmsd,volume,bias";
  }

  my @selarr=split(/,/,$sellist);

  for (my $i=$from; $i<=$rex->{trun} && $i<=$to; $i+=$step) {
    my $c=$rex->getClientData(undef,$i);
    my @cdat=();
    my $scid=undef;
    foreach my $cid (@{$rex->{clientid}}) {
      if ($c->{$cid}->{cond}->{inx} == $condinx) {
	$scid=$cid;
        push(@cdat,$c->{$cid});
      }
    }
    
    die "cannot find data record"
      if (!defined $cdat[0]);

    &showdata($rex,$i,$scid,\@cdat,\@selarr,$files,$dir);
  }
} elsif ($rank) {
  my %sumrank;
  my %bestrank;
  my $nrank=0;
  $rankfrom=1 if (!defined $rankfrom || $rankfrom==0);
  $rankfrom=$rex->{trun}-$ranklast if ($rankfrom<0);
  $rankfrom=1 if ($rankfrom<=0);
  $rankto=$rex->{trun} if (!defined $rankto || $rankto>$rex->{trun});
  for (my $i=$rankfrom; $i<=$rankto; $i++) {
    my $c=$rex->getClientData(undef,$i);
    foreach my $cid (@{$rex->{clientid}}) {
      $sumrank{$cid}+=$c->{$cid}->{cond}->{inx}+1;
      if ($c->{$cid}->{cond}->{inx} == 0) {
        $bestrank{$cid}++;
      } 
    }
    $nrank++;
  }
  for (my $ci=0; $ci<=$#{$rex->{clientid}}; $ci++) {
     my $cid=$rex->{clientid}->[$ci];
     printf "%-2d %s %1.2f %1.2f\n",
     $ci+1,$cid,$sumrank{$cid}/$nrank,$bestrank{$cid}/$nrank*100.0;
  }
}

sub showdata {
  my $rex=shift;
  my $run=shift;
  my $cid=shift;
  my $cdatarr=shift;
  my $sel=shift;
  my $files=shift;
  my $dir=shift;

  my @fcid=split(/:/,$cid);

  for (my $icid=0; $icid<=$#fcid; $icid++) {
   my $cidv=$fcid[$icid];
   my $cdat=$cdatarr->[$icid];
   if (defined $files && $files) {
    printf "%s/%s/%s/%s/final.pdb ",$dir,$cidv,lc $rex->getMode($run),&GenUtil::dataDir($rex->getRun($run));
   } else {
    for (my $i=0; $i<=$#{$sel}; $i++) {
      my $item=$sel->[$i];
      $item=lc $item;
      printf "%-4d ",$run                     if ($item eq "trun");
      printf "%-4s ",$rex->getMode($run)      if ($item eq "mode");
      printf "%-4d ",$rex->getRun($run)       if ($item eq "run");
      printf "%-6s ",$cidv                    if ($item eq "client" || $item eq "clientid");
      printf "%-4d ",$cdat->{cond}->{inx}     if ($item eq "cond");
      printf "%-4d ",$cdat->{cond}->{inx}+1   if ($item eq "condinx");
      printf "%-12.4f ",$cdat->{ener}        if ($item eq "ener" || $item eq "energy");
      printf "%-9.4f ",$cdat->{cond}->{temp} if ($item eq "temp" || $item eq "temperature");
      printf "%-9.4f ",$cdat->{rmsd}         if ($item eq "rmsd");
      printf "%-9.4f ",$cdat->{volume}       if ($item eq "volume");
      
      if ($item eq "bias") {
	for (my $j=0; $j<$rex->nUmbrellas(); $j++) {
	  printf "%-9.4f %-9.4f ",
	  $cdat->{cond}->{bias}->[$j]->{target},$cdat->{val}->[$j];
	}
      } elsif ($item eq "biasval") {
	for (my $j=0; $j<$rex->nUmbrellas(); $j++) {
	  printf "%-9.4f ",$cdat->{val}->[$j];
	}
      } elsif ($item =~ /target([0-9]+)/) {
	printf "%-9.4f ",$cdat->{cond}->{bias}->[$1-1]->{target};
      } elsif ($item =~ /val([0-9]+)/) {
	printf "%-9.4f ",$cdat->{val}->[$1-1];
      }
    }
   }
  }
  printf "\n";
}
