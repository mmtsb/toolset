#!/usr/bin/env perl
#
# calculate missing properties in ensemble 
#
# http://mmtsb.scripps.edu/doc/calcprop.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   calcprop.pl [options] tag\n";
  printf STDERR "options: [-dir datadir]\n";
  printf STDERR "         [-natpdb file] [-chain id]\n";
  printf STDERR "         [-l min:max[=...]]\n";
  printf STDERR "         [PARALLELoptions]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use Sys::Hostname;

use GenUtil;
use Ensemble;
use JobServer;
use Client;
use JobClient;

my $tag;
my $dir=".";
my $mp=0;
my $keepmpdir=0;

my $srec={};
my $cpus=1;
my $natpdb;

my $hostfile;
my $fraglist;
my $chain;

my $jobrank;
my $saveid;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-rserv") {
    shift @ARGV;
    my $rfile=shift @ARGV;
    if (&GenUtil::checkFile($rfile)) {
      open RINP,"$rfile";
      my $line=<RINP>;
      chomp $line;
      ($srec->{name},$srec->{port},$srec->{id})=
      split(/:/,$line);
      close RINP;
    }
  } elsif ($ARGV[0] eq "-cpus") {
    shift @ARGV;
    $cpus=shift @ARGV;
  } elsif ($ARGV[0] eq "-mp") {
    shift @ARGV;
    $mp=1;
  } elsif ($ARGV[0] eq "-keepmpdir") {
    shift @ARGV;
    $keepmpdir=1;
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-natpdb") {
    shift @ARGV;
    $natpdb=shift @ARGV;
  } elsif ($ARGV[0] eq "-chain") {
    shift @ARGV;
    $chain=shift @ARGV;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fraglist=shift @ARGV;
  } elsif ($ARGV[0] eq "-hosts") {
    shift @ARGV;
    $hostfile=shift @ARGV;
  } elsif ($ARGV[0] =~/^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } elsif ($ARGV[0] eq "-saveid") {
    shift @ARGV;
    $saveid=shift @ARGV;
  } elsif ($ARGV[0] eq "-jobenv") {
    shift @ARGV;
    $jobrank=$ENV{shift @ARGV};
    $saveid="save.id" unless (defined $saveid);
  } else {
    $tag=shift @ARGV
      if (!defined $tag);
  }    
}

&usage() if (!defined $tag);

if (defined $jobrank && $jobrank>0) {
  do {
    sleep 10;
    if (&GenUtil::checkFile($saveid)) {
      open RINP,"$saveid";
      my $line=<RINP>;
      chomp $line;
      ($srec->{name},$srec->{port},$srec->{id})=
	split(/:/,$line);
      close RINP;
    }
  } while(!defined $srec->{name});
}

my $mpClient;
if ($mp && defined $srec->{name}) {
  $mpClient=Client->new("mp",$srec);
  &GenUtil::makeDir($dir);
  $mpClient->getFile("$dir/ens.cfg");
  $mpClient->getFile("$dir/$tag$Ensemble::OptionFileSuffix");
}

my $ens=Ensemble->new($tag,$dir);

$ens->set(natpdb=>$natpdb, fraglist=>$fraglist, chain=>$chain);

$ens->save()
  if (!defined $srec->{name});

if (($cpus>1 || defined $jobrank) && !defined $srec->{name}) {
  my $jlist=$ens->jobList(undef,undef,"rmsdall");

  die "nothing to do"
    if (!defined $jlist || $#{$jlist}<0);
  
  my $jobServer=JobServer->new($jlist,$ens,100);
  ($srec->{port},$srec->{id},$srec->{pid})=$jobServer->run(4000);
  $srec->{name}=hostname;

  sleep 5;
  
  if (defined $saveid) {
    open OUT,">$saveid";
    printf OUT "%s:%s:%s\n",
    $srec->{name},$srec->{port},$srec->{id};
    close OUT;
    system "chmod 600 $saveid";
  }
}

if (defined $srec->{name}) {
  my @pidlist;

  if (defined $hostfile) {
    my $hostlist=&GenUtil::readHostFile($hostfile);
    my $cpusleft=$cpus;
    my $workdir=$ENV{PWD};
    chomp $workdir;

     for (my $i=0; $i<=$#{$hostlist} && $cpusleft>0; $i++) {
       my $pid;
       my $topt=" -dir $dir";

       if ($mp) {
	 $topt.=" -mp";
       } else {
	 $hostlist->[$i]->{localdir}=$workdir;
       }

       $topt.=" $tag";

       ($pid,$cpusleft)=&GenUtil::submitRemote($hostlist->[$i],$srec,$cpusleft,
					       (split(/\//,$0))[-1],$topt,$mp,!$keepmpdir);
      push(@pidlist,$pid);
    }
  } else {
    if ($mp && defined $ens->{par}->{natpdb} && !-r $ens->{par}->{natpdb}) {
      $mpClient->getFile($ens->{par}->{natpdb},"$ens->{dir}/nat.pdb");
      $ens->set(natpdb=>"$ens->{dir}/nat.pdb");
    }

    for (my $i=0; $i<$cpus; $i++) {
      my $pid=fork();
      if (!$pid) {
	my $jobClient=JobClient->new($srec);
	$jobClient->initialize();
	$jobClient->establishConnection();

	my $job;
	my $lastprop;
	while (($job=$jobClient->nextJob($lastprop))!~/^NO/) {
	  if ($mp) {
	    my $d=$ens->{dir}."/".&GenUtil::dataDir($job);
	    &GenUtil::makeDir($d);
#	    $jobClient->getFile("$d/$ens->{tag}.elog");
	    $jobClient->getFile("$d/$ens->{tag}.pdb");
	  }
	  $ens->getStructInfo($job);
	  $lastprop=$ens->getPropString($job);
	}
	$jobClient->finish();
	exit 0;
      } else {
	push(@pidlist,$pid);
      }
    }
  }

  foreach my $p ( @pidlist ) {
    waitpid($p,0);
  }

  waitpid($srec->{pid},0) 
    if (defined $srec->{pid});
} else {
  $ens->update();
  $ens->save();
}

