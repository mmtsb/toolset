#!/usr/bin/env perl

# runs minimizations with CHARMM for ensembles
#
# http://mmtsb.scripps.edu/doc/ensmin.pl.html
# 2001, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   enseval.pl [options] tag\n";
  printf STDERR "options: [-par CHARMMparams]\n";
  printf STDERR "         [-set proptag=value,...]\n";
  printf STDERR "         [-overwrite]\n";
  printf STDERR "         [-run [from:]to]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-update frq]\n";
  printf STDERR "         [PARALLELoptions]\n";
  printf STDERR "         [-log file]\n";
  printf STDERR "         [-custom file]\n";
  printf STDERR "         [-charmmexec charmmexec]\n";
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

use Sys::Hostname;

use GenUtil;
use Molecule;
use CHARMM;
use Ensemble;
use JobServer;
use JobClient;
use Client;

my $intag;
my $dir=".";
my $srec={};
my $cpus=1;
my $hostfile;

my $overwrite=0;
my ($from,$to);

my $logfile;

my $customfile;

my $charmmexec=undef;

my %par=();

my $update=10;

my $mp=0;
my $keepmpdir=0;

my $jobrank;
my $saveid;

my @slist;

my $setstr="";
my $parstr="";

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-set") {
    shift @ARGV;
    $setstr=shift @ARGV;
    @slist=split(/,/,$setstr);
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-update") {
    shift @ARGV;
    $update=shift @ARGV;
  } elsif ($ARGV[0] eq "-run") {
    shift @ARGV;
    ($from,$to)=split(/:/,shift @ARGV);
    if (!defined $to) {
      $to=$from;
      $from=1;
    }
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
  } elsif ($ARGV[0] eq "-jobenv") {
    shift @ARGV;
    $jobrank=$ENV{shift @ARGV};
    $saveid="save.id" unless (defined $saveid);
  } elsif ($ARGV[0] eq "-hosts") {
    shift @ARGV;
    $hostfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-mp") {
    shift @ARGV;
    $mp=1;
  } elsif ($ARGV[0] eq "-keepmpdir") {
    shift @ARGV;
    $keepmpdir=1;
  } elsif ($ARGV[0] eq "-saveid") {
    shift @ARGV;
    $saveid=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-overwrite") {
    shift @ARGV;
    $overwrite=1;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    $parstr=shift @ARGV;
    &GenUtil::parsePar(\%par,$parstr);
  } elsif ($ARGV[0] eq "-charmmexec") {
    shift @ARGV;
    $charmmexec=shift @ARGV;
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    $intag=shift @ARGV if (!defined $intag);
  }
}

&usage() if (!defined $intag);

if (defined $charmmexec) {
  $ENV{CHARMMEXEC}=$charmmexec;
  $CHARMM::exec=$charmmexec;
}

push(@slist,"etotal=total") if ($#slist<0);

my $setlist=();
foreach my $s ( @slist ) {
  my @f=split(/=/,$s);
  my $rec={};
  $rec->{ptag}=$f[0];
  $rec->{prop}=$f[1];
  push (@{$setlist},$rec);
}

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
}

my $ens=Ensemble->new($intag,$dir);

my $jtag=($overwrite)?"none":$setlist->[0]->{ptag};

if (($cpus>1 || defined $jobrank) && !defined $srec->{name}) {
  my $jlist=$ens->jobList($from,$to,$jtag);

  die "nothing to do"
    if (!defined $jlist || $#{$jlist}<0);

  my $jobServer=JobServer->new($jlist,$ens,$update);
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
      $topt.=" -log $logfile-$i" if (defined $logfile);
      $topt.=" -par $parstr" if (defined $parstr);
      $topt.=" -set $setstr" if (defined $setstr);
      $topt.=" -overwrite" if ($overwrite);
      $topt.=" -custom $customfile" if (defined $customfile);
      $topt.=" -charmmexec $charmmexec" if (defined $charmmexec);

      if ($mp) {
	$topt.=" -mp";
      } else {
	$hostlist->[$i]->{localdir}=$workdir;
      }

      $topt.=" $intag";

      ($pid,$cpusleft)=&GenUtil::submitRemote($hostlist->[$i],$srec,$cpusleft,
					      (split(/\//,$0))[-1],$topt,$mp,!$keepmpdir);
      push(@pidlist,$pid);
    }
  } else {
    for (my $i=0; $i<$cpus; $i++) {
      my $pid=fork();
      
      if (!$pid) {
	my $jobClient=JobClient->new($srec);
	$jobClient->initialize();
	$jobClient->establishConnection();

	my $charmm=&CHARMM::new( (defined $logfile)?"$i-$logfile":undef );
	$charmm->loadParameters(%par);
  
#        if (defined $customfile && &GenUtil::checkFile($customfile)) {
#          my $custom=&GenUtil::readData(&GenUtil::getInputFile($customfile));
#          $charmm->stream($custom);
#        }
      
	my $job;
	my $lastprop;
	my $sendfiles;
	while (($job=$jobClient->nextJob($lastprop,$sendfiles))!~/^NO/) {
	  my $datadir;
	  if ($mp) {
	    $datadir=$ens->{dir}."/".&GenUtil::dataDir($job);
	    &GenUtil::makeDir($datadir);
	    $jobClient->getFile("$datadir/$intag.pdb");
	  }

	  &dojob($job,$ens,$charmm,$setlist,$customfile);

	  $lastprop=$ens->getPropString($job);
	}
	$jobClient->finish();
	
	$charmm->finish();
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
  my $jlist=$ens->jobList($from,$to,$jtag);

  die "nothing to do"
    if (!defined $jlist || $#{$jlist}<0);

  my $charmm=&CHARMM::new($logfile);
  $charmm->loadParameters(%par);
#  if (defined $customfile && &GenUtil::checkFile($customfile)) {
#     my $custom=&GenUtil::readData(&GenUtil::getInputFile($customfile));
#     $charmm->stream($custom);
#  }
  my $inx=0;
  foreach my $i ( @{$jlist} ) {
    &dojob($i,$ens,$charmm,$setlist,$customfile);
    $ens->save() if (++$inx%$update==0);
  }
  $charmm->finish();   
  $ens->save();
}

exit 0;

## dojob ######

sub dojob {
  my $job=shift;
  my $ens=shift;
  my $charmm=shift;
  my $slist=shift;
  my $customfile=shift;

  my $datadir=$ens->{dir}."/".&GenUtil::dataDir($job);

  my $inpfile=$datadir."/".$ens->{tag}.".pdb";

  if (!$charmm->{_setup}) {
    $charmm->setupFromPDB($inpfile);
    if (defined $customfile && &GenUtil::checkFile($customfile)) {
       my $custom=&GenUtil::readData(&GenUtil::getInputFile($customfile));
       $charmm->stream($custom);
    }
  } else {
    $charmm->initCoordinates();
    $charmm->readFromPDB($inpfile);
    $charmm->clearEnergy();
  }

  $charmm->setupEnergy();


  my $ener=$charmm->getEnergy();

  foreach my $s ( @{$slist} ) {
    my $o=$s->{prop};
    my $sum=0.0;
    my @addlist=($o=~/([+-]*)([^+-]+)/g);
    while (@addlist) {
      my $sgn=shift @addlist;
      my $a=shift @addlist;
      my $mult=1.0;
      if (defined $ener->{$a}) {
	if ($sgn eq "-") {
	  $mult=-1;
	}
	$sum+=$ener->{$a}*$mult;
      }
    }
    $ens->setProp($s->{ptag},$job,$sum);
  }
}
