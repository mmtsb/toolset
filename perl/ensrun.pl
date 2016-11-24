#!/usr/bin/env perl

# runs user command for ensemble structures
#
# http://mmtsb.scripps.edu/doc/ensrun.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   ensrun.pl [options] tag [command]\n";
  printf STDERR "options: [-new tag]\n";
  printf STDERR "         [-set prop[:index][,prop[:index]]]\n";
  printf STDERR "         [-overwrite]\n";
  printf STDERR "         [-noinp]\n";
  printf STDERR "         [-natpdb pdbFile]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-update frq]\n";
  printf STDERR "         [-[no]compress]\n";
  printf STDERR "         [-run [from:]to]\n";
  printf STDERR "         [PARALLELoptions]\n";
  printf STDERR "         [-cmd file]\n";
  printf STDERR "         [-log file]\n";
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
use Ensemble;
use JobServer;
use JobClient;
use Client;

my $intag;
my $dir=".";
my $srec={};
my $cpus=1;
my $hostfile;
my $logfile;

my $update=10;

my $useinp=1;

my ($from,$to);

my $natpdb;
my $newtag;
my $setprop=();

my $setarg;

my $cmdfile;
my $cmd;

my $mp=0;
my $keepmpdir=0;

my $overwrite=0;

my $jobrank;
my $saveid;
my $compress;

while ($#ARGV>=0 && !defined $cmd) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
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
  } elsif ($ARGV[0] eq "-new") {
    shift @ARGV;
    $newtag=shift @ARGV;
  } elsif ($ARGV[0] eq "-nocompress") {
    shift @ARGV;
    $compress=0;
  } elsif ($ARGV[0] eq "-compress") {
    shift @ARGV;
    $compress=1;
  } elsif ($ARGV[0] eq "-set") {
    shift @ARGV;
    $setarg=shift @ARGV;
  } elsif ($ARGV[0] eq "-overwrite") {
    shift @ARGV;
    $overwrite=1;
  } elsif ($ARGV[0] eq "-natpdb") {
    shift @ARGV;
    $natpdb=shift @ARGV;
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logfile=shift @ARGV;
    &GenUtil::setLogFile($logfile);
  } elsif ($ARGV[0] eq "-noinp") {
    shift @ARGV;
    $useinp=0;
  } elsif ($ARGV[0] =~ /^-/) {
    die "unknown option $ARGV[0]";
  } else {
    if (!defined $intag) {
      $intag=shift @ARGV;
    } else {
      $cmd=join(" ",@ARGV);
    }
  }
}

&usage() if (!defined $intag);

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

my @cmdlist;

my $tag=(defined $newtag)?$newtag:$intag;

my $mpClient;
if ($mp && defined $srec->{name}) {
  $mpClient=Client->new("mp",$srec);
  &GenUtil::makeDir($dir);
  $mpClient->getFile("$dir/ens.cfg");
  $mpClient->getFile("$dir/$tag$Ensemble::OptionFileSuffix");
}

my $ens=Ensemble->new($tag,$dir);

$ens->set(natpdb=>$natpdb);
$ens->set(compress=>$compress) if (defined $compress);
$ens->save()
  if (!defined $srec->{name});

if (defined $cmd) {
  push (@cmdlist,$cmd);
  $cmdfile=undef;
} else {
  if ($mp && defined $srec->{name} && 
      defined $cmdfile && !&GenUtil::checkFile($cmdfile)) {
    $mpClient->getFile($cmdfile,"$dir/tmpcmd-$$");
    $cmdfile="$dir/tmpcmd-$$";
  }
  my $cmdhandle=&GenUtil::getInputFile($cmdfile);
  while(<$cmdhandle>) {
    chomp;
    push(@cmdlist,$_);
  }
  undef $cmdhandle;
}

my $lookuptag="none";

if (defined $setarg) {
  my $tinx=0;
  foreach my $f ( split(/,/,$setarg) ) {
    my @t=split(/:/,$f);
    my $rec={};
    $rec->{tag}=$t[0];
    $rec->{index}=(defined $t[1])?$t[1]:$tinx+1;
    $tinx=$rec->{index};
    push(@{$setprop},$rec);
  }
  $lookuptag=$setprop->[0]->{tag} if (!$overwrite);
}

if (($cpus>1 || defined $jobrank) && !defined $srec->{name}) {
  my $jlist=$ens->jobList($from,$to,$lookuptag);

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
    if (!defined $cmdfile) {
      $cmdfile="tmpcmd-$$";
      my $cfile=&GenUtil::getOutputFile($cmdfile);
      foreach my $c ( @cmdlist ) {
	print $cfile $c,"\n";
      }
      close $cfile;
    }

    my $hostlist=&GenUtil::readHostFile($hostfile);
    my $cpusleft=$cpus;
    my $workdir=$ENV{PWD};
    chomp $workdir;

    for (my $i=0; $i<=$#{$hostlist} && $cpusleft>0; $i++) {
      my $pid;
      my $topt=" -dir $dir";
      $topt.=" -new $newtag" if (defined $newtag);
      $topt.=" -set $setarg" if (defined $setarg);
      $topt.=" -cmd $cmdfile";
      $topt.=" -log $logfile" if (defined $logfile);
      $topt.=" -noinp" if (!$useinp);
      
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
    if ($mp) {
      if (defined $ens->{par}->{natpdb} && !-r $ens->{par}->{natpdb}) {
	$mpClient->getFile($ens->{par}->{natpdb},"$ens->{dir}/nat.pdb");
	$ens->set(natpdb=>"$ens->{dir}/nat.pdb");
      }
    }

    for (my $i=0; $i<$cpus; $i++) {
      my $pid=fork();
      
      if (!$pid) {
	my $jobClient=JobClient->new($srec);
	$jobClient->initialize();
	$jobClient->establishConnection();
	
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

	  &dojob($job,$ens,$intag,\@cmdlist,$setprop,$useinp);

	  if ($mp && defined $newtag) {
	    $sendfiles=();
	    push(@{$sendfiles},
		 { local  => "$datadir/$newtag.pdb",
		   remote => "$datadir/$newtag.pdb" });
	  }

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

  &GenUtil::remove($cmdfile)
    if (defined $cmdfile && $cmdfile=~/tmpcmd/);

  waitpid($srec->{pid},0) 
    if (defined $srec->{pid});
} else {
  my $jlist=$ens->jobList($from,$to,$lookuptag);

  die "nothing to do"
    if (!defined $jlist || $#{$jlist}<0);

  my $cnt=0;
  foreach my $i ( @{$jlist} ) {
    &dojob($i,$ens,$intag,\@cmdlist,$setprop,$useinp);
    $ens->saveProp()
      if (++$cnt%$update==0);
  }
  $ens->save();
}

exit 0;

## dojob ######

sub dojob {
  my $job=shift;
  my $ens=shift;
  my $intag=shift;
  my $cmdlist=shift;
  my $prop=shift;
  my $useinp=shift;

  my $datadir=$ens->{dir}."/".&GenUtil::dataDir($job);
  my $inpfile=$intag.".pdb";

  my $ret;
  for (my $i=0; $i<=$#{$cmdlist}; $i++) {
    my $cmd=$cmdlist->[$i];
    if ($i==0 && $useinp) {
      if (-r "$datadir/$inpfile.gz") {
	$cmd="cd $datadir; gunzip -c $inpfile.gz | $cmd";
      } else {
	$cmd="cd $datadir; $cmd < $inpfile";
      }
    } else {
      $cmd="cd $datadir; $cmd";
    }

    if ($i==$#{$cmdlist}) {
      if ($intag ne $ens->{tag}) {
	$cmd="$cmd > $ens->{tag}.pdb";
	system "$cmd";
	$ens->update($job);
	$ens->cleanUp($job);
      } elsif (defined $prop) {
	$ret=`$cmd`;      
	chomp $ret;
	$ret=~s/^ +//;
	my @rf=split(/[ \n\t]+/,$ret);
	foreach my $p ( @{$prop} ) {
	  $ens->setProp($p->{tag},$job,$rf[$p->{index}-1]);
	}
      } else {
	system $cmd;
      }
    } else {
      system $cmd;
    }
  }
}




