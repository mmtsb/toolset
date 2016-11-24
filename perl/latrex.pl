#!/usr/bin/env perl
#
# Lattice replica exchange simulations
#
# http://mmtsb.scripps.edu/doc/latrex.pl.html
# 2001, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   latrex.pl [options] [files]\n";
  printf STDERR "options: [-n runs]\n";
  printf STDERR "         [-par initruns=value,equilruns=value,\n";
  printf STDERR "               [no]save,savebestfreq=value,archive,\n";
  printf STDERR "               [no]rebuild,ensmode=add|replace,\n";
  printf STDERR "               natpdb=file,seq=file]\n"; 
  printf STDERR "         [-temp nwin:min:max]\n";
  printf STDERR "         [-condfile file]\n";
  printf STDERR "         [-input pdb|chain]\n";
  printf STDERR "         [-f listfile]\n";
  printf STDERR "         [-simpar ncycle=val,icycle=val,\n";
  printf STDERR "                  stiff=val,short=val,central=val,kdcore=val]\n";
  printf STDERR "         [-simopt limforce=value,gridsize=value]\n";
  printf STDERR "         [-l refPDB min:max[=min:max ...]]\n";
  printf STDERR "         [-d force res1:res2[=res1:res2...]]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-ens tag] [-ensdir dir]\n";
  printf STDERR "         [PARALLELoptions]\n";
  printf STDERR "         [-log file]\n";
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

use ReXServer;
use LatSimReXServer;
use Client;
use ReXClient;
use GenUtil;
use Sequence;
use Ensemble;
use LatEnsemble;
use Molecule;
use SICHO;
use MONSSTER;

my %par=();
my %defpar = (
   initruns     => 0,
   equilruns    => 0,
   save         => 1,
   savebestfreq => -1,
   templist     => undef
);

my $nwindows;
my $mintemp;
my $maxtemp;

my $nruns=100;
my $condfile;

my $listfile;
my @initfile;

my $srec={};

my $inputmode="pdb";

my $sport;
my $sname;
my $sid;

my $from;
my $to;
my $cpus;

my $hostfile;
my $mp=0;
my $keepmpdir=0;

my %simopt=();
my %defsimopt = (
 limforce  => 1.0,
 gridsize  => 100
);

my %simpar=();
my %defsimpar = ( 
 ncycle  =>  5,
 icycle  =>  20
);

my ($drestforce, $drestlist);

my $serverlog;

my $saveid;

my $prun=0;

my $dir=".";

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-n") {
    shift @ARGV;
    $nruns=shift @ARGV;
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-ens") {
    shift @ARGV;
    $par{enstag}=shift @ARGV;
  } elsif ($ARGV[0] eq "-ensdir") {
    shift @ARGV;
    $par{ensdir}=shift @ARGV;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      if (defined $val) {
	$par{$key}=$val;
      } else {
	if ($key=~/^no(.+)$/) {
	  $par{$1}=0;
	} else {
	  $par{$key}=1;
	}
      }
    }
  } elsif ($ARGV[0] eq "-temp") {
    shift @ARGV;
    ($nwindows,$mintemp,$maxtemp)=split(/:/,shift @ARGV);
  } elsif ($ARGV[0] eq "-condfile") {
    shift @ARGV;
    $condfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $serverlog=shift @ARGV;
  } elsif ($ARGV[0] eq "-f") {
    shift @ARGV;
    $listfile=shift @ARGV;
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
  } elsif ($ARGV[0] eq "-jobs") {
    shift @ARGV;
    ($from,$to)=split(/:/,shift @ARGV);
    $to=$from if (!defined $to);
  } elsif ($ARGV[0] eq "-jobenv") {
    shift @ARGV;
    $from=$to=$ENV{shift @ARGV}+1;
    $cpus=1;
    $prun=1;
    $saveid="save.id" unless (defined $saveid);
  } elsif ($ARGV[0] eq "-cpus") {
    shift @ARGV;
    $cpus=shift @ARGV;
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
  } elsif ($ARGV[0] eq "-simpar") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      $simpar{$key}=(defined $val)?$val:1;
    }
  } elsif ($ARGV[0] eq "-simopt") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      if (defined $val) {
	$simopt{$key}=$val;
      } else {
	if ($key=~/^no(.+)$/) {
	  $simopt{$1}=0;
	} else {
	  $simopt{$key}=1;
	}
      }
    }
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $par{fragref}=shift @ARGV;
    $par{fraglist}=shift @ARGV;
  } elsif ($ARGV[0] eq "-input") {
    shift @ARGV;
    $inputmode=shift @ARGV;

  } elsif ($ARGV[0] eq "-d") {
    shift @ARGV;
    $drestforce=shift @ARGV;
    $drestlist=shift @ARGV;
  } elsif ($ARGV[0]=~/^-/) {
    printf STDERR "unknown option %s\n",shift @ARGV;
    &usage();
  } else {
    push(@initfile,shift @ARGV);
  }
}

if ($prun && $from>1) {
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
} else {
  if (defined $serverlog) {
    &GenUtil::makeDir("$dir");
    &GenUtil::setLogFile("$dir/$serverlog");
  }
}

my $mpClient;

if ($mp && defined $srec->{name}) {
  $mpClient=Client->new("mp",$srec);
  &GenUtil::makeDir($dir);
  $mpClient->getFile("$dir/rexclient.cfg");
  $mpClient->getFile("$dir/rexclient.options");
}

my $rex=Ensemble->new("rexclient",$dir,"rexclient.cfg");
$rex->set(fragref=>$par{fragref},fraglist=>$par{fraglist},runs=>$nruns,seq=>$par{seq});
$rex->setOption(%simopt);
$rex->setOption(drestforce=>$drestforce, drestlist=>$drestlist);
$rex->setPar(%simpar);

foreach my $p ( keys %defsimpar ) {
  if (!defined $rex->getPar() || !defined $rex->getPar()->{$p}) {
    $rex->setPar($p=>$defsimpar{$p});
  }
}
  
foreach my $p ( keys %defsimopt ) {
  if (!defined $rex->{opt} || !defined $rex->{opt}->{$p}) {
    $rex->setOption($p=>$defsimopt{$p});
  }
}

if (!defined $srec->{name}) {
  &GenUtil::log("rexserver","starting");
  
  if (defined $listfile && &GenUtil::checkFile($listfile)) {
    my $lfile=&GenUtil::getInputFile($listfile);
    while (<$lfile>) {
      chomp;
      push(@initfile,$_);
    }
    close $lfile;
  }
  
  my $server=LatSimReXServer->new($rex->{par}->{runs},$dir,\@initfile,%par);
  
  if (defined $server->{par}->{enstag}) {
    my $ens=LatEnsemble->new($server->{par}->{enstag},$server->{par}->{ensdir});
    $ens->set(fragref=>$rex->{par}->{fragref},
	      fraglist=>$rex->{par}->{fraglist},
	      seq=>$rex->{par}->{seq});
    $server->setEnsemble($ens);
    $server->{par}->{savebestfreq}=1 if ($server->{par}->{savebestfreq}<=0);
  }
  
  foreach my $p ( keys %defpar ) {
    if (!defined $server->{par} || !defined $server->{par}->{$p}) {
      $server->set($p=>$defpar{$p});
    }
  }
  
  if (defined $nwindows && !defined $condfile) {
    my $tlist=&ReXServer::expTempList($nwindows,$mintemp,$maxtemp);
    my $out=&GenUtil::getOutputFile("tmp-cond$$");
    foreach my $t ( @{$tlist} ) {
      printf $out "%f\n",$t;
    }
    close $out;
    $server->setup("tmp-cond$$");
    &GenUtil::remove("tmp-cond$$");
  } else {
    $server->setup($condfile);
  }

  $server->save();
  
  ($srec->{port},$srec->{id},$srec->{pid})=
    $server->run(4100);
  
  $srec->{name}=hostname;

  $nwindows=$server->nWindows();

  $rex->save();
  
  sleep 5;

  if (defined $saveid) {
    open OUT,">$saveid";
    printf OUT "%s:%s:%s\n",
    $srec->{name},$srec->{port},$srec->{id};
    close OUT;
    system "chmod 600 $saveid";
  }
  
}

my @joblist;
if (defined $from && defined $to) {
  for (my $i=$from; $i<=$to; $i++) {
    push(@joblist,sprintf("lat%d",$i));
  }
} elsif (defined $nwindows) {
  for (my $i=1; $i<=$nwindows; $i++) {
    push(@joblist,sprintf("lat%d",$i));
  }
}

die "no jobs to run" 
  if ($#joblist<0);


my @pidlist;

if (defined $hostfile) {
  my $hostlist=&GenUtil::readHostFile($hostfile);
  my $jobsdone=0;
  my $workdir=$ENV{PWD};
  chomp $workdir;

  my $xtrajobs;
  my $xtra=0;
  if (defined $cpus && $cpus<$#joblist+1) {
    $xtra=(($#joblist+1)-$cpus);
    $xtrajobs=int($xtra/$cpus);
  }

  for (my $i=0; $i<=$#{$hostlist} && $jobsdone<$#joblist+1; $i++) {
    my $maxcpus=&GenUtil::remoteCPUs($hostlist->[$i]);
    my $jfrom=$jobsdone+1;
    my $jto=$jobsdone+$maxcpus;

    $jto+=$xtrajobs*$maxcpus if ($xtra>0);
    $xtra-=$xtrajobs*$maxcpus;
    
    $jto=$#joblist+1 if ($jto>$#joblist+1);

    my $topt="";
    $topt.=" -jobs $jfrom:$jto";
    $topt.=" -log $serverlog-$i" if (defined $serverlog);
    $topt.=" -dir $dir";
    $topt.=" -input $inputmode";

    if ($mp) {
      $topt.=" -mp";
    } else {
      $hostlist->[$i]->{localdir}=$workdir;
    }

    my ($pid,$cpusleft)=&GenUtil::submitRemote($hostlist->[$i],$srec,$maxcpus,
				       (split(/\//,$0))[-1],"$topt",$mp,!$keepmpdir);
#($jto-$jfrom+1),
	
    $jobsdone=$jto;
    push(@pidlist,$pid);
  }
} else {
  my $hostid=sprintf("%s.%s",hostname,$$);

  for (my $i=0; $i<=$#joblist; $i++) {
    my $pid=fork();
    if (!$pid) {
      my $rexClient=ReXClient->new($joblist[$i],$srec);
      die "client cannot connect to server"
	if (!defined $rexClient);
      $rexClient->initialize($hostid,(defined $cpus)?$cpus:$#joblist+1);
      &doJob($mp,$dir,$joblist[$i],$rex,$rexClient,$inputmode);
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

exit 0;

### doJob ######

sub doJob {
  my $mp=shift;
  my $dir=shift;
  my $job=shift;
  my $rex=shift;
  my $rexClient=shift;
  my $inputmode=shift;

  &GenUtil::makeDir("$dir/$job");

  my $initfile=$rexClient->initFile();
  $initfile=undef if ($initfile eq "");

  if ($mp) {
    if (defined $initfile) {
      $rexClient->getFile($initfile,"$dir/$job/init");
      $initfile="$dir/$job/init";
    }
    if (!-r $rex->{par}->{seq}) {
      $rexClient->getFile($rex->{par}->{seq},"$dir/$job/seq");
      $rex->set(seq=>"$dir/$job/seq");
    }
  }

  my $inpchain=SICHO::new(gridsize=>$rex->{opt}->{gridsize});
  if (defined $initfile || defined $rex->getFragList()) {
    if ($inputmode eq "pdb") {
      if (defined $initfile) {
	my $mol=Molecule::new($initfile);
	$mol->selectChain("");
	$mol->center() if (!defined $rex->getFragList());
	$inpchain->genMONSSTERFromAllAtom($mol,fraglist=>$rex->getFragList());
      } else {
	if ($mp && defined $rex->{par}->{fragref} && !-r $rex->{par}->{fragref}) {
	  $rexClient->getFile($rex->{par}->{fragref},"$dir/$job/fragref.pdb");
	  $rex->set(fragref=>"$dir/$job/fragref.pdb");
	}
	$inpchain->genMONSSTERFromAllAtom($rex->getFragRef(),fraglist=>$rex->getFragList());
      }
    } else {
      $inpchain->readChain($initfile);
    }
  } else {
    $inpchain->genRandomMONSSTER($#{$rex->getSeq()->{sequence}}+1);
  }
  my $nchain=$#{$inpchain->{sidechain}}+1;
  my $nres=$#{$rex->getSeq()->{sequence}}+1;

  die "sequence and chain lengths do not match (seq: $nres, chain: $nchain)"
    if ($nres != $nchain-2);

  my $monsster=MONSSTER::new($rex->getSeq());

  $monsster->setDirectory("$dir/$job");


  $monsster->setParameters(%{$rex->getPar()});

  if (defined $rex->getFragList()) {
    $rex->getSeq()->setValidResidues($rex->getFragList(),1);
    $monsster->
      setPositionalRestraints(&GenUtil::gradForceList($rex->getSeq()->listFromValid(),
						      $rex->{opt}->{limforce}));
  }

  $monsster->setDistanceRestraints($rex->{opt}->{drestforce},
				   &GenUtil::fragListFromOption($rex->{opt}->{drestlist})) 
  if (defined $rex->{opt}->{drestforce} && defined $rex->{opt}->{drestlist});

  open OUT,">$dir/$job/monsster.job";
  printf OUT "%d\n%d\n%s\n%s\n%s\n%s\n%s\n",
  $nres+2,$mp,$job,$rexClient->{serverName},
  $rexClient->{serverPort},$rexClient->{serverID},
  "$dir/$job";
  close OUT;

  $monsster->run($inpchain);
}
