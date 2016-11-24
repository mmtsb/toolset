#!/usr/bin/env perl
#
# All-atom replica exchange simulations with NAMD
#
# 2009, Michael Feig

sub usage {
  printf STDERR "usage:   aarexNAMD.pl [options] [coorfiles]\n";
  printf STDERR "options: [-n runs]\n";
  printf STDERR "         [-par initruns=value,equilruns=value,\n";
  printf STDERR "               [no]save,savebestfreq=value,archive\n";
  printf STDERR "               ensmode=add|replace,\n";
  printf STDERR "               psf=file,pdb=file]\n";
  printf STDERR "         [-temp nwin:min:max]\n";
  printf STDERR "         [-condfile file]\n";
  printf STDERR "         [-f listfile]\n";
  printf STDERR "         [-mdpar NAMDparams]\n";
  printf STDERR "         [-mdopt [no]trajout,[no]restout]\n";
  printf STDERR "         [-opt optionsfile]\n";
  printf STDERR "         [-custom file]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-ens tag] [-ensdir dir]\n";
  printf STDERR "         [PARALLELoptions]\n";
  printf STDERR "         [-log file] [-namdlog file]\n";
  printf STDERR "         [-namddir name] [-namdcpus num]\n";
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

use Server;
use ReXServer;
use NAMDReXServer;
use Client;
use ReXClient;
use GenUtil;
use Molecule;
use Ensemble;
use CHARMM;
use NAMD;

my %par=();
my %defpar = (
	      initruns     => 0,
	      equilruns    => 0,
	      save         => 1,
	      savebestfreq => -1,
);

my $condfile;
my $nwindows;
my $mintemp;
my $maxtemp;

my $nruns=100;

my $listfile;
my @initfile;

my $srec={};

my $from;
my $to;
my $cpus;

my $hostfile;

my $saveid;

my $clogfile;

my $cons;

my $namddir="/apps/namd/namd2.6x64";
my $namdcpus=2;

my %mdopt=();
my %defmdopt = (
		limforce  => 1.0,
		limsel    => "cab",
		trajout   => 0,
		restout   => 1
);

my %mdpar=();
my %defmdpar = ( 
		dynoutfrq =>  undef,
		dynens    =>  "NPT",
		shake     =>  1,
                dyntemp   =>  298,
		mmtsbcycle => 500
); 

my $serverlog;

my $dir=".";

my $optfile;

my $customfile;

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
  } elsif ($ARGV[0] eq "-opt") {
    shift @ARGV;
    $optfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
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
    my $rservfile=shift @ARGV;
    if (&GenUtil::checkFile($rservfile)) {
      open RINP,"$rservfile";
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
    $saveid="save.id" unless (defined $saveid);
  } elsif ($ARGV[0] eq "-cpus") {
    shift @ARGV;
    $cpus=shift @ARGV;
  } elsif ($ARGV[0] eq "-hosts") {
    shift @ARGV;
    $hostfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-saveid") {
    shift @ARGV;
    $saveid=shift @ARGV;
  } elsif ($ARGV[0] eq "-mdpar") {
    shift @ARGV;
    &GenUtil::parsePar(\%mdpar,shift @ARGV);
  } elsif ($ARGV[0] eq "-mdopt") {
    shift @ARGV;
    &GenUtil::parsePar(\%mdopt,shift @ARGV);
  } elsif ($ARGV[0] eq "-namdlog") {
    shift @ARGV;
    $clogfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-namddir") {
    shift @ARGV;
    $namddir=shift @ARGV;
  } elsif ($ARGV[0] eq "-namdcpus") {
    shift @ARGV;
    $namdcpus=shift @ARGV;
  } elsif ($ARGV[0]=~/^-/) {
    printf STDERR "unknown option %s\n",shift @ARGV;
    &usage();
  } else { 
    push(@initfile,shift @ARGV);
  }
}

if (defined $serverlog) {
  &GenUtil::makeDir("$dir");
  &GenUtil::setLogFile("$dir/$serverlog");
}

my $rex=Ensemble->new("rexclient",$dir,"rexclient.cfg");
if (defined $optfile) {
  foreach my $o ( split(/:/,$optfile) ) {
    $rex->readOptions($o);
  }
}

$rex->set(runs=>$nruns,pdb=>$par{pdb},psf=>$par{psf});
$rex->setOption(%mdopt);

if (defined $cons) {
  my $consstr="$cons->{sel},$cons->{type},$cons->{list},$cons->{reffile}";
  $rex->setOption(cons=>$consstr);
}

$rex->setPar(%mdpar);

foreach my $p ( keys %defmdpar ) {
  if (!defined $rex->getPar() || !defined $rex->getPar()->{$p}) {
    $rex->setPar($p=>$defmdpar{$p});
  }
}

$rex->setPar(dynrestfrq=>$rex->getPar()->{mmtsbcycle});
$rex->setPar(dynoutfrq=>$rex->getPar()->{mmtsbcycle}) unless (defined $rex->getPar()->{dynoutfrq});
$rex->setPar(dynsteps=>($rex->{par}->{runs}+1)*$rex->getPar()->{mmtsbcycle});

foreach my $p ( keys %defmdopt ) {
  if (!defined $rex->{opt} || !defined $rex->{opt}->{$p}) {
    $rex->setOption($p=>$defmdopt{$p});
  }
}

my $custom={};
&GenUtil::readCustomFile($custom,"$dir/rexclient.custom");

$custom->{prod}=&GenUtil::readData(&GenUtil::getInputFile($customfile))
  if (&GenUtil::checkFile($customfile));

&GenUtil::writeCustomFile($custom,"$dir/rexclient.custom")
  if (defined $custom);

if (!defined $srec->{name}) {
  &GenUtil::log("rexserver","starting");

  if (defined $listfile && &GenUtil::checkFile($listfile)) {
    my $lfile=&GenUtil::getInputFile($listfile);
    while (<$lfile>) {
      chomp;
      push(@initfile,$_);
      &GenUtil::log("rexserver","adding init file $_");
    }
    close $lfile;
  }

  if ($#initfile<0 && defined $rex->{par}->{pdb}) {
    push(@initfile,$rex->{par}->{pdb});
  }

  if (!defined $rex->{par}->{pdb} && $#initfile>=0) {
    $rex->{par}->{pdb}=$initfile[0];
  }

  my $server=NAMDReXServer->new($rex->{par}->{runs},$dir,\@initfile,%par);

  if (defined $server->{par}->{enstag}) {
    my $ens=Ensemble->new($server->{par}->{enstag},$server->{par}->{ensdir});
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
    push(@joblist,sprintf("aa%d",$i));
  }
} elsif (defined $nwindows) {
  for (my $i=1; $i<=$nwindows; $i++) {
    push(@joblist,sprintf("aa%d",$i));
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
    $topt.=" -namdlog $clogfile" if (defined $clogfile);
    $topt.=" -dir $dir";      
    $topt.=" -namddir $namddir";
    $topt.=" -namdcpus $namdcpus";
    $topt.=" -n $rex->{par}->{runs}";

    $hostlist->[$i]->{localdir}=$workdir;

    my ($pid,$cpusleft)=&GenUtil::submitRemote($hostlist->[$i],$srec,$maxcpus,
			                       (split(/\//,$0))[-1],"$topt",0,1);
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
      die "client cannot connect to server" if (!defined $rexClient);
      my $lasttemp=$rexClient->initialize($hostid,(defined $cpus)?$cpus:$#joblist+1);
#      $rexClient->establishConnection();
      my $thostfile="/tmp/$$-$joblist[$i]-hosts";

      open OUT,">$thostfile";
      printf OUT "group main\n";
      for (my $i=0; $i<$namdcpus; $i++) {
	printf OUT "host %s\n",hostname;
      }
      close OUT;

      $ENV{"CONV_RSH"}="ssh" if (!defined $ENV{"CONV_RSH"});
      $NAMD::exec=sprintf("%s/charmrun ++nodelist %s +p%d %s/namd2 +giga",$namddir,$thostfile,$namdcpus,$namddir);

      &doJob($dir,$joblist[$i],$rex,$rexClient,
             $lasttemp,$clogfile);
      if (-r $thostfile) {
        &GenUtil::remove($thostfile);
      }
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
  my $dir=shift;
  my $job=shift;
  my $rex=shift;
  my $rexClient=shift;
  my $lasttemp=shift;
  my $clogfile=shift;

  &GenUtil::makeDir("$dir/$job");

  my $initfile=$rexClient->initFile();

  if ($initfile eq "") {
    $rexClient->terminateServer();
    die "need initial structure for all-atom replica exchange";
  }

  my $namd=&NAMD::new();
  $namd->setParameter(%{$rex->getPar()});
  $namd->setParameter(mmtsbserver=>$rexClient->{serverName});
  $namd->setParameter(mmtsbport=>$rexClient->{serverPort});
  $namd->setParameter(mmtsbsid=>$rexClient->{serverID});
  $namd->setParameter(mmtsbjobid=>$job);

  $namd->fixParameters();
  $namd->runDynamics((defined $clogfile)?"$dir/$job/$clogfile":undef,undef,$rex->{par}->{psf},$initfile,undef,undef,undef,"$dir/$job/$NAMDReXServer::namdout","$dir/$job/$NAMDReXServer::namdrest",undef,0,$custom->{prod},undef);
}
