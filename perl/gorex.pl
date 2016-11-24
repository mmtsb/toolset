#!/usr/bin/env perl
#
# Note: mp mode is currently not usable, because Server.pm won't accept .crd
#    files for file transfers (line 285 and subsequent of Server.pm)
#
#

# Go model replica exchange simulations
#
# http://mmtsb.scripps.edu/doc/gorex.pl.html
# 2001, Michael Feig, Brooks group, TSRI
# 2001, John Karanicolas, Brooks group, TSRI

sub usage {
    printf STDERR "usage:   gorex.pl [options] PREFIX\n";
    printf STDERR "options: [-n runs]\n";
    printf STDERR "         [-f listfile]\n";
    printf STDERR "         [-dir workdir]\n";
    printf STDERR "         [-par initruns=value,equilruns=value,\n";
    printf STDERR "               [no]save,savebestfreq=value]\n";
    printf STDERR "         [-temp nwin:min:max]\n";
    printf STDERR "         [-condfile file]\n";
    printf STDERR "         [-ens tag] [-ensdir dir]\n";
    printf STDERR "         [PARALLELoptions]\n";
    printf STDERR "         [-log file] [-elog file] [-charmmlog file]\n";
    printf STDERR "         [-opt optionsfile]\n";
    printf STDERR "         [-custom setup|pre|post[:init|equi|prod] file]\n";
    printf STDERR "         [-mdpar steps=val,fbeta=val,timestep=val,outfrq=val,dynupdnb=val,dynupdimg=val]\n";
    printf STDERR "         [-mdopt [no]trajout,[no]enerout]\n";
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
use Client;
use ReXClient;
use GenUtil;
use Ensemble;
use CHARMM;

my %par=();
my %defpar = (
   initruns     => 10,
   equilruns    => 250,
   save         => 1,
   savebestfreq => -1,
   archive      => 1
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
my $mp=0;
my $keepmpdir=0;

my $saveid;

my $elogfile;
my $clogfile;

my $prefix;

my %mdopt=();
my %defmdopt = (
 trajout   => 1,
);

my %mdpar=();
my %defmdpar = ( 
 timestep =>  undef,
 fbeta    =>  undef,
 steps    =>  20000,
 outfrq   =>  500,
 dynupdnb =>  0,
 dynupdimg =>  0
); 

my $serverlog;

my $dir=".";

my $optfile;

my %customfile;

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
	&GenUtil::makeDir("$dir");
	&GenUtil::setLogFile("$dir/$serverlog");
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
    } elsif ($ARGV[0] eq "-mdpar") {
	shift @ARGV;
	foreach my $p ( split(/,/,shift @ARGV) ) {
	    my ($key,$val)=split(/=/,$p);
	    $mdpar{$key}=(defined $val)?$val:1;
	}
    } elsif ($ARGV[0] eq "-mdopt") {
	shift @ARGV;
	foreach my $p ( split(/,/,shift @ARGV) ) {
	    my ($key,$val)=split(/=/,$p);
	    if (defined $val) {
		$mdopt{$key}=$val;
	    } else {
		if ($key=~/^no(.+)$/) {
		    $mdopt{$1}=0;
		} else {
		    $mdopt{$key}=1;
		}
	    }
	}
    } elsif ($ARGV[0] eq "-elog") {
	shift @ARGV;
	$elogfile=shift @ARGV;
    } elsif ($ARGV[0] eq "-charmmlog") {
	shift @ARGV;
	$clogfile=shift @ARGV;
    } elsif ($ARGV[0] eq "-custom") {
	shift @ARGV;
	my $ctag=shift @ARGV;
	$customfile{$ctag}=shift @ARGV;
    } elsif ($ARGV[0]=~/^-/) {
	printf STDERR "unknown option %s\n",shift @ARGV;
	&usage();
    } else {
	$prefix=shift @ARGV;
	my $pdbfile="GO_" . $prefix . ".pdb";
	push(@initfile,$pdbfile);
    }
}

die "Prefix must be specified for Go models"
    if (! defined $prefix);
die "Timestep must be specified for Go models"
    if (! defined $mdpar{timestep});
die "Fbeta must be specified for Go models"
    if (! defined $mdpar{fbeta});


my $mpClient;

if ($mp && defined $srec->{name}) {
    $mpClient=Client->new("mp",$srec);
    &GenUtil::makeDir($dir);
    $mpClient->getFile("$dir/rexclient.cfg");
    $mpClient->getFile("$dir/rexclient.options");
    $mpClient->getFile("$dir/rexclient.custom");
}

my $rex=Ensemble->new("rexclient",$dir,"rexclient.cfg");
if (defined $optfile) {
    foreach my $o ( split(/:/,$optfile) ) {
	$rex->readOptions($o);
    }
}
$rex->set(runs=>$nruns,natpdb=>$par{natpdb});
$rex->setOption(%mdopt);
$rex->setPar(%mdpar);

foreach my $p ( keys %defmdpar ) {
    if (!defined $rex->getPar() || !defined $rex->getPar()->{$p}) {
	$rex->setPar($p=>$defmdpar{$p});
    }
}

foreach my $p ( keys %defmdopt ) {
    if (!defined $rex->{opt} || !defined $rex->{opt}->{$p}) {
	$rex->setOption($p=>$defmdopt{$p});
    }
}

my $custom={};
&GenUtil::readCustomFile($custom,"$dir/rexclient.custom");

foreach my $c ( keys %customfile ) {
    $custom->{$c}=&GenUtil::readData(&GenUtil::getInputFile($customfile{$c}))
	if (&GenUtil::checkFile($customfile{$c}));
}    

&GenUtil::writeCustomFile($custom,"$dir/rexclient.custom")
    if (defined $custom);


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
    
    my $server=ReXServer->new($rex->{par}->{runs},$dir,\@initfile,%par);
    
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
	    printf $out "%1.4f\n",$t;
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
  
    if (defined $saveid) {
	open OUT,">$saveid";
	printf OUT "host: %s, port: %s, id: %s\n",
	$srec->{name},$srec->{port},$srec->{id};
	close OUT;
	chmod 0600, $saveid;
    }
    
    $nwindows=$server->nWindows();
    
    $rex->save();
    
    sleep 10;
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
  my $workdir=`pwd`;
  chomp $workdir;

  for (my $i=0; $i<=$#{$hostlist} && $jobsdone<$#joblist+1; $i++) {
    my $maxcpus=&GenUtil::remoteCPUs($hostlist->[$i]);
    my $jfrom=$jobsdone+1;
    my $jto=$jobsdone+$maxcpus;
    $jto=$#joblist+1 if ($jto>$#joblist+1);

    my $topt="";
    $topt.=" -jobs $jfrom:$jto";
    $topt.=" -log $serverlog-$i" if (defined $serverlog);
    $topt.=" -charmmlog $clogfile" if (defined $clogfile);
    $topt.=" -elog $elogfile" if (defined $elogfile);
    $topt.=" -dir $dir";

    if ($mp) {
      $topt.=" -mp";
    } else {
      $hostlist->[$i]->{localdir}=$workdir;
    }

    $topt .= " -mdpar fbeta=" . $mdpar{fbeta} . ",timestep=" . $mdpar{timestep};
    $topt .= " " . $prefix;
    my $cmd=$0;
    if ($ENV{'CHARMMEXEC'} eq "") {
	$cmd="setenv CHARMMEXEC ~/bin/go_charmm; " . $0;
#    $exec=$ENV{'CHARMMEXEC'};
    } #else {
#    $exec=&GenUtil::findExecutable("charmm");
#    }
    my ($pid,$cpusleft)=&GenUtil::submitRemote($hostlist->[$i],$srec,($jto-$jfrom+1),
					       $cmd,"$topt",$mp,!$keepmpdir);
    $jobsdone=$jto;
    push(@pidlist,$pid);
  }
} else {
  my $hostid=sprintf("%s.%s",hostname,$$);

  if ($mp) {
    if (defined $rex->{par}->{natpdb} && !-r $rex->{par}->{natpdb}) {
      $mpClient->getFile($rex->{par}->{natpdb},"$rex->{dir}/nat.pdb");
      $rex->set(natpdb=>"$rex->{dir}/nat.pdb");
    }
  }

  for (my $i=0; $i<=$#joblist; $i++) {
    my $pid=fork();
    if (!$pid) {
      my $rexClient=ReXClient->new($joblist[$i],$srec);
      die "client cannot connect to server"
	if (!defined $rexClient);

      my $lasttemp=$rexClient->initialize($hostid,(defined $cpus)?$cpus:$#joblist+1);
      $rexClient->establishConnection();

      &doJob($mp,$dir,$joblist[$i],$rex,$rexClient,
	     $lasttemp,$clogfile,$elogfile);
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
  my $lasttemp=shift;
  my $clogfile=shift;
  my $elogfile=shift;

  my $last={};
  $last->{temp}=$lasttemp;
  $last->{ener}="N/A";
  $last->{rmsd}="N/A";
  $last->{val}=();
  
  &GenUtil::makeDir("$dir/$job");

  my $initfile=$rexClient->initFile();

  if ($initfile eq "") {
    $rexClient->terminateServer();
    die "need initial structure for Go model replica exchange";
  }

  if ($mp) {
    $rexClient->getFile($initfile,"$dir/$job/init.pdb");
    $initfile="$dir/$job/init.pdb";
    $rexClient->getFile("$dir/$job/$ReXServer::restartfile") 
  }

  my $bias=$rexClient->biasInfo();

  my $charmm=&CHARMM::new((defined $clogfile)?"$dir/$job/$clogfile":undef);

  $charmm->verbose("eten on");
  $charmm->verbose("stream \"GO_" . $prefix . ".top\"");
  $charmm->verbose("stream \"GO_" . $prefix . ".param\"");
  $charmm->verbose("stream \"GO_" . $prefix . ".seq\"");
  $charmm->{molecule}=1;

  $charmm->setEnergyLogFile("$dir/$job/$elogfile") 
    if (defined $elogfile);

  $charmm->verbose("open unit 1 read form name \"" . $initfile . "\"");
  $charmm->verbose("read coor pdb unit 1");
  $charmm->verbose("close unit 1");

  $charmm->verbose("update cutnb 999.0 ctofnb 995.0 ctonnb 990.0");

  $charmm->shake(shake=>1,shakemode=>"all");

  $charmm->stream($custom->{setup}) 
      if (defined $custom && exists $custom->{setup});
  
  my $flist=();
  do {
    my $par={};
    $par->{ener}=$last->{ener};
    $par->{rmsd}=sprintf("%1.5lf",$last->{rmsd}) if ($last->{rmsd} ne "N/A");

    if ($#{$last->{val}}>=0) {
      $par->{val}=join(":",@{$last->{val}});
    }

    my $ret=$rexClient->nextCycle($flist,%{$par});

    if (!defined $ret) {

      $charmm->finish();
      return;
    }
    
    $charmm->stream($custom->{pre}) 
	if (defined $custom && exists $custom->{pre});
    $charmm->stream($custom->{"pre:init"}) 
	if (defined $custom && exists $custom->{"pre:init"} && $ret->{mode} eq "INIT");
    $charmm->stream($custom->{"pre:equi"}) 
	if (defined $custom && exists $custom->{"pre:equi"} && $ret->{mode} eq "EQUI");
    $charmm->stream($custom->{"pre:prod"}) 
	if (defined $custom && exists $custom->{"pre:prod"} && $ret->{mode} eq "PROD");
    
    my $charmmtrajfile=($rex->{opt}->{trajout})?"$dir/$job/traj.crd":undef;

    if ($ret->{mode} ne "INIT") {
      for (my $i=0; $i<=$#{$bias}; $i++) {
	$charmm->applyBias(%{$bias->[$i]},%{$ret->{bias}->[$i]});
      }
    }

    $charmm->runDynamics(undef, undef, $charmmtrajfile, undef, lang=>1,
                         dyntstep=>$mdpar{timestep}, langfbeta=>$mdpar{fbeta},
			 dynupdnb=>$rex->getPar()->{dynupdnb}, 
			 dynupdimg=>$rex->getPar()->{dynupdimg},
                         dynsteps=>$rex->getPar()->{steps}, dyntemp=>$ret->{temp}, 
			 dynoutfrq=>((defined $rex->getPar()->{outfrq})?
				$rex->getPar()->{outfrq}:$rex->getPar()->{steps}));

    my $tag=(lc $ret->{mode})."/".&GenUtil::dataDir($ret->{run});
    $tag .= ":".$ret->{temp};
    for (my $i=0; $i<=$#{$bias}; $i++) {
      $tag .= ":".$ret->{bias}->[$i]->{target};
    }
    $charmm->logMDEnergy($tag);

    for (my $i=0; $i<=$#{$bias}; $i++) {
      $charmm->clearBias(%{$bias->[$i]})
	if ($ret->{mode} ne "INIT");
      $last->{val}->[$i]=$charmm->getBiasVal(%{$bias->[$i]},%{$ret->{bias}->[$i]})+0.0;
    }

    $charmm->stream($custom->{post}) 
	if (defined $custom && exists $custom->{pre});
    $charmm->stream($custom->{"post:init"}) 
	if (defined $custom && exists $custom->{"post:init"} && $ret->{mode} eq "INIT");
    $charmm->stream($custom->{"post:equi"}) 
	if (defined $custom && exists $custom->{"post:equi"} && $ret->{mode} eq "EQUI");
    $charmm->stream($custom->{"post:prod"}) 
	if (defined $custom && exists $custom->{"post:prod"} && $ret->{mode} eq "PROD");

    $last->{ener}=$charmm->getEnergy()->{total}+0.0;
    $last->{temp}=$ret->{temp};

    $charmm->writePDB("$dir/$job/$ReXServer::outpdb");

    system "mv $dir/$job/new.restart $dir/$job/$ReXServer::restartfile"
      if (-r "$dir/$job/new.restart");

    if ($mp) {
      if (!defined $flist || $#{$flist}<0) {
	if (-r "$dir/$job/$ReXServer::restartfile") {
	  my $frec={};
	  $frec->{local}="$dir/$job/$ReXServer::restartfile";
	  $frec->{remote}=$frec->{local};
	  push(@{$flist},$frec);
	}
 
	my $frec={};
	$frec->{local}="$dir/$job/$ReXServer::outpdb";
	$frec->{remote}=$frec->{local};
	push(@{$flist},$frec);
 
#        if (defined $charmmtrajfile) {
#	  my $frec={};
#	  $frec->{local}="$charmmtrajfile";
#	  $frec->{remote}=$frec->{local};
#	  push(@{$flist},$frec);
#	}
      }
    }

    if (defined $rex->getNatPDB()) {
      if (!defined $rex->{obj}->{analyze}) {
	$rex->{obj}->{analyze}=Analyze::new($rex->getNatPDB());
      }

      my $exclmode=(defined $rex->{par}->{exclmode})?$rex->{par}->{exclmode}:1;
      $last->{rmsd}=undef;
    }

  } while(1);
}
