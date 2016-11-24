#!/usr/bin/env perl
#
# All-atom replica exchange simulations
#
# http://mmtsb.scripps.edu/doc/aarexAmber.pl.html
# 2001, Michael Feig, Brooks group, TSRI
# 2001, John Karanicolas, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   aarexAmber.pl [options] [files]\n";
  printf STDERR "options: [-n runs]\n";
  printf STDERR "         [-par initruns=value,equilruns=value,\n";
  printf STDERR "               [no]save,savebestfreq=value,archive\n";
  printf STDERR "               ensmode=add|replace,natpdb=file,partop=file\n";
  printf STDERR "         [-temp nwin:min:max]\n";
  printf STDERR "         [-condfile file]\n";
  printf STDERR "         [-f listfile]\n";
  printf STDERR "         [-mdpar AmberParams]\n";
  printf STDERR "         [-mdopt [no]trajout,[no]conslim,\n";
  printf STDERR "                 limforce=value,limsel=ca|cb|cab|heavy],\n";
  printf STDERR "                 [no]translate,ambercoor\n";
  printf STDERR "         [-l refPDB min:max[=min:max ...]]\n";
  printf STDERR "         [-cons [ca|cb|cab|heavy] ref|self min:max[_force][=...]]\n";
  printf STDERR "         [-opt optionsfile]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-ens tag] [-ensdir dir]\n";
  printf STDERR "         [PARALLELoptions]\n";
  printf STDERR "         [-log file] [-elog file] [-amberlog file]\n";
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
use AmberReXServer;
use Client;
use ReXClient;
use GenUtil;
use Molecule;
use Ensemble;
use Amber;

my %par=();
my %defpar = (
   initruns     => 0,
   equilruns    => 0,
   save         => 1,
   savebestfreq => -1,
   partop       => undef
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

my $prun=0;

my $saveid;

my $elogfile;
my $clogfile;

my $cons;

my %mdopt=();
my %defmdopt = (
 limforce  => 1.0,
 limsel    => "cab",
 conslim   => 0,
 trajout   => 0,
 translate => 1,
 ambercoor => 0
);

my %mdpar=();
my %defmdpar = ( 
 param     =>  22,
 gb        =>  1,
 cutoff    =>  16.0,
 dynsteps  =>  1000,
 rex       =>  1
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
  } elsif ($ARGV[0] eq "-mdpar") {
    shift @ARGV;
    &GenUtil::parsePar(\%mdpar,shift @ARGV);
  } elsif ($ARGV[0] eq "-mdopt") {
    shift @ARGV;
    &GenUtil::parsePar(\%mdopt,shift @ARGV);
  } elsif ($ARGV[0] eq "-elog") {
    shift @ARGV;
    $elogfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-amberlog") {
    shift @ARGV;
    $clogfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $par{fragref}=shift @ARGV;
    $par{fraglist}=shift @ARGV;
  } elsif ($ARGV[0] eq "-cons") {
    shift @ARGV;
    $cons={};
    if ($ARGV[0] =~ /^(ca|cb|cab|heavy)$/) {
      $cons->{sel}=shift @ARGV;
    } else {
      $cons->{sel}="heavy";
    }
    if ($ARGV[0] eq "self") {
      shift @ARGV;
      $cons->{type}="self";
    } else {
      $cons->{type}="ref";
      $cons->{reffile}=shift @ARGV;
    }
    $cons->{list}=shift @ARGV;
    $cons->{exclmode}=0;
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    my $ctag=shift @ARGV;
    $customfile{$ctag}=shift @ARGV;
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

$mdopt{conslim}=1 
  if (!defined $mdopt{conslim} && (defined $mdopt{limsel} || defined $mdopt{limforce}));

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

$rex->set(fragref=>$par{fragref},fraglist=>$par{fraglist},
	  runs=>$nruns,natpdb=>$par{natpdb},partop=>$par{partop});
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

$rex->{opt}->{conslim}=1
  if (!defined $rex->{opt}->{conslim} && defined $rex->getFragList() && 
      !defined $rex->{opt}->{cons});

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
      &GenUtil::log("rexserver","adding init file $_");
    }
    close $lfile;
  }
  
  my $server=AmberReXServer->new($rex->{par}->{runs},$dir,\@initfile,%par);
  
  if (defined $server->{par}->{enstag}) {
    my $ens=Ensemble->new($server->{par}->{enstag},$server->{par}->{ensdir});
    $ens->set(fraglist=>$rex->{par}->{fraglist});
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

if (defined $clogfile) {
  foreach my $j (@joblist) {
    if (&GenUtil::checkFile("$dir/$j/$clogfile")) {
    DONE:
      for (my $i=1; $i<=99999; $i++) {
	if (!&GenUtil::checkFile("$dir/$j/$clogfile.$i")) {
	  system "mv $dir/$j/$clogfile $dir/$j/$clogfile.$i";
	  last DONE;
	}
      }
    }
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

  for (my $i=0; $i<=$#{$hostlist} && $jobsdone<$#joblist+1; $i++) {
    my $maxcpus=&GenUtil::remoteCPUs($hostlist->[$i]);
    my $jfrom=$jobsdone+1;
    my $jto=$jobsdone+$maxcpus;
    $jto=$#joblist+1 if ($jto>$#joblist+1);

    my $topt="";
    $topt.=" -jobs $jfrom:$jto";
    $topt.=" -log $serverlog-$i" if (defined $serverlog);
    $topt.=" -amberlog $clogfile" if (defined $clogfile);
    $topt.=" -elog $elogfile" if (defined $elogfile);
    $topt.=" -dir $dir";      

    if ($mp) {
      $topt.=" -mp";
    } else {
      $hostlist->[$i]->{localdir}=$workdir;
    }

    my ($pid,$cpusleft)=&GenUtil::submitRemote($hostlist->[$i],$srec,($jto-$jfrom+1),
					       (split(/\//,$0))[-1],"$topt",$mp,!$keepmpdir,
					       $clogfile,$dir);
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
#      $rexClient->establishConnection();

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
  $last->{boxsize}="N/A";
  $last->{volume}="N/A";
  $last->{val}=();
  
  &GenUtil::makeDir("$dir/$job");

  my $partop=$rex->{par}->{partop};

  my $initfile=$rexClient->initFile();

  if ($initfile eq "") {
    $rexClient->terminateServer();
    die "need initial structure for all-atom replica exchange";
  }

  if ($mp) {
    if (!&GenUtil::checkFile($initfile)) {
      $rexClient->getFile($initfile,"$dir/$job/init.pdb");
      $initfile="$dir/$job/init.pdb";
    }
    $rexClient->getFile("$dir/$job/$ReXServer::restartfile");
    if (!&GenUtil::checkFile($partop)) {
      $rexClient->getFile($partop,"$dir/$job/partop");
      $partop="$dir/$job/partop";
    }
  }

  die "need partop file" if (!defined $partop || !&GenUtil::checkFile($partop));

  my $bias=$rexClient->biasInfo();

  my $allcons=();

#  if (defined $rex->{opt}->{cons}) {
#    my @cstr=split(/,/,$rex->{opt}->{cons});
#    my $c={ sel=>$cstr[0],
#	    type=>$cstr[1],
#	    list=>&GenUtil::fragListFromOption($cstr[2]),
#	    reffile=>$cstr[3],
#	    exclmode=>$cstr[4] };
#    push(@{$allcons},$c);
#  } 

#  if (defined $rex->getFragList() && $rex->{opt}->{conslim}) {
#    my $c={ sel=>(defined $rex->{opt}->{limsel})?$rex->{opt}->{limsel}:"heavy" ,
#	    force=>(defined $rex->{opt}->{limforce})?$rex->{opt}->{limforce}:1.0,
#	    type=>"ref",
#	    reffile=>$rex->{par}->{fragref},
#	    list=>$rex->getFragList(),
#	    exclmode=>1 
#          };
#
#    push(@{$allcons},$c);
#  }

#  if ($mp) {
#    my $nref=1;
#    foreach my $ac ( @{$allcons} ) {
#      if (defined $ac->{reffile} && !-r $ac->{reffile}) {
#	$rexClient->getFile($ac->{reffile},"$dir/$job/ref$nref");
#	$ac->{reffile}="$dir/$job/ref$nref";
#	$nref++;
#      }
#    }
#  }

  my $ambinpcmd=lc "$dir/$job/$$.inp.cmd";

  my $ambinpcoor=lc "$dir/$job/$$.inpcoor";
  my $restart;
  if (&GenUtil::checkFile("$dir/$job/$ReXServer::restartfile")) {
    system "cp $dir/$job/$ReXServer::restartfile $ambinpcoor";
    $restart=1;
  } else {
    if ($rex->{opt}->{ambercoor}) {
      system "cp $initfile $ambinpcoor";
    } else {
      my $mol=Molecule::new();
      $mol->readPDB($initfile);
      $mol->fixHistidine($rex->{opt}->{hsd},$rex->{opt}->{hse});
      $mol->translate("AMBER") if ($rex->{opt}->{translate});
      $mol->writeAmber($ambinpcoor);
    }
    $restart=0;
  }

  my $amber=Amber::new();

  my $amberinp=new IO::File;
  $amberinp->open(">$ambinpcmd");
  $rex->setPar(dyntemp=>$lasttemp) if (defined $lasttemp && $lasttemp>0);
  if (!defined $rex->getPar()->{rexmode}) {
    if ($#{$bias}<0) {
      $rex->setPar(rexmode=>1);
    } elsif ($bias->[0]->{type} eq "lambda") {
      $rex->setPar(rexmode=>2);
    } else {
      die "unknown bias type";
    }
  }

  $amber->genInputDynamics($amberinp,($restart),($rex->{opt}->{trajout}),%{$rex->getPar()});

  my $ambrefcoor="$dir/$job/$$.ref.x";

  if ($#{$allcons}>=0) {
#    my $refmol=Molecule::new();
#    $refmol->readAmber($partop,$ambinpcoor);
#    for (my $ic=0; $ic<=$#{$mol->{chain}}; $ic++) {
#      $refmol->numberReset($map[$ic],$mol->{chain}->[$ic]->{id});
#    }
#    
#    if ($cons->{type} eq "ref" && defined $cons->{reffile}) {
#      $refmol->zapCoordinates();
#      $refmol->fillCoorFromPDB($cons->{reffile});
#    }
#    $refmol->writeAmber($ambrefcoor);
#    
#    if (defined $cons->{exclmode} && $cons->{exclmode}) {
#      $refmol->setValidResidues($cons->{list},1);
#      $cons->{list}=&GenUtil::gradForceList($refmol,$cons->{force});
#    }	
#    $amber->genInputRestraints($amberinp,$cons);
  }

  undef $amberinp;

  my $rexcontrol="$dir/$job/amberrex.job";
  open OUT,">$dir/$job/amberrex.job";
  printf OUT "%d\n%d\n%s\n%s\n%s\n%s\n%s\n",
  0,$mp,$job,$rexClient->{serverName},
  $rexClient->{serverPort},$rexClient->{serverID},
  "$dir/$job";
  close OUT;

  my $amboutlog=(defined $clogfile)?"$dir/$job/$clogfile":"$dir/$job/$$.log";
  my $amboutcoor="$dir/$job/$ReXServer::restartfile";

  $amber->runSander(input=>$ambinpcmd, partop=>$partop,inpcoor=>$ambinpcoor, 
		    outcoor=>$amboutcoor,log=>$amboutlog,restraint=>($#{$allcons}>=0)?$ambrefcoor:undef,
		    trajout=>($rex->{opt}->{trajout})?"$dir/$job/$AmberReXServer::ambertrajout":undef,
		    rexcontrol=>$rexcontrol);
}  

