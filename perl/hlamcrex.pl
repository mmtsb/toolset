#!/usr/bin/env perl
#
# Lattice/all-atom replica exchange simulations
#
# http://mmtsb.scripps.edu/doc/hlamcrex.pl.html
# 2001, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   hlamcrex.pl [options] [files]\n";
  printf STDERR "options: [-n runs]\n";
  printf STDERR "         [-f listfile]\n";
  printf STDERR "         [-par initruns=value,equilruns=value,\n";
  printf STDERR "               [no]save,savebestfreq=value,\n";
  printf STDERR "               natpdb=file,seq=file]\n";
  printf STDERR "         [-mcpar mcruns=value,mode=coupled|free,\n";
  printf STDERR "                 [no]keepchain,[no]eval,\n";
  printf STDERR "                 caforce=value,hmcmforce=force]\n";
  printf STDERR "         [-latpar icycle=val,ncycle=value,\n";
  printf STDERR "                  stiff=val,short=val,central=val,kdcore=val]\n";
  printf STDERR "         [-latopt gridsize=value]\n";
  printf STDERR "         [-aapar  CHARMMparams]\n";
  printf STDERR "         [-evalpar CHARMMparams]\n";
  printf STDERR "         [-temp nwin:min:max]\n";
  printf STDERR "         [-ltemp min:max]\n";
  printf STDERR "         [-condfile file]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-ens tag] [-ensdir dir]\n";
  printf STDERR "         [PARALLELoptions]\n";
  printf STDERR "         [-log file] [-aalog file] [-evallog file]\n";
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
use LMCReXServer;
use Client;
use ReXClient;
use GenUtil;
use Molecule;
use Ensemble;
use CHARMM;
use Sequence;
use LatEnsemble;
use SICHO;
use MONSSTER;

my %par=();
my %defpar = (
   initruns     => 0,
   equilruns    => 0,
   save         => 1,
   savebestfreq => -1
);

my %mcpar=();
my %defmcpar = (
  caforce   => 0.5,
  hmcmforce => 0.2,
  mcruns    => 10,
  mode      => "free",
  keepchain => 1,		 
  eval      => 1
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
my $alogfile;

my %aapar=();
my %defaapar = ( 
 sdsteps   =>  50,		 
 minsteps  =>  100,
 dyntemp   =>  298,
 dynsteps  =>  0,
 dynoutfrq =>  20
); 

my %evalpar=();
my %defevalpar = (); 

my %latopt=();
my %deflatopt = (
 gridsize  => 100
);

my %latpar=();
my %deflatpar = ( 
 icycle  =>  20,
 ncycle  =>  1
);

my $serverlog;

my $dir=".";

my $lmint;
my $lmaxt;

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
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-mcpar") {
    shift @ARGV;
    &GenUtil::parsePar(\%mcpar,shift @ARGV);
  } elsif ($ARGV[0] eq "-ltemp") {
    shift @ARGV;
    ($lmint,$lmaxt)=split(/:/,shift @ARGV);
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
  } elsif ($ARGV[0] eq "-aapar") {
    shift @ARGV;
    &GenUtil::parsePar(\%aapar,shift @ARGV);
  } elsif ($ARGV[0] eq "-evalpar") {
    shift @ARGV;
    &GenUtil::parsePar(\%evalpar,shift @ARGV);
  } elsif ($ARGV[0] eq "-latopt") {
    shift @ARGV;
    &GenUtil::parsePar(\%latopt,shift @ARGV);
  } elsif ($ARGV[0] eq "-latpar") {
    shift @ARGV;
    &GenUtil::parsePar(\%latpar,shift @ARGV);
  } elsif ($ARGV[0] eq "-aalog") {
    shift @ARGV;
    $alogfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-evallog") {
    shift @ARGV;
    $elogfile=shift @ARGV;
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
  $mpClient->getFile("$dir/aarexclient.cfg");
  $mpClient->getFile("$dir/aarexclient.options");
  $mpClient->getFile("$dir/evalrexclient.cfg");
  $mpClient->getFile("$dir/evalrexclient.options");
  $mpClient->getFile("$dir/latrexclient.cfg");
  $mpClient->getFile("$dir/latrexclient.options");
}

my $aarex=Ensemble->new("aarexclient",$dir,"aarexclient.cfg");
$aarex->set(runs=>$nruns,natpdb=>$par{natpdb});
$aarex->setPar(%aapar);
$aarex->setDefPar(%defaapar);
$aarex->setOption(%mcpar);
$aarex->setDefOption(%defmcpar);

my $latrex=Ensemble->new("latrexclient",$dir,"latrexclient.cfg");
$latrex->set(runs=>$nruns,seq=>$par{seq});
$latrex->setOption(%latopt);
$latrex->setDefOption(%deflatopt);
$latrex->setPar(%latpar);
$latrex->setDefPar(%deflatpar);

my $evalrex=Ensemble->new("evalrexclient",$dir,"evalrexclient.cfg");
$evalrex->setPar(%evalpar);
$evalrex->setDefPar(%defevalpar);



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

  if (defined $nwindows && defined $lmint && defined $lmaxt) {
    my $ltlist=&ReXServer::expTempList($nwindows,$lmint,$lmaxt);
    $par{ltemplist}=join(":",@{$ltlist});
  }
  
  my $server=LMCReXServer->new($latrex->{par}->{runs},$dir,\@initfile,%par);
  
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
    $server->run(4200);
  
  $srec->{name}=hostname;

  $nwindows=$server->nWindows();

  $aarex->save();
  $evalrex->save();
  $latrex->save();
  
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
    push(@joblist,sprintf("ls%d",$i));
  }
} elsif (defined $nwindows) {
  for (my $i=1; $i<=$nwindows; $i++) {
    push(@joblist,sprintf("ls%d",$i));
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
    $topt.=" -aalog $alogfile" if (defined $alogfile);
    $topt.=" -evallog $elogfile" if (defined $elogfile);
    $topt.=" -dir $dir";

    if ($mp) {
      $topt.=" -mp";
    } else {
      $hostlist->[$i]->{localdir}=$workdir;
    }

    my ($pid,$cpusleft)=&GenUtil::submitRemote($hostlist->[$i],$srec,($jto-$jfrom+1),
					       (split(/\//,$0))[-1],"$topt",$mp,!$keepmpdir);
    $jobsdone=$jto;
    push(@pidlist,$pid);
  }
} else {
  my $hostid=sprintf("%s.%s",hostname,$$);

  if ($mp) {
    if (defined $aarex->{par}->{natpdb} && !-r $aarex->{par}->{natpdb}) {
      $mpClient->getFile($aarex->{par}->{natpdb},"$aarex->{dir}/nat.pdb");
      $aarex->set(natpdb=>"$aarex->{dir}/nat.pdb");
    }
    if (defined $latrex->{par}->{seq} && !-r $latrex->{par}->{seq}) {
      $mpClient->getFile($latrex->{par}->{seq},"$latrex->{dir}/seq");
      $latrex->set(seq=>"$latrex->{dir}/seq");
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

      &doJob($mp,$dir,$joblist[$i],$aarex,$evalrex,$latrex,$rexClient,
	     $lasttemp,$alogfile,$elogfile);
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
  my $aarex=shift;
  my $evalrex=shift;
  my $latrex=shift;
  my $rexClient=shift;
  my $lasttemp=shift;
  my $alogfile=shift;
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
    die "need initial structure for all-atom replica exchange";
  }

  if ($mp) {
    $rexClient->getFile($initfile,"$dir/$job/init.pdb");
    $initfile="$dir/$job/init.pdb";
  }

  my $mcharmm=&CHARMM::new((defined $alogfile)?"$dir/$job/$alogfile":undef);
  $mcharmm->loadParameters(%{$aarex->getPar()});
  $mcharmm->setupFromPDB($initfile);
  $mcharmm->setupEnergy();

  my $echarmm;
  
  if ($aarex->{opt}->{"eval"} != 0 ) {
    $echarmm=&CHARMM::new((defined $elogfile)?"$dir/$job/$elogfile":undef);
    $echarmm->loadParameters(%{$evalrex->getPar()});
    $echarmm->setupFromPDB($initfile);
    $echarmm->setupEnergy();
  }

  my $first=1;
  my $outmol=Molecule::new();

  my $monsster=MONSSTER::new($latrex->getSeq());
  $monsster->setDirectory("$dir/$job");
  $monsster->setParameters(%{$latrex->getPar()});
  $monsster->setParameters(tsteps=>1);

  my $inpchain=SICHO::new(gridsize=>$latrex->{opt}->{gridsize});
  my $initmol=Molecule::new($initfile);
  $initmol->center();
  $inpchain->genMONSSTERFromAllAtom($initmol);

  my $cons=();

  if ($aarex->{opt}->{mode} eq "coupled" && $aarex->{opt}->{hmcmforce}>0.0) {
    my $c={};
    $c->{type}="hmcm";
    $c->{list}=&GenUtil::fragListFromOption(sprintf("-9999:9999_%1.1f",$aarex->{opt}->{hmcmforce}));
    $c->{reffile}="$dir/$MONSSTER::fileName{finalchain}";
    push(@{$cons},$c);
  } elsif ($aarex->{opt}->{mode} eq "free" && $aarex->{opt}->{caforce}>0.0) {
    my $c={};
    $c->{sel}="ca";
    $c->{type}="self";
    my $caforce=$aarex->{par}->{caforce};
    $c->{list}=&GenUtil::fragListFromOption(sprintf("-9999:9999_%1.1f",$aarex->{opt}->{caforce}));
    $c->{exclmode}=0;
    push(@{$cons},$c);
  }

  my $flist=();
  do {
    my $par={};
    $par->{ener}=$last->{ener};
    $par->{rmsd}=sprintf("%1.5f",$last->{rmsd}) if ($last->{rmsd} ne "N/A");

    if ($#{$last->{val}}>=0) {
      $par->{val}=join(":",@{$last->{val}});
    }

    my $ret=$rexClient->nextCycle($flist,%{$par});

    if (!defined $ret) {
      $mcharmm->finish();
      $echarmm->finish()
	if (defined $echarmm);
      return;
    }

    my $ltemp=$ret->{bias}->[0]->{ltemp};

    $monsster->setParameters(temp=>$ltemp);

    my $nmcruns=0;

#    printf STDERR "client $$ running lattice temp %f, MC temp: %f, run: %d\n",
#    $ltemp,$ret->{temp},$ret->{run};

    do {
      $monsster->run($inpchain);

      $mcharmm->clearEnergy() if (!$first);
      $mcharmm->clearRestraints() if (!$first);
      $mcharmm->clearShake() if (!$first);

      my $outchain;
      
      if ($aarex->{opt}->{mode} eq "free") {
	$outchain=SICHO::new();
	$outchain->readChain("$dir/$job/$MONSSTER::fileName{finalchain}");

	my $rebmol=Molecule::new();
	$rebmol->rebuildFromSICHO($latrex->getSeq(),$outchain);
	$rebmol->writePDB("$dir/$job/rebuilt.pdb",translate=>"CHARMM22");

	$mcharmm->initCoordinates();
	$mcharmm->readFromPDB("$dir/$job/rebuilt.pdb");
      }

      $mcharmm->setupEnergy();

      if ($aarex->{opt}->{mode} eq "free") {
	$mcharmm->setupRestraints(1.0,$cons);
      } else {
	$mcharmm->hmcmRestraint($cons->[0]);
      }    
      
      $mcharmm->minimizeSD() if ($mcharmm->{par}->{sdsteps}>0);
      $mcharmm->minimize()   if ($mcharmm->{par}->{minsteps}>0);

      my $me=$mcharmm->_processEneOutput()
	if ($mcharmm->{par}->{minsteps}+$mcharmm->{par}->{sdsteps}>0);

      if ($mcharmm->{par}->{dynsteps}>0 && 
	  (!defined $me ||
	   ($me->[$#{$me}]->{total}<0.0 && $mcharmm->{_lastOutput}!~/far from minimum/))) {
	$mcharmm->shake();
	$mcharmm->runDynamics(undef,undef,undef,undef,
			      dyntwin=>0.02*$mcharmm->{par}->{dyntemp}); 
      }

      my $chmoutpdb="$dir/$job/chmout.pdb";
      $mcharmm->writePDB($chmoutpdb);
	
      my $etrial;
      if ($aarex->{opt}->{"eval"}!=0) {
	$echarmm->clearEnergy();
	$echarmm->initCoordinates();
	$echarmm->readFromPDB("$dir/$job/chmout.pdb");
	$echarmm->setupEnergy();
	$etrial=$echarmm->getEnergy()->{total}+0.0;
      } else {
	my $e=$mcharmm->getEnergy();
	$etrial=$e->{total}-$e->{constr};
      }
 
      if ($first || $etrial<$last->{ener} ||
	  rand()<exp(-($etrial-$last->{ener})/($ret->{temp}))) {
	$last->{ener}=$etrial;
	$outmol->readPDB($chmoutpdb,
			 translate=>&CHARMM::getConvType($mcharmm->{par}->{param}),chainfromseg=>1);
	&GenUtil::remove($chmoutpdb);

	if ($aarex->{opt}->{keepchain}) {
	  $inpchain=$outchain;
	} else {
	  $inpchain=SICHO::new(gridsize=>$latrex->{opt}->{gridsize});
	  $outmol->center();
	  $inpchain->genMONSSTERFromAllAtom($outmol);
	  if ($#{$inpchain->{sidechain}}<0) {
	    my $tf=sprintf "tmp-%d.pdb",$job;
	    printf STDERR "no chain! writing out $tf\n";
	    $outmol->setSSBonds($mcharmm->{molecule}->getSSBonds());
	    $outmol->writePDB($tf);
	  }
	}
      }
      $first=0;
    } while (++$nmcruns<$aarex->{opt}->{mcruns});

    $outmol->setSSBonds($mcharmm->{molecule}->getSSBonds());
    $outmol->writePDB("$dir/$job/$ReXServer::outpdb",translate=>"CHARMM22");

    if ($mp) {
      if (!defined $flist || $#{$flist}<0) {
	my $frec={};
	$frec->{local}="$dir/$job/$ReXServer::outpdb";
	$frec->{remote}=$frec->{local};
	push(@{$flist},$frec);
      }
    }

    if (defined $aarex->getNatPDB()) {
      if (!defined $aarex->{obj}->{analyze}) {
	$aarex->{obj}->{analyze}=Analyze::new($aarex->getNatPDB());
      }
      
      $aarex->{obj}->{analyze}->lsqfit($outmol,"cab",0);
      my $rmsd=$aarex->{obj}->{analyze}->rmsd($outmol,0);
      $last->{rmsd}=$rmsd->{CA};
    }

  } while(1);
}
