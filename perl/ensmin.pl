#!/usr/bin/env perl

# runs minimizations with CHARMM for ensembles
#
# http://mmtsb.scripps.edu/doc/ensmin.pl.html
# 2001, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   ensmin.pl [options] intag outtag\n";
  printf STDERR "options: [-par CHARMMparams]\n";
  printf STDERR "         [-opt file[:file]]\n";
  printf STDERR "         [-l refPDB min:max[=min:max ...]]\n";
  printf STDERR "         [-[no]conslim] [-limforce value]\n";
  printf STDERR "         [-limsel=ca|cb|cab|heavy]\n";
  printf STDERR "         [-cons [ca|cb|cab|heavy] self|refpdb min:max[_force][=...]]\n";
  printf STDERR "         [-run [from:]to]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-natpdb pdbFile]\n";
  printf STDERR "         [-[no]compress]\n";
  printf STDERR "         [-update frq]\n";
  printf STDERR "         [PARALLELoptions]\n";
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
use CHARMM;
use Ensemble;
use JobServer;
use JobClient;
use Client;

my ($intag,$outtag);
my $dir=".";
my $srec={};
my $cpus=1;
my $hostfile;

my $natpdb;
my ($from,$to);

my ($fraglist,$fragref);

my $logfile;
my %par=();
my $cons;

my %paropt=();

my %defopt = (
 limforce  => 1.0,
 limsel    => "cab",
 conslim   => 0
);

my %defpar = (); 

my $optfile;

my $mp=0;
my $keepmpdir=0;

my $jobrank;
my $saveid;

my $update=10;

my $compress;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
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
  } elsif ($ARGV[0] eq "-update") {
    shift @ARGV;
    $update=shift @ARGV;
  } elsif ($ARGV[0] eq "-opt") {
    shift @ARGV;
    $optfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-natpdb") {
    shift @ARGV;
    $natpdb=shift @ARGV;
  } elsif ($ARGV[0] eq "-nocompress") {
    shift @ARGV;
    $compress=0;
  } elsif ($ARGV[0] eq "-compress") {
    shift @ARGV;
    $compress=1;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fragref=shift @ARGV;
    $fraglist=shift @ARGV;
  } elsif ($ARGV[0] eq "-mp") {
    shift @ARGV;
    $mp=1;
  } elsif ($ARGV[0] eq "-keepmpdir") {
    shift @ARGV;
    $keepmpdir=1;
  } elsif ($ARGV[0] eq "-saveid") {
    shift @ARGV;
    $saveid=shift @ARGV;
  } elsif ($ARGV[0] eq "-conslim") {
    shift @ARGV;
    $paropt{conslim}=1;
  } elsif ($ARGV[0] eq "-noconslim") {
    shift @ARGV;
    $paropt{conslim}=0;
  } elsif ($ARGV[0] eq "-limforce") {
    shift @ARGV;
    $paropt{limforce}=shift @ARGV;
  } elsif ($ARGV[0] eq "-limsel") {
    shift @ARGV;
    $paropt{limsel}=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logfile=shift @ARGV;
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
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    if (!defined $intag) {
      $intag=shift @ARGV;
    } elsif (!defined $outtag) {
      $outtag=shift @ARGV;
    }
  }
}

&usage() if (!defined $intag || !defined $outtag);

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

$paropt{conslim}=1
  if (!defined $paropt{conslim} && (defined $paropt{limsel} || defined $paropt{limforce}));

my $mpClient;
if ($mp && defined $srec->{name}) {
  $mpClient=Client->new("mp",$srec);
  &GenUtil::makeDir($dir);
  $mpClient->getFile("$dir/ens.cfg");
  $mpClient->getFile("$dir/$outtag$Ensemble::OptionFileSuffix");
}

my $ens=Ensemble->new($outtag,$dir);

if (defined $optfile) {
  foreach my $o ( split(/:/,$optfile) ) {
    $ens->readOptions($o);
  }
}

$ens->set(natpdb=>$natpdb, fragref=>$fragref, fraglist=>$fraglist);
$ens->set(compress=>$compress) if (defined $compress);

$ens->setOption(%paropt);

if (defined $cons) {
  my $consstr="$cons->{sel},$cons->{type},$cons->{list},$cons->{reffile},$cons->{exclmode}";
  $ens->setOption(cons=>$consstr);
}

$ens->setPar(%par);

foreach my $p ( keys %defpar ) {
  if (!defined $ens->getPar() || !defined $ens->getPar()->{$p}) {
    $ens->setPar($p=>$defpar{$p});
  }
}

$ens->{opt}->{conslim}=1
  if (!defined $ens->{opt}->{conslim} && defined $ens->getFragList() && 
      !defined $ens->{opt}->{cons});

foreach my $p ( keys %defopt ) {
  if (!defined $ens->{opt} || !defined $ens->{opt}->{$p}) {
    $ens->setOption($p=>$defopt{$p});
  }
}

$ens->save()
  if (!defined $srec->{name});

my $allcons=();

if (defined $ens->{opt}->{cons}) {
  my @cstr=split(/,/,$ens->{opt}->{cons});
  my $c={ sel=>$cstr[0],
	  type=>$cstr[1],
	  list=>&GenUtil::fragListFromOption($cstr[2]),
          reffile=>$cstr[3],
          exclmode=>$cstr[4] };
  push(@{$allcons},$c);
}

if (defined $ens->getFragList() && $ens->{opt}->{conslim}) {
  my $c={ sel=>(defined $ens->{opt}->{limsel})?$ens->{opt}->{limsel}:"heavy" ,
	  force=>(defined $ens->{opt}->{limforce})?$ens->{opt}->{limforce}:1.0,
	  type=>"ref",
	  reffile=>$ens->{par}->{fragref},
	  list=>$ens->getFragList(),
          exclmode=>1 };
  push(@{$allcons},$c);
}

if (($cpus>1 || defined $jobrank) && !defined $srec->{name}) {
  my $jlist=$ens->jobList($from,$to,"etot");

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

      if ($mp) {
	$topt.=" -mp";
      } else {
	$hostlist->[$i]->{localdir}=$workdir;
      }

      $topt.=" $intag $outtag";

      ($pid,$cpusleft)=&GenUtil::submitRemote($hostlist->[$i],$srec,$cpusleft,
					      (split(/\//,$0))[-1],$topt,$mp,!$keepmpdir);
      push(@pidlist,$pid);
    }
  } else {
    if ($mp) {
      my $nref=1;
      foreach my $ac ( @{$allcons} ) {
	if (defined $ac->{reffile} && !-r $ac->{reffile}) {
	  $mpClient->getFile($ac->{reffile},"$ens->{dir}/ref$nref");
	  $ac->{reffile}="$ens->{dir}/ref$nref";
	  $nref++;
	}
      }

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

	my $charmm=&CHARMM::new( (defined $logfile)?"$i-$logfile":undef );
	$charmm->loadParameters(%{$ens->getPar()});
      
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

	  &dojob($job,$ens,$intag,$charmm,$allcons);

	  if ($mp) {
	    $sendfiles=();
	    push(@{$sendfiles},
		 { local  => "$datadir/$outtag.pdb",
		   remote => "$datadir/$outtag.pdb" });
	    push(@{$sendfiles},
		 { local  => "$datadir/$outtag.elog",
		   remote => "$datadir/$outtag.elog" });
	  }

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
  my $jlist=$ens->jobList($from,$to,"etot");

  die "nothing to do"
    if (!defined $jlist || $#{$jlist}<0);

  my $inx=0;
  my $charmm=&CHARMM::new($logfile);
  $charmm->loadParameters(%{$ens->getPar()});
  foreach my $i ( @{$jlist} ) {
    &dojob($i,$ens,$intag,$charmm,$allcons);
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
  my $intag=shift;
  my $charmm=shift;
  my $cons=shift;

  my $par=$ens->getPar();
  
  my $datadir=$ens->{dir}."/".&GenUtil::dataDir($job);

  my $inpfile=$datadir."/".$intag.".pdb";

  if (!$charmm->{_setup}) {
    $charmm->setupFromPDB($inpfile);
  } else {
    $charmm->clearRestraints() if ($#{$cons}>=0);
    $charmm->initCoordinates();
    $charmm->readFromPDB($inpfile);
    $charmm->clearEnergy();
    $charmm->clearShake();
  }

  $charmm->setupEnergy();
  $charmm->shake();

  $charmm->setEnergyLogFile($datadir."/".$ens->{tag}.".elog");

  $charmm->setupRestraints(1.0,$cons)
    if ($#{$cons}>=0);

  if ($charmm->{par}->{sdsteps}>0) {
    $charmm->minimizeSD();
    $charmm->logEnergy("SD");
  }

  $charmm->minimize();
  $charmm->logEnergy("Min final");

  my $chmoutpdb=$datadir."/".$intag.".chmout.pdb";
  $charmm->writePDB($chmoutpdb);

  my $outmol=Molecule::new();
  $outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),chainfromseg=>1);
  $outmol->setSSBonds($charmm->{molecule}->getSSBonds());
  $outmol->writePDB($datadir."/$ens->{tag}.pdb",translate=>"CHARMM22");

  &GenUtil::remove($chmoutpdb);

  $ens->update($job);
  $ens->cleanUp($job);
}
