#!/usr/bin/env perl

# runs lattice simulations and rebuilds all-atom structures
# for ensembles
#
# http://mmtsb.scripps.edu/doc/enslatsim.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   enslatsim.pl [options]\n";
  printf STDERR "options: [-seq seqFile]\n";
  printf STDERR "         [-rnd | -chain file | -pdb file]\n";
  printf STDERR "         [-sa temp] [-const temp]\n";
  printf STDERR "         [-par tsteps=val,ncycle=val,icycle=val,\n";
  printf STDERR "               stiff=val,short=val,central=val,kdcore=val]\n";
  printf STDERR "         [-g gridsize] [-limforce value]\n";
  printf STDERR "         [-d force res1:res2[=res1:res2 ...]]\n";
  printf STDERR "         [-l refPDB min:max[=min:max ...]]\n";
  printf STDERR "         [-[no]compress]\n";
  printf STDERR "         [-natpdb pdbFile]\n";
  printf STDERR "         [-opt file[:file]]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-run [from:]to]\n";
  printf STDERR "         [-keeptraj]\n";
  printf STDERR "         [PARALLELoptions]\n";
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
use Sequence;
use Molecule;
use SICHO;
use MONSSTER;
use LatEnsemble;
use JobServer;
use JobClient;
use Client;

my $tag="lat";
my $dir=".";
my $srec={};
my $cpus=1;
my $hostfile;

my $natpdb;
my $seq;
my ($from,$to);
my ($chainfile, $pdbfile, $rndchain);
my $gridsize;
my ($fraglist,$fragref,$limforce);
my $limforce=1.0;
my ($drestforce, $drestlist);
my %simpar=();
my $keeptra;

my $optfile;

my $mp=0;
my $keepmpdir=0;

my $jobrank;
my $saveid;

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
  } elsif ($ARGV[0] eq "-keeptraj") {
    shift @ARGV;
    $keeptra=1;
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
  } elsif ($ARGV[0] eq "-saveid") {
    shift @ARGV;
    $saveid=shift @ARGV;
  } elsif ($ARGV[0] eq "-hosts") {
    shift @ARGV;
    $hostfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-mp") {
    shift @ARGV;
    $mp=1;
  } elsif ($ARGV[0] eq "-keepmpdir") {
    shift @ARGV;
    $keepmpdir=1;
  } elsif ($ARGV[0] eq "-opt") {
    shift @ARGV;
    $optfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-seq") {
    shift @ARGV;
    $seq=shift @ARGV;
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
  } elsif ($ARGV[0] eq "-chain") {
    shift @ARGV;
    $chainfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-rnd") {
    shift @ARGV;
    $rndchain=1;
  } elsif ($ARGV[0] eq "-pdb") {
    shift @ARGV;
    $pdbfile=shift @ARGV;

  } elsif ($ARGV[0] eq "-g") {
    shift @ARGV;
    $gridsize=shift @ARGV;
  } elsif ($ARGV[0] eq "-limforce") {
    shift @ARGV;
    $limforce=shift @ARGV;
  } elsif ($ARGV[0] eq "-d") {
    shift @ARGV;
    $drestforce=shift @ARGV;
    $drestlist=shift @ARGV;
  } elsif ($ARGV[0] eq "-sa") {
    shift @ARGV;
    $simpar{temp}=(shift @ARGV).":1.0"
  } elsif ($ARGV[0] eq "-const") {
    shift @ARGV;
    $simpar{temp}=(shift @ARGV);
    $simpar{temp}.=":$simpar{temp}";
    $simpar{central}=0.5;
    $simpar{tsteps}=1;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      $simpar{$key}=$val;
    }
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  }
}

&usage() if (!defined $rndchain && !defined $chainfile && !defined $pdbfile);

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

my $ens=LatEnsemble->new($tag,$dir);

if (defined $optfile) {
  foreach my $o ( split(/:/,$optfile) ) {
    $ens->readOptions($o);
  }
}

$ens->set(seq=>$seq, natpdb=>$natpdb, 
	  fragref=>$fragref, fraglist=>$fraglist);
$ens->set(compress=>$compress) if (defined $compress);
$ens->setOption(gridsize=>$gridsize, limforce=>$limforce, 
		drestforce=>$drestforce, drestlist=>$drestlist);
$ens->setPar(%simpar);

$ens->save() 
  if (!defined $srec->{name});

my $inpchain;
my $monsster;
if (!defined $hostfile) {
  if ($mp && defined $srec->{name}) {
    if (defined $chainfile && !&GenUtil::checkFile($chainfile)) {
      $mpClient->getFile($chainfile,"inpchain$$");
      $chainfile="inpchain$$";
    }

    if (defined $pdbfile && !&GenUtil::checkFile($pdbfile)) {
      $mpClient->getFile($pdbfile,"inppdb$$");
      $pdbfile="inppdb$$";
    }

    if (defined $ens->{par}->{fragref} && !&GenUtil::checkFile($ens->{par}->{fragref})) {
      $mpClient->getFile($ens->{par}->{fragref},"$ens->{dir}/fragref.pdb");
      $ens->set(fragref=>"$ens->{dir}/fragref.pdb");
    }
    
    if (defined $ens->{par}->{seq} && !&GenUtil::checkFile($ens->{par}->{seq})) {
      $mpClient->getFile($ens->{par}->{seq},"$ens->{dir}/seq");
      $ens->set(seq=>"$ens->{dir}/seq");
    }
    
    if (defined $ens->{par}->{natpdb} && !&GenUtil::checkFile($ens->{par}->{natpdb})) {
      $mpClient->getFile($ens->{par}->{natpdb},"$ens->{dir}/nat.pdb");
      $ens->set(natpdb=>"$ens->{dir}/nat.pdb");
    }
  }

  $inpchain=SICHO::new(gridsize=>$ens->{opt}->{gridsize});
  if (defined $chainfile) {
    $inpchain->readChain($chainfile);
  } elsif (defined $pdbfile) {
    my $mol=Molecule::new($pdbfile);
    $inpchain->genMONSSTERFromAllAtom($mol,fraglist=>$ens->getFragList());
  }
  my $nchain=$#{$inpchain->{sidechain}}+1;
  
  die "sequence file is missing"
    if (!defined $ens->getSeq());
  
  my $nres=$#{$ens->getSeq()->{sequence}}+1;
  
  die "sequence and chain lengths do not match (seq: $nres, chain: $nchain)"
    if ($nres != $nchain-2 && !$rndchain);
  
  $monsster=MONSSTER::new($ens->getSeq());
  $monsster->setParameters(%{$ens->getPar()});
  
  if (defined $ens->getFragList()) {
    $ens->getSeq()->setValidResidues($ens->getFragList(),1);
    $monsster->
      setPositionalRestraints(&GenUtil::gradForceList($ens->getSeq()->listFromValid(),
						      $ens->{opt}->{limforce}));
  }
  
  $monsster->setDistanceRestraints(
	$ens->{opt}->{drestforce},
	&GenUtil::fragListFromOption($ens->{opt}->{drestlist}))
    if (defined $ens->{opt}->{drestforce} && 
	defined $ens->{opt}->{drestlist});
}

if (($cpus>1 || defined $jobrank) && !defined $srec->{name}) {
  my $jlist=$ens->jobList($from,$to,"etot");

  die "nothing to do"
    if (!defined $jlist || $#{$jlist}<0);

  my $jobServer=JobServer->new($jlist,$ens,10);
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
      $topt.=" -chain $chainfile" if (defined $chainfile);
      $topt.=" -rnd" if ($rndchain);
      $topt.=" -pdb $pdbfile" if (defined $pdbfile);

      if ($mp) {
	$topt.=" -mp";
      } else {
	$hostlist->[$i]->{localdir}=$workdir;
      }
      
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
	
	my $job;
	my $lastprop;
	my $sendfiles;
	while (($job=$jobClient->nextJob($lastprop,$sendfiles))!~/^NO/) {
	  &dojob($job,$ens,$rndchain,$inpchain,$monsster,$keeptra);
	  
	  if ($mp) {
	    my $datadir=$ens->{dir}."/".&GenUtil::dataDir($job);
	    $jobClient->makeDir("$datadir");
	    $sendfiles=();
	    push(@{$sendfiles},
		 { local  => "$datadir/$ens->{tag}.pdb",
		   remote => "$datadir/$ens->{tag}.pdb" });

	    push(@{$sendfiles},
		 { local  => "$datadir/$MONSSTER::fileName{finalchain}",
		   remote => "$datadir/$MONSSTER::fileName{finalchain}" });

	    push(@{$sendfiles},
		 { local  => "$datadir/$MONSSTER::fileName{output}",
		   remote => "$datadir/$MONSSTER::fileName{output}" });
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

  waitpid($srec->{pid},0) 
    if (defined $srec->{pid});
} else {
  my $jlist=$ens->jobList($from,$to,"etot");

  die "nothing to do"
    if (!defined $jlist || $#{$jlist}<0);

  foreach my $i ( @{$jlist} ) {
    &dojob($i,$ens,$rndchain,$inpchain,$monsster,$keeptra);
  }
  $ens->save();
}

exit 0;

## dojob ######

sub dojob {
  my $job=shift;
  my $ens=shift;
  my $rndchain=shift;
  my $inpchain=shift;
  my $monsster=shift;
  my $keepTra=shift;

  if (defined $rndchain && $rndchain) {
    if (defined $ens->getFragList()) {
      $inpchain->
	genMONSSTERFromAllAtom($ens->getFragRef(),
			       fraglist=>$ens->getFragList());
    } else {
      $inpchain->genRandomMONSSTER($#{$ens->getSeq()->{sequence}}+1);
    }
  }

  my $datadir=$ens->{dir}."/".&GenUtil::dataDir($job);

  $monsster->setDirectory($datadir);
  $monsster->run($inpchain);

  my $rebmol=Molecule::new();

  $rebmol->rebuildFromSICHO($ens->getSeq(),$monsster->finalChain(),
			    $ens->{par}->{fraglist},$ens->{par}->{fragref});

  my $rebpdb=$datadir."/$ens->{tag}.pdb";
  $rebmol->writePDB($rebpdb,translate=>"CHARMM22",ssbond=>0);

  $ens->update($job);
  $ens->cleanUp($job,$keepTra);
}
   


