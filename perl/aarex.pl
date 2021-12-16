#!/usr/bin/env perl
#
# All-atom replica exchange simulations
#
# http://mmtsb.scripps.edu/doc/aarex.pl.html
# 2001, Michael Feig, Brooks group, TSRI
# 2001, John Karanicolas, Brooks group, TSRI
# PHMD added by Jana Khandogin, 2005

sub usage {
  printf STDERR "usage:   aarex.pl [options] [files]\n";
  printf STDERR "options: [-n runs]\n";
  printf STDERR "         [-par initruns=value,equilruns=value,\n";
  printf STDERR "               [no]save,savebestfreq=value,archive\n";
  printf STDERR "               ensmode=add|replace,natpdb=file,psf=file\n";
  printf STDERR "               arcmode=add|replace,[no]removecons]\n";
  printf STDERR "         [-temp nwin:min:max]\n";
  printf STDERR "         [-condfile file]\n";
  printf STDERR "         [-f listfile]\n";
  printf STDERR "         [-mdpar CHARMMparams]\n";
  printf STDERR "         [-mdopt [no]trajout,[no]restout,[no]conslim,\n";
  printf STDERR "                 limforce=value,limsel=ca|cb|cab|heavy]\n";
  printf STDERR "         [-l refPDB min:max[=min:max ...]]\n";
  printf STDERR "         [-cons [ca|cb|cab|heavy] ref|self min:max[_force][=...]]\n";
  printf STDERR "         [-opt optionsfile]\n";
  printf STDERR "         [-custom setup|pre|post[:init|equi|prod] file]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-ens tag] [-ensdir dir]\n";
  printf STDERR "         [PARALLELoptions]\n";
  printf STDERR "         [-mpirun n hosts mpirunexec]\n";
  printf STDERR "         [-ibrun n mpirunexec]\n";
  printf STDERR "         [-charmmexec charmmexec]\n";
  printf STDERR "         [-openmpi n openmpirunexec]\n";
  printf STDERR "         [-gpu n]\n";
  printf STDERR "         [-log file] [-elog file] [-charmmlog file]\n";
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
use Molecule;
use Ensemble;
use CHARMM;

my %par=();
my %defpar = (
	      initruns     => 0,
	      equilruns    => 0,
	      save         => 1,
	      savebestfreq => -1,
              removecons   => 0,
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

my $nmpi=0;
my $mpihosts=undef;
my $mpiexec=undef;

my $nib=0;

my $nopenmpi=0;
my $openmpiexec=undef;

my $gpu=undef;

my $charmmexec=undef;

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
		restout   => 1
);

my %mdpar=();
my %defmdpar = ( 
		param     =>  22,
		gb        =>  1,
		cutnb     =>  20.0,
		cutoff    =>  18.0,
		cuton     =>  16.0,
		dynsteps  =>  1000,
		dyneqfrq  =>  50,
		shake     =>  1
); 

my $serverlog;

my $dir=".";

my $optfile;

my $prunnodes;

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
  } elsif ($ARGV[0] eq "-charmmexec") {
    shift @ARGV;
    $charmmexec=shift @ARGV;
  } elsif ($ARGV[0] eq "-gpu") {
    shift @ARGV;
    $gpu=shift @ARGV;
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
  } elsif ($ARGV[0] eq "-mpirun") {
    shift @ARGV;
    $nmpi=shift @ARGV;
    $mpihosts=shift @ARGV;
    $mpiexec=shift @ARGV;
  } elsif ($ARGV[0] eq "-ibrun") {
    shift @ARGV;
    $nib=shift @ARGV;
    $mpiexec=shift @ARGV;
  } elsif ($ARGV[0] eq "-openmpi"){
    shift @ARGV;
    $nopenmpi=shift @ARGV;
    $openmpiexec=shift @ARGV;
  } elsif ($ARGV[0] eq "-charmmexec"){
    shift @ARGV;
    $charmmexec=shift @ARGV;
  } elsif ($ARGV[0] eq "-prunmulti") {
    shift @ARGV;
    $from=$to=$ENV{RMS_RANK}+1;
    $prunnodes=shift @ARGV;
    $cpus=1;
    $prun=1;
    $saveid="save.id" unless (defined $saveid);
  } elsif ($ARGV[0] eq "-mdpar") {
    shift @ARGV;
    &GenUtil::parsePar(\%mdpar,shift @ARGV);
  } elsif ($ARGV[0] eq "-mdopt") {
    shift @ARGV;
    &GenUtil::parsePar(\%mdopt,shift @ARGV);
  } elsif ($ARGV[0] eq "-elog") {
    shift @ARGV;
    $elogfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-charmmlog") {
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
      ($srec->{name},$srec->{port},$srec->{id})=split(/:/,$line);
      close RINP;
    }
  } while(!defined $srec->{name});
} else {
  if (defined $serverlog) {
    &GenUtil::makeDir("$dir");
    &GenUtil::setLogFile("$dir/$serverlog");
  }
}

my @mpih=();
if ($nmpi>0 && -r $mpihosts) {
  open INP,"$mpihosts";
  while (<INP>) {
    chomp;
    push(@mpih,$_);  
  }
  close INP;
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
          runs=>$nruns,natpdb=>$par{natpdb},psf=>$par{psf},
          removecons=>$par{removecons});
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

  my $server=ReXServer->new($rex->{par}->{runs},$dir,\@initfile,%par);

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

die "no jobs to run" 
  if ($#joblist<0);

if (defined $charmmexec) {
  $ENV{CHARMMEXEC}=$charmmexec;
  $CHARMM::exec=$charmmexec;
}

if (defined $gpu) {
  $ENV{OPENMM_DEVICE}=$gpu;
}

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
    $topt.=" -charmmlog $clogfile" if (defined $clogfile);
    $topt.=" -elog $elogfile" if (defined $elogfile);
    $topt.=" -dir $dir";      
    $topt.=" -openmpi $nopenmpi $openmpiexec" if (defined $nopenmpi && $nopenmpi>0);
    $topt.=" -charmmexec $charmmexec" if (defined $charmmexec);
    $topt.=" -gpu $hostlist->[$i]->{gpu}" if (defined $hostlist->[$i]->{gpu});

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

  if ($mp) {
    if (defined $rex->{par}->{natpdb} && !-r $rex->{par}->{natpdb}) {
      $mpClient->getFile($rex->{par}->{natpdb},"$rex->{dir}/nat.pdb");
      $rex->set(natpdb=>"$rex->{dir}/nat.pdb");
    }
    if (defined $rex->{par}->{psf} && !-r $rex->{par}->{psf}) {
      $mpClient->getFile($rex->{par}->{psf},"$rex->{dir}/struct.psf");
      $rex->set(psf=>"$rex->{dir}/struct.psf");
    }
  }

  for (my $i=0; $i<=$#joblist; $i++) {
    my $pid=fork();
    if (!$pid) {
      my $rexClient=ReXClient->new($joblist[$i],$srec);
      die "client cannot connect to server" if (!defined $rexClient);
      my $lasttemp=$rexClient->initialize($hostid,(defined $cpus)?$cpus:$#joblist+1);
      $rexClient->establishConnection();
      my $pnodes;
      my $pbase;
      if (defined $prun && $prun && defined $prunnodes) {
        $pnodes=$prunnodes;
        if ($joblist[$i] =~/aa([0-9]+)/) {
          my $tnode=int($1*$prunnodes/4+0.999)-1;
          $pbase=`/usr/local/bin/offset2base $tnode`;
          chomp $pbase;
        }
        printf STDERR "prun cpus: $pnodes, base: $pbase\n";
      }  
      my $thostfile="/tmp/$$-$joblist[$i]-hosts";
      my $mpicmd;
      if ($nmpi>0 && $#mpih>=0) {
         open OUT,">$thostfile";
         for (my $j=0; $j<$nmpi; $j++) {
           my $jinx=$j+$i*$nmpi;
           if ($jinx>$#mpih) {
             printf STDERR "not enough hosts in MPI host file\n";
           } else {
             printf OUT "%s\n",$mpih[$jinx];
           }
         }
         close OUT;
         $mpicmd=sprintf("%s -np %d -hostfile %s %s",$mpiexec,$nmpi,$thostfile,$ENV{CHARMMEXEC});
      } elsif ($nopenmpi>0) {
	$mpicmd=sprintf("%s -np %d %s",$openmpiexec,$nopenmpi,$charmmexec);
      } elsif ($nib>0) {
        my $offset=$i*$nib;
        $mpicmd=sprintf("ibrun -n %d -o %d %s",$nib,$offset,$mpiexec);
        printf STDERR "mpicmd: %s\n",$mpicmd; 
      }
#      printf STDERR "mpicmd: %s\n",$mpicmd; 
      &doJob($mp,$dir,$joblist[$i],$rex,$rexClient,
             $lasttemp,$clogfile,$elogfile,$pnodes/4,$pnodes,$pbase,$mpicmd);
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
  my $mp=shift;
  my $dir=shift;
  my $job=shift;
  my $rex=shift;
  my $rexClient=shift;
  my $lasttemp=shift;
  my $clogfile=shift;
  my $elogfile=shift;
  my $prunnodes=shift;
  my $pruncpus=shift;
  my $prunbase=shift;
  my $mpicmd=shift;
  my $gpu=shift;

  my $dowham=1;
  my $whamopen=0;

  my $consbias=0;

  my $last={};
  $last->{temp}=$lasttemp;
  $last->{ener}="N/A";
  $last->{rmsd}="N/A";
  $last->{boxsize}="N/A";
  $last->{volume}="N/A";
  $last->{val}=();
  $last->{valx}=(); #hmcm
  $last->{valy}=(); #hmcm
  $last->{valz}=(); #hmcm
  
  &GenUtil::makeDir("$dir/$job");

  my $initfile=$rexClient->initFile();

  if ($initfile eq "") {
    $rexClient->terminateServer();
    die "need initial structure for all-atom replica exchange";
  }

  if ($mp) {
    if (!&GenUtil::checkFile($initfile)) {
      if (defined $rex->{par}->{psf}) {
	$rexClient->getFile($initfile,"$dir/$job/init.crd");
	$initfile="$dir/$job/init.crd";
      } else {
	$rexClient->getFile($initfile,"$dir/$job/init.pdb");
	$initfile="$dir/$job/init.pdb";
      }
    }
    $rexClient->getFile("$dir/$job/$ReXServer::restartfile") 
  }

  my $bias=$rexClient->biasInfo();

  my $allcons=();

  if (defined $rex->{opt}->{cons}) {
    my @cstr=split(/,/,$rex->{opt}->{cons});
    my $c={ sel=>$cstr[0],
	    type=>$cstr[1],
	    list=>&GenUtil::fragListFromOption($cstr[2]),
	    reffile=>$cstr[3],
	    exclmode=>$cstr[4] };
    push(@{$allcons},$c);
  } 

  if (defined $rex->getFragList() && $rex->{opt}->{conslim}) {
    my $c={ sel=>(defined $rex->{opt}->{limsel})?$rex->{opt}->{limsel}:"heavy" ,
	    force=>(defined $rex->{opt}->{limforce})?$rex->{opt}->{limforce}:1.0,
	    type=>"ref",
	    reffile=>$rex->{par}->{fragref},
	    list=>$rex->getFragList(),
	    exclmode=>1 
          };

    push(@{$allcons},$c);
  }

  if ($mp) {
    my $nref=1;
    foreach my $ac ( @{$allcons} ) {
      if (defined $ac->{reffile} && !-r $ac->{reffile}) {
	$rexClient->getFile($ac->{reffile},"$dir/$job/ref$nref");
	$ac->{reffile}="$dir/$job/ref$nref";
	$nref++;
      }
    }
  }

  my $charmm=&CHARMM::new((defined $clogfile)?"$dir/$job/$clogfile":undef,undef,$prunnodes,$pruncpus,$prunbase,$mpicmd);
  $charmm->setEnergyLogFile("$dir/$job/$elogfile") 
    if (defined $elogfile);
  
  $charmm->loadParameters(%{$rex->getPar()});

  if (defined $rex->{par}->{psf}) {
    $charmm->setupFromPSF($rex->{par}->{psf},$initfile);
  } else {
    $charmm->setupFromPDB($initfile);
  }

  if ($charmm->{par}->{qfep}) {
    $charmm->setupPert($charmm->{par}->{qfepsel},(defined $charmm->{par}->{qfepfix} && $charmm->{par}->{qfepfix}==1)?$charmm->{par}->{qfepsel}:undef);
  }
## pHMD
  if (defined $charmm->{par}->{phmdpar}) {
      if( -r "$dir/$job/final.lamb" ) {
	  system("cp $dir/$job/final.lamb $dir/$job/final.lamb.last");
      }
      $charmm->writePHMD("$dir/$job/final.lamb");
  }
###

  $charmm->setupEnergy();
  
  if ($#{$allcons}>=0) {
    $charmm->setupRestraints(1.0,$allcons);
  }

  $charmm->noeRestraints();

  $charmm->shake();

  $charmm->stream($custom->{setup}) 
    if (defined $custom && exists $custom->{setup});
  
  my $outmol=Molecule::new();

  my $flist=();
  do {
    my $par={};
    $par->{ener}=$last->{ener};
    $par->{rmsd}=sprintf("%1.5f",$last->{rmsd}) if ($last->{rmsd} ne "N/A");
    $par->{boxsize}=$last->{boxsize};
    $par->{volume}=$last->{volume};

    if ($#{$last->{val}}>=0) {
      $par->{val}=join(":",@{$last->{val}});
    }
    if ($#{$last->{valx}}>=0) { #hmcm
      $par->{valx}=join(":",@{$last->{valx}});
    }
    if ($#{$last->{valy}}>=0) { #hmcm
      $par->{valy}=join(":",@{$last->{valy}});
    }
    if ($#{$last->{valz}}>=0) { #hmcm
      $par->{valz}=join(":",@{$last->{valz}});
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
    
    my $charmmtrajfile=($rex->{opt}->{trajout})?"$dir/$job/$ReXServer::trajout":undef;
    my $charmmrestfile=undef;
    my $charmminprest=undef;

    my $simulationtemp=$ret->{temp};
    
    if ($rex->{opt}->{restout}) {
      $charmmrestfile="$dir/$job/new.restart";
      if (-r "$dir/$job/$ReXServer::restartfile" && 
	  ((defined $last->{temp} && $simulationtemp == $last->{temp}))) {
	$charmminprest="$dir/$job/$ReXServer::restartfile";
      }
    }
    
    my $lambdacmd;
    if ($ret->{mode} ne "INIT") {
      for (my $i=0; $i<=$#{$bias}; $i++) {
	if ($bias->[$i]->{type} eq "lambda") {
	  $charmm->_sendCommand("LAMB $ret->{bias}->[$i]->{target}");
	  $charmm->_sendCommand("WLAM 89");
	  $lambdacmd="";
#	  $lambdacmd="pstart 1 pstop $charmm->{par}->{dynsteps} lstart 0.0 lstop 1.0 lambda $ret->{bias}->[$i]->{target} pwind";
	  if ($dowham) {
#	    $lambdacmd.=" wham 89";
	    $charmm->_sendCommand("open write card unit 89 name $dir/$job/wham.out");
	    $whamopen=1;
	  }
	} else {
	  $charmm->applyBias(%{$bias->[$i]},%{$ret->{bias}->[$i]});
	}
      }
    }

    $charmm->setParameter(dyntrfrq=>$charmm->{par}->{dynsteps})
      if (($charmm->{explicitWater}) && ($ret->{mode} eq "PROD"));

#    my $dextra=$charmm->{par}->{dynextra};
#    if (defined $dextra && defined $lambdacmd) {
#      $dextra.=" ".$lambdacmd;
#    } elsif (defined $lambdacmd) {
#      $dextra=$lambdacmd;
#    }

## pHMD
    if (defined $charmm->{par}->{phmdpar}) {
       $charmm->appendPHMD("$dir/$job/final.lamb");
    }

##
    $charmm->runDynamics($charmminprest,$charmmrestfile,$charmmtrajfile,undef,
			 dyntemp=>$simulationtemp, dyntwin=>0.02*$simulationtemp,
                         dynextra=>$lambdacmd);

    $last->{boxsize}=$charmm->{mdener}->[$#{$charmm->{mdener}}]->{boxsize};
    $last->{volume}=$charmm->{mdener}->[$#{$charmm->{mdener}}]->{volume};
    $last->{pressure}=$charmm->{par}->{dynpress};

    my $tag=(lc $ret->{mode})."/".&GenUtil::dataDir($ret->{run});
    $tag .= ":".$ret->{temp};
    for (my $i=0; $i<=$#{$bias}; $i++) {
      $tag .= ":".$ret->{bias}->[$i]->{target};
    }
    $charmm->logMDEnergy($tag);

    my $lambdabias=0;
    for (my $i=0; $i<=$#{$bias}; $i++) {
      $charmm->clearBias(%{$bias->[$i]})
	if ($ret->{mode} ne "INIT");
      if ($bias->[$i]->{type} eq "cons") {
	$last->{val}->[$i]=$ret->{bias}->[$i]->{force};
	$consbias=1;
      } elsif ($bias->[$i]->{type} eq "lambda") {
#	$charmm->_sendCommand("energy lstart 0.0 lstop 1.0 lambda 0.0");
	$charmm->_sendCommand("WLAM -1");
	$charmm->_sendCommand("LAMB 0.0");
	$charmm->_sendCommand("ENERGY");
	my $te=$charmm->_processEneOutput();
	$last->{ener}=$te->[0]->{total};
#	$charmm->_sendCommand("energy lstart 0.0 lstop 1.0 lambda 1.0");
	$charmm->_sendCommand("WLAM -1");
	$charmm->_sendCommand("LAMB 1.0");
	$charmm->_sendCommand("ENERGY");
	$te=$charmm->_processEneOutput();
	$last->{val}->[$i]=$te->[0]->{total};
	$lambdabias=1;
      } elsif ($bias->[$i]->{type} eq "hmcm"){
	($last->{valx}->[$i],$last->{valy}->[$i],$last->{valz}->[$i],$last->{val}->[$i])=($charmm->getBiasVal(%{$bias->[$i]},%{$ret->{bias}->[$i]}));
      } else {
	$last->{val}->[$i]=$charmm->getBiasVal(%{$bias->[$i]},%{$ret->{bias}->[$i]})+0.0;
      }
    }

    if ($whamopen) {
      $charmm->_sendCommand("close unit 89");
      $whamopen=0;
    }

    $charmm->stream($custom->{post}) 
      if (defined $custom && exists $custom->{pre});
    $charmm->stream($custom->{"post:init"}) 
      if (defined $custom && exists $custom->{"post:init"} && $ret->{mode} eq "INIT");
    $charmm->stream($custom->{"post:equi"}) 
      if (defined $custom && exists $custom->{"post:equi"} && $ret->{mode} eq "EQUI");
    $charmm->stream($custom->{"post:prod"}) 
      if (defined $custom && exists $custom->{"post:prod"} && $ret->{mode} eq "PROD");

    if (!$lambdabias) {
      if (defined $rex->{opt}->{average} && $rex->{opt}->{average}) {
	my $nn=0;
	my $sum=0;
	for (my $i=$#{$charmm->{mdener}}/2; $i<=$#{$charmm->{mdener}}; $i++) {
	  $nn++;
	  my $val=$charmm->{mdener}->[$i]->{pot};
	  if ($consbias) {
	    $val-=$charmm->{mdener}->[$i]->{constr};
	  }
	  $sum+=$val;
	}
	$last->{ener}=$sum/$nn;
      } else {
	my $ener=$charmm->getEnergy();
	
	$last->{ener}=$ener->{total}+0.0;  
	if ($consbias) {
	  $last->{ener}-=$ener->{constr};
	} elsif (defined $rex->{par}->{removecons} && $rex->{par}->{removecons}) {
          $last->{ener}-=$ener->{constr};
          $last->{ener}-=$ener->{resd};
          $last->{ener}-=$ener->{cdih};
          $last->{ener}-=$ener->{geo};
        }
      }
    }

    $last->{temp}=$simulationtemp;

    my $chmoutpdb="$dir/$job/chmout.pdb";
    $charmm->writePDB($chmoutpdb);
    $outmol->readPDB($chmoutpdb,translate=>$rex->getPar()->{param}==19?"CHARMM19":"",
		     chainfromseg=>(defined $rex->{par}->{psf})?0:1);
    $outmol->setSSBonds($charmm->{molecule}->getSSBonds());
    $outmol->writePDB("$dir/$job/$ReXServer::outpdb",translate=>"CHARMM22");
    &GenUtil::remove($chmoutpdb);

    $charmm->writeCRD("$dir/$job/$ReXServer::outcrd")
      if (defined $rex->{par}->{psf});

    system "mv $dir/$job/new.restart $dir/$job/$ReXServer::restartfile"
      if (-r "$dir/$job/new.restart");
## pHMD
    $charmm->closePHMD()
	if (defined $charmm->{par}->{phmdpar});
##
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

	if (defined $rex->{par}->{psf}) {
	  my $frec={};
	  $frec->{local}="$dir/$job/$ReXServer::outcrd";
	  $frec->{remote}=$frec->{local};
	  push(@{$flist},$frec);
	}
## pHMD
	if (defined $charmm->{par}->{phmdpar}) {
	  my $frec={};
	  $frec->{local}="$dir/$job/$ReXServer::outphmd";
	  $frec->{remote}=$frec->{local};
	  push(@{$flist},$frec);
	}
##
      }
    }

    if (defined $rex->getNatPDB()) {
      if (!defined $rex->{obj}->{analyze}) {
	$rex->{obj}->{analyze}=Analyze::new($rex->getNatPDB());
      }
      
      my $exclmode=(defined $rex->{par}->{exclmode})?$rex->{par}->{exclmode}:1;
      $outmol->setValidResidues($rex->getFragList(),$exclmode)
	if (defined $rex->getFragList());
      $rex->{obj}->{analyze}->lsqfit($outmol,"cab",0);
      $outmol->setValidResidues($rex->getFragList(),0)
	if (defined $rex->getFragList() && $exclmode);
      my $rmsd=$rex->{obj}->{analyze}->rmsd($outmol,0,($exclmode)?0:1);
      $last->{rmsd}=$rmsd->{CA};
    }
  } while(1);
}
