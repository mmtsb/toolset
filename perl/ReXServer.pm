# ReXServer package
# provide server for replica exchange sampling
#
# http://mmtsb.scripps.edu/doc/ReXServer.pm.html
# 2001, Michael Feig, Brooks group, TSRI
# 2001, John Karanicolas, Brooks group, TSRI
#
# derived from SimData, Server

package ReXServer;

require 5.004;

use strict;

use Fcntl;
use IO::File;
use IO::Handle;
use IO::Socket;
use Net::hostent;
use Sys::Hostname;

use GenUtil;
use Server;
use Molecule;
use Ensemble;
use Analyze;
use SimData;

## data: par -> { initruns equilruns save render
## data:          savebestfreq savedatafreq saverestart ensmode }
## exchange simulation parameters

## data: runsleft
## number of exchange steps left to do

## data: runscompleted
## number of exchange steps completed

## bias: biasdata[] -> { type sel etc. }
## data for each umbrella

## data: cond[] -> { temp bias[] -> { force target ref etc. } 
## data:             inx clientid nlist nindex }
## temperature and biasing potential parameters for each replica

## data: initstruct[]
## initial structures

## data: trun
## current run index

## data: status -> { clientwaiting accepted rejected }
## server status information

## data: clientid[]
## list of client ids

## data: clientdata[] -> { id } 
## data:    { ener rmsd volume exchanged cond val[] }
## client data for each replica exchange cycle

## data: min -> { id } -> { temp [target[]] }
## overall minimum window visited

## data: max -> { id } -> { temp [target[]] }
## overall maximum window visited

## data: host -> { hostid } -> { nrun, maxrun, wait[] }
## host information for scheduling client jobs

## data: ens
## ensemble object for saving conformations

## data: analyze 
## Analyze object for calculating RMSD values
## during the simulation

# use vars qw ( @ISA $kB $outpdb $restartfile $trajout $outcrd $defaultCondFile );
# Added for pHMD
use vars qw ( @ISA $kB $outpdb $outpdb2 $outpdb3 $restartfile $trajout $outcrd $outcrd2 $outcrd3 $defaultCondFile $outphmd );

@ISA = ( "Server" );

BEGIN {
  $kB=0.001987;
  $outpdb="final.pdb";
  $outpdb2="final2.pdb";
  $outpdb3="final3.pdb";
  $outcrd="final.crd";
  $outcrd2="final2.crd";
  $outcrd3="final3.crd";
  $restartfile="restart";
  $trajout="traj.crd";
  $defaultCondFile="rexserver.cond";
  $outphmd="final.lamb";
}

## constructor: new(runs,dir[,initfiles[,parameters])
## creates a new ReXServer object. The number of requested replica 
## exchange steps is required as argument, other parameters may be 
## given as additional arguments.

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $runs=shift;
  my $dir=shift;
  my $initfiles=shift;
  my $self=$class->SUPER::new($dir,"rexserver.cfg");

  $self->{runs}=$runs;

  $self->{clientdata}=();
  $self->{clientid}=();
  $self->{min}={};
  $self->{max}={};

  $self->{host}={};

  $self->{status}={};
  $self->{status}->{clientwaiting}=0;
  $self->{status}->{accepted}=0;
  $self->{status}->{rejected}=0;
  $self->{trun}=0;
  $self->{biasdata}=();
  $self->{cond}=();
  $self->{initstruct}=$initfiles;
  $self->{arcrun}=0;
  
  $self->set(@_);

  $self->{_boltzFactor}=$kB;
  $self->{_pressure}=undef;

  return $self;
}

sub setup {
  my $self=shift;
  my $condfile=shift;

  $self->readCondFile($condfile);
  $self->setupCond();
  $self->setupPar();
  $self->setupStruct();
}

sub setupPar {
  my $self=shift;

  my %default = ( initruns      => 10,
		  equilruns     => 40,
                  savebestfreq  => -1,
                  xchangefreq   => 1,
		  savedatafreq  => 5,
                  saverestart   => 0,
                  save          => 0,
                  archive       => 0,
                  render        => 0,
                  ensmode       => "replace",
		  arcmode       => "add" );

  foreach my $n ( keys %default ) {
    $self->{par}->{$n}=$default{$n} 
      if (!defined $self->{par}->{$n});
  }
}

sub setupStruct {
  my $self=shift;

  if (defined $self->{initstruct} && $#{$self->{initstruct}}>=0) {
    for (my $i=0; $i<$self->nWindows(); $i++) {
      if (!defined $self->{initstruct}->[$i]) {
	$self->{initstruct}->[$i]=$self->{initstruct}->[$i-1];
      }
    }
  }
}

sub writeCondFile {
  my $self=shift;
  my $condfile=shift;

  $condfile=$self->{dir}."/$defaultCondFile" if (!defined $condfile);

  my $out=&GenUtil::getOutputFile($condfile);
  for (my $i=0; $i<$self->nUmbrellas(); $i++) {
    my $b=$self->{biasdata}->[$i];
    my $opt=join(",",map "$_=$b->{$_}",grep($_ ne "type",keys %{$b}));
    if ($opt ne "") {
      print $out "bias $b->{type} $opt\n";
    } else {
      print $out "bias $b->{type}\n";
    }
  }
  for (my $i=0; $i<$self->nWindows(); $i++) {
    my $c=$self->{cond}->[$i];
    my @opt=();
    for (my $i=0; $i<$self->nUmbrellas(); $i++) {
      push(@opt,join(",",map "$_=$c->{bias}->[$i]->{$_}",keys %{$c->{bias}->[$i]}));
    }
    if ($#opt>=0) {
      print $out "cond $c->{temp} ".join(" ",@opt)."\n";
    } else {
      print $out "cond $c->{temp}\n";
    }
  }
  close $out;
}

sub readCondFile {
  my $self=shift;
  my $condfile=shift;

  $condfile=$self->{dir}."/$defaultCondFile" if (!defined $condfile);

  die "cannot read condition file" if (!&GenUtil::checkFile($condfile));

  $self->{biasdata}=();
  $self->{cond}=();

  my $inp=&GenUtil::getInputFile($condfile);
  while(<$inp>) {
    if (!/^\#/) {
      chomp;
      s/^[ \t]+//;
      if ($_ ne "") {
	my @f=split(/[ \t]+/);
	my $key=$f[0];
	if (lc($key) eq "bias") {
	  my $brec={};
	  %{$brec}=split(/[,=]/,$f[2])
	    if (defined $f[2]);
	  $brec->{type}=lc $f[1];
	  push(@{$self->{biasdata}},$brec);
	} else {
	  shift @f if (lc($key) eq "cond");
	  die "syntax error" if ($f[0]!~/^[0-9\.]+$/);
	  my $condrec={ temp=>shift @f, bias=>() };
	  foreach my $b ( @f ) {
	    my $crec={};
# ensure lower case
	    %{$crec}=split(/[,=]/,$b);
	    push(@{$condrec->{bias}},$crec);
	  }
	  push(@{$self->{cond}},$condrec);
	}
      }
    }
  }
  close $inp;
}

sub setupCond {
  my $self=shift;

  die "no conditions available"
    if ($#{$self->{cond}}<0);


  for (my $i=0; $i <= $#{$self->{cond}}; $i++) {
    my $sc=$self->{cond}->[$i];
    $sc->{inx}=$i;
    $sc->{clientid}=undef;
    $sc->{nlist}=();
    $sc->{nindex}=0;
  }
  $self->getNeighbours();
  $self->readNeighbourData();
}

## method: getNeighbours()
## generates a list of neighbouring conditions,
## to be used in exchanges.
## Neighbours are defined as replicas which are
## identical in all except one conditions, and
## no other replica has an intermediate value in
## this condition.

sub getNeighbours {
  my $self=shift;

  foreach my $c ( @{$self->{cond}} ) {
    $c->{nlist}=();

    my @l=();
    foreach my $cc ( @{$self->{cond}} ) {
      my $ok=1;
      for (my $j=0; $j<$self->nUmbrellas() && $ok; $j++) {
	$ok=0 
	  if ($cc->{bias}->[$j]->{target} ne $c->{bias}->[$j]->{target});
      }
      push(@l,$cc) if ($ok);
    }
    if ($#l>0) {
      my @sl=sort {$a->{temp} <=> $b->{temp} } @l;
      for (my $is=0; $is<=$#sl; $is++) {
	if ($sl[$is]->{inx} == $c->{inx} && $is>0) {
	  my $rec;
	  $rec->{nid}=$sl[$is-1]->{inx};
	  $rec->{acc}=0;
	  $rec->{tot}=0;
	  push(@{$c->{nlist}},$rec);
	}
	if ($sl[$is]->{inx} == $c->{inx} && $is<$#sl) {
	  my $rec;
	  $rec->{nid}=$sl[$is+1]->{inx};
	  $rec->{acc}=0;
	  $rec->{tot}=0;
	  push(@{$c->{nlist}},$rec);
	}
      }
    }

    for (my $i=0; $i<$self->nUmbrellas(); $i++) {
      my @l=();
      foreach my $cc ( @{$self->{cond}} ) {
	my $ok=($c->{temp} eq $cc->{temp});
	for (my $j=0; $j<$self->nUmbrellas() && $ok; $j++) {
	  $ok=0 
	    if ($j!=$i && $cc->{bias}->[$j]->{target} ne $c->{bias}->[$j]->{target});
	}
	push(@l,$cc) if ($ok);
      }
      if ($#l>0) {
	my @sl=sort {$a->{bias}->[$i]->{target} <=> 
		     $b->{bias}->[$i]->{target}} @l;
	for (my $is=0; $is<=$#sl; $is++) {
	  if ($sl[$is]->{inx} == $c->{inx} && $is>0) {
	    my $rec;
	    $rec->{nid}=$sl[$is-1]->{inx};
	    $rec->{acc}=0;
	    $rec->{tot}=0;
	    push(@{$c->{nlist}},$rec)
	  }

	  if ($sl[$is]->{inx} == $c->{inx} && $is<$#sl) {
	    my $rec;
	    $rec->{nid}=$sl[$is+1]->{inx};
	    $rec->{acc}=0;
	    $rec->{tot}=0;
	    push(@{$c->{nlist}},$rec)
	  }
	}
      }
    }
  }
}

## method: setEnsemble(ens)
## provides an ensemble object for saving conformations
## from the replica exchange simulation 

sub setEnsemble {
  my $self=shift;
  my $ens=shift;

  $self->{ens}=$ens;
  $self->{ens}->readFileList();
}

## method: num = nWindows()
## returns the number of windows from the temperature array

sub nWindows {
  my $self=shift;
  return $#{$self->{cond}}+1;
}

## method: num = nUmbrellas()
## returns the number of umbrellas being used (not counting temperature)

sub nUmbrellas {
  my $self=shift;
  return $#{$self->{biasdata}}+1;
}

sub getClientData {
  my $self=shift;
  my $clientid=shift;
  my $runnum=shift;

  $runnum=$self->{trun} if (!defined $runnum);

  return undef if ($runnum<=0);

  return $self->{clientdata}->[$runnum-1] 
    if (! defined $clientid);

  return $self->{clientdata}->[$runnum-1]->{$clientid}
}

sub getMode {
  my $self=shift;
  my $trun=shift;

  $trun=$self->{trun} if (!defined $trun);
  
  if ($trun<=$self->{par}->{initruns}) {
    return "INIT";
  } elsif ($trun<=$self->{par}->{initruns}+$self->{par}->{equilruns}) {
    return "EQUI";
  } else {
    return "PROD";
  }
}

sub getRun {
  my $self=shift;
  my $trun=shift;

  $trun=$self->{trun} if (!defined $trun);

  if ($trun<=$self->{par}->{initruns}) {
    return $trun;
  } elsif ($trun<=$self->{par}->{initruns}+$self->{par}->{equilruns}) {
    return $trun-$self->{par}->{initruns};
  } else {
    return $trun-$self->{par}->{initruns}-$self->{par}->{equilruns};
  }
}

sub getInitialStructure {
  my $self=shift;
  my $clientid=shift;

  if (defined $self->{initstruct}) {
    for (my $i=0; $i<=$#{$self->{clientid}}; $i++) {
      if ($clientid eq $self->{clientid}->[$i]) {
	return $self->{initstruct}->[$i]
	  if (defined $self->{initstruct}->[$i]);
      }
    }
  }

  return "$self->{dir}/$clientid/$outcrd"
    if (&GenUtil::checkFile("$self->{dir}/$clientid/$outcrd"));

  return "$self->{dir}/$clientid/$outpdb"
    if (&GenUtil::checkFile("$self->{dir}/$clientid/$outpdb"));

  return "$self->{dir}/$clientid/$outphmd"
    if (&GenUtil::checkFile("$self->{dir}/$clientid/$outphmd"));

  return undef;
}

sub getInitialCond {
  my $self=shift;
  my $clientid=shift;

  for (my $i=0; $i<=$#{$self->{clientid}}; $i++) {
    if ($clientid eq $self->{clientid}->[$i]) {
      return $self->{cond}->[$i];
    }
  }
  return undef;
}

## method: writeData([file])
## writes status information for last completed exchange
## cycle

sub writeData {
  my $self=shift;
  my $fname=shift;

  return if (!defined $self->{clientid} || $#{$self->{clientid}}<0);

  $fname=$self->{dir}."/rexserver.data" if (!defined $fname);

  my $handle=&GenUtil::getOutputFile($fname);

  print $handle "CLIENTID ".join(" ",@{$self->{clientid}})."\n";
  
  printf $handle "STATUS %d %d\n",
  $self->{status}->{accepted},$self->{status}->{rejected};

# append
# improve efficiency
  print $handle "MIN temp ",join(" ",map $self->{min}->{$_}->{temp},@{$self->{clientid}}),"\n";
  print $handle "MAX temp ",join(" ",map $self->{max}->{$_}->{temp},@{$self->{clientid}}),"\n";

  for (my $i=0; $i<$self->nUmbrellas(); $i++) {
    print $handle "MIN $i ",join(" ",map $self->{min}->{$_}->{target}->[$i],@{$self->{clientid}}),"\n";
    print $handle "MAX $i ",join(" ",map $self->{max}->{$_}->{target}->[$i],@{$self->{clientid}}),"\n";
  }

  for (my $i=1; $i<=$self->{trun}; $i++) {
    printf $handle "RUN %d %s %d\n",$i,$self->getMode($i),$self->getRun($i);
    my $cdat=$self->getClientData(undef,$i);
    print $handle "ENER ",join(" ",map $cdat->{$_}->{ener},@{$self->{clientid}}),"\n";
    print $handle "RMSD ",join(" ",map $cdat->{$_}->{rmsd},@{$self->{clientid}}),"\n";
    print $handle "VOLUME ",join(" ",map $cdat->{$_}->{volume},@{$self->{clientid}}),"\n";
    print $handle "COND ",join(" ",map $cdat->{$_}->{cond}->{inx},@{$self->{clientid}}),"\n";
    for (my $j=0; $j < $self->nUmbrellas(); $j++) {
      print $handle "VAL $j ",join(" ",map $cdat->{$_}->{val}->[$j],@{$self->{clientid}}),"\n";
    }
  }

  close $handle;
}

## method: readData([file])
## reads status information from last exchange cycle for
## restarting replica exchange simulations

sub readData {
  my $self=shift;
  my $fname=shift;
  
  $fname=$self->{dir}."/rexserver.data" if (!defined $fname);

  return if (!&GenUtil::checkFile($fname));

  my $handle=&GenUtil::getInputFile($fname);

  $self->{trun}=0;

  my $cdat=undef;

  while(<$handle>) {
    chomp;
    my @f=split(/ +/);
    my $key=shift @f;

    if ($key eq "CLIENTID") {
      @{$self->{clientid}}=@f;
    } elsif ($key eq "STATUS") {
      $self->{status}->{accepted}=$f[0];
      $self->{status}->{rejected}=$f[1];
    } elsif ($key eq "MIN" || $key eq "MAX") {
      my $subkey=lc shift @f;
      foreach my $c (@{$self->{clientid}}) {
	if ($subkey eq "temp") {
	  $self->{lc $key}->{$c}->{temp}=shift @f;
	} else {
	  $self->{lc $key}->{$c}->{target}->[$subkey]=shift @f;
	}
      }
    } elsif ($key eq "RUN") {
      $self->{clientdata}->[$self->{trun}]={};
      $self->{trun}++;
      $cdat=$self->getClientData();
      $self->{arcrun}=1 if (!$self->{arcrun});
    } elsif (($key eq "ENER") || ($key eq "RMSD") || ($key eq "VOLUME")) {
      foreach my $c (@{$self->{clientid}}) {
	$cdat->{$c}->{lc $key}=shift @f;
      }
    } elsif ($key eq "COND") {
      foreach my $c (@{$self->{clientid}}) {
	$cdat->{$c}->{cond}=$self->{cond}->[shift @f];
      }
    } elsif ($key eq "VAL") {
      my $inx=shift @f;
      foreach my $c (@{$self->{clientid}}) {
	$cdat->{$c}->{val}->[$inx]=shift @f;
      }
    }      
  }	  
  close $handle;
}


## method: writeNeighbourData([file])
## writes the neighbour indices for last completed exchange
## cycle and acceptance ratios for each neighbour pair

sub writeNeighbourData {
  my $self=shift;
  my $fname=shift;

  return if (!defined $self->{cond} || $#{$self->{cond}}<0);

  $fname=$self->{dir}."/rexserver.ninx" if (!defined $fname);

  my $handle=&GenUtil::getOutputFile($fname);
  print $handle join(" ",map $_->{nindex},@{$self->{cond}}),"\n";

  for (my $inx=0; $inx < $#{$self->{cond}}; $inx++) {
    my $ci=$self->{cond}->[$inx];
    for (my $jnx=$inx+1; $jnx <= $#{$self->{cond}}; $jnx++) {
      my $cj=$self->{cond}->[$jnx];
      my $acc=0;
      my $tot=0;
      foreach my $k (@{$self->{cond}->[$inx]->{nlist}}) {
        if ($k->{nid} == $jnx) {
	  $tot+=$k->{tot};
	  $acc+=$k->{acc};
        }
      }
      foreach my $k (@{$self->{cond}->[$jnx]->{nlist}}) {
        if ($k->{nid} == $inx) {
	  $tot+=$k->{tot};
	  $acc+=$k->{acc};
        }
      }

      printf $handle "%3d %3d %5d %5d %6.2f\n", $inx, $jnx, $acc, $tot, (100*$acc/$tot)
	if ($tot > 0);
    }
  }

  close $handle;
}

## method: readNeighbourData([file])
## reads status information from last exchange cycle for
## restarting replica exchange simulations

sub readNeighbourData {
  my $self=shift;
  my $fname=shift;
  
  return if (!defined $self->{cond});

  $fname=$self->{dir}."/rexserver.ninx" if (!defined $fname);

  return if (!&GenUtil::checkFile($fname));

  my $handle=&GenUtil::getInputFile($fname);
  my $line=<$handle>;
  chomp $line;
  my @f=split(/ +/,$line);
  for (my $i=0; $i<=$#f; $i++) {
    $self->{cond}->[$i]->{nindex}=$f[$i];
  }

  while (<$handle>) {
    chomp();
    my ($inx,$jnx,$acc,$tot,$junk)=split(' ');

    foreach my $k (@{$self->{cond}->[$inx]->{nlist}}) {
      if ($k->{nid} == $jnx) {
	$k->{tot}=$tot;
	$k->{acc}=$acc;
      }
    }
  }
  close $handle;

  return 1;
}

sub serverSetup {
  my $self=shift;

  $self->readData();

  $self->{runscompleted}=0;
  $self->{runsleft}=($self->{runs}>0)?$self->{runs}:
    ($self->{par}->{initruns}+$self->{par}->{equilruns}-$self->{runs}-$self->{trun});
}

sub maxClients {
  my $self=shift;
  return $self->nWindows();
}

sub renderImages {
  my $self=shift;
  foreach my $cid ( @{$self->{clientid}} ) {
    if (&GenUtil::checkFile("$self->{dir}/$cid/$outpdb")) {
      system "cd $self->{dir}/$cid; cp $outpdb fit.pdb;  molauto fit.pdb | molscript -s -size 110 110 -r | render -jpeg > final.jpeg 2> /dev/null";
    }
  }
}

sub serverTask {
  my $self=shift;
  my $ret=0;

  if ($self->{status}->{clientwaiting} >= $self->nWindows()) {
    $self->renderImages() if ($self->{par}->{render});
    $self->pushHTML();
    $ret=$self->exchange();
  } 
  return $ret;
}

sub serverFinish {
  my $self=shift;

  $self->writeData();
  $self->SUPER::serverFinish();
}

sub processCommand {
  my $self=shift;
  my $client=shift;
  my $cmd=shift;
  my $htmlfile=shift;
  my @arg=@_;

  if ($cmd eq "IFILE") {
    my $ifile=$self->getInitialStructure($arg[0]);
    &GenUtil::writeToSocket($client,"$ifile\n");
    &GenUtil::log("ReXServer::processCommand","IFILE $arg[0] -> $ifile");
  } elsif ($cmd eq "BIAS") {
    my @l=();
    for (my $i=0; $i<$self->nUmbrellas(); $i++) {
      push(@l,join(",",map "$_=$self->{biasdata}->[$i]->{$_}",keys %{$self->{biasdata}->[$i]}));
    }
    &GenUtil::writeToSocket($client,join(" ",@l)."\n");
    &GenUtil::log("ReXServer::processCommand","BIAS $arg[0] -> ".join(" ",@l));
  } else {
    return $self->SUPER::processCommand($client,$cmd,$htmlfile,@arg);
  }

  close $client;
  return 0;
}

sub processClientResponse {
  my $self=shift;
  my $clientid=shift;
  my $data=shift;

  my $inp;
  if (defined $data) {
    $data=~s/^ +//;
    chomp $data;
    %{$inp}=split(/ +/,$data);
  }

  $self->prepareNew($clientid);

  my $hostid=$self->{client}->{$clientid}->{host};
  die "host id undefined" if (!defined $hostid);

  my $host=$self->{host}->{$hostid};
  die "host undefined" if (!defined $host);
    
  $host->{nrun}-- if ($host->{nrun}>0);
  $self->giveClearance(pop(@{$host->{wait}}))
    if ($#{$host->{wait}}>=0);

  if ($self->{trun}>0) {
    my $cl = $self->getClientData($clientid);
    $cl->{ener}=$inp->{ener}
      if (defined $inp->{ener} && $inp->{ener} ne "N/A");

    $cl->{rmsd}=$inp->{rmsd}
      if (defined $inp->{rmsd} && $inp->{rmsd} ne "N/A");

    $cl->{volume}=$inp->{volume}
      if (defined $inp->{volume} && $inp->{volume} ne "N/A");

    $self->{_pressure}=$inp->{pressure}
      if (defined $inp->{pressure} && $inp->{pressure} ne "N/A");

    @{$cl->{val}}=split(/:/,$inp->{val})
      if (defined $inp->{val});

    @{$cl->{valx}}=split(/:/,$inp->{valx})
      if (defined $inp->{valx}); #hmcm

    @{$cl->{valy}}=split(/:/,$inp->{valy})
      if (defined $inp->{valy}); #hmcm

    @{$cl->{valz}}=split(/:/,$inp->{valz})
      if (defined $inp->{valz}); #hmcm

    $self->postProcess($clientid);
  }

  $self->{status}->{clientwaiting}++;
}

sub getInfo {
  my $self=shift;
  my $subcmd=shift;

  if ($subcmd eq "SUMMARY") {
    return $self->getSummary();
  } elsif ($subcmd eq "CLIENT") {
    return $self->getClientInfo();
  } elsif ($subcmd eq "HOST") {
    return $self->getHostInfo();
  } elsif ($subcmd eq "ID") {
    return $self->getIDInfo();
  }
}

sub initClient {
  my $self=shift;
  my $clientid=shift;
  my $hostid=shift;
  my $maxcpus=shift;
  
  if (!defined $self->{host}->{$hostid}) {
    $self->{host}->{$hostid}={};
    $self->{host}->{$hostid}->{id}=$hostid;
    $self->{host}->{$hostid}->{nrun}=0;
    $self->{host}->{$hostid}->{maxrun}=$maxcpus;
    $self->{host}->{$hostid}->{wait}=();
  }

  $self->{client}->{$clientid}->{host}=$hostid;

  if ($#{$self->{clientid}}+1==$self->nWindows()) {

    my $cdat=$self->getClientData($clientid);
    return "$cdat->{cond}->{temp}";
  } else {
    &GenUtil::makeDir("$self->{dir}/$clientid");
    push(@{$self->{clientid}},$clientid);
    return "";
  }
}

sub giveClearance {
  my $self=shift;
  my $clientid=shift;

  my $c=$self->getClientData($clientid);
  my $s=sprintf("%s %d %d %f",$self->getMode(),$self->getRun(),$self->{trun},$c->{cond}->{temp});
  for (my $i=0; $i<$self->nUmbrellas(); $i++) {
    $s.=" ".join(",",map "$_=$c->{cond}->{bias}->[$i]->{$_}",keys %{$c->{cond}->{bias}->[$i]});
  }

  $self->respondToClient($clientid,$s);
}

sub prepareNew {
  my $self=shift;
  my $clientid=shift;

  $self->{clientdata}->[$self->{trun}]={}
    if (!defined $self->{clientdata}->[$self->{trun}]);

  my $pdat=$self->getClientData($clientid);
  my $newrec={};

  if (defined $pdat->{cond}) {
    $newrec->{cond}=$pdat->{cond};
  } else {
    $newrec->{cond}=$self->getInitialCond($clientid);
  }

  $newrec->{cond}->{clientid}=$clientid;

  $newrec->{ener}="N/A";
  $newrec->{rmsd}="N/A";
  $newrec->{volume}="N/A";
  $newrec->{exchanged}=0;
  $newrec->{val}=();
  for (my $i=0; $i < $self->nUmbrellas(); $i++) {
    $newrec->{val}->[$i]=0.0;
  }

  $self->{clientdata}->[$self->{trun}]->{$clientid}=$newrec;

}

sub testLimits {
  my $self=shift;

  foreach my $cid ( @{$self->{clientid}} ) {
    my $c=$self->getClientData($cid);
    $self->{min}->{$cid}->{temp}=$c->{cond}->{temp}
    if (!defined $self->{min}->{$cid}->{temp} || 
	$c->{cond}->{temp}<$self->{min}->{$cid}->{temp});
    $self->{max}->{$cid}->{temp}=$c->{cond}->{temp}
    if (!defined $self->{max}->{$cid}->{temp} || 
	$c->{cond}->{temp}>$self->{max}->{$cid}->{temp});

    for (my $i=0; $i<$self->nUmbrellas(); $i++) {
      $self->{min}->{$cid}->{target}=() 
	if (!defined $self->{min}->{$cid}->{target});
      $self->{min}->{$cid}->{target}->[$i]=$c->{cond}->{bias}->[$i]->{target}
      if (!defined $self->{min}->{$cid}->{target}->[$i] || 
	  $c->{cond}->{bias}->[$i]->{target}<$self->{min}->{$cid}->{target}->[$i]);

      $self->{max}->{$cid}->{target}=() 
	if (!defined $self->{max}->{$cid}->{target});
      $self->{max}->{$cid}->{target}->[$i]=$c->{cond}->{bias}->[$i]->{target}
      if (!defined $self->{max}->{$cid}->{target}->[$i] || 
	  $c->{cond}->{bias}->[$i]->{target}>$self->{max}->{$cid}->{target}->[$i]);
    }
  }
}

sub exchange {
  my $self=shift;

  if ($self->{trun}>0) {

    $self->testLimits();

    if ($self->{par}->{savedatafreq}>0 && 
	($self->{trun} % $self->{par}->{savedatafreq}) == 0) {
      $self->writeData();
      $self->writeNeighbourData();
      $self->backupClientFiles();
    }

    $self->saveBest() 
      if ($self->getMode() eq "PROD" && 
	  $self->{par}->{savebestfreq}>0 && 
	  ($self->getRun()%$self->{par}->{savebestfreq})==0);
  }

  return 1 if (--$self->{runsleft}<0);

  $self->{runscompleted}++;

  $self->{trun}++;
  $self->{arcrun}++;

  $self->xchangeStep()
    if ($self->getMode() eq "PROD" && 
	(($self->getRun()%$self->{par}->{xchangefreq})==0));

  foreach my $cid ( @{$self->{clientid}} ) {
    my $hostid=$self->{client}->{$cid}->{host};
    die "host id undefined" if (!defined $hostid);

    my $host=$self->{host}->{$hostid};
    die "host undefined" if (!defined $host);

    if (!defined $host->{maxrun} || $host->{nrun}<$host->{maxrun}) {
      $self->giveClearance($cid);
      $host->{nrun}++;
    } else {
      push(@{$host->{wait}},$cid);
    }
  }

  $self->{status}->{clientwaiting}=0;

  return 0;
}

sub backupClientFiles {
  my $self=shift;

  foreach my $cid ( @{$self->{clientid}} ) {
    my $d="$self->{dir}/$cid";
    system "cp $d/$outpdb $d/backup.pdb"          if (-r "$d/$outpdb");
    system "cp $d/$restartfile $d/backup.restart" if (-r "$d/$restartfile");
  }
}

sub postProcess {
  my $self=shift;
  my $clientid=shift;
  my $d="$self->{dir}/$clientid";

  if ($self->{par}->{save} && $self->{trun}>0) {
    if ($self->{par}->{archive}) {
      if (-r "$d/$outpdb") {
	my $all;
	my $pdb=&GenUtil::getInputFile("$d/$outpdb");
	while (<$pdb>) {
	  if (/ATOM/) {
	    $all.=substr($_,30,24);
	  }
	}
	undef $pdb;
	if ($self->{par}->{arcmode} eq "replace"){
	    if ($self->{arcrun}==1){
		unlink("$d/".(lc $self->getMode()).".coor.archive");
	    }
	    my $saw=$self->getRun();
	    print "$saw\n";
	&GenUtil::writeArchiveFile("$d/".(lc $self->getMode()).".coor.archive",
				   $all,$self->{arcrun},1);
	}
	else{
	&GenUtil::writeArchiveFile("$d/".(lc $self->getMode()).".coor.archive",
				   $all,$self->getRun(),1);
    }
      }
      if (-r "$d/$restartfile" && $self->{par}->{saverestart}) {
	&GenUtil::archiveFile("$d/".(lc $self->getMode()).".restart.archive",
			      "$d/$restartfile",$self->getRun());
      }
      if (-r "$d/$trajout") {
	&GenUtil::archiveFile("$d/".(lc $self->getMode()).".traj.archive",
			      "$d/$trajout",$self->getRun());
      }
    } else {
      my $tag=(lc $self->getMode())."/".&GenUtil::dataDir($self->getRun());
      my $savedir="$d/$tag";
      &GenUtil::makeDir($savedir);
      if (-r "$d/$outpdb") {
	system "cp $d/$outpdb $savedir";
	&GenUtil::compress("$savedir/$outpdb");
      }
      if (-r "$d/$outphmd") {
	system "cp $d/$outphmd $savedir";
	&GenUtil::compress("$savedir/$outphmd");
      }
      if (-r "$d/wham.out") {
	system "cp $d/wham.out $savedir";
	&GenUtil::compress("$savedir/wham.out");
      }
      if (-r "$d/$restartfile" && $self->{par}->{saverestart}) {
	system "cp $d/$restartfile $savedir";
	&GenUtil::compress("$savedir/$restartfile");
      }

      if (-r "$d/$trajout") {
	system "mv $d/$trajout $savedir";
	&GenUtil::compress("$savedir/$trajout");
      }
    }
  }

  if (defined $self->getNatPDB() && 
      !defined $self->getClientData($clientid)->{rmsd}) {
    if (!defined $self->{obj}->{analyze}) {
      $self->{obj}->{analyze}=Analyze::new($self->getNatPDB());
    }

    if (-r "$d/$outpdb") {
      my $cmp=Molecule::new("$d/$outpdb");
      my $exclmode=(defined $self->{par}->{exclmode})?$self->{par}->{exclmode}:1;
      $cmp->setValidResidues($self->getFragList(),$exclmode)
	if (defined $self->getFragList());
      $self->{obj}->{analyze}->lsqfit($cmp,"cab",0);
      $cmp->setValidResidues($self->getFragList(),0)
	if (defined $self->getFragList() && $exclmode);
      my $rmsd=$self->{obj}->{analyze}->rmsd($cmp,0,($exclmode)?0:1);
      $self->getClientData($clientid)->{rmsd}=$rmsd->{CA};
    }
  }
}

sub save {
  my $self=shift;

  $self->writeCondFile();
  $self->SUPER::save();
}

sub saveBest {
  my $self=shift;
  my $cdat=$self->getClientData();

  $self->saveReplica($self->{cond}->[0]->{clientid});
}

sub saveReplica {
  my $self=shift;
  my $clientid=shift;

  if (defined $self->{ens}) {
    my $at=($self->{par}->{ensmode} eq "add")?$self->{ens}->{par}->{runs}+1:$self->getRun();
    my $fname=sprintf("%s/%s",$self->getRun(),$clientid);
    $self->{ens}->setFileList($at,$fname);
    $self->ensembleCheckin($clientid,$at);
    $self->{ens}->save();
  }
}

sub ensembleCheckin {
  my $self=shift;
  my $clientid=shift;
  my $at=shift;

  my $mol=Molecule::new("$self->{dir}/$clientid/$outpdb");
  $self->{ens}->checkinPDB($at,$mol,$self->getClientData($clientid)->{ener});
}

sub setPairs {
  my $self=shift;
  my $parr=shift;
  
  my $xp=();

  foreach my $p ( split(/::/,$parr) ) {
    my @f=split(/:/,$p);
    my $prec={};
    $prec->{i}=$self->{cond}->[$f[0]]->{clientid};
    $prec->{j}=$self->{cond}->[$f[1]]->{clientid};
    push(@{$xp},$prec);
  }

  return $xp;
}

sub xchangeStep {
  my $self=shift;

  &GenUtil::log("ReXServer::xchangeStep","start exchange");

  return if ($self->{trun}<=1);

  my $xpairs;

  if (defined $self->{par}->{evenpairs} && $self->{trun} % 2 == 0) {
    $xpairs=$self->setPairs($self->{par}->{evenpairs});
  } elsif (defined $self->{par}->{oddpairs} && $self->{trun} % 2 == 1) {
    $xpairs=$self->setPairs($self->{par}->{oddpairs});
  } else {
    $xpairs=$self->xchangePairs();
  }

  foreach my $p ( @{$xpairs} ) {

    my $pci=$self->getClientData($p->{i},$self->{trun}-1);
    my $pcj=$self->getClientData($p->{j},$self->{trun}-1);

    &GenUtil::log("ReXServer::xchangeStep",
       "testing conds $p->{i} $p->{j} $pci->{cond}->{inx} $pcj->{cond}->{inx}");

    my $bm=1.0/($self->{_boltzFactor}*$pci->{cond}->{temp});
    my $bn=1.0/($self->{_boltzFactor}*$pcj->{cond}->{temp});

    my $emi=$pci->{ener};
    my $eni=$emi;
    my $emj=$pcj->{ener};
    my $enj=$emj;
    
    for (my $k=0; $k < $self->nUmbrellas(); $k++) {
      if ($self->{biasdata}->[$k]->{type} eq "cons") {
	if (abs($pci->{cond}->{temp}-$pcj->{cond}->{temp})<0.1) {
	  if ($pci->{cond}->{bias}->[$k]->{target}<$pcj->{cond}->{bias}->[$k]->{target}) {
	    $enj=0.0;
	    $eni=0.0;
	  } else {
	    $emj=0.0;
	    $emi=0.0;
	  }
	}
      } elsif ($self->{biasdata}->[$k]->{type} eq "lambda") {
	&GenUtil::log("ReXServer::xchangeStep",
		      "lambda $k $pci->{cond}->{bias}->[$k]->{target}:$pci->{val}->[$k] $pcj->{cond}->{bias}->[$k]->{target}:$pcj->{val}->[$k]");

        my $tvi=$pci->{val}->[$k];
        my $tvj=$pcj->{val}->[$k];
        my @fvi=split(/=/,$tvi);
        my @fvj=split(/=/,$tvj);

        my $vi=$fvi[0];
        my $vj=$fvj[0];

        my $cli=($#fvi>0)?$fvi[1]:0;
        my $clj=($#fvj>0)?$fvj[1]:0;

        my $str=sprintf("%f %f | %f %f\n",$vi,$cli,$vj,$clj);
	&GenUtil::log("ReXServer::xchangeStep",$str);

	my $bi = $pci->{cond}->{bias}->[$k]->{target};
	my $bj = $pcj->{cond}->{bias}->[$k]->{target};

	$emi = (1.0-$bi)*$emi + $bi*$vi; # i at lambdai
        $emj = (1.0-$bi)*$emj + $bi*$vj; # j at lambdai
	$eni = (1.0-$bj)*$eni + $bj*$vi; # i at lambdaj
        $enj = (1.0-$bj)*$enj + $bj*$vj; # j at lambdaj

        $emi += $cli if ($bi>0.00000001 && $bi<0.999999999);
        $emj += $clj if ($bi>0.00000001 && $bi<0.999999999);
        $eni += $cli if ($bj>0.00000001 && $bj<0.999999999);
        $enj += $clj if ($bj>0.00000001 && $bj<0.999999999);
      } elsif ($self->{biasdata}->[$k]->{type} eq "hmcm" || $self->{biasdata}->[$k]->{type} eq "hmcr") {
	&GenUtil::log("ReXServer::xchangeStep",
		      "umbrella $k ($pci->{valx}->[$k],$pci->{valy}->[$k],$pci->{valz}->[$k]):($pci->{cond}->{bias}->[$k]->{refx},$pci->{cond}->{bias}->[$k]->{refy},$pci->{cond}->{bias}->[$k]->{refz}):$pci->{cond}->{bias}->[$k]->{target}:$pci->{cond}->{bias}->[$k]->{force}:$pci->{val}->[$k] ($pcj->{valx}->[$k],$pcj->{valy}->[$k],$pcj->{valz}->[$k]):($pcj->{cond}->{bias}->[$k]->{refx},$pcj->{cond}->{bias}->[$k]->{refy},$pcj->{cond}->{bias}->[$k]->{refz}):$pcj->{cond}->{bias}->[$k]->{target}:$pcj->{cond}->{bias}->[$k]->{force}:$pcj->{val}->[$k]");
	
	my $dxii= $pci->{valx}->[$k]-$pci->{cond}->{bias}->[$k]->{refx};
	my $dyii= $pci->{valy}->[$k]-$pci->{cond}->{bias}->[$k]->{refy};
	my $dzii= $pci->{valz}->[$k]-$pci->{cond}->{bias}->[$k]->{refz};

	my $dxjj= $pcj->{valx}->[$k]-$pcj->{cond}->{bias}->[$k]->{refx};
	my $dyjj= $pcj->{valy}->[$k]-$pcj->{cond}->{bias}->[$k]->{refy};
	my $dzjj= $pcj->{valz}->[$k]-$pcj->{cond}->{bias}->[$k]->{refz};

	my $dxij= $pci->{valx}->[$k]-$pcj->{cond}->{bias}->[$k]->{refx};
	my $dyij= $pci->{valy}->[$k]-$pcj->{cond}->{bias}->[$k]->{refy};
	my $dzij= $pci->{valz}->[$k]-$pcj->{cond}->{bias}->[$k]->{refz};

	my $dxji= $pcj->{valx}->[$k]-$pci->{cond}->{bias}->[$k]->{refx};
	my $dyji= $pcj->{valy}->[$k]-$pci->{cond}->{bias}->[$k]->{refy};
	my $dzji= $pcj->{valz}->[$k]-$pci->{cond}->{bias}->[$k]->{refz};

	my $bi = $pci->{cond}->{bias}->[$k]->{target};
	my $bj = $pcj->{cond}->{bias}->[$k]->{target};
	my $fi = $pci->{cond}->{bias}->[$k]->{force};
	my $fj = $pcj->{cond}->{bias}->[$k]->{force};
	my $vii = sqrt($dxii*$dxii+$dyii*$dyii+$dzii*$dzii);
	my $vjj = sqrt($dxjj*$dxjj+$dyjj*$dyjj+$dzjj*$dzjj);
	my $vij = sqrt($dxij*$dxij+$dyij*$dyij+$dzij*$dzij);
	my $vji = sqrt($dxji*$dxji+$dyji*$dyji+$dzji*$dzji);

	$emi += 0.5*$fi * ($vii - $bi) * ($vii - $bi);
	$emj += 0.5*$fi * ($vji - $bi) * ($vji - $bi);
	$eni += 0.5*$fj * ($vij - $bj) * ($vij - $bj);
	$enj += 0.5*$fj * ($vjj - $bj) * ($vjj - $bj);
      } elsif ($self->{biasdata}->[$k]->{type} eq "path") {
	&GenUtil::log("ReXServer::xchangeStep",
		      "umbrella $k $pci->{cond}->{bias}->[$k]->{target}:$pci->{cond}->{bias}->[$k]->{force}:$pci->{val}->[$k] $pcj->{cond}->{bias}->[$k]->{target}:$pcj->{cond}->{bias}->[$k]->{force}:$pcj->{val}->[$k]");
	my $bi = $pci->{cond}->{bias}->[$k]->{target};
	my $bj = $pcj->{cond}->{bias}->[$k]->{target};
	my $fi = $pci->{cond}->{bias}->[$k]->{force};
	my $fj = $pcj->{cond}->{bias}->[$k]->{force};
	my $vi = $pci->{val}->[$k];
	my $vj = $pcj->{val}->[$k];

	#Note that force constant is not halved!

	$emi += $fi * ($vi - $bi) * ($vi - $bi);
	$emj += $fi * ($vj - $bi) * ($vj - $bi);
	$eni += $fj * ($vi - $bj) * ($vi - $bj);
	$enj += $fj * ($vj - $bj) * ($vj - $bj);
      } elsif ($self->{biasdata}->[$k]->{type} eq "dihe") {
	&GenUtil::log("ReXServer::xchangeStep",
		      "umbrella $k $pci->{cond}->{bias}->[$k]->{target}:$pci->{cond}->{bias}->[$k]->{force}:$pci->{val}->[$k] $pcj->{cond}->{bias}->[$k]->{target}:$pcj->{cond}->{bias}->[$k]->{force}:$pcj->{val}->[$k]");
	my $bi = $pci->{cond}->{bias}->[$k]->{target};
	my $bj = $pcj->{cond}->{bias}->[$k]->{target};
	#Convert force constant from kcal/mol/rad/rad to kcal/mol/deg/deg
	my $pi=4*atan2(1,1);
	my $convf=($pi/180)*($pi/180);
	my $fi = $convf*$pci->{cond}->{bias}->[$k]->{force};
	my $fj = $convf*$pcj->{cond}->{bias}->[$k]->{force};
	my $vi = $pci->{val}->[$k];
	my $vj = $pcj->{val}->[$k];

  my $dvb=abs($vi - $bi);
  $dvb -= 360 if ($dvb>180);
  #$emi += 0.5*$fi * ($vi - $bi) * ($vi - $bi);
  $emi += 0.5*$fi * ($dvb) * ($dvb);
  
  $dvb=abs($vj - $bi);
  $dvb -= 360 if ($dvb>180);
  #$emj += 0.5*$fi * ($vj - $bi) * ($vj - $bi);
  $emj += 0.5*$fi * ($dvb) * ($dvb);
  
  $dvb=abs($vi - $bj);
  $dvb -= 360 if ($dvb>180);
  #$eni += 0.5*$fj * ($vi - $bj) * ($vi - $bj);
  $eni += 0.5*$fj * ($dvb) * ($dvb);
  
  $dvb=abs($vi - $bj);
  $dvb -= 360 if ($dvb>180);
  #$enj += 0.5*$fj * ($vj - $bj) * ($vj - $bj);
  $enj += 0.5*$fj * ($dvb) * ($dvb);
  
      } else {
	&GenUtil::log("ReXServer::xchangeStep",
		      "umbrella $k $pci->{cond}->{bias}->[$k]->{target}:$pci->{cond}->{bias}->[$k]->{force}:$pci->{val}->[$k] $pcj->{cond}->{bias}->[$k]->{target}:$pcj->{cond}->{bias}->[$k]->{force}:$pcj->{val}->[$k]");
	
	my $bi = $pci->{cond}->{bias}->[$k]->{target};
	my $bj = $pcj->{cond}->{bias}->[$k]->{target};
	my $fi = $pci->{cond}->{bias}->[$k]->{force};
	my $fj = $pcj->{cond}->{bias}->[$k]->{force};
	my $vi = $pci->{val}->[$k];
	my $vj = $pcj->{val}->[$k];

	$emi += 0.5*$fi * ($vi - $bi) * ($vi - $bi);
	$emj += 0.5*$fi * ($vj - $bi) * ($vj - $bi);
	$eni += 0.5*$fj * ($vi - $bj) * ($vi - $bj);
	$enj += 0.5*$fj * ($vj - $bj) * ($vj - $bj);
      }
    }

    my $delta=$bm*($emj-$emi) - $bn*($enj-$eni);

    $delta+=$self->{_pressure}*1.4583963E-5*($bm-$bn)*
              ($pcj->{volume}-$pci->{volume})
      if ((defined $self->{_pressure}) && (defined $pci->{volume}) &&
          (defined $pcj->{volume}));

    &GenUtil::log("ReXServer::xchangeStep","delta: $delta, emj/emi: $emj/$emi, enj/eni: $enj/$eni, bm: $bm, bn: $bn");

    my $ci=$self->getClientData($p->{i});
    my $cj=$self->getClientData($p->{j});

    my $nptr;
    foreach my $n ( @{$ci->{cond}->{nlist}} ) {
      $nptr=$n
	if ($n->{nid} == $cj->{cond}->{inx});
    }
    $nptr->{tot}=0 if (! defined $nptr->{tot});
    $nptr->{acc}=0 if (! defined $nptr->{acc});
    $nptr->{tot}++;

    if ($delta<=0.0 || exp(-$delta)>rand) {
      &GenUtil::log("ReXServer::xchangeStep","exchange accepted");
      $self->{status}->{accepted}++;
      $nptr->{acc}++;

      my $t;
      $ci->{cond}=$pcj->{cond};
      $cj->{cond}=$pci->{cond};

      $ci->{cond}->{clientid}=$p->{i};
      $cj->{cond}->{clientid}=$p->{j};

      $ci->{exchanged}=1;
      $cj->{exchanged}=1;

    } else {
      &GenUtil::log("ReXServer::xchangeStep","exchange rejected");
      $self->{status}->{rejected}++;
    }
  }
}

sub xchangePairs {
  my $self=shift;

  my $xp=();

  my %trnd;

  my @unused=();
  for (my $i=0; $i <= $#{$self->{cond}}; $i++) {
    $unused[$i]=1;
    $trnd{sprintf("%1.4f",rand())}=$i 
      if ($#{$self->{cond}->[$i]->{nlist}}>=1);
  }

#  my $inx=int(rand($#unused+1));
#  while ($#{$self->{cond}->[$inx]->{nlist}} < 1) {
#    $inx=int(rand($#unused+1));
#  }

#  my $inx=0;
#  my $startinx=$inx;

#  my $outflag=1;
#  while ($outflag) {

  foreach my $ki ( sort keys %trnd ) {
    my $inx=$trnd{$ki};
    if ($unused[$inx]) {
      if (defined $self->{cond}->[$inx]->{nlist}) {
        my $Ninx=$self->{cond}->[$inx]->{nindex};
        my $startNinx=$Ninx;

        my $inflag=1;
        while ($inflag) {
	  my $jnx=$self->{cond}->[$inx]->{nlist}->[$Ninx]->{nid};
	  if ($unused[$jnx]) {
	    my $prec={};
	    $prec->{i}=$self->{cond}->[$inx]->{clientid};
	    $prec->{j}=$self->{cond}->[$jnx]->{clientid};
	    push(@{$xp},$prec);
            $unused[$jnx]=0;
	    $inflag=0;
          } else {
	  }
	  $Ninx++;
	  $Ninx=0 if ($Ninx == $#{$self->{cond}->[$inx]->{nlist}}+1);
	  $inflag=0 if ($Ninx == $startNinx);
        }
        $self->{cond}->[$inx]->{nindex}=$Ninx;
        $unused[$inx]=0;
      }
    }
#    $inx++;
#    $inx=0 if ($inx == $#unused+1);
#    $outflag = 0 if ($inx == $startinx);
  }

  return $xp;
}

sub getRMSD {
  my $self=shift;
  my $reffile=shift;
  my $cmpfile=shift;

  my $exclmode=(defined $self->{par}->{exclmode})?$self->{par}->{exclmode}:1;
  my $refmol=Molecule::new($reffile);
  my $ana=Analyze::new($refmol);
  my $cmpmol=Molecule::new($cmpfile);
  
  $cmpmol->setValidResidues($self->getFragList(),$exclmode)
    if (defined $self->getFragList());
  $ana->lsqfit($cmpmol,"cab",0);
  $cmpmol->setValidResidues($self->getFragList(),0)
    if (defined $self->getFragList() && $exclmode);

  return $ana->rmsd($cmpmol,0,($exclmode)?0:1);
}

sub getSummary {
  my $self=shift;

  return sprintf("%d %d %d %s %d %d %d\n",
  $self->nWindows(),$self->{status}->{clientwaiting},
  $self->{trun},$self->getMode(),$self->getRun(),
  $self->{status}->{accepted},$self->{status}->{rejected});
}
  
sub getClientInfo {
  my $self=shift;

  my $str="";
  if ($self->{trun}>1) {
    foreach my $cid ( sort @{$self->{clientid}} ) {
      my $dat=$self->getClientData($cid,$self->{trun}-1);
      $str.=
	sprintf("%s %s %f %f %f %f %f",
		$cid,$self->{client}->{$cid}->{host},
		$dat->{ener},$dat->{rmsd},$dat->{volume},$dat->{cond}->{temp},
	        $self->{min}->{$cid}->{temp},$self->{max}->{$cid}->{temp});
      for (my $i=0; $i<$self->nUmbrellas(); $i++) {
	$str.=sprintf(" %s:%f:%f %f %f %f",
		      $self->{biasdata}->[$i]->{type},
		      $dat->{bias}->[$i]->{force},$dat->{bias}->[$i]->{target},
		      $dat->{val}->[$i],
		      $self->{min}->{$cid}->{target}->[$i],
		      $self->{max}->{$cid}->{target}->[$i]);
      }      
      $str .= "\n";
    }
    return $str;
  } else {
    return "\n";
  }
}

sub getIDInfo {
  my $self=shift;

  return join(" ",@{$self->{clientid}})."\n";
}

sub getHostInfo {
  my $self=shift;

  return "" if (!defined $self->{host});

  my $buffer="";
  foreach my $h ( sort keys (%{$self->{host}})) {
    my $hrec=$self->{host}->{$h};
    $buffer.=sprintf("%s %d %d %d\n",
    $h,$hrec->{nrun},$hrec->{maxrun},$#{$hrec->{wait}}+1);
  }
  return $buffer;
}

sub getHTMLPage {
  my $self=shift;
  my $request=shift;

  my $run=$self->{trun}-1;

  my $page=();

  push(@{$page},"<HTML>");
  push(@{$page},"<META HTTP-EQUIV=\"Refresh\" CONTENT=\"30\">") 
    if ($request=~/pull/);
#  if ($request=~/viz/) {
#  push(@{$page},"<HEAD><script language=JavaScript>");
#  foreach my $cid ( @{$self->{clientid}} ) {
#    push(@{$page},"image$cid = new Image(); image$cid = \"$self->{dir}/$cid/final.jpeg?$self->{serverID}\";");
#  }
#  push(@{$page},"</script></head>");
#  } else {
  push(@{$page},"<HEAD></HEAD>");
#  }
  if ($request=~/viz/) {
    push(@{$page},"<BODY BGCOLOR=#404040>");
  } else {
    push(@{$page},"<BODY BGCOLOR=#101010>");
  }
  push(@{$page},"<CENTER>");
  push(@{$page},"<TABLE>");
  push(@{$page},"<TR><TD COLSPAN=2>");
  push(@{$page},"<FONT FACE=sans-serif,ARIAL,HELVETICA COLOR=#5050FF SIZE=3><B>MMTSB Replica Exchange Monitor</B></FONT>");
  push(@{$page},"</TD></TR>");

  if (defined $self->{clientid} && $#{$self->{clientid}}>=0 && $run>=1) {
    push(@{$page},"<TR><TD>&nbsp;</TD></TR>");
    push(@{$page},"<TR><TD ALIGN=left>");
    push(@{$page},
	 sprintf("<FONT FACE=sans-serif,ARIAL,HELVETICA COLOR=white SIZE=3><B>%s %d</B></FONT>",
		 $self->getMode($run),$self->getRun($run)));
    push(@{$page},"</TD>");

    my @tlist;
    foreach my $cid ( @{$self->{clientid}} ) {
      my $cdat=$self->getClientData($cid,$run);
      my $trec={};
      $trec->{cid}=$cid;
      $trec->{cdat}=$cdat;
      push(@tlist,$trec);
    }

    if ($request=~/viz/) {
      push(@{$page},"</TR>");
      push(@{$page},"<TR><TD><TABLE CELLSPACING=3 CELLPADDING=3>");
      push(@{$page},"<TR>");
      my $n=0;
      foreach my $crec ( sort { $a->{cdat}->{cond}->{temp} <=> 
                                $b->{cdat}->{cond}->{temp} } @tlist ) {

#      foreach my $crec ( @tlist ) {
        my $cid=$crec->{cid};
        my $cdat=$crec->{cdat};
        my $info=sprintf("%s%s</FONT> %s%1.1f</FONT> %s%1.2f</FONT>",
			  &GenUtil::HTMLFont(2,"white"),$cid,
                          &GenUtil::HTMLFont(2,"#5050FF"),$cdat->{cond}->{temp},
			  &GenUtil::HTMLFont(2,"#FF5050"),$cdat->{rmsd});
      
        push(@{$page},"<TD>$info<BR><IMG SRC=\"$cid/final.jpeg?$self->{serverID}\"></TD>");
        if (++$n%6==0 && $n<=$#{$self->{clientid}}) {
         push(@{$page},"</TR><TR>");
        }
      } 
      push(@{$page},"</TR>");
    } else {
    push(@{$page},"<TD ALIGN=right>");
    if ($self->{status}->{rejected}+$self->{status}->{accepted}>0) {
      push(@{$page},
	   sprintf("<FONT FACE=COURIER COLOR=#FF50FF SIZE=2>%1d/%1d %1.2f</FONT>",
		   $self->{status}->{accepted},$self->{status}->{rejected},
		   $self->{status}->{accepted}/($self->{status}->{accepted}+
						$self->{status}->{rejected})));
    }
    push(@{$page},"</TD>");
    push(@{$page},"</TR>");
    push(@{$page},"<TR><TD COLSPAN=2>");
    push(@{$page},"<TABLE CELLSPACING=0 CELLPADDING=0>");

    foreach my $crec ( sort { $a->{cdat}->{cond}->{temp} <=> 
			      $b->{cdat}->{cond}->{temp} } @tlist ) {
      my $cid=$crec->{cid};
      my $cdat=$crec->{cdat};
      push(@{$page},"<TR>");
      my $color=($cdat->{exchanged}?"#FF5050":"#50FF50");
      push(@{$page},sprintf("<TD WIDTH=35>%s<B>%s</B></FONT></TD>",
			    &GenUtil::HTMLFont(2,"white"),$cid));
      push(@{$page},sprintf("<TD WIDTH=60 ALIGN=right>%s<B>%1.2f</B></FONT></TD>",
			    &GenUtil::HTMLFont(2,$color),$cdat->{cond}->{temp}));
      push(@{$page},sprintf("<TD WIDTH=60 ALIGN=right>%s%1.2f</FONT></TD>",
			    &GenUtil::HTMLFont(2,"white"),$self->{min}->{$cid}->{temp}));
      push(@{$page},sprintf("<TD WIDTH=60 ALIGN=right>%s%1.2f</FONT></TD>",
			    &GenUtil::HTMLFont(2,"white"),$self->{max}->{$cid}->{temp}));
      push(@{$page},sprintf("<TD WIDTH=85 ALIGN=right>%s<B>%1.2f</B></FONT></TD>",
			    &GenUtil::HTMLFont(2,$color),$cdat->{ener}));
      push(@{$page},sprintf("<TD WIDTH=60 ALIGN=right>%s<B>%1.2f</B></FONT></TD>",
			    &GenUtil::HTMLFont(2,"#FFFF50"),$cdat->{rmsd}));
      push(@{$page},sprintf("<TD WIDTH=60 ALIGN=right>%s<B>%1.2f</B></FONT></TD>",
			    &GenUtil::HTMLFont(2,"#FFFF50"),$cdat->{volume}));

      for (my $j=0; $j<$self->nUmbrellas(); $j++) {
	push(@{$page},sprintf("<TD WIDTH=60 ALIGN=right>%s<B>%1.2f</B></FONT></TD>",
			      &GenUtil::HTMLFont(2,$color),$cdat->{val}->[$j]));
	push(@{$page},sprintf("<TD WIDTH=60 ALIGN=right>%s<B>%1.2f</B></FONT></TD>",
			      &GenUtil::HTMLFont(2,$color),$cdat->{cond}->{bias}->[$j]->{target}));
	push(@{$page},sprintf("<TD WIDTH=60 ALIGN=right>%s%1.2f</FONT></TD>",
			      &GenUtil::HTMLFont(2,"white"),$self->{min}->{$cid}->{target}->[$j]));
	push(@{$page},sprintf("<TD WIDTH=60 ALIGN=right>%s%1.2f</FONT></TD>",
			      &GenUtil::HTMLFont(2,"white"),$self->{max}->{$cid}->{target}->[$j]));
      }

      push(@{$page},"</TR>");
    }
    }
    push(@{$page},"</TABLE>");
  } else {
    push(@{$page},"<TR><TD>");
    push(@{$page},"<FONT FACE=sans-serif,ARIAL,HELVETICA COLOR=white SIZE=3><B>server is starting up ...</B></FONT>");
  }
  push(@{$page},"</TD></TR>");
  push(@{$page},"</TABLE>");
  push(@{$page},"</CENTER>");
  push(@{$page},"</BODY>");
  push(@{$page},"</HTML>");
 
  return $page;
}


## function: expTempList(nwin,mintemp,maxtemp)
## generates an exponentially spaced temperature list
##  from <mark>mintemp</mark> to <mark>maxtemp</mark>.

sub expTempList {
  my $nwindows=shift;
  my $mintemp=shift;
  my $maxtemp=shift;

  my @arr=();

  if (defined $nwindows) {
    $mintemp=298.0 if (!defined $mintemp);
    $maxtemp=500.0 if (!defined $maxtemp);
    
    my $f=exp(log($maxtemp/$mintemp)/($nwindows-1));
    my $currT=$mintemp;
    for (my $i=0; $i<$nwindows; $i++) {
      push(@arr,sprintf("%1.4f",$currT));
      $currT*=$f;
    }
  }
  return \@arr;
}

1;
