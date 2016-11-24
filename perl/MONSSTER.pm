# MONSSTER package
# setup and run MONSSTER program
#
# http://mmtsb.scripps.edu/doc/MONSSTER.pm.html
# 2000, Michael Feig, Brooks group, TSRI

package MONSSTER;

require 5.004;

use strict;

use IPC::Open2;

use GenUtil;
use Sequence;
use SICHO;

## data: sequence
## Sequence object used to generate input for MONSSTER

## data: workdir
## Working directory for running MONSSTER

## data: par -> { random, ncycle, icycle, tsteps,
## data:          resmin, resmax, temp, softcore, 
## data:          central, stiff, pair, kdcore,
## data:          hbond, short, burial,
## data:          multibody, threebody,
## data:          rgyr, krgyr, rho, krho }
## MONSSTER simulation parameters

## data: posrest[] -> { from, to, force }
## list of positional restraints

## data: distrest[] -> { from, to }
## list of distance restraints

## data: distrestforce
## restraint force for distance restraints

## data: output -> { final  -> { short, long, total, central },
## data:             list[] -> { cycle[] -> { r2, s2, ener, ecent }, 
## data:                         step, temp } }
## simulation data extracted from MONSSTER output

use vars qw ( %fileName );

BEGIN {
  %fileName = (
     datadir    => "monsster.datadir",
     input      => "monsster.input",
     output     => "monsster.output",
     sequence   => "monsster.seq",
     initchain  => "monsster.init.chain",
     finalchain => "monsster.final.chain",
     trajectory => "monsster.tra",
     restraints => "monsster.restraints"
  );

}

## constructor: new([sequence]);
## creates a new MONSSTER object. A Sequence object
## may be given as argument

sub new {
  my $self={};

  $self->{sequence}=shift;
  $self->{distrestforce}=1.0;
  $self->{workdir}=".";

  my %defaultpar= (
    random     => -1,
    ncycle     => 20,
    icycle     => 50,
    tsteps     => 10,
    resmin     => 0,
    resmax     => 0,
    temp       => "2.5:1.0",
    softcore   => 4.0,
    central    => 0.25,
    stiff      => 1.25,
    pair       => 1.0,
    kdcore     => 0.35,
    hbond      => -0.3,
    short      => 0.5,
    burial     => 1.0,
    multibody  => 0.5,
    threebody  => 0.25,
    rgyr       => 0.0,
    krgyr      => 0.0,
    rho        => 0.0,
    krho       => 0.0
  );

  $self->{par}=\%defaultpar;
  $self->{output}->{_readData}=1;

  bless $self;
  return $self;


}

## method: setDirectory(directory)
## set working directory

sub setDirectory {
  my $self=shift;
  my $newdir=shift;

  if (defined $newdir && $newdir ne $self->{workdir}) {
    &GenUtil::makeDir($newdir)
      if (!-d $newdir);
    $self->{workdir}=$newdir;
    $self->{output}->{_readData}=1;
  }
}

## method: $status = complete([requiretrajectory])
## determines files from a MONSSTER simulation run are
## complete.
## If the argument is set to a non-zero value it
## requires the trajectory file to be present in addition
## to the final chain and output file.

sub complete {
  my $self=shift;
  my $requiretra=shift;

  my $d=$self->{workdir};

  $requiretra=0 if (!defined $requiretra);

  return (&GenUtil::checkFile("$d/$fileName{output}") &&
	  &GenUtil::checkFile("$d/$fileName{finalchain}") &&
	  (!$requiretra || &GenUtil::checkFile("$d/$fileName{trajectory}")));
}

## method: setParameter(parameters)
## sets MONSSTER simulation parameters that
## are expected as hash-type key=>value pair arguments.

sub setParameters {
  my $self=shift;
  my %par=@_;

  foreach my $n ( keys %par ) {
    if (defined $self->{par}->{$n}) {
      $self->{par}->{$n}=$par{$n};
    }
  }
}

## method: setPositionalRestraints(list)
## sets the positional restraint list

sub setPositionalRestraints {
  my $self=shift;
  my $posrest=shift;

  @{$self->{posrest}}=@{$posrest};

  $self->_translateList("posrest");
}

## method: setDistanceRestraints(force,list)
## sets the distance restraint list and force

sub setDistanceRestraints {
  my $self=shift;
  my $force=shift;
  my $distrest=shift;

  @{$self->{distrest}}=@{$distrest};
  $self->{distrestforce}=$force;

  $self->_translateList("distrest");
}

sub _translateList {
  my $self=shift;
  my $listname=shift;
  
  return if (!defined $self->{sequence});

  my $offset=$self->{sequence}->{sequence}->[0]->{index}-1;

  return if ($offset==0);

  foreach my $l ( @{$self->{$listname}} ) {
    $l->{from}-=$offset;
    $l->{to}-=$offset;
  }
}

## method: setup(input)
## sets up all necessary files to run MONSSTER.
## An input SICHO lattice chain is required as argument.

sub setup {
  my $self=shift;
  my $input=shift;

  die "sequence not defined" 
    if (!defined $self->{sequence} && !ref $self->{sequence});

  $input->writeChain("$self->{workdir}/$fileName{initchain}");
  $self->{sequence}->writeMONSSTER("$self->{workdir}/$fileName{sequence}");
  $self->writeInput("$self->{workdir}/$fileName{input}");

  if (defined $self->{posrest}) {
    $self->writePosRestraints("$self->{workdir}/$fileName{restraints}");
  } elsif (-r "$self->{workdir}/$fileName{restraints}") {
    unlink "$self->{workdir}/$fileName{restraints}";
  }

  my $mmtsbdata=&GenUtil::findDataDir();
  die "cannot find data directory for MONSSTER" 
    if (!-d $mmtsbdata);

  open DDIROUT,">$self->{workdir}/$fileName{datadir}";
  print DDIROUT "$mmtsbdata\n";
  close DDIROUT;
}

## method: energy(input)
## sets up and runs mener to get MONSSTER energies
## An input SICHO lattice chain is required as argument.

sub energy {
  my $self=shift;
  my $input=shift;

  die "sequence not defined" 
    if (!defined $self->{sequence} && !ref $self->{sequence});

  my $menerbin=&GenUtil::findExecutable("mener");
  die "cannot find MONSSTER energy binary"
    if (!defined $menerbin);

  my $mmtsbdata=&GenUtil::findDataDir();
  die "cannot find data directory for MONSSTER" 
    if (!-d $mmtsbdata);

  local(*READ,*WRITE);
  my $pid=open2(*READ,*WRITE,"$menerbin");

  my @t=split(/:/,$self->{par}->{temp});

  printf WRITE "%s\n",$mmtsbdata;
  printf WRITE "%f %f %f\n",
  $self->{par}->{softcore},$self->{par}->{central},$self->{par}->{stiff};
  printf WRITE "%f %f %f\n",
  $self->{par}->{pair},$self->{par}->{kdcore},$self->{par}->{hbond};
  printf WRITE "%f %f %f\n",
  $self->{par}->{short},$self->{par}->{burial},$self->{par}->{multibody};
  printf WRITE "%f %f\n",
  $self->{par}->{threebody},$t[$#t];

  printf WRITE "%d\n",$#{$self->{sequence}->{sequence}}+1;

  $self->{sequence}->writeMONSSTER(\*WRITE);

  $input->writeChain(\*WRITE);
  
  close WRITE;

  my $ret=<READ>;
  close READ;
  
  waitpid($pid,0);

  $ret=~s/^ +//;
  return split(/ +/,$ret);
}

## method: run(input)
## sets up and runs a MONSSTER simulation
## An input SICHO lattice chain is required as argument.

sub run {
  my $self=shift;
  my $input=shift;

  my $monssterbin=&GenUtil::findExecutable("mfold");
  die "cannot find MONSSTER"
    if (!defined $monssterbin);

  $self->setup($input);

  system "cd $self->{workdir}; $monssterbin";

  $self->{output}->{_readData}=1;
}

## method: initChain()
## returns the input chain from a MONSSTER run as SICHO object

sub initChain {
  my $self=shift;
  
  my $chain=&SICHO::new();
  $chain->readChain("$self->{workdir}/$fileName{initchain}");
  return $chain;
}

## method: finalChain()
## returns the final chain from a MONSSTER run as SICHO object

sub finalChain {
  my $self=shift;
  
  my $chain=&SICHO::new();
  $chain->readChain("$self->{workdir}/$fileName{finalchain}");
  return $chain;
}

## method: trajFrame(timestep,cycle)
## reads a MONSSTER trajectory file and returns the frame
## corresponding to the time step and cycle arguments

sub trajFrame {
  my $self=shift;
  my $tstep=shift;
  my $cycle=shift;

  my $chain=&SICHO::new();
  $chain->readMONSSTERTraj("$self->{workdir}/$fileName{trajectory}",$tstep,$cycle);
  return $chain;
}

## method: finalEnergy()
## returns the final energy from the MONSSTER output

sub finalEnergy { 
  my $self=shift;
  
  $self->_processOutput() if ($self->{output}->{_readData});

  return $self->{output}->{final};
}

## method: nSteps()
## returns the number of time steps from the MONSSTER output

sub nSteps {
  my $self=shift;

  $self->_processOutput() if ($self->{output}->{_readData});

  return $#{$self->{output}->{list}}+1;
}

## method: nCycle()
## returns the number of cycles at each time step from the
## MONSSTER output

sub nCycle {
  my $self=shift;
  
  $self->_processOutput() if ($self->{output}->{_readData});
  
  return $#{$self->{output}->{list}->[0]->{cycle}}+1;
}

## method: $step = stepInfo(index)
## returns a data structure with MONSSTER output for a 
## given time step

sub stepInfo {
  my $self=shift;
  my $step=(shift)-1;

  die "step index out of range" 
    if ($step<0 || $step>$#{$self->{output}->{list}});
  
  $self->_processOutput() if ($self->{output}->{_readData});

  return $self->{output}->{list}->[$step];
}

## method: $cycle = cycleInfo(index)
## returns a data structure with MONSSTER output for a
## given cycle

sub cycleInfo {
  my $self=shift;
  my $step=(shift)-1;
  my $cycle=(shift)-1;

  die "step index out of range" 
    if ($step<0 || $step>$#{$self->{output}->{list}} ||
	$cycle<0 || $cycle>$#{$self->{output}->{list}->[$step]->{cycle}});

  $self->_processOutput() if ($self->{output}->{_readData});

  return $self->{output}->{list}->[$step]->{cycle}->[$cycle];
}

## method: writeInput(file)
## writes a MONSSTER input parameter file from
## previously set parameters.

sub writeInput {
  my $self=shift;
  my $monssterinp=&GenUtil::getOutputFile(shift);

  my $random;
  if ($self->{par}->{random}<=0) {
    srand();
    $random=int(rand(1000000));
  } else {
    $random=$self->{par}->{random};
  }

  printf $monssterinp "%d %d %d %d\n",
  $random, $self->{par}->{ncycle}, 
  $self->{par}->{icycle}, $self->{par}->{tsteps};

  printf $monssterinp "%d %d\n",
  $self->{par}->{resmin}, $self->{par}->{resmax};

  printf $monssterinp "%1.6f %1.6f %1.6f %1.6f\n",
  $self->{par}->{rgyr},$self->{par}->{krgyr},
  $self->{par}->{rho},$self->{par}->{krho};
  
  my @t=split(/:/,$self->{par}->{temp});
  printf $monssterinp "%1.2f %1.2f %1.3f %1.3f\n",
  $t[0], ((defined $t[1])?$t[1]:$t[0]), 
  $self->{par}->{softcore}, $self->{par}->{central};

  printf $monssterinp "%1.3f %1.3f %1.3f %1.3f %1.3f\n",
  $self->{par}->{stiff}, $self->{par}->{pair}, $self->{par}->{kdcore}, 
  $self->{par}->{hbond}, $self->{par}->{short};

  printf $monssterinp "%1.3f %1.3f %1.3f\n",
  $self->{par}->{burial}, $self->{par}->{multibody}, $self->{par}->{threebody};

  if (defined $self->{distrest}) {
    printf $monssterinp "%d %f\n",
      $#{$self->{distrest}}+1,$self->{distrestforce};
    foreach my $n (@{$self->{distrest}}) {
      printf $monssterinp "%d %d\n",$n->{from},$n->{to};
    }
  } else {
    print $monssterinp "0 0.0\n";
  }
  print $monssterinp "0\n0\n0\n0\n";  
  undef $monssterinp;
}

## method: writePosRestraints(file)
## writes a MONSSTER restraint file from a
## previously set positional restraint list

sub writePosRestraints {
  my $self=shift;
  my $restrfile=&GenUtil::getOutputFile(shift);
  
  my $cnt=0;
  foreach my $n (@{$self->{posrest}}) {
    $cnt+=($n->{to}-$n->{from}+1);
  }
  
  printf $restrfile "%5d\n",$cnt;
  foreach my $n (@{$self->{posrest}}) {
    for (my $i=$n->{from}; $i<=$n->{to}; $i++) {
      printf $restrfile "%5d%9.5f\n",$i,$n->{force};
    }
  }
  undef $restrfile;
}

sub _processOutput {
  my $self=shift;
  
  my $fname="$self->{workdir}/$fileName{output}";

  my $outinp=&GenUtil::getInputFile($fname);

  $self->{output}->{list}=();
  $self->{output}->{final}={};

  my $curr;

  while(<$outinp>) {
    if (/step +([0-9]+) +Temperature = +([0-9\.])/) {
      push (@{$self->{output}->{list}},$curr) if (defined $curr);
      $curr={};
      $curr->{step}=$1;
      $curr->{temp}=$2;
      $curr->{cycle}=();
    } elsif (defined $curr && /^( +[0-9\.\-]+){9}/) {
      chomp;
      my @f=split(/ +/);
      my $rec={};
      $rec->{r2}=$f[2];
      $rec->{s2}=$f[3];
      $rec->{ener}=$f[8];
      $rec->{ecent}=$f[9];
      push (@{$curr->{cycle}},$rec);
    } elsif (/TEST OF EXCLUDED VOLUME/) {
      push (@{$self->{output}->{list}},$curr) if (defined $curr);
      undef $curr;
    } elsif (/ +short range +([0-9\.\-]+)/) {
      $self->{output}->{final}->{short}=$1;
    } elsif (/ +long range +([0-9\.\-]+)/) {
      $self->{output}->{final}->{long}=$1;
    } elsif (/ +FINAL ENERGY = +([0-9\-\.]+) +([0-9\-\.]+)/) {
      $self->{output}->{final}->{total}=$1;
      $self->{output}->{final}->{central}=$2;
    }   
  }
  undef $outinp;
  
  $self->{output}->{_readData}=0;
}


1;
