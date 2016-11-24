# Amber package
# run Amber program
#
# http://mmtsb.scripps.edu/doc/Amber.pm.html
# 2000, Michael Feig, Brooks group, TSRI

package Amber;

require 5.004;

use strict;

use IO::Handle;
use IO::File;
use Fcntl;

use GenUtil;
use SICHO;
use Molecule;

## data: par -> { itersd, iterconj, etol, step,
## data:          gb, gbsa, cutoff, epsilon, dielec,
## data:          restrain }
## Amber minimization parameters

use vars qw ( $leapexec $sanderexec );

BEGIN {
  if (!exists $ENV{'AMBERHOME'} || !-d $ENV{'AMBERHOME'}) {
    die "Please set environment variable AMBERHOME";
  } else {
    $leapexec=(exists $ENV{'LEAPEXEC'}) ? 
      $ENV{'LEAPEXEC'} : "$ENV{'AMBERHOME'}/exe/teLeap";
    $sanderexec=(exists $ENV{'SANDEREXEC'}) ?
      $ENV{'SANDEREXEC'} : "$ENV{'AMBERHOME'}/exe/sander";
  }
}

## constructor: new()
## creates a new Amber object.

sub new {
  my $self={};
  
  die "leap executable not available"
    if (!-x $leapexec);
  die "sander executable not available"
    if (!-x $sanderexec);

  my %parhash = (
   param     => "ff94",   # force field: ff94, ff96, ff86, ff98, ff99, ff02, ff02EP, gaff
   sdsteps   => 50,       # number of minimization steps for steepest descent min.
   minsteps  => 300,      # number of minimization steps for ABNR min.
   minetol   => 1E-2,     # energy gradient tolerance for min.
   minstepsz => 0.005,    # initial min. step size
   minoutfrq => 10,       # minimization output frequency
   gb        => undef,    # use GB potential: tc6, tc, obc, jsb, mgb
   gbsa      => 0,        # include surface area term in GB potential
   gbeps     => 80,       # external dielectric for GB
   cutoff    => 14,       # cutoff for non-bonded interactions
   cut       => 1,        # use cutoff
   epsilon   => 1.0,      # dielectric constant
   dielec    => 'CDIE',   # dielectric: cdie, rdie
   dynsteps  => 100,      # dynamics steps
   dyntstep  => 0.002,    # dynamics step size
   dyntemp   => 300,      # dynamics temperature
   dynitemp  => undef,    # dynamics initial temperature
   dynseed   => undef,    # random number generator seed 
   dynoutfrq => 10,       # minimization output frequency
   restrain  => 0,        # restraints from restraint file
   rex       => 0,        # flag for running replica exchange
   lambda    => undef,    # lambda value for thermodynamic integration
   saltcon   => undef,    # salt concentration for use with GB
   rexmode   => undef,    # 1: temperature, 2: lambda
   dyntautp  => 1.0,      # temperature coupling 
   hsd       => "",
   hse       => ""
         );			 	 

  $self->{par}=\%parhash;
 
  bless $self;
  return $self;
}

## method: generateTopCoor(molecule,topfile,coorfile,[parameters])
## generates Amber topology and coordinate files from
## a given Molecule object

sub generateTopCoor {
  my $self=shift;
  my $mol=shift;
  my $top=shift;
  my $coor=shift;
  
  $self->_getpar(@_);
 
  my $tmppdb="tmp$$.pdb";
  $mol->writePDB($tmppdb); #,"NOH");  #,"AMBER");

  my $leapinp="tmp$$.leap.inp";

  my $param=$self->{par}->{param};

  die "cannot find force field input script for $param" 
    if (!-r "$ENV{AMBERHOME}/dat/leap/cmd/leaprc.$param");
  
  open FFINP,"$ENV{AMBERHOME}/dat/leap/cmd/leaprc.$param";

  open LEAPINP,">$leapinp";
  while (<FFINP>) {
    print LEAPINP;
  }
  close FFINP;

  my $gb=$self->{par}->{gb};
  if (defined $gb) {
    if ($gb eq "tc6") {
      print LEAPINP "set default PBradii amber6\n";
    } elsif ($gb eq "tc" || $gb == 1) {
      print LEAPINP "set default PBradii mbondi\n";
    } elsif ($gb eq "obc" || $gb eq "obc2") {
      print LEAPINP "set default PBradii bondi\n";
    } elsif ($gb eq "jsb") {
      print LEAPINP "set default PBradii gbjsb\n";
    } elsif ($gb eq "mgb") {
      print LEAPINP "set default PBradii mgbjsb\n";
    } else {
      die "unknown GB mode $gb\n";
    }
  }

  print LEAPINP "pdb = loadpdb $tmppdb\n";
  print LEAPINP "saveamberparm pdb $top $coor\n";
  print LEAPINP "quit\n";
  close LEAPINP;
  
  system "$leapexec -I$ENV{AMBERHOME}/dat -I$ENV{AMBERHOME}/dat/leap/lib -I$ENV{AMBERHOME}/dat/leap/cmd -I$ENV{AMBERHOME}/dat/leap/parm -f $leapinp > leap.out 2>&1";

  &GenUtil::remove($leapinp);
  &GenUtil::remove("leap.log");
  &GenUtil::remove($tmppdb);
}

## method: runSander(files) e,topfile,coorfile,outfile,logfile,restfile,
## method:           restartinfile,restartoutfile,trajinfile,trajoutfile)
## runs Amber to perform a minimization or dynamics run. Required input are
## a command input file, a topology file, and a coordinate file. A restraint
## file is optional. The final
## conformation after the simulation run is completed is written to 
## the output file given, program output is written to the log file.
## The methods <mark>generateTopCoor</mark>, <mark>genInputMinimize</mark>, and
## <mark>genInputRestraints</mark> are available to create the proper Amber 
## files expected by this method.

sub runSander {
  my $self=shift;
  my %file=@_;

#  my $inputfile=shift;
#  my $topfile=shift;
#  my $coorfile=shift;
#  my $outfile=shift;
#  my $logfile=shift;
#  my $restfile=shift;
#  my $restartinfile=shift;
#  my $restartoutfile=shift;
#  my $trajoutfile=shift;

  die "cannot read input file" unless (defined $file{input} && -r $file{input});
  die "cannot read parameter file" unless (defined $file{partop} && -r $file{partop});
  die "cannot read input coordinates" unless (defined $file{inpcoor} && -r $file{inpcoor});

  my $inffile="$$.mdinf";

  my $fileopt=sprintf("-i %s -p %s -c %s -O -o %s -r %s -inf %s",
		      $file{input},$file{partop},$file{inpcoor},$file{log},$file{outcoor},$inffile);

  if (defined $file{restraint}) {
    die "cannot read restraint file" unless (-r $file{restraint});
    $fileopt.=" -ref $file{restraint}";
  }

  if (defined $file{trajout}) {
    $fileopt.=" -x $file{trajout}";
  }

  if (defined $file{rexcontrol}) {
    $fileopt.=" -mmtsb $file{rexcontrol}";
  }

#  printf STDERR "running $sanderexec $fileopt\n";
#  system "$sanderexec $fileopt >& err$$.out";
  system "$sanderexec $fileopt 2> /dev/null";
  &GenUtil::remove($inffile);
#  sleep 3600;
}

## method: genInputDynamics(file[,parameters])
## generates an Amber command input file for minimization
## runs. A file name for the generated file is required as
## argument. In addition hash-type key=>value arguments
## may be given to set or modify minimization parameters.

sub genInputDynamics {
  my $self=shift;
  my $file=&GenUtil::getOutputFile(shift);
  my $restart=shift;
  my $trajout=shift;

  $self->_getpar(@_);
  
  my %gbmode = ( tc6 => 1,
		 tc  => 1,
		 obc => 2,
		 jsb => 3,
		 mgb => 4,
                 obc2 => 5,
                 1    => 1);

  print $file  " dynamics\n";
  print $file  " &cntrl\n";
  printf $file "  imin=0, ntmin=1, ntpr=%d,\n",
  $self->{par}->{dynoutfrq},$self->{par}->{itersd},$self->{par}->{iterconj}+$self->{par}->{itersd};
  if ($restart) {
    printf $file "  ntx=5, irest=1,\n";
  } else {
    printf $file "  ntx=1,\n";
  }
  if ($trajout) {
    printf $file "  ntwx=%d,\n",
      $self->{par}->{dynoutfrq};
  }

  printf $file "  cut=%1.1f,\n",($self->{par}->{cut})?$self->{par}->{cutoff}:9999.0;

  if (defined $self->{par}->{gb}) {
    printf $file "  igb=%d, gbsa=%d,\n",$gbmode{$self->{par}->{gb}},($self->{par}->{gbsa}==0)?0:1;
    printf $file "  surften=%f,\n",$self->{par}->{gbsa} 
      if ($self->{par}->{gbsa}!=0 && $self->{par}->{gbsa}!=1);
    printf $file "  extdiel=%f,\n",$self->{par}->{gbeps};
    printf $file "  saltcon=%f,\n",$self->{par}->{saltcon}
      if (defined $self->{par}->{saltcon} && $self->{par}->{saltcon}>0);
  } else {
    printf $file "  igb=0, dielc=%f,\n",$self->{par}->{epsilon};
  }

  $self->{par}->{dynitemp}=$self->{par}->{dyntemp}
    if (!defined $self->{par}->{dynitemp});
  $self->{par}->{dynseed}=2*int(rand 10000000)+1
    if (!defined $self->{par}->{dynseed});

  if ($self->{par}->{rex}) {
    printf $file  "  nstlim=999999999, dt=%f, temp0=%f, tempi=%f, ig=%d, \n",
      $self->{par}->{dyntstep},$self->{par}->{dyntemp},
	$self->{par}->{dynitemp}, $self->{par}->{dynseed};
    printf $file  "  mmtsb_switch=%d, mmtsb_iterations=%d, \n",
      $self->{par}->{rexmode},$self->{par}->{dynsteps};
  } else {
    printf $file  "  nstlim=%d, dt=%f, temp0=%f, tempi=%f, ig=%d, \n",
      $self->{par}->{dynsteps}, $self->{par}->{dyntstep},$self->{par}->{dyntemp},
	$self->{par}->{dynitemp}, $self->{par}->{dynseed};
  }

  if (defined $self->{par}->{lambda}) {
    printf $file "  icfe=1, clambda=%f,\n",$self->{par}->{lambda};
  }

  printf $file  "  ntt=1, ntc=2, ntf=2, tautp=%f\n",$self->{par}->{dyntautp};
  printf $file  "  ntb=0, ntr=%d\n",$self->{par}->{restrain}?1:0;
  print  $file  " &end\n";

  undef $file;
}

## method: genInputMinimize(file[,parameters])
## generates an Amber command input file for minimization
## runs. A file name for the generated file is required as
## argument. In addition hash-type key=>value arguments
## may be given to set or modify minimization parameters.

sub genInputMinimize {
  my $self=shift;
  my $file=&GenUtil::getOutputFile(shift);
  $self->_getpar(@_);
  
  my %gbmode = ( tc6 => 1,
		 tc  => 1,
		 obc => 2,
		 jsb => 3,
		 mgb => 4,
                 obc2 => 5,
                 1    => 1);

  print $file  " minimization\n";
  print $file  " &cntrl\n";
  printf $file "  imin=1, ntmin=1, ntpr=%d, ncyc=%d, maxcyc=%d,\n",
  $self->{par}->{minoutfrq},$self->{par}->{sdsteps},$self->{par}->{minsteps}+$self->{par}->{sdsteps};
  printf $file "  drms=%f, dx0=%f,\n",
  $self->{par}->{minetol},$self->{par}->{minstepsz};
  printf $file "  cut=%1.1f,\n",($self->{par}->{cut})?$self->{par}->{cutoff}:9999.0;

  if (defined $self->{par}->{gb}) {
    printf $file "  igb=%d, gbsa=%d,\n",$gbmode{$self->{par}->{gb}},($self->{par}->{gbsa}==0)?0:1;
    printf $file "  surften=%f,\n",$self->{par}->{gbsa} 
      if ($self->{par}->{gbsa}!=0 && $self->{par}->{gbsa}!=1);
    printf $file "  extdiel=%f,\n",$self->{par}->{gbeps};
    printf $file "  saltcon=%f,\n",$self->{par}->{saltcon}
      if (defined $self->{par}->{saltcon} && $self->{par}->{saltcon}>0);
  } else {
    printf $file "  igb=0, dielc=%f,\n",$self->{par}->{epsilon};
  }

  printf $file  "  ntx=1, ntb=0, ntr=%d\n",$self->{par}->{restrain}?1:0;
  print  $file  " &end\n";

  undef $file;
}


## method: genInputRestraints(file,cons)
## generates an Amber restraint input file from a
## <mark>cons</mark> data structure (see: <docmark>CHARMM.pm</docmark>).

sub genInputRestraints {
  my $self=shift;
  my $file=&GenUtil::getOutputFile(shift);
  my $cons=shift;

  my $selpat=&_getSel($cons->{sel});

  my $n=0;
  foreach my $c (@{$cons->{list}}) {
    printf $file "restraint %d\n",++$n;
    printf $file "%f\n",$c->{force};
    print  $file "FIND\n";
    printf $file "%s\n",$selpat;
    print  $file "SEARCH\n";
    printf $file "RES %d %d\n",$c->{from},$c->{to};
    print  $file "END\n";
  }
  print $file "END\n";

  undef $file;
}

## method: getEnergy(logfile)
## extracts energy information from Amber output log file

sub getEnergy {
  my $self=shift;
  my $logfile=&GenUtil::getInputFile(shift);

  my $erec={};

  my $on=0;
 DONE:
  while (<$logfile>) {
    $on=1 if (/FINAL RESULTS/);
    if ($on && /^[ \t]+NSTEP[ \t]+ENERGY[ \t]+RMS[ \t]+GMAX[ \t]+NAME[ \t]+NUMBER/) {
      my ($step,$total,$bond,$angle,$dihed,$vdwaals,$eel);
      my ($vdw14,$eel14,$constr,$egb,$esurf);
      $egb=$esurf=0.0;
      my $line=<$logfile>;

      ($step=substr($line,0,7))=~s/ +//g;
      ($total=substr($line,7,17))=~s/ +//g;

      <$logfile>;

      $line=<$logfile>;
      ($bond=substr($line,11,13))=~s/ +//g;
      ($angle=substr($line,36,13))=~s/ +//g;
      ($dihed=substr($line,64,13))=~s/ +//g;

      $line=<$logfile>;
      ($vdwaals=substr($line,11,13))=~s/ +//g;
      ($eel=substr($line,36,13))=~s/ +//g;
      ($egb=substr($line,64,13))=~s/ +//g if ($line =~ /EGB/);

      $line=<$logfile>;
      ($vdw14=substr($line,11,13))=~s/ +//g;
      ($eel14=substr($line,36,13))=~s/ +//g;
      ($constr=substr($line,64,13))=~s/ +//g;

      $line=<$logfile>;
      ($esurf=substr($line,11,13))=~s/ +//g if ($line =~/ESURF/);

      $erec->{total}=$total;
      $erec->{bonds}=$bond;
      $erec->{angles}=$angle;
      $erec->{dihedrals}=$dihed;
      $erec->{vdwaals}=$vdwaals+$vdw14;
      $erec->{vdw14}=$vdw14;
      $erec->{elec}=$eel+$eel14;
      $erec->{elec14}=$eel14;
      $erec->{gb}=$egb;
      $erec->{sasa}=$esurf;
      $erec->{constr}=$constr;
      last DONE;
    }
  }
  
  $logfile->close();

  return $erec;
}

## method: logEnergy(logfile,energylogfile)
## produces an energy log file from the Amber output log file

sub logEnergy {
  my $self=shift;
  my $logfile=&GenUtil::getInputFile(shift);
  my $efile=&GenUtil::getOutputFile(shift);

 ELOGIO:
  while (<$logfile>) {
    last ELOGIO if (/FINAL RESULTS/);
    if (/^[ \t]+NSTEP[ \t]+ENERGY[ \t]+RMS[ \t]+GMAX[ \t]+NAME[ \t]+NUMBER/) {
      my ($step,$total,$bond,$angle,$dihed,$vdwaals,$eel);
      my ($vdw14,$eel14,$constr,$egb,$esurf);
      my $line=<$logfile>;

      ($step=substr($line,0,7))=~s/ +//g;
      ($total=substr($line,7,17))=~s/ +//g;

      <$logfile>;

      $line=<$logfile>;
      ($bond=substr($line,11,13))=~s/ +//g;
      ($angle=substr($line,36,13))=~s/ +//g;
      ($dihed=substr($line,64,13))=~s/ +//g;

      $line=<$logfile>;
      ($vdwaals=substr($line,11,13))=~s/ +//g;
      ($eel=substr($line,36,13))=~s/ +//g;
      ($egb=substr($line,64,13))=~s/ +//g if ($line =~ /EGB/);

      $line=<$logfile>;
      ($vdw14=substr($line,11,13))=~s/ +//g;
      ($eel14=substr($line,36,13))=~s/ +//g;
      ($constr=substr($line,64,13))=~s/ +//g;
      
      if (defined $egb) {
	$line=<$logfile>;
	($esurf=substr($line,11,13))=~s/ +//g;
      } else {
	$egb=$esurf=0.0;
      }

      my $tag=($step<=$self->{par}->{itersd})?"SD":"CONJ";
      my $istep=($step<=$self->{par}->{itersd})?
	$step:$step-$self->{par}->{itersd};

      printf $efile "%-15s %4d %10.2f %10.2f %10.2f %10.2f %10.2f\n",
      $tag,$istep,$total,$vdwaals+$vdw14,$eel+$eel14,$egb+$esurf,$constr;
    }
  }
  
  $efile->close();
  $logfile->close();
}

### _getSel ######

sub _getSel {
  my $a=lc shift;
  
  if ($a eq "ca") {
    return "CA * * *";
  } elsif ($a eq "cb") {
    return "CB * * *";
  } elsif ($a eq "cab") {
    return "CA * * *\nCB * * *";
  } elsif ($a eq "heavy") {
    return "* N3 * *\n* CT * *\n* N2 * *\n* CA * *\n* C * *\n* O * *\n* N * *\n* SH * *\n* S * *\n* O2 * *\n* CC * *\n* NA * *\n* CR * *\n* CV * *\n* OH * *\n* NB * *\n* CW * *\n* CB * *\n* CN * *";
  } else {
    return "* * * *";
  }
}

### _getpar ######

sub _getpar {
  my $self=shift;
  my %arg=@_;
  foreach my $n ( keys %arg ) {
    if (exists $self->{par}->{$n}) {
      $self->{par}->{$n}=$arg{$n};
    } else {
      printf STDERR "Unknown Amber parameter $n will be ignored!\n";
    }
  }
}

1;
