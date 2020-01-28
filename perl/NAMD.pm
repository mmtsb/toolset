# NAMD package
# interact with NAMD program
#
# 2007, Michael Feig, Michigan State University

package NAMD;

require 5.004;

use strict;

use IO::Handle;
use IO::File;
use IPC::Open2;
use Sys::Hostname;
use Fcntl;

use GenUtil;
use Molecule;

## data: tag 
## temporary tag for input/output files

## data: par -> { ... } 
## NAMD parameters used in various commands

## data: mdener[] -> { step, time, temp, total, kine, 
## data:               pot, vdwaals, elec, gb, asp, constr }
## energy data extracted from the last molecular dynamics run

## data: molecule
## Molecule object for current molecule structure

use vars qw ( $exec $datadir );

BEGIN {
  if ($ENV{'NAMDEXEC'} ne "") {
    $exec=$ENV{'NAMDEXEC'};
  } else {
    $exec=&GenUtil::findExecutable("namd2");
  }
  
  if ($ENV{'CHARMMDATA'} ne "") {
    $datadir=$ENV{'CHARMMDATA'};
  } elsif ($datadir eq "" && $ENV{'MMTSBDIR'} ne "") {
    $datadir=$ENV{'MMTSBDIR'}."/data";
  } else {
    $datadir="./";
  }
}

## constructor: new()
## creates a new NAMD object. A log file for 
## NAMD output and a command log file for commands
## sent to NAMD may be given as an option.

sub new {
  my $self={};

  srand(time ^ ($$ + ($$ << 15)));

  bless $self;

  my @ef=split(/\s+/,$exec);
  die "cannot execute binary" if ((!-x $ef[0]) && ($exec !~ "mpi") && ($exec !~ "charmrun") );

  $self->{tag}="namd-".hostname."-$$";

  my %parhash = ( 
   bomlev    => undef,          # set bomb level
   prnlev    => undef,          # set print level

   param     => 36,             # parameter set: 19 (CHARMM19), 22 (CHARMM22), opls, a94 (amber 94), eef1, mmff
   xpar      => undef,          # extra parameter files
   xtop      => undef,
   patch     => undef,
   cmap      => 1,              # use CMAP correction
   geomvdw   => undef,          # geometric combination of sigmas for OPLS L-J terms
   scale14   => 1.0,            # 1-4 scaling

   epsilon   => 1.0,            # dielectric constant for non-bonded interactions
   cutnb     => undef,          # cutoff for non-bonded list generation
   cutoff    => undef,          # cutoff for non-bonded interactions
   cuton     => undef,          # onset of switching function for non-bonded interaction
   cut       => 1,              # use cutoff
   trunc     => undef,          # shift/fshift/switch/fswitch

   shake     => 0,              # use shake
   shaketol  => 1.0E-8,         # shake tolerance

   periodic  => 1,              # use periodic boundaries
   ewald     => 1,              # use Ewald for explicit solvent
   npme      => undef,          #64,             # number of grid points for PME grid
   npmex     => undef,          # number of grid points in x direction
   npmey     => undef,          # number of grid points in y direction
   npmez     => undef,          # number of grid points in z direction
   pmegrid   => undef,          # pme grid spacing
   margin    => undef,          # patch margin 
   ldb       => 1,              # load balancing

   boxshape  => undef,          # box shape for periodic boundaries
   boxsize   => undef,          # box size
   boxx      => undef,          # box x dimension
   boxy      => undef,          # box y dimension
   boxz      => undef,          # box z dimension

   dynens    => 'NVT',          # dynamics ensemble 'NVT', 'NPT', 'NPH', 'NVE'
   dynber    => 0,              # use Berendsen thermostat
   tcouplefile => undef,        # file with Berendsen time constant
   dyntstep  => 0.002,          # dynamics time step in picoseconds
   dynsteps  => 100,            # dynamics steps
   dyntemp   => 298,            # dynamics temperature
   dynseed   => undef,          # dynamics random seed
   dyntrfrq  => undef,          # trans/rot removal frequency
   dynremovecom => undef,       # remove center of mass motion when restarting
   dyneqfrq  => 0,              # velocity reassignment frequency
   dynscalefrq => 0,            # temperature rescaling frequency
   dynoutfrq => 500,            # frequency of energy/trajectory output
   dynouttime => 1000,          # frequency of writing timing information
   dynrestfrq=> 1000,           # frequency of saving restart files
   dyneoutfrq=> undef,          # frequency of writing out energy information
   dynupdnb  => 10,             # update of non-bonded list
   dynpress  => 1.00,           # the pressure (in atm) for NPT simulations
   dynpperiod => 200.0,         # gamma value for NPT simulations
   dynpdecay  => 100.0,         # Langevin piston mass for NPT simulations
   dynnpat    => 0,             # constant surface area for membrane simulations
  
   lang      => 1,              # flag for running Langevin dynamics
   langBAOAB => 0,              # BAOAB integrator from Leimkuhler et al.
   langfbeta => 1.0,            # friction coefficient for Langevin dynamics

   tmdrmsd   => undef,          # set initial TMD RMSD to turn it on
   tmdrmsdto => undef,          # set final TMD RMSD, by default equal to initial RMSD
   tmdk      => 100,            # force constant
   tmdoutfrq => undef,          # TMD output frequency
   tmdref    => undef,          # TMD reference file

   tmdrmsd2   => undef,          # set initial TMD RMSD to turn it on
   tmdrmsd2to => undef,          # set final TMD RMSD, by default equal to initial RMSD
   tmdk2      => 100,            # force constant
   tmdref2    => undef,          # TMD reference file

   mmtsbserver => undef,          # MMSTB server host name
   mmtsbport   => undef,          # MMTSB server port number
   mmtsbsid    => undef,          # MMTSB server ID
   mmtsbjobid  => undef,          # MMTSB job ID
   mmtsbcycle  => undef,          # replica exchange cycle length

   cons      => 0,              # set to 1 and provide reference file through -consref

   colvars     => undef,        # turn on and provide config file for colvars

   extrabonds  => undef,        # turn on and provide extra bonds input file

   minsteps  => undef,          # minimize
   mintinystep => 1.0E-6,
   minbabystep => 1.0E-2,
   minlinegoal => 1.0E-4,
  );

  $self->{par}=\%parhash;
  return $self;
}

DESTROY {
  my $self=shift;
}

### setParameter ######

sub setParameter {
  my $self=shift;
  $self->_getpar(@_);
}

### fixParameters #####

sub fixParameters {
  my $self=shift;

  if (defined $self->{par}->{boxsize}) { 
    $self->{par}->{boxshape}="cubic";
    $self->{par}->{boxx}=$self->{par}->{boxsize};
    $self->{par}->{boxy}=$self->{par}->{boxsize};
    $self->{par}->{boxz}=$self->{par}->{boxsize};
  } elsif (defined $self->{par}->{boxx} && defined $self->{par}->{boxy} && defined $self->{par}->{boxz}) {
    $self->{par}->{boxshape}="ortho" if (!defined $self->{par}->{boxshape});
  }

  if ($self->{par}->{cut} && !defined $self->{par}->{cutoff}) {
    $self->{par}->{cutoff}=10;
  } elsif (!$self->{par}->{cut}) {
    $self->{par}->{cutoff}=990;
  }
  if (!defined $self->{par}->{cuton}) {
    $self->{par}->{cuton}=$self->{par}->{cutoff}-1.5;
  }
  if (!defined $self->{par}->{cutnb}) {
    $self->{par}->{cutnb}=$self->{par}->{cutoff}+2.5;
  }

  $self->{par}->{dynseed}=2*int(rand 10000000)+1
    if (!defined $self->{par}->{dynseed});

  if (defined $self->{par}->{npme}) {
    $self->{par}->{npmex}=$self->{par}->{npme};
    $self->{par}->{npmey}=$self->{par}->{npme};
    $self->{par}->{npmez}=$self->{par}->{npme};
  }

  if (!defined $self->{par}->{pmegrid}) {
    if (defined $self->{par}->{boxx} && defined $self->{par}->{boxy} && defined $self->{par}->{boxz}) {
      if (!defined $self->{par}->{npmex}) {
	foreach my $n ( qw ( 16 18 24 32 36 48 54 64 72 96 108 128 144 162 192 216 256 288 324 384 432 486 512 576 648 768 864 972 1024 ) ) {
	  if ($n>$self->{par}->{boxx}*0.8 && !defined $self->{par}->{npmex})  {
	    $self->{par}->{npmex}=$n;
	  }
	}
      }
      
      if (!defined $self->{par}->{npmey}) {
	foreach my $n ( qw ( 16 18 24 32 36 48 54 64 72 96 108 128 144 162 192 216 256 288 324 384 432 486 512 576 648 768 864 972 1024 ) ) {
	  if ($n>$self->{par}->{boxy}*0.8 && !defined $self->{par}->{npmey})  {
	    $self->{par}->{npmey}=$n;
	  }
	}
      }
      
      if (!defined $self->{par}->{npmez}) {
	foreach my $n ( qw ( 16 18 24 32 36 48 54 64 72 96 108 128 144 162 192 216 256 288 324 384 432 486 512 576 648 768 864 972 1024 ) ) {
	  if ($n>$self->{par}->{boxz}*0.8 && !defined $self->{par}->{npmez})  {
	    $self->{par}->{npmez}=$n;
	  }
	}
      }
    } elsif (!defined $self->{par}->{npmex} || 
             !defined $self->{par}->{npmey} ||
             !defined $self->{par}->{npmez}) {
      $self->{par}->{pmegrid}=1.0;
    }
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
      printf STDERR "Unknown NAMD parameter $n will be ignored!\n";
    }
  }
}

## method: runDynamics($log,$cmd,options ...)

sub runDynamics {
  my $self=shift;
  my $log=shift;
  my $cmd=shift;
  my $psffile=shift;
  my $pdbfile=shift;
  my $coorfile=shift;
  my $extfile=shift;
  my $trajout=shift;
  my $outname=shift;
  my $restname=shift;
  my $enerout=shift;
  my $firststep=shift;
  my $customfile=shift;
  my $consfile=shift;
  
  my $cfgname=(defined $cmd)?$cmd:$self->{tag}.".inp";
  $self->writeConfigFile($cfgname,$psffile,$pdbfile,$coorfile,$extfile,$trajout,$outname,$restname,$firststep,$customfile,$consfile);
  
  my $logout;
  if (defined $log) {
    $logout=&GenUtil::getOutputFile($log);
    $logout->autoflush(1);
  }

  my $eout;
  if (defined $enerout) {
    $eout=&GenUtil::getOutputFile($enerout);
    $eout->autoflush(1);
  }

  open INP,"$exec $cfgname |";
  my @etags;
  while (<INP>) {
    print $logout $_ if (defined $log);

    if (defined $enerout) {
      if (/^ETITLE:/) {
	chomp;
	@etags=split(/\s+/);
      } elsif (/^ENERGY:/) {
	chomp;
	my @evals=split(/\s+/);
	my $e={};
	for (my $i=0; $i<=$#evals; $i++) {
	  $e->{$etags[$i]}=$evals[$i];
	}

	$eout->printf(
	 "%-9d %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f\n",
	 $e->{TS},($e->{TS}+$firststep)*$self->{par}->{dyntstep},
         $e->{TEMP},$e->{TOTAL},$e->{KINETIC},$e->{TOTAL}-$e->{KINETIC},
	 $e->{VDW},$e->{ELECT},$e->{BOND},$e->{ANGLE},$e->{DIHED},$e->{IMPRP},$e->{MISC},$e->{VOLUME},$e->{PRESSURE});      
      }
    }
  }
  close INP;

  close $logout if (defined $log);

  close $eout if (defined $enerout);

  &GenUtil::remove($cfgname) if (!defined $cmd);
}

## method: writeConfigFile(fname,psffile,pdbfile[,coorfile[,extfile[,trajout[,outname[,restname[,firststep]]]]]])
## writes NAMD config (input) file

sub writeConfigFile {
  my $self=shift;
  my $fname=shift;
  my $psffile=shift;
  my $pdbfile=shift;
  my $coorfile=shift;
  my $extfile=shift;
  my $trajout=shift;
  my $outname=shift;
  my $restname=shift;
  my $firststep=shift;
  my $customfile=shift;
  my $consfile=shift;

  $firststep=0 if (!defined $firststep);

  die "need PSF file\n" if (!defined $psffile || !-r $psffile);
  die "need PDB file\n" if (!defined $pdbfile || !-r $pdbfile);

  my $fout=&GenUtil::getOutputFile($fname);

  printf $fout "structure                %s\n",$psffile;
  printf $fout "coordinates              %s\n",$pdbfile;
  printf $fout "\n";
  printf $fout "outputName               %s\n",$outname;
  printf $fout "binaryoutput             yes\n";
  printf $fout "binaryrestart            yes\n";
  printf $fout "restartname              %s\n",$outname;
  printf $fout "restartfreq              %d\n",$self->{par}->{dynrestfrq};
  printf $fout "xstfreq                  %d\n";$self->{par}->{dynrestfrq};
  printf $fout "DCDfile                  %s\n",$trajout if (defined $trajout);
  printf $fout "DCDfreq                  %d\n",$self->{par}->{dynoutfrq} if (defined $trajout);
  printf $fout "outputEnergies           %d\n",
    (defined ($self->{par}->{dyneoutfrq}))?$self->{par}->{dyneoutfrq}:$self->{par}->{dynoutfrq};
  printf $fout "outputTiming             %d\n",$self->{par}->{dynouttime};
  printf $fout "\n";

  if (defined $restname && -r $restname.".coor" && -r $restname.".vel" && -r $restname.".xsc") {
    printf $fout "binCoordinates           %s.coor\n",$restname if (-r $restname.".coor");
    printf $fout "binVelocities            %s.vel\n",$restname  if (-r $restname.".vel");
    printf $fout "extendedSystem           %s.xsc\n",$restname  if (-r $restname.".xsc");
    printf $fout "firsttimestep            %d\n",$firststep;
  } else {
    printf $fout "temperature              %f\n",$self->{par}->{dyntemp};
    printf $fout "binCoordinates           %s\n",$coorfile if (defined $coorfile && -r $coorfile);

    if (defined $extfile && -r $extfile) {
      printf $fout "extendedSystem           %s\n",$extfile;
    } else {
      if ($self->{par}->{periodic}) {
	if ($self->{par}->{boxshape} eq "cubic" || $self->{par}->{boxshape} eq "ortho") {  
	  printf $fout "cellBasisVector1      %f 0 0\n",$self->{par}->{boxx};
	  printf $fout "cellBasisVector2      0 %f 0\n",$self->{par}->{boxy};
	  printf $fout "cellBasisVector3      0 0 %f\n",$self->{par}->{boxz};
	  printf $fout "cellOrigin            0 0 0\n";
	}
      }
    }
  }

  if (defined $self->{par}->{colvars}) {
    if (-r $self->{par}->{colvars}) {
      printf $fout "colvars         on\n";
      printf $fout "colvarsConfig   %s\n",$self->{par}->{colvars};
      printf $fout "colvarsInput    %s.colvars.state\n",$restname if (defined $restname && -r $restname.".colvars.state");
    } else {
      printf STDERR "cannot read colvars config file >%s<\n",$self->{par}->{colvars};
    }
  }

  printf $fout "\n";

  printf $fout "paraTypeCharmm           on\n";
  
  my @parfiles=$self->getParameterFiles();
  foreach my $p ( @parfiles ) { 
    printf $fout "parameters               %s\n",$p;
  }

  printf $fout "vdwGeometricSigma        on\n" if (defined $self->{par}->{geomvdw} && $self->{par}->{geomvdw});
  
  printf $fout "exclude                  scaled1-4\n";
  printf $fout "1-4scaling               %f\n",$self->{par}->{scale14};
  printf $fout "switching                on\n";

  printf $fout "cutoff                   %f\n",$self->{par}->{cutoff};
  printf $fout "switchdist               %f\n",$self->{par}->{cuton};
  printf $fout "pairlistdist             %f\n",$self->{par}->{cutnb};

  printf $fout "\n";

  printf $fout "stepspercycle            %d\n",$self->{par}->{dynupdnb};

  printf $fout "timestep                 %f\n",$self->{par}->{dyntstep}*1000.0;

  printf $fout "seed                     %d\n",$self->{par}->{dynseed};

  printf $fout "\n";

  if ($self->{par}->{shake}) {
    printf $fout "rigidBonds               all\n";
    printf $fout "rigidTolerance           %e\n",$self->{par}->{shaketol};
  } else {
    printf $fout "rigidBonds               none\n";
  }

  printf $fout "\n";

  printf $fout "nonbondedFreq            1\n";
  printf $fout "fullElectFrequency       1\n";

  printf $fout "\n";

  printf $fout "wrapWater                on\n";
  printf $fout "wrapAll                  on\n";

  printf $fout "margin                   %f\n",$self->{par}->{margin}
    if (defined $self->{par}->{margin});

  printf $fout "\n";

  if ($self->{par}->{ewald}) {
    printf $fout "PME                      yes\n";
    if (defined $self->{par}->{pmegrid}) {
      printf $fout "PMEGridSpacing           %f\n",$self->{par}->{pmegrid};
    } elsif (defined $self->{par}->{npmex} && defined $self->{par}->{npmey} && defined $self->{par}->{npmez}) {
      printf $fout "PMEGridSizeX             %d\n",$self->{par}->{npmex};
      printf $fout "PMEGridSizeY             %d\n",$self->{par}->{npmey};
      printf $fout "PMEGridSizeZ             %d\n",$self->{par}->{npmez};
    }
  } else {
    printf $fout "PME                      no\n";
  }

  printf $fout "\n";

  printf $fout "useGroupPressure         yes\n";
  if ($self->{par}->{dynnpat}) {
    printf $fout "useFlexibleCell          yes\n";
    printf $fout "useConstantArea          yes\n";
    printf $fout "COMmotion                yes\n";
    printf $fout "wrapNearest              on\n";
  } else {
    printf $fout "useFlexibleCell          no\n";
    printf $fout "useConstantArea          no\n";
    printf $fout "wrapNearest              off\n";
  }
  printf $fout "useConstantRatio         no\n";

  printf $fout "fixedAtoms               off\n";

  if (defined $self->{par}->{dynremovecom} && !$self->{par}->{dynremovecom}) {
    printf $fout "COMmotion                yes\n";
  } else {
    printf $fout "COMmotion               no\n";
  } 
  printf $fout "zeroMomentum             yes\n" if ($self->{par}->{dyntrfrq}>0);

  printf $fout "\n";

  if ((uc $self->{par}->{dynens}) eq "NVT" || (uc $self->{par}->{dynens}) eq "NPT") {
    if ($self->{par}->{lang}) {
      printf $fout "langevin                 on\n";
      printf $fout "langevinTemp             %f\n",$self->{par}->{dyntemp};
      printf $fout "langevinDamping          %f\n",$self->{par}->{langfbeta};
      printf $fout "langevinBAOAB            on\n" if ($self->{par}->{langBAOAB});
    } elsif ($self->{par}->{dynber}) {
      printf $fout "tCouple                  on\n";
      printf $fout "tCoupleTemp              %f\n",$self->{par}->{dyntemp};
      if (defined $self->{par}->{tcouplefile} && -r $self->{par}->{tcouplefile}) {
        printf $fout "tCoupleFile              %s\n",$self->{par}->{tcouplefile};
      } else {
        printf $fout "tCoupleFile              %s\n",$pdbfile;
      }
      printf $fout "tCoupleCol               O\n";
    } elsif ($self->{par}->{dyneqfrq}>0) {
      printf $fout "reassignFreq             %d\n",$self->{par}->{dyneqfrq};
      printf $fout "reassignTemp             %f\n",$self->{par}->{dyntemp};
    } elsif ($self->{par}->{dynscalefrq}>0) {
      printf $fout "rescaleFreq              %d\n",$self->{par}->{dynscalefrq};
      printf $fout "rescaleTemp              %f\n",$self->{par}->{dyntemp};
    }
  }

  printf $fout "\n";

  if ((uc $self->{par}->{dynens} eq "NPT") || (uc $self->{par}->{dynens}) eq "NPH") {
    printf $fout "LangevinPiston           on\n";
    printf $fout "LangevinPistonTarget     %f\n",$self->{par}->{dynpress};
    printf $fout "LangevinPistonPeriod     %f\n",$self->{par}->{dynpperiod};
    printf $fout "LangevinPistonDecay      %f\n",$self->{par}->{dynpdecay};
    printf $fout "LangevinPistonTemp       %f\n",$self->{par}->{dyntemp};
  }

  printf $fout "\n";

  if (defined $self->{par}->{tmdrmsd}) {
    printf $fout "\n";
    printf $fout "TMD                     on\n";
    printf $fout "TMDk                    %f\n",$self->{par}->{tmdk};
    printf $fout "TMDOutputFreq           %d\n",
      (defined $self->{par}->{tmdoutfrq})?$self->{par}->{tmdoutfrq}:$self->{par}->{dynoutfrq};
    printf $fout "TMDFile                 %s\n",
      (defined $self->{par}->{tmdref} && -r $self->{par}->{tmdref})?$self->{par}->{tmdref}:$pdbfile;
    printf $fout "TMDFirstStep            0\n";
    printf $fout "TMDLastStep             %d\n",$self->{par}->{dynsteps};
    printf $fout "TMDInitialRMSD          %f\n",$self->{par}->{tmdrmsd};
    printf $fout "TMDFinalRMSD            %f\n",
      (defined $self->{par}->{tmdrmsdto})?$self->{par}->{tmdrmsdto}:$self->{par}->{tmdrmsd};
  }

  if (defined $self->{par}->{tmdrmsd2}) {
    printf $fout "\n";
    printf $fout "TMD2                    on\n";
    printf $fout "TMDk2                   %f\n",$self->{par}->{tmdk2};
    printf $fout "TMDFile2                %s\n",
      (defined $self->{par}->{tmdref2} && -r $self->{par}->{tmdref2})?$self->{par}->{tmdref2}:$pdbfile;
    printf $fout "TMDInitialRMSD2         %f\n",$self->{par}->{tmdrmsd2};
    printf $fout "TMDFinalRMSD2           %f\n",
      (defined $self->{par}->{tmdrmsd2to})?$self->{par}->{tmdrmsd2to}:$self->{par}->{tmdrmsd2};
  }


  if (defined ($self->{par}->{mmtsbserver})) {
    printf $fout "MMTSB                   On\n";
    printf $fout "MMTSBHost               %s\n",$self->{par}->{mmtsbserver};
    printf $fout "MMTSBPort               %s\n",$self->{par}->{mmtsbport};
    printf $fout "MMTSBServerId           %s\n",$self->{par}->{mmtsbsid};
    printf $fout "MMTSBJobId              %s\n",$self->{par}->{mmtsbjobid};
    printf $fout "MMTSBCycle              %d\n",$self->{par}->{mmtsbcycle};
  }

  if (defined $self->{par}->{cons} && $self->{par}->{cons}) {
    printf $fout "\n";
    printf $fout "constraints              on\n";
    printf $fout "consexp                  2\n";
    printf $fout "consref                  %s\n",
      (defined $consfile && -r $consfile)?$consfile:$pdbfile;
    printf $fout "conskfile                %s\n",
      (defined $consfile && -r $consfile)?$consfile:$pdbfile;
    printf $fout "conskcol                 B\n";
  }

  if (defined $self->{par}->{extrabonds} && -r $self->{par}->{extrabonds}) {
    printf $fout "extraBonds on\n";
    printf $fout "extraBondsFile %s\n",$self->{par}->{extrabonds};
  }
 
  if (!$self->{par}->{ldb}) {
#    printf $fout "LdbStrategy other\n";
    printf $fout "LdbPeriod 10000000000\n";
    printf $fout "FirstLdbStep 100000000000\n";
  }

#  printf $fout "printBadContacts    True\n";
  if (defined $customfile && -r $customfile) {
    open CINP,"$customfile";
    while (<CINP>) {
      print $fout $_;
    }
    close CINP;
  }

  printf $fout "\n";

  if (defined $self->{par}->{minsteps} && $self->{par}->{minsteps}>0) {
#    printf $fout "minimization             on\n";
    printf $fout "minTinyStep              %e\n",$self->{par}->{mintinystep};
    printf $fout "minBabyStep              %e\n",$self->{par}->{minbabystep};
    printf $fout "minLineGoal              %e\n",$self->{par}->{minlinegoal};
    printf $fout "minimize                 %d\n",$self->{par}->{minsteps};
    printf $fout "reinitvels               %f\n",$self->{par}->{dyntemp};
    printf $fout "run                      %d\n",$self->{par}->{dynsteps};
  } else {
    printf $fout "numsteps                 %d\n",$self->{par}->{dynsteps};
  }

  printf $fout "\n";

  close $fout;
  undef $fout;
}

## method: getParameterFiles 
## return the parameter files to be loaded

sub getParameterFiles {
  my $self=shift;
  
  my @para=();

  if ($self->{par}->{param} eq "19") {
    push(@para,"$datadir/param19.inp");
  } elsif ($self->{par}->{param} eq "22") {
    if ($self->{par}->{cmap} && 
        -r "$datadir/par_all22_prot_cmap.inp") {
      push(@para,"$datadir/par_all22_prot_cmap.inp");
    } else {
      push(@para,"$datadir/par_all22_prot.inp");
    }
  } elsif ($self->{par}->{param} eq "opls") {
    push(@para,"$datadir/opls.prm");
    $self->{par}->{geomvdw}=1;
  } elsif ($self->{par}->{param} eq "a94") {
    push(@para,"$datadir/cornell_all.prm");
    $self->{par}->{scale14}=0.83;
  } elsif ($self->{par}->{param} eq "lpdb") {
    push(@para,"$datadir/parm.prm-lpdb");
  } elsif ($self->{par}->{param} eq "simdb") {
    push(@para,"$datadir/simdb.prm");
  } elsif ($self->{par}->{param} eq "27") {
    if ($self->{par}->{cmap} && 
        -r "$datadir/par_all27_prot_na_cmap.prm") {
      push(@para,"$datadir/par_all27_prot_na_cmap.prm");
    } else {
      push(@para,"$datadir/par_all27_prot_na.prm");
    }
  } elsif ($self->{par}->{param} eq "36") {
    push(@para,"$datadir/par_all36_prot.prm");
    push(@para,"$datadir/par_all36_na.prm");
  }
  
  if (defined $self->{par}->{xpar}) {
    foreach my $n ( split(/:/,$self->{par}->{xpar}) ) {
      $n="$datadir/$n" if (!-r $n);
      push(@para,$n) if (-r $n);
    }
  }

  die "no parameter file available\n" if ($#para<0);

  return @para;
}

1;
