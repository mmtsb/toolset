# Ensemble package
# manage ensemble properties for lattice simulation/minimization runs
#
# http://mmtsb.scripps.edu/doc/Ensemble.pm.html
# 2000, Michael Feig, Brooks group, TSRI
#
# derived from SimData

package Ensemble;

require 5.004;

use strict;

use Fcntl ':flock';

use SimData;

use GenUtil;
use Molecule;
use Analyze;
use CHARMM;

## data: tag
## structure tag to identify a set of structures in an ensemble

## data: prop -> { etot[], eelec[], evdw[], egb[], esasa[], easp[],
## data:           rgyr[], rmsdpol[], rmsdhyd[], rmsdchg[],
## data:           rmsdall[], rmsdca[], rmsdback[], 
## data:           pphi[], ppsi[], rmsdphi[], rmsdpsi[],
## data:           pchi1[], rmsdchi1[], cont[], rho[],
## data:           gdtts[] }
## ensemble strucuture property values. Shown are commonly 
## used property tags but
## actually available tags depend on data availability.

## data: opt -> { ... }
## application defined options associated with the ensemble tag

## data: optpar -> { ... }
## application defined parameters associated with the ensemble tag.
## parameters are intended to be used to separate CHARMM and MONSSTER
## simulation parameters from other more general options

## data: filelist[]
## list of original file names if external files were checked in

use vars qw ( @ISA $PropFileSuffix $OptionFileSuffix );

@ISA = ( "SimData" );

BEGIN {
  $PropFileSuffix=".prop.dat";
  $OptionFileSuffix=".options";
}

## constructor: new(tag,dir[,configfile])
## creates a new Ensemble object. Required arguments are a tag and
## the ensemble directory. The ensemble configuration file as well 
## as option and property files associated with the tag are read
## automatically.

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;

  my $tag=shift;

  my $self=$class->SUPER::new(@_);

  $self->{tag}=$tag;

  $self->set(runs=>0) 
    if (!defined $self->{par}->{runs});

  $self->set(compress=>1)
    if (!defined $self->{par}->{compress});

  $self->readProp();
  $self->readOptions();

  undef $self->{_save};
  undef $self->{_saveoptions};
  undef $self->{_savefilelist};

  return $self
}

## method: save()
## save configuration data, options, parameters, properties
## and the file list if available and if values have been
## modified or added since the last time they were written.

sub save {
  my $self=shift;

  $self->SUPER::save();  
  $self->saveOptions();
  $self->saveProp();
  $self->saveFileList();
}

## method: saveProp()
## saves property values if changed since the last write

sub saveProp {
  my $self=shift;

  $self->writeProp() 
    if (defined $self->{_save});
}

## method: saveOptions()
## saves option/parameter values associated with the 
## current tag if changed since the last write

sub saveOptions {
  my $self=shift;
  
  $self->writeOptions()
    if (defined $self->{_saveoptions});
}

## method: saveFileList()
## saves the file list if changed since the last write

sub saveFileList {
  my $self=shift;

  $self->writeFileList()
    if (defined $self->{_savefilelist});
}

## method: readProp([file])
## reads property values from a file for the current tag.
## A default name is used if no file name is given.

sub readProp {
  my $self=shift;
  my $fname=shift;

  $fname="$self->{dir}/$self->{tag}".$PropFileSuffix
    if (!defined $fname);

  return if (!-r $fname);
  my $handle=&GenUtil::getInputFile($fname);

  my $count=0;
  my $maxruns=-1;

  my $header=<$handle>;
  chomp $header;
  my @ptag=();
  if ($header=~/^\#/) {
    @ptag=split(/ +/,$header);
    shift @ptag;
  }
  
  while (<$handle>) {
    if (!/^\#/ && !/^[ \t\n\r]+$/) {
      chomp;
      my @s=split(/[ \t]+/);
      $count++;
      $maxruns=$s[0];
      for (my $i=0; $i<=$#ptag; $i++) {
	$self->{prop}->{$ptag[$i]}->[$s[0]]=$s[$i+1]
	  if ($s[$i+1] ne "N/A");
      }
    }
  }

  close $handle;

  $self->set(runs=>$maxruns)
    if ($maxruns>$self->{par}->{runs});

  return $count;
}

## method: writeProp([file])
## writes property values to a file for the current tag.
## A default name is used if no file name is given.

sub writeProp {
  my $self=shift;
  my $fname=shift;

  $fname="$self->{dir}/$self->{tag}".$PropFileSuffix
    if (!defined $fname);

  return if (!defined $self->{prop});

  open OUT,">$fname";
  flock(OUT,LOCK_EX);

  my @ptag=();
  foreach my $k ( sort keys %{$self->{prop}} ) {
    push (@ptag,$k) if (defined $self->{prop}->{$k});
  }

  print OUT "# ",join(" ",@ptag),"\n";
  print OUT "# property file\n";
  my $now=localtime();
  printf OUT "# automatically generated on: %s\n",$now;

  for (my $i=1; $i<=$self->{par}->{runs}; $i++) {
    printf OUT "%d",$i;
    foreach my $t ( @ptag ) {
      if (defined $self->{prop}->{$t}->[$i]) {
	printf OUT " %f",$self->{prop}->{$t}->[$i];
      } else {
	print OUT " N/A";
      }
    }
    print OUT "\n";
  }

  flock(OUT,LOCK_UN);
  close OUT;

  undef $self->{_save};
}

## method: readOptions([file])
## reads options and parameters from a file for the current tag.
## A default name is used if no file name is given

sub readOptions {
  my $self=shift;
  my $fname=shift;
  
  $fname="$self->{dir}/$self->{tag}".$OptionFileSuffix
    if (!defined $fname);

  return if (!-r $fname);

  open OPTIN,"$fname";

  while (<OPTIN>) {
    if (!/^\#/ && !/^[ \t\n\r]+$/) {
      s/^[ \t]+//;
      my @f=split(/[ \t\n]+/);
      my $type=shift @f;
      my $tag=shift @f;
      my $value=join(" ",@f);

      if ($type eq "opt") {
	$self->setOption($tag => $value);
      } elsif ($type eq "par") {
	$self->setPar($tag => $value);
      } elsif ($type eq "par2") {
	$self->setPar2($tag => $value);
      } elsif ($type eq "par3") {
	$self->setPar3($tag => $value);
      }
    }
  }

  close OPTIN;
}

## method: writeOptions([file])
## writes options and parameters to a file for the current tag.
## A default name is used if no file name is given.

sub writeOptions {
  my $self=shift;
  my $fname=shift;

  return if (!defined $self->{opt} && !defined $self->{optpar} && !defined $self->{optpar2} && !defined $self->{optpar3});

  $fname="$self->{dir}/$self->{tag}".$OptionFileSuffix
    if (!defined $fname);
  
  open OPTOUT,">$fname";
  flock(OPTOUT,LOCK_EX);
  
  printf OPTOUT "# $self->{tag} options file\n";
  my $now=localtime();
  printf OPTOUT "# automatically generated on: %s\n",$now;
    
  foreach my $n ( sort keys %{$self->{opt}} ) {
    printf OPTOUT "opt $n %s\n",$self->{opt}->{$n};
  }

  foreach my $n ( sort keys %{$self->{optpar}} ) {
    printf OPTOUT "par $n %s\n",$self->{optpar}->{$n};
  }

  foreach my $n ( sort keys %{$self->{optpar2}} ) {
    printf OPTOUT "par2 $n %s\n",$self->{optpar2}->{$n};
  }
  
  foreach my $n ( sort keys %{$self->{optpar3}} ) {
    printf OPTOUT "par3 $n %s\n",$self->{optpar3}->{$n};
  }
  
  flock(OPTOUT,LOCK_UN);
  close OPTOUT;

  undef $self->{_saveoptions};
}

## method: readFileList([file])
## reads the file list from a file for the current tag.
## A default name is used if no file name is given.

sub readFileList {
  my $self=shift;
  my $fname=shift;
  
  $fname="$self->{dir}/ens.files"
    if (!defined $fname);

  return if (!-r $fname);

  $self->{filelist}=();

  open FINP,"$fname";
  while (<FINP>) {
    chomp;
    my ($n,$name)=split(/ +/);
    $self->{filelist}->[$n]=$name;
  }
  close FINP;
}

## method: writeFileList([file])
## writes the file list to a file for the current tag.
## A default name is used if no file name is given.

sub writeFileList {
  my $self=shift;
  my $fname=shift;
  
  $fname="$self->{dir}/ens.files"
    if (!defined $fname);

  return if (!defined $self->{filelist} || $#{$self->{filelist}}<0);

  open FOUT,">$fname";
  for (my $i=1; $i<=$#{$self->{filelist}}; $i++) {
    printf FOUT "%d %s\n",$i,$self->{filelist}->[$i];
  }
  close FOUT;
  undef $self->{_savefilelist};
}

## method: setPar(parameters)
## set parameters associated with the current tag.
## Parameters are distinguished from options to separate 
## simulation parameters for CHARMM and MONSSTER from other
## general options that may be set with <mark>setOption</mark>.
## (use <mark>set</mark> inherited from <docmark>SimData.pm</docmark>
## to set common ensemble parameters like the fragment list or the 
## native PDB structure).

sub setPar {
  my $self=shift;
  my %par=@_;

  foreach my $n ( keys %par ) {
    $self->{optpar}->{$n}=$par{$n} if (defined $par{$n});
  }
  $self->{_saveoptions}=1;
}

sub setPar2 {
  my $self=shift;
  my %par=@_;

  foreach my $n ( keys %par ) {
    $self->{optpar2}->{$n}=$par{$n} if (defined $par{$n});
  }
  $self->{_saveoptions}=1;
}

sub setPar3 {
  my $self=shift;
  my %par=@_;

  foreach my $n ( keys %par ) {
    $self->{optpar3}->{$n}=$par{$n} if (defined $par{$n});
  }
  $self->{_saveoptions}=1;
}

## method: setDefPar(parameters)
## set parameters that have not been defined yet from 
## default parameters given as arguments

sub setDefPar {
  my $self=shift;
  my %par=@_;

  foreach my $n ( keys %par ) {
    $self->{optpar}->{$n}=$par{$n} 
      if (defined $par{$n} && 
	  (!defined $self->{optpar} || !defined $self->{optpar}->{$n}));
  }
  $self->{_saveoptions}=1;
}

sub setDefPar2 {
  my $self=shift;
  my %par=@_;

  foreach my $n ( keys %par ) {
    $self->{optpar2}->{$n}=$par{$n} 
      if (defined $par{$n} && 
	  (!defined $self->{optpar2} || !defined $self->{optpar2}->{$n}));
  }
  $self->{_saveoptions}=1;
}

sub setDefPar3 {
  my $self=shift;
  my %par=@_;

  foreach my $n ( keys %par ) {
    $self->{optpar3}->{$n}=$par{$n} 
      if (defined $par{$n} && 
	  (!defined $self->{optpar3} || !defined $self->{optpar3}->{$n}));
  }
  $self->{_saveoptions}=1;
}

## method: $par = getPar()
## returns a hash reference to the list of parameters associated 
## with the current tag

sub getPar {
  my $self=shift;
  $self->{optpar}={} if (!defined $self->{optpar});
  return $self->{optpar};
}

sub getPar2 {
  my $self=shift;
  $self->{optpar2}={} if (!defined $self->{optpar2});
  return $self->{optpar2};
}

sub getPar3 {
  my $self=shift;
  $self->{optpar3}={} if (!defined $self->{optpar3});
  return $self->{optpar3};
}

## method: setOption(parameters)
## set options associated with the current tag

sub setOption {
  my $self=shift;
  my %par=@_;

  foreach my $n ( keys %par ) {
    $self->{opt}->{$n}=$par{$n} if (defined $par{$n});
  }
  $self->{_saveoptions}=1;
}

## method: setDefOption(parameters)
## sets default options associated with the current tag
## if not defined previously

sub setDefOption {
  my $self=shift;
  my %par=@_;

  foreach my $n ( keys %par ) {
    $self->{opt}->{$n}=$par{$n} 
      if (defined $par{$n} && (!defined $self->{opt} || !defined $self->{opt}->{$n}));
  }
  $self->{_saveoptions}=1;
}

## method: setFileList(index,name)
## set a file name for the given index. The file list
## is maintained to keep original external file names 
## that have been checked in with <docmark>checkin.pl</mark>
## or <docmark>evaluate.pl</mark>.

sub setFileList {
  my $self=shift;
  my $inx=shift;
  my $name=shift;

  $self->{filelist}=() 
    if (!defined $self->{filelist});
  $self->{filelist}->[$inx]=$name;
  $self->{_savefilelist}=1;
}

## method: getEnergy(index)
## parses the energy log file for the given index, that has been
## generated during a CHARMM minimization or simulation run.
## The final energies are used to update the ensemble property values 
## accordingly.

sub getEnergy {
  my $self=shift;
  my $inx=shift;

  my $d=$self->{dir}."/".&GenUtil::dataDir($inx);
  if (-d $d) {
    $self->set(runs=>$inx) if ($inx>$self->{par}->{runs});

    if (!defined $self->{prop}->{etot}->[$inx]) {
      my $elog=$self->{tag}.".elog";
      if (&GenUtil::checkFile("$d/$elog")) {
	my ($ener,$sasa)=CHARMM::readEnergy("$d/$elog");
	my $lastener=$ener->[$#{$ener}];

	$self->{prop}->{etot}->[$inx]=$lastener->{total};
	$self->{prop}->{eelec}->[$inx]=$lastener->{elec};
	$self->{prop}->{evdw}->[$inx]=$lastener->{vdwaals};
	$self->{prop}->{easp}->[$inx]=$lastener->{asp};
	$self->{prop}->{egb}->[$inx]=$lastener->{gb};
	$self->{prop}->{esasa}->[$inx]=$sasa;

      } elsif (&GenUtil::checkFile("$d/$self->{tag}.pdb")) {
	$self->{prop}->{etot}->[$inx]=0.0;
      }
    }
    $self->{_save}=1;
  }
}

## method: getStructInfo(index)
## analyzes the PDB structure for the given index and
## updates the corresponding property values accordingly.

sub getStructInfo {
  my $self=shift;
  my $inx=shift;

  my $d=$self->{dir}."/".&GenUtil::dataDir($inx);
  if (-d $d) {
    $self->set(runs=>$inx) if ($inx>$self->{par}->{runs});

    if (!defined $self->{prop}->{rgyr}->[$inx] ||
	!defined $self->{prop}->{rmsdpol}->[$inx] ||
	!defined $self->{prop}->{rmsdhyd}->[$inx] ||
	!defined $self->{prop}->{rmsdchg}->[$inx] ||
	(defined $self->getNatPDB() && 
	 (!defined $self->{prop}->{pphi}->[$inx] ||
	  !defined $self->{prop}->{pchi1}->[$inx] ||
	  !defined $self->{prop}->{cont}->[$inx] ||
	  !defined $self->{prop}->{gdtts}->[$inx] ||
	  !defined $self->{prop}->{rmsdall}->[$inx]))) {

      my $pdb=$self->{tag}.".pdb";    
      if (&GenUtil::checkFile("$d/$pdb")) {
	my $mol=Molecule::new("$d/$pdb");
	
	if (defined $mol) {
	  $mol->selectChain($self->{par}->{chain})
	    if (defined $self->{par}->{chain});
	  
	  my $dchain=$mol->activeChains()->[0];
	  
	  if (defined $dchain && $#{$dchain->{res}}>0) {
	    $self->{prop}->{rgyr}->[$inx]=
	      &Analyze::radiusOfGyration($mol)
		if (!defined $self->{prop}->{rgyr}->[$inx]);
	    
	    $self->{prop}->{rmsdpol}->[$inx]=
	      &Analyze::radialRMSD($mol,"pol")
		if (!defined $self->{prop}->{rmsdpol}->[$inx]);
	    
	    $self->{prop}->{rmsdhyd}->[$inx]=
	      &Analyze::radialRMSD($mol,"hyd")
		if (!defined $self->{prop}->{rmsdhyd}->[$inx]);
	    
	    $self->{prop}->{rmsdchg}->[$inx]=
	      &Analyze::radialRMSD($mol,"chg")
		if (!defined $self->{prop}->{rmsdchg}->[$inx]);
	    
	    if (defined $self->getNatPDB()) {
	      my $analyze=$self->getAnalyze();
	      my $exclmode=(defined $self->{par}->{exclmode})?$self->{par}->{exclmode}:1;
	      
	      die "cannot get analysis object" if (!defined $analyze);
	      
	      $mol->setValidResidues($self->getFragList(),0)
		if (defined $self->getFragList());
	      
	      if (!defined $self->{prop}->{pphi}->[$inx]) {
		my ($phi,$psi,$cphi,$cpsi)=$analyze->phipsiRMSD($mol);
		$self->{prop}->{pphi}->[$inx]=$cphi;
		$self->{prop}->{ppsi}->[$inx]=$cpsi;
		$self->{prop}->{rmsdphi}->[$inx]=$phi;
		$self->{prop}->{rmsdpsi}->[$inx]=$psi;
	      }
	      
	      if (!defined $self->{prop}->{pchi1}->[$inx]) {
		my ($chi1,$cchi1)=$analyze->chi1RMSD($mol);
		$self->{prop}->{pchi1}->[$inx]=$cchi1;
		$self->{prop}->{rmsdchi1}->[$inx]=$chi1;
	      }
	      
	      if (!defined $self->{prop}->{cont}->[$inx]) {
		my ($ncont,$pcont,$rho)=$analyze->contacts($mol);
		$self->{prop}->{cont}->[$inx]=$pcont;
		$self->{prop}->{rho}->[$inx]=$rho;
	      }
	  
	      if (!defined $self->{prop}->{gdtts}->[$inx]) {
#		$self->{prop}->{gdtts}->[$inx]=$analyze->gdtts($mol);
	      }
	      
	      if (!defined $self->{prop}->{rmsdall}->[$inx]) {
		$mol->setValidResidues($self->getFragList(),$exclmode)
		  if (defined $self->getFragList() && $exclmode);
		$analyze->lsqfit($mol,"cab",0);
		$mol->setValidResidues($self->getFragList(),0)
		  if (defined $self->getFragList() && $exclmode);
		my $rmsd=$analyze->rmsd($mol,0,($exclmode)?0:1);
		
		$self->{prop}->{rmsdall}->[$inx]=$rmsd->{all};
		$self->{prop}->{rmsdca}->[$inx]=$rmsd->{CA};
		$self->{prop}->{rmsdback}->[$inx]=$rmsd->{back};
	      }
	    }
	  }
	}
      }
    }
    $self->{_save}=1;
  }
}

## method: cleanUp([index])
## cleans up the ensemble directory for the given
## index or for all indices if no index is given.
## At this time 'clean up' means compressing 
## the PDB and energy log files.

sub cleanUp {
  my $self=shift;
  my $inx=shift;

  my $from=(defined $inx)?$inx:1;
  my $to=(defined $inx)?$inx:$self->{par}->{runs};

  for (my $i=$from; $i<=$to; $i++) {
    my $dir=$self->{dir}."/".&GenUtil::dataDir($i);
    &GenUtil::compress("$dir/$self->{tag}.elog") if ($self->{par}->{compress});
    &GenUtil::compress("$dir/$self->{tag}.pdb")  if ($self->{par}->{compress});
  }
}

## method: checkinPDB(index,mol[,ener])
## checks in an all-atom structure given as a Molecule object <mark>mol</mark>
## under the given index. An energy value for the total energy
## may be given if available to set the corresponding property
## value.

sub checkinPDB {
  my $self=shift;
  my $inx=shift;
  my $mol=shift;
  my $ener=shift;
  my $mode=shift;
  my $links=shift;
  
  $mode="CHARMM22" unless (defined $mode);

#  die "non continuous index $inx"
#    if ($inx>$self->{par}->{runs}+1);

  $self->set(runs=>$inx)
    if ($inx>$self->{par}->{runs});

  my $dir=$self->{dir}."/".&GenUtil::dataDir($inx);
  &GenUtil::makeDir($dir);
  if (defined $links && $links) {
    my $extension="pdb";
    if ($mol=~/.+\.([^\.]+)$/) {
      $extension=$1;
    }
    system("cd $dir; ln -s $mol $self->{tag}.$extension");
  } else {
    $mol->writePDB("$dir/$self->{tag}.pdb",translate=>$mode);
    &GenUtil::compress("$dir/$self->{tag}.pdb") if ($self->{par}->{compress});
  }
  $self->setProp("etot",$inx,$ener) if (defined $ener);
}

## method: update([index])
## updates missing energy and structural properties for the given
## index or all available indices if no index is given.
## Depending on the number and size of structures this may
## be relatively compute intensive.

sub update {
  my $self=shift;
  my $inx=shift;

  my $from=(defined $inx)?$inx:1;
  my $to=(defined $inx)?$inx:$self->{par}->{runs};

  for (my $i=$from; $i<=$to; $i++) {
    $self->getEnergy($i);
    $self->getStructInfo($i);
  }

  $self->{_save}=1;
}

## method: $list = jobList([from[,to[,propertytag]]])
## creates a job list from missing property values.
## By default the property <mark>etot</mark> is used
## but other ones can be specified. <mark>none</mark>
## will always give a list of all run indices.<BR>
## The <mark>from</mark> and <mark>to</mark> arguments
## may be used to limit the range of indices for
## which the job list is generated.

sub jobList {
  my $self=shift;
  my $from=shift;
  my $to=shift;
  my $ptag=shift;
  
  $ptag="etot" if (!defined $ptag);

  $from=1 if (!defined $from);
  $to=$self->{par}->{runs} if (!defined $to);

  my $list=();
  for (my $i=$from; $i<=$to; $i++) {
    push(@{$list},$i)
      if ($ptag eq "none" || 
	  !defined $self->{prop}->{$ptag} || 
	  !defined $self->{prop}->{$ptag}->[$i]);
  }
  return $list;
}

## method: $list = fileList([from[,to]])
## returns a list of PDB file names for ensemble 
## structures in the given index range

sub fileList {
  my $self=shift;
  my $from=shift;
  my $to=shift;

  my $list=();

  if (defined $from && ref $from) {
    foreach my $i ( @{$from} ) {
      my $t=$self->_addFile($i);
      push(@{$list},$t) 
	if (defined $t);
    }
  } else {
    $from=1 if (!defined $from);
    $to=$self->{par}->{runs} if (!defined $to);

    for (my $i=$from; $i<=$to; $i++) {
      my $t=$self->_addFile($i);
      push(@{$list},$t) 
	if (defined $t);
    }
  }
  return $list;
}  

sub _addFile {
  my $self=shift;
  my $inx=shift;

  my $d=$self->{dir}."/".&GenUtil::dataDir($inx);    
  my $fname="$d/$self->{tag}.pdb";
  if (&GenUtil::checkFile($fname,1)) {
    my $frec={};
    $frec->{name}=$fname;
    $frec->{inx}=$inx;
    return $frec;
  }
  return undef;
}

## method: $analyze = getAnalyze()
## creates an Analyze object from the native PDB structure
## set for the ensemble previously (if available).

sub getAnalyze {
  my $self=shift;
  
  if (!defined $self->{obj}->{analyze} && defined $self->getNatPDB()) {
    $self->{obj}->{analyze}=Analyze::new($self->getNatPDB());
  }

  return $self->{obj}->{analyze};
}

## method: $list = getPropTags()
## returns a list of all available property tags

sub getPropTags {
  my $self=shift;
  my @ptag=sort keys %{$self->{prop}};
  return \@ptag;
}

## method: $string = getPropString(index)
## returns all properties for a given run index in a single
## string. This is used for interprocess communication in 
## parallel jobs. The string representation is parsed
## back into the property data structure with <mark>setPropsFromString</mark>.

sub getPropString {
  my $self=shift;
  my $inx=shift;

  my @ps;
  foreach my $p ( sort keys %{$self->{prop}} ) {
    if (defined $self->{prop}->{$p}->[$inx]) {
      push(@ps,sprintf("%s=%s",$p,$self->{prop}->{$p}->[$inx]));
    }
  }
  return join(",",@ps);
}

## method: setPropsFromString(index,string)
## updates the property data structure for the given index 
## with values extracted from a string representation
## as generated with <mark>getPropString</mark>.

sub setPropsFromString {
  my $self=shift;
  my $inx=shift;
  my $str=shift;

  foreach my $s ( split(/,/,$str) ) {
    my @f=split(/=/,$s);
    $self->{prop}->{$f[0]}->[$inx]=$f[1]
      if ($f[1] ne "N/A");
  }

  $self->set(runs=>$inx)
    if ($inx>$self->{par}->{runs});

  $self->{_save}=1;
}

## method: removeProp(proptag) 
## removes all data for the given property

sub removeProp {
  my $self=shift;
  my $tag=shift;

  $self->{prop}->{$tag}=undef;

  $self->{_save}=1;
}

## method: setProp(proptag,index,value)
## sets an individual property value. The property tag,
## run index, and value are expected as arguments.

sub setProp {
  my $self=shift;
  my $tag=shift;
  my $inx=shift;
  my $val=shift;

  $self->{prop}->{$tag}=() 
    if (!defined $self->{prop}->{$tag});
  
  $self->{prop}->{$tag}->[$inx]=($val eq "N/A")?undef:$val
    if ($inx>0 && $inx<=$self->{par}->{runs});

  $self->{_save}=1;
}

## method: $list = getSortedList(tag,list)
## sorts the list given as argument according to
## values in the given property tag

sub getSortedList {
  my $self=shift;
  my $tag=shift;
  my $ilist=shift;
  my @xtags=@_;

  my $plist=$self->getPropList($tag,$ilist);
  
  for (my $i=0; $i<=$#{$plist}; $i++) {
    $ilist->[$i]->{val}=$plist->[$i]->{val};
  }

  if ($#xtags>=0) {
    foreach my $x ( @xtags ) {
      my $xlist=$self->getPropList($x,$ilist);
      for (my $i=0; $i<=$#{$xlist}; $i++) {
	$ilist->[$i]->{xtag}=() if (!exists $ilist->[$i]->{xtag});
	push(@{$ilist->[$i]->{xtag}},$xlist->[$i]->{val});
      }
    }
  }
  
  my @olist=sort { $a->{val} <=> $b->{val} } @{$ilist};
  return \@olist;
}

## method: $list = getList(tag,list)
## returns the list given as argument according to
## values in the given property tag

sub getList {
  my $self=shift;
  my $tag=shift;
  my $ilist=shift;
  my @xtags=@_;

  my $plist=$self->getPropList($tag,$ilist);
  
  for (my $i=0; $i<=$#{$plist}; $i++) {
    $ilist->[$i]->{val}=$plist->[$i]->{val};
  }

  if ($#xtags>=0) {
    foreach my $x ( @xtags ) {
      my $xlist=$self->getPropList($x,$ilist);
      for (my $i=0; $i<=$#{$xlist}; $i++) {
	$ilist->[$i]->{xtag}=() if (!exists $ilist->[$i]->{xtag});
	push(@{$ilist->[$i]->{xtag}},$xlist->[$i]->{val});
      }
    }
  }
  
  return $ilist;
}

## method: $list = getPropList(tag[,sellist])
## returns the list of values for a given property tag.
## A selection list may be given to request values
## only for a subset of all available runs.

sub getPropList {
  my $self=shift;
  my $tag=shift;
  my $sellist=shift;

  my $prop=();

  my @taglist=($tag=~/([+-]*)([^+-]+)/g);
  my @plist=();
  while (@taglist) {
    my $sgn=shift @taglist;
    my $tt=shift @taglist;
    my @ft=split(/\*/,$tt);
    my $t;
    my $fac;
    if ($#ft==0) {
	$t=$tt;
        $fac=1.0;
    } else {
        $fac=$ft[0];
        $t=$ft[1];
    }
    if (defined $self->{prop}->{$t}) {
      my $prec={};
      if ($sgn eq "-") {
	$prec->{sign}=-1;
      } else {
	$prec->{sign}=+1;
      }
      $prec->{factor}=$fac;

      $prec->{prop}=$self->{prop}->{$t};

      push(@plist,$prec);
    }
  }

  return undef 
    if ($#plist<0);

  if (!defined $sellist) {
    for (my $i=1; $i<=$self->{par}->{runs}; $i++) {
      my $sum=0;
      my $have=1;
      foreach my $p ( @plist) {
	if (defined $p->{prop}->[$i]) {
	  $sum+=$p->{sign}*$p->{prop}->[$i]*$p->{factor};
	} else {
	  $have=0;
	}
      }
      if ($have) {
	my $rec={};
	$rec->{inx}=$i;
	$rec->{val}=$sum;
	push(@{$prop},$rec);
      } else {
        my $rec={};
        $rec->{inx}=$i;
        $rec->{val}=undef;
        push(@{$prop},$rec);
      }
    }
  } else {
    foreach my $s ( @{$sellist} ) { 
      my $i=$s->{inx};
      my $sum=0;
      my $have=1;
      foreach my $p ( @plist) {
	if (defined $p->{prop}->[$i]) {
	  $sum+=$p->{sign}*$p->{prop}->[$i]*$p->{factor};
	} else {
	  $have=0;
	}
      }
      if ($have) {
	my $rec={};
	$rec->{inx}=$i;
	$rec->{val}=$sum;
	push(@{$prop},$rec);
      } else {
        my $rec={};
        $rec->{inx}=$i;
        $rec->{val}=undef;
        push(@{$prop},$rec);
      }
    }
  }    
  
  return $prop;
}

1;

