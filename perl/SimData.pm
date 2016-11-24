# SimData package
# keep common data for lattice simulation/minimization ensembles
#
# http://mmtsb.scripps.edu/doc/SimData.pm.html
# 2000, Michael Feig, Brooks group, TSRI

package SimData;

require 5.004;

use strict;

use GenUtil;
use Molecule;
use MONSSTER;

use Fcntl ':flock';

## data: dir
## ensemble working directory

## data: par -> { fraglist fragref hsdlist hselist 
## data:          seq natpdb }
## parameters in text string form that are used to create
## corresponding objects when requested

## data: obj -> { fraglist fragref seq natpdb }
## object references for objects created from the corresponding 
## parameter values

## constructor: new(dir)
## creates a new SimData object with the working directory
## given as argument. 
## The default configuration file is read automatically
## if present.

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;

  my $self={};
  bless($self,$class);

  my $dir=shift;

  $dir="." if (!defined $dir);
  &GenUtil::makeDir($dir) if (!-d $dir);
  $self->{dir}=$dir;

  $self->{configfile}=shift;
  $self->{configfile}="ens.cfg"
    if (!defined $self->{configfile});

  $self->readConfig();

  undef $self->{_saveconfig};

  return $self;
}

## method: set(parameters)
## sets new parameters or modifies existing ones

sub set {
  my $self=shift;
    
  my %par=@_;

  foreach my $n ( keys %par ) {
    if (defined $par{$n}) {
      $self->{par}->{$n}=$par{$n};

      undef $self->{obj}->{fraglist} if ($n eq "fraglist");
      undef $self->{obj}->{fragref}  if ($n eq "fragref");
      undef $self->{obj}->{fragref}  if ($n eq "hsdlist");
      undef $self->{obj}->{fragref}  if ($n eq "hselist");
      undef $self->{obj}->{seq}      if ($n eq "seq");
      undef $self->{obj}->{natpdb}   if ($n eq "natpdb");

      $self->{_saveconfig}=1;
    }
  }
}

## method: save()
## saves all parameters to the default configuration file
## if new parameters have been set or changed

sub save {
  my $self=shift;

  $self->writeConfig()
    if (defined $self->{_saveconfig});
}

## method: readConfig([file])
## reads a configuration file. If no argument is
## given it reads the default configuration file.

sub readConfig {
  my $self=shift;
  my $fname=shift;
  
  $fname="$self->{dir}/$self->{configfile}" 
    if (!defined $fname);

  return if (!-r $fname);

  open CFGIN,"$fname";

  while (<CFGIN>) {
    if (!/^\#/ && !/^[ \t\n\r]+$/) {
      s/^[ \t]+//;
      my @f=split(/[ \t\n]+/);
      $self->set($f[0] => $f[1]);
    }
  }

  close CFGIN;
}

## method: writeConfig([file])
## writes all parameters to a configuration file. If no
## argument is given the parameters are written to the
## default configuration file

sub writeConfig {
  my $self=shift;
  my $fname=shift;
  
  $fname="$self->{dir}/$self->{configfile}"
    if (!defined $fname);

  open CFGOUT,">$fname";
  flock(CFGOUT,LOCK_EX);

  printf CFGOUT "# config file\n";
  my $now=localtime();
  printf CFGOUT "# automatically generated on: %s\n",$now;
  
  foreach my $n ( sort keys %{$self->{par}} ) {
    printf CFGOUT "$n %s\n",$self->{par}->{$n};
  }

  flock(CFGOUT,LOCK_UN);
  close CFGOUT;

  undef $self->{_saveconfig};
}

## method: $fraglist = getFragList()
## returns a fragment list object generated
## from the fragment list parameter <mark>fraglist</mark> in string
## format

sub getFragList {
  my $self=shift;

  if (!defined $self->{obj}->{fraglist} && 
      defined $self->{par}->{fraglist}) {
    $self->{obj}->{fraglist}=
      &GenUtil::fragListFromOption($self->{par}->{fraglist});
  }
  
  return $self->{obj}->{fraglist};
}

## method: $mol = getFragRef()
## returns a Molecule object created from the <mark>fragref</mark> 
## parameter containing
## a PDB file name for the reference template structure used in 
## fragment/loop modeling

sub getFragRef {
  my $self=shift;

  if (!defined $self->{obj}->{fragref} && 
      &GenUtil::checkFile($self->{par}->{fragref})) {
    $self->{obj}->{fragref}=
      Molecule::new($self->{par}->{fragref});
    $self->{obj}->{fragref}->
      fixHistidine($self->{par}->{hsdlist},$self->{par}->{hselist});
  }

  return $self->{obj}->{fragref};
}  

## method: $seq = getSeq()
## returns a Sequence object created from the parameter <mark>seq</mark>
## describing a MONSSTER sequence file

sub getSeq {
  my $self=shift;

  if (!defined $self->{obj}->{seq} &&
      &GenUtil::checkFile($self->{par}->{seq})) {
    $self->{obj}->{seq}=Sequence::new();
    $self->{obj}->{seq}->readMONSSTER($self->{par}->{seq});
  }

  return $self->{obj}->{seq};
}

## method: $mol = getNatPDB()
## returns a Molecule object created from the <mark>natpdb</mark> parameter 
## containing a file name for the native structure.

sub getNatPDB {
  my $self=shift;
  
  if (!defined $self->{obj}->{natpdb} &&
      &GenUtil::checkFile($self->{par}->{natpdb})) {
    $self->{obj}->{natpdb}=
      Molecule::new($self->{par}->{natpdb});
  }

  return $self->{obj}->{natpdb};
}

1;
