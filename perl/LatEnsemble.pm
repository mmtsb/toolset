# LatEnsemble package
# manage ensemble properties for lattice simulation
#
# http://mmtsb.scripps.edu/doc/LatEnsemble.pm.html
# 2000, Michael Feig, Brooks group, TSRI
#
# derived from Ensemble package

package LatEnsemble;

require 5.004;

use strict;

use Fcntl ':flock';

use Ensemble;

use GenUtil;
use Molecule;
use Analyze;
use CHARMM;
use MONSSTER;

use vars qw ( @ISA );

@ISA = ( "Ensemble" );

## constructor: new(tag,dir)
## creates a new LatEnsemble object. The required tag and directory 
## arguments are passed on to the constructor of the Ensemble package.

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $self=$class->SUPER::new(@_);
  return $self
}

## method: getEnergy(index)
## replaces the <mark>getEnergy</mark> method in
## <docmark>Ensembe.pm</docmark> to extract energy
## information from MONSSTER output.

sub getEnergy {
  my $self=shift;
  my $inx=shift;

  my $d=$self->{dir}."/".&GenUtil::dataDir($inx);

  if (-d $d) {
    $self->set(runs=>$inx) if ($inx>$self->{par}->{runs});

    my $monsster=MONSSTER::new();
    $monsster->setDirectory($d);

    if (!defined $self->{prop}->{etot}->[$inx]) {
      if (&GenUtil::checkFile("$d/$MONSSTER::fileName{output}")) {
	$self->{prop}->{etot}->[$inx]=$monsster->finalEnergy()->{total};
      } elsif (&GenUtil::checkFile("$d/$self->{tag}.pdb")) {
	$self->{prop}->{etot}->[$inx]=0.0;
      }
    }
    $self->{_save}=1;
  }
}

## method: checkinSICHO(index,chain,energy)
## checks in a SICHO lattice chain under the given index and rebuilds
## automatically an all-atom structure. An energy value for the total energy
## may be given if available to set the corresponding property
## value.

sub checkinSICHO {
  my $self=shift;
  my $inx=shift;
  my $sicho=shift;
  my $ener=shift;
  
  die "non continuous index $inx"
    if ($inx>$self->{par}->{runs}+1);

  $self->set(runs=>$inx)
    if ($inx>$self->{par}->{runs});

  my $dir=$self->{dir}."/".&GenUtil::dataDir($inx);

  &GenUtil::makeDir($dir);

  $sicho->writeChain("$dir/$MONSSTER::fileName{finalchain}");
  $self->setProp("etot",$inx,$ener) if (defined $ener);

  my $rebmol=Molecule::new();
  $rebmol->rebuildFromSICHO($self->getSeq(),$sicho);
#			    $self->{par}->{fraglist},$self->{par}->{fragref});
  $rebmol->writePDB("$dir/$self->{tag}.pdb",translate=>"CHARMM22");

  &GenUtil::compress("$dir/$MONSSTER::fileName{finalchain}") if ($self->{par}->{compress});  
  &GenUtil::compress("$dir/$self->{tag}.pdb") if ($self->{par}->{compress});
}

## method: cleanUp([index])
## replaces the <mark>cleanUp</mark> in <docmark>Ensemble.pm</docmark> 
## to remove and compress MONSSTER files in the ensemble directory for the given
## index or for all indices if no index is given.

sub cleanUp {
  my $self=shift;
  my $inx=shift;
  my $keepTra=shift;

  my $from=(defined $inx)?$inx:1;
  my $to=(defined $inx)?$inx:$self->{par}->{runs};

  for (my $i=$from; $i<=$to; $i++) {
    my $dir=$self->{dir}."/".&GenUtil::dataDir($i);

    &GenUtil::remove("$dir/$MONSSTER::fileName{datadir}");
    &GenUtil::remove("$dir/$MONSSTER::fileName{sequence}");
    &GenUtil::remove("$dir/$MONSSTER::fileName{trajectory}")
      if (!defined $keepTra || !$keepTra);
    &GenUtil::remove("$dir/$MONSSTER::fileName{restraints}");
    &GenUtil::remove("$dir/$MONSSTER::fileName{initchain}");
    &GenUtil::remove("$dir/$MONSSTER::fileName{input}");

    &GenUtil::compress("$dir/$MONSSTER::fileName{output}");
    &GenUtil::compress("$dir/$MONSSTER::fileName{finalchain}");
    &GenUtil::compress("$dir/$MONSSTER::fileName{trajectory}");
  }
  $self->SUPER::cleanUp($inx);
}

1;

