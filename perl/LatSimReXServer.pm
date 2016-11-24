# LatSimReXServer package
# provide server for lattice-based replica exchange sampling
# this class inherits from ReXServer
#
# http://mmtsb.scripps.edu/doc/LatSimReXServer.pm.html
# 2001, Michael Feig, Brooks group, TSRI

package LatSimReXServer;

require 5.004;

use lib "$ENV{MMTSBDIR}/perl";
use strict;

use ReXServer;
use SICHO;
use GenUtil;
use Molecule;
use MONSSTER;

use vars qw (@ISA);

@ISA=( "ReXServer" );

## data: par -> { rebuild }
## flag for rebuilding all-atom structures

sub setupPar {
  my $self=shift;

  $self->SUPER::setupPar(@_);

  $self->{par}->{rebuild}=0
    if (!defined $self->{par}->{rebuild});

  $self->{par}->{ensmode}="add"
    if (!defined $self->{par}->{ensmode});

  $self->{_boltzFactor}=1.0;
}

sub ensembleCheckin {
  my $self=shift;
  my $clientid=shift;
  my $at=shift;

  my $chain=SICHO::new();
  $chain->readChain("$self->{dir}/$clientid/$MONSSTER::fileName{finalchain}");
  $self->{ens}->checkinSICHO($at,$chain,$self->getClientData($clientid)->{ener});
}

sub postProcess {
  my $self=shift;
  my $clientid=shift;
  my $d="$self->{dir}/$clientid";

  if ($self->{par}->{save} && $self->{trun}>0) {
    if ($self->{par}->{archive}) {
      &GenUtil::archiveFile("$d/".(lc $self->getMode()).".finalchain.archive",
			    "$d/$MONSSTER::fileName{finalchain}",$self->getRun());
    } else {
      my $tag=(lc $self->getMode())."/".&GenUtil::dataDir($self->getRun());
      my $savedir="$d/$tag";
      &GenUtil::makeDir($savedir);
      system "cp $d/$MONSSTER::fileName{finalchain} $savedir" 
	if (-r "$d/$MONSSTER::fileName{finalchain}");
      &GenUtil::compress("$savedir/$MONSSTER::fileName{finalchain}");
    }

    if ($self->{par}->{rebuild}) {
      my $chain=SICHO::new();
      $chain->readChain("$d/$MONSSTER::fileName{finalchain}");

      if (!defined $self->getSeq()) {
	$self->set(seq=>"$d/$MONSSTER::fileName{sequence}");
      }

      my $rebmol=Molecule::new();
      $rebmol->rebuildFromSICHO($self->getSeq(),$chain,
				$self->{par}->{fraglist},$self->{par}->{fragref});

      $rebmol->writePDB("$d/$ReXServer::outpdb",translate=>"CHARMM22");
    }
  }
  $self->SUPER::postProcess($clientid);
}

1;
