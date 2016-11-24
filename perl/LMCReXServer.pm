# LMCReXServer package
# provide server for lattice-based replica exchange sampling
# this class inherits from ReXServer
#
# http://mmtsb.scripps.edu/doc/LMCReXServer.pm.html
# 2001, Michael Feig, Brooks group, TSRI

package LMCReXServer;

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

sub setup {
  my $self=shift;
  my $condfile=shift;

  $self->SUPER::setup($condfile);

  die "need lattice temperature list" 
    if (!defined $self->{par}->{ltemplist});

  my @tarr=split(/:/,$self->{par}->{ltemplist});

  die "number of lattice temperatures does not match the number of replicas"
    if ($#tarr+1 != $self->nWindows());

  for (my $il=0; $il<=$#tarr; $il++) {
    $self->{cond}->[$il]->{ltemp}=$tarr[$il];
  }
}

sub giveClearance {
  my $self=shift;
  my $clientid=shift;

  my $c=$self->getClientData($clientid);
  my $s=sprintf("%s %d %d %f ltemp=%f",
		$self->getMode(),$self->getRun(),$self->{trun},
		$c->{cond}->{temp},$c->{cond}->{ltemp});

  $self->respondToClient($clientid,$s);
}

1;
