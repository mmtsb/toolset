# AmberReXServer package
# provide server for lattice-based replica exchange sampling
# this class inherits from ReXServer
#
# http://mmtsb.scripps.edu/doc/AmberReXServer.pm.html
# 2003, Michael Feig, Brooks group, TSRI

package AmberReXServer;

require 5.004;

use lib "$ENV{MMTSBDIR}/perl";
use strict;

use ReXServer;
use GenUtil;
use Molecule;

use vars qw (@ISA $ambertrajout );

@ISA=( "ReXServer" );

BEGIN {
  $ambertrajout="traj.x";
}

## data: par -> { partop }
## flag for rebuilding all-atom structures

sub setupPar {
  my $self=shift;

  $self->SUPER::setupPar(@_);

  die "need partop file" if (!&GenUtil::checkFile($self->{par}->{partop}));

  $self->{par}->{ensmode}="add"
    if (!defined $self->{par}->{ensmode});
}

sub postProcess {
  my $self=shift;
  my $clientid=shift;
  my $d="$self->{dir}/$clientid";

  if ($self->{par}->{save} && $self->{trun}>0) {
    if ($self->{par}->{archive}) {
      if (-r "$d/$ambertrajout") {
	&GenUtil::archiveFile("$d/".(lc $self->getMode()).".trajx.archive",
			      "$d/$ambertrajout",$self->getRun());
      }
    } else {
      my $tag=(lc $self->getMode())."/".&GenUtil::dataDir($self->getRun());
      my $savedir="$d/$tag";
      &GenUtil::makeDir($savedir);
 
      if (-r "$d/$ambertrajout") {
	system "mv $d/$ambertrajout $savedir";
	&GenUtil::compress("$savedir/$ambertrajout");
      }
    }
  }
  my $mol=Molecule::new();
  $mol->readAmber($self->{par}->{partop},"$d/$ReXServer::restartfile");
  $mol->writePDB("$d/$ReXServer::outpdb",translate=>"CHARMM22");

  $self->SUPER::postProcess($clientid);
}

1;
