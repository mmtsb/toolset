# NAMDReXServer package
# provide server for NAMD replica exchange sampling (requires custom NAMD executable)
# this class inherits from ReXServer
#
# 2009, Michael Feig, MSU

package NAMDReXServer;

require 5.004;

use lib "$ENV{MMTSBDIR}/perl";
use strict;

use ReXServer;
use GenUtil;
use Molecule;

use vars qw (@ISA $namdout $namdrest $namdcoorout $namdvelout $namdxscout $lastnamdcoorout $lastnamdvelout $lastnamdxscout );

@ISA=( "ReXServer" );

BEGIN {
  $namdout="final";
  $namdrest="last";
  $namdcoorout=$namdout.".coor";
  $namdvelout=$namdout.".vel";
  $namdxscout=$namdout.".xsc";
  $lastnamdcoorout=$namdrest.".coor";
  $lastnamdvelout=$namdrest.".vel";
  $lastnamdxscout=$namdrest.".xsc";
}

## data: par -> { pdb }
## flag for rebuilding all-atom structures

sub setupPar {
  my $self=shift;

  $self->SUPER::setupPar(@_);

  die "need pdb" if (!&GenUtil::checkFile($self->{par}->{pdb}));
  die "need psf" if (!&GenUtil::checkFile($self->{par}->{psf}));

  $self->{par}->{ensmode}="add"
    if (!defined $self->{par}->{ensmode});
}

sub postProcess {
  my $self=shift;
  my $clientid=shift;
  my $d="$self->{dir}/$clientid";

  if (-r "$d/$namdcoorout") {
    system "cp $d/$namdcoorout $d/$lastnamdcoorout";

    my $mol=Molecule::new();
    $mol->readPDB($self->{par}->{pdb});

    my $coorfile=&GenUtil::getInputFile("$d/$namdcoorout");
    binmode $coorfile;
    
    my $buffer="";
    read($coorfile,$buffer,4);
    my $natoms=unpack("L",$buffer);
    
    my $nn=0;
    foreach my $c ( @{$mol->{chain}}) {
      foreach my $a ( @{$c->{atom}}) {
	read($coorfile,$buffer,24);
	my ($xval,$yval,$zval)=unpack("d*",$buffer);
	$a->{xcoor}=$xval;
	$a->{ycoor}=$yval;
	$a->{zcoor}=$zval;
	$nn++;
      }
    }
    if ($nn != $natoms) {
      printf STDERR "%d atoms read from coordinate file, %d atoms expected from header\n",$nn,$natoms;
    }
    
    close $coorfile;
    undef $coorfile;

    $mol->writePDB("$d/$ReXServer::outpdb",translate=>"CHARMM22");
  }

  if (-r "$d/$namdvelout") {
    system "cp $d/$namdvelout $d/$lastnamdvelout";
  }

  if (-r "$d/$namdxscout") {
    system "cp $d/$namdxscout $d/$lastnamdxscout";
  }

  if ($self->{par}->{save} && $self->{trun}>0) {
    if ($self->{par}->{archive}) {
      if (-r "$d/$namdcoorout") {
	&GenUtil::archiveFile("$d/".(lc $self->getMode()).".ncoor.archive",
			      "$d/$namdcoorout",$self->getRun());
	system "rm $d/$namdcoorout";
      }
      if (-r "$d/$namdvelout") {
	&GenUtil::archiveFile("$d/".(lc $self->getMode()).".vel.archive",
			      "$d/$namdvelout",$self->getRun());
	system "rm $d/$namdvelout";
      }
      if (-r "$d/$namdxscout") {
	open TOUT,">$d/$namdxscout.reg";
        my $inp=&GenUtil::getInputFile("$d/$namdxscout");
        while (<$inp>) {
	  if (/^#/) {
	    print TOUT;
          } else {
	    chomp;
	    my $buf=$_;
	    my $add="";
	    for (my $i=0; $i<120-length($buf); $i++) {
	      $add.=" ";
	    }
	    print TOUT $buf.$add."\n";
	  }
        }
	undef $inp;
	close TOUT;

	&GenUtil::archiveFile("$d/".(lc $self->getMode()).".xsc.archive",
			      "$d/$namdxscout.reg",$self->getRun());
	system "rm $d/$namdxscout";
	system "rm $d/$namdxscout.reg";
      }
    } else {
      my $tag=(lc $self->getMode())."/".&GenUtil::dataDir($self->getRun());
      my $savedir="$d/$tag";
      &GenUtil::makeDir($savedir);
 
      if (-r "$d/$namdcoorout") {
	system "mv $d/$namdcoorout $savedir";
	&GenUtil::compress("$savedir/$namdcoorout");
      }
      if (-r "$d/$namdvelout") {
	system "mv $d/$namdvelout $savedir";
	&GenUtil::compress("$savedir/$namdvelout");
      }
      if (-r "$d/$namdxscout") {
	system "mv $d/$namdxscout $savedir";
	&GenUtil::compress("$savedir/$namdxscout");
      }
    }
  }


  $self->SUPER::postProcess($clientid);
}

1;
