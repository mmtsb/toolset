# ReXClient package
# connect to replica exchange server
#
# http://mmtsb.scripps.edu/doc/ReXClient.pm.html
# 2001, Michael Feig, Brooks group, TSRI
#
# derived from Client

package ReXClient;

require 5.004;

use strict;

use IO::Socket;
use Net::hostent;

use GenUtil;
use Client;

use vars qw ( @ISA );

@ISA = ( "Client" );

## constructor: new(clientid,serverrec)
## creates a new ReXClient object for connecting
## a replica exchange client to a replica exchange server

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  return $class->SUPER::new(@_);
}

## method: $ret = initFile() 
## obtains the file name for the initial structure from 
## the client

sub initFile {
  my $self=shift;
  
  my $ret=&GenUtil::sendCommand(
   &GenUtil::connectToServer($self->{serverName},$self->{serverPort}),
   "IFILE $self->{serverID} $self->{id}");

  die "error connecting to server $self->{serverName} at port $self->{serverPort}"
    unless (defined $ret);

  return $ret;
}

## method: $ret = biasInfo()
## obtains bias information from the server

sub biasInfo {
  my $self=shift;

  my $ret=&GenUtil::sendCommand(
   &GenUtil::connectToServer($self->{serverName},$self->{serverPort}),
   "BIAS $self->{serverID} $self->{id}");

  die "error connecting to server $self->{serverName} at port $self->{serverPort}"
    unless (defined $ret);

  chomp $ret;
  my $retrec=();
  my @f=split(/[ \t]+/,$ret);
  foreach my $tf (@f) {
    my @sf=split(/[=,]/,$tf);
    my $trec={};
    %{$trec}=@sf;
    push(@{$retrec},$trec);
  }
  return $retrec;
}

## method: $ret = nextCycle(filelist,values)
## sends the energy, etc. from the last simulation cycle
## to the server and then waits for a new temperature
## to be sent back for running the next cycle.

sub nextCycle {
  my $self=shift;
  my $flist=shift;
  my %par=@_;

  $par{ener}=0.0 if (!defined $par{ener});

  my $buffer="";
  foreach my $kpar ( keys %par ) {
    $buffer.="$kpar $par{$kpar} "
      if (defined $par{$kpar});
  }

  my $ret=$self->nextJob($buffer,$flist);

  if (defined $ret) {
    chomp $ret;
    my $retrec={};
    my @f=split(/[ \t]+/,$ret);
    $retrec->{mode}=uc shift @f;
    $retrec->{run}=shift @f;
    $retrec->{trun}=shift @f;
    $retrec->{temp}=shift @f;
    $retrec->{bias}=();
    foreach my $tf (@f) {
      my @sf=split(/[=,]/,$tf);
      my $trec={};
      %{$trec}=@sf;
      push(@{$retrec->{bias}},$trec);
    }
    return $retrec;
  } else {
    return undef;
  }
}



1;
