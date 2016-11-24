# JobClient package
# connect to job server
#
# http://mmtsb.scripps.edu/doc/JobClient.pm.html
# 2001, Michael Feig, Brooks group, TSRI
#
# derived from Client

package JobClient;

require 5.004;

use strict;

use IO::Socket;
use Net::hostent;
use Sys::Hostname;

use GenUtil;
use Client;

use vars qw ( @ISA );

@ISA = ( "Client" );

## constructor: new(serverrec)
## creates a new JobClient object for connecting
## a job client to a job server. Server information
## in the <mark>serverrec</mark> argument is 
## passed to the constructor in <docmark>Client.pm</docmark>.

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $id=sprintf("%s%d",hostname,$$);
  return $class->SUPER::new($id,@_);
}

## method: makeDir(dir)
## requests a directory to be created unless
## it already exists

sub makeDir {
  my $self=shift;
  my $dir=shift;

  &GenUtil::sendCommand(
    &GenUtil::connectToServer($self->{serverName},$self->{serverPort}),
    ".MKDIR $self->{serverID} $dir");
}

1;
