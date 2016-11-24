# Client package
# basic client functionality for connecting to server
#
# http://mmtsb.scripps.edu/doc/Client.pm.html
# 2001, Michael Feig, Brooks group, TSRI

package Client;

require 5.004;

use strict;

use IO::Socket;
use Net::hostent;

use GenUtil;

## constructor: new(clientid,serverrec)
## creates a new Client object for connecting
## to a server from <docmark>Server.pm</docmark>. 
## Required arguments are a unique client ID and
## a reference to a data structure with the following
## elements: <mark>name</mark> -> server host name, 
## <mark>port</mark> -> server 
## TCP/IP port, and <mark>id</mark> -> server ID.

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;

  my $self={};

  $self->{id}=shift;

  my $srec=shift;
  $self->{serverName}=$srec->{name};
  $self->{serverPort}=$srec->{port};
  $self->{serverID}=$srec->{id};
  $self->{connection}=undef;

  bless($self,$class);
  return $self;
}

## method: ret = initialize(hostid,maxcpus)
## registers the client with the server. A host ID
## and maximum number of CPUs on the given host
## are required as arguments for client scheduling
## by the server. The method returns data sent
## back from the server that may be used for
## restarting a previous simulation

sub initialize {
  my $self=shift;
  my $hostid=shift;
  my $maxcpus=shift;

  my $cmd="INIT $self->{serverID} $self->{id}";
  $cmd.=" $hostid" if (defined $hostid);
  $cmd.=" $maxcpus" if (defined $maxcpus);

  my $ret=&GenUtil::sendCommand(&GenUtil::connectToServer($self->{serverName},$self->{serverPort}),$cmd,0);
  die "error initializing connection to server $self->{serverName} at port $self->{serverPort}"
    unless (defined $ret);
  return ($ret=~/^ +/)?undef:$ret;
}

## method: establishConnection()
## establishes a persistent socket connection with the server.
## <mark>initialize</mark> needs to be called first before
## this method can be used.

sub establishConnection {
  my $self=shift;
  
  $self->{connection}=
    &GenUtil::connectToServer($self->{serverName},$self->{serverPort});

  my $cmd="CONNECT $self->{serverID} $self->{id}";

  my $ret=&GenUtil::sendCommand($self->{connection},$cmd,1);
  die "error establishing connection to server $self->{serverName} at port $self->{serverPort}"
    unless (defined $ret && $ret=~/connected/);
}

## method: $ret = nextJob(data,filelist)
## sends data and files to the server and waits
## for information about the next job to be
## run by the client. The return value is 
## undefined if no more jobs are available.
## This method requires that
## a persistent socket connection with the server
## is established with <mark>establishConnection</mark>.

sub nextJob {
  my $self=shift;
  my $data=shift;
  my $flist=shift;

  die "no connection available"
    if (!defined $self->{connection});

  if (defined $flist) {
    foreach my $f ( @{$flist} ) {
      if (&GenUtil::checkFile($f->{local})) {
	my $buffer='';
	&GenUtil::writeToSocket($self->{connection},"FILE $f->{remote}\n");
	my $ret=&GenUtil::readFromSocket($self->{connection},"[\r\n]\$");
	if (defined $ret && $ret =~ /accepted/) {
	  my $inp=&GenUtil::getInputFile($f->{local});
	  read $inp,$buffer,1000000000;
	  close $inp;
	  &GenUtil::sendFixedData($self->{connection},$buffer);
	}
      }
    }
  }

  &GenUtil::writeToSocket($self->{connection},"DATA $data\n");

  my $ret=&GenUtil::readFromSocket($self->{connection},"[\r\n]\$");
  chomp $ret;
  $ret=~s/^ +//g;
  
  if ($ret=~/^DONE/) {
    close $self->{connection};
    $self->{connection}=undef;
    return undef;
  } 

  return $ret;
}

## method: terminateServer()
## sends the <mark>KILLIT</mark> command to terminate the server
## Under normal circumstances the server should terminate
## on its own and close client connections once all jobs are done.
## This function is intended to handle client-side exceptions that 
## might leave the server running indefinitely.

sub terminateServer {
  my $self=shift;

  my $ret=&GenUtil::sendCommand(
   &GenUtil::connectToServer($self->{serverName},$self->{serverPort}),
   "KILLIT $self->{serverID}");
}

## method: finish()
## closes the persistent socket connection with the server

sub finish {
  my $self=shift;
  
  die "no connection available"
    if (!defined $self->{connection});

  &GenUtil::writeToSocket($self->{connection},"DONE\n");
  sleep 1;
  close $self->{connection};
  $self->{connection}=undef;
}

## method: getFile(remotefile[, localfile])
## retrieves a file from the server. The file names for
## the remote and local files are expected as arguments.

sub getFile {
  my $self=shift;
  my $remote=shift;
  my $locfile=shift;

  $locfile=$remote if (!defined $locfile);

  my $sock=&GenUtil::connectToServer($self->{serverName},$self->{serverPort});
  die "cannot connect to server"
    if (!defined $sock);

  &GenUtil::writeToSocket($sock,".GET $self->{serverID} $remote\n");
  my $data=&GenUtil::recvFixedData($sock);
  close $sock;

  if (defined $data) {
    if ($data!~/^\#=EMPTY/) {
      my $lfile=&GenUtil::getOutputFile($locfile);
      print $lfile $data;
      close $lfile;
    }
  } else {
    printf STDERR "ERROR: socket connection lost while receiving file $locfile\n";
  }
}

## method: $ret = getInfo(name)
## obtains status information from the server according
## to the given argument.

sub getInfo {
  my $self=shift;
  my $name=shift;

  my $ret=&GenUtil::sendCommand(
     &GenUtil::connectToServer($self->{serverName},$self->{serverPort}),
     "INFO $self->{serverID} $name");

  return $ret;
}

1;
