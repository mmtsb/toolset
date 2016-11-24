# Server package
# provide base server functionality
#
# http://mmtsb.scripps.edu/doc/Server.pm.html
# 2001, Michael Feig, Brooks group, TSRI

package Server;

require 5.004;

use strict;

use Fcntl;
use IO::File;
use IO::Handle;
use IO::Socket;
use Net::hostent;
use Sys::Hostname;

use GenUtil;
use SimData;
use MONSSTER;
use ReXServer;
use AmberReXServer;

use vars qw ( @ISA );

@ISA = ( "SimData" );

$SIG{PIPE}=sub{};

## data: serverPort
## server TCP/IP port address

## data: serverID
## server ID

## data: serverSocket
## server socket handle

## data: connected
## number of clients currently connected

## data: client -> { id } -> { socket host }
## information about connected clients

## data: html[] -> { client request alive }
## HTML clients expecting server push

## constructor: new([arguments])
## creates a new Server object. If arguments are given they
## are passed on to the parent class <docmark>SimData.pm</docmark>.

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;

  my $self=$class->SUPER::new(@_);

  $self->{client}={};
  $self->{connected}=0;

  return $self;
}

## method: (port,id,pid) = run(baseport[,presetid])
## starts the server as a separate process and
## returns the assigned TCP/IP port, internal ID, and 
## process ID. It is expected that the caller
## uses the process ID to wait for the server
## to finish. An address given as
## <mark>baseport</mark> argument can be used
## to set the beginning of the TCP/IP port range
## used for the server.

sub run {
  my $self=shift;
  my $baseport=shift;
  my $presetid=shift;

  $baseport=4000 if (!defined $baseport);

  &GenUtil::log("Server::run","starting");

  $self->{serverPort}=&GenUtil::getServerPort($baseport);
  $self->{serverID}=(defined $presetid)?$presetid:&GenUtil::getServerID();
  $self->{html}=();

  my $pid=fork();
  if (!$pid) {
    $self->{serverSocket}=&GenUtil::getServerSocket($self->{serverPort});

    $self->serverSetup();

    $self->{readselect}='';
    vec($self->{readselect},fileno($self->{serverSocket}),1)=1;

    my $quit;

    do {
      my $rout;

      my $nfound=select($rout=$self->{readselect},undef,undef,180);

      $quit=0;

      if (!$nfound) {
	foreach my $cid ( keys %{$self->{client}} ) {
	  if (!defined $self->{client}->{$cid}->{socket}) {
	    printf STDERR "Client $cid appears to be dead.\n";
	    #printf STDERR "Client $cid appears to be dead. Server is terminating.\n";
	    #$self->serverFinish();
	    #exit 0;
	  }
	}
      } else {
	if (vec($rout,fileno($self->{serverSocket}),1)) {
	  my $client=$self->{serverSocket}->accept();
	  fcntl($client,F_SETFL,O_NONBLOCK);
	  my $cmdline=&GenUtil::readFromSocket($client,"[\n\r]\$");
	  if (defined $cmdline) {
	    my ($cmd,@cmdargs)=split(/[ \r\n]+/,$cmdline);
	    $cmd=uc $cmd;

	    if ($cmd eq "TEST") {
	      &GenUtil::writeToSocket($client,"alive and kicking!\n");
	    } else {
	      my $htmlfile;
	      my $sid=shift @cmdargs;
	      if (defined $sid && $sid =~ /\/(.*)\?(.*)/) {
		$sid=$2;
		$htmlfile=$1;
	      }
	      if (defined $sid && $sid eq "$self->{serverID}") {
		$quit=$self->processCommand($client,$cmd,$htmlfile,@cmdargs);
	      } else {
		&GenUtil::writeToSocket($client,"ERROR: Invalid ID\n");
	      }
	    }
	  } else {
	    close $client;
	    &GenUtil::log("Server::run",
			  "socket connection lost while reading command line",1);
	  }
	}

	if (defined $self->{client}) {
	SERVICEDONE:
	  foreach my $cid ( keys %{$self->{client}} ) {
	    if (defined $self->{client}->{$cid}->{socket} && 
		vec($rout,fileno($self->{client}->{$cid}->{socket}),1)) {
	      $quit=$self->serviceClient($cid);
	    }
	    last SERVICEDONE if ($quit);
	  }
	}

	$quit=$self->serverTask()
	  if (!$quit);
      }
    } while (!$quit);

    $self->serverFinish();
    exit 0;
  } else {
    return ($self->{serverPort},$self->{serverID},$pid);
  }
}

sub processCommand {
  my $self=shift;
  my $client=shift;
  my $cmd=shift;
  my $htmlfile=shift;
  my @arg=@_;

  &GenUtil::log("Server::processCommand","$cmd ".join(" ",@arg));

  my $closeclient=1;
  if ($cmd eq "KILLIT") {
    &GenUtil::writeToSocket($client,"exiting ...\n");
    return 1;
  } elsif ($cmd eq "INIT") {
    my $id=shift @arg;
    $self->{client}->{$id}={};
    &GenUtil::writeToSocket($client,$self->initClient($id,@arg)."\n");    
  } elsif ($cmd eq "CONNECT") {
    if (++$self->{connected}>$self->maxClients()) {
      &GenUtil::writeToSocket($client,"ERROR: maximum number of clients reached\n");
    } else {
      my $id=shift @arg;
      if (!defined $self->{client}->{$id}) {
	&GenUtil::writeToSocket($client,"ERROR: need to initialize client first\n");
      } elsif (defined $self->{client}->{$id}->{socket}) {
	&GenUtil::writeToSocket($client,"ERROR: client $id has open socket already\n");
      } else {
	$self->{client}->{$id}->{socket}=$client;
	my $hostinfo = gethostbyaddr($client->peeraddr);
	$self->{client}->{$id}->{hostname}=$hostinfo->name || $client->peerhost;
	&GenUtil::log("Server::processCommand",
		      "client $id connected from $self->{client}->{$id}->{hostname}");
	vec($self->{readselect},fileno($client),1)=1;
	&GenUtil::writeToSocket($client,"$id connected\n");
	$closeclient=0;
      }
    }
  } elsif ($cmd eq ".GET") {
    my $localfile=$arg[0];
    my $data="";
    if (&GenUtil::safeToSend($localfile)) {
      my $lfile=&GenUtil::getInputFile($localfile);
      read $lfile,$data,1000000000;
      close $lfile;
    } else {
      &GenUtil::log("Server::processCommand",
		    "request for file $localfile refused",1)
	if (&GenUtil::checkFile($localfile));
      $data="#=EMPTY";
    } 
    &GenUtil::sendFixedData($client,$data);
  } elsif ($cmd eq "INFO") {
    &GenUtil::writeToSocket($client,$self->getInfo(uc $arg[0],$arg[1]));
  } elsif ($cmd eq "GET" && defined $htmlfile && $arg[0] =~ /^HTTP/ ) {
    my $httpheader=join(" ",@arg);
    $httpheader.=&GenUtil::readFromSocket($client,"\n *[\r]*\n")
      if ($httpheader !~ /User-Agent:/);
    my $explorer=($httpheader=~/MSIE/);

    &GenUtil::writeToSocket($client,"HTTP/1.0  200  OK\r\n");
    &GenUtil::writeToSocket($client,"Server: MMTSB\r\n");
    &GenUtil::writeToSocket($client,"Expires: now\r\n");

    if ($htmlfile=~/aa[0-9]+\/final\.jpeg/) {
      if (-r "$self->{dir}/$htmlfile") {
      my $buffer="";
      $buffer.="Content-type: image/jpeg\r\nExpires: now\r\n\r\n";
#$nowstring = strftime "%a, %e %b %Y %H:%M:%S -0700", localtime;
#my $header="Content-type: application/binary\nExpires: $nowstring\n\n";
      &GenUtil::writeToSocket($client,$buffer);

      open INP,"$self->{dir}/$htmlfile";
      $buffer="";
      my $nread=0;
      do {
        $nread=sysread(INP,$buffer,100000);
        &GenUtil::writeToSocket($client,$buffer)
          if ($nread>0);
      } while (defined $nread && $nread>0);
      }
    } elsif ($htmlfile eq "viz") {
      if ($explorer) {
	&GenUtil::writeToSocket($client,"Content-type: text/html\r\n\r\n");
        &GenUtil::writeToSocket($client,join("\r\n",@{$self->getHTMLPage("pullviz")})."\r\n");
      } else {
	&GenUtil::writeToSocket($client,"Content-type: multipart/x-mixed-replace;boundary=---981251251---\r\n");
	&GenUtil::writeToSocket($client,"\r\n---981251251---\r\n");
	
	my $rec={};
	$rec->{client}=$client;
	$rec->{request}=$htmlfile;
	$rec->{alive}=1;
	push(@{$self->{html}},$rec);
	$self->sendHTML($rec);
	$closeclient=0;
      }
    } elsif ($htmlfile eq "data") {
      if ($explorer) {
	&GenUtil::writeToSocket($client,"Content-type: text/html\r\n\r\n");
        &GenUtil::writeToSocket($client,join("\r\n",@{$self->getHTMLPage("pull")})."\r\n");
      } else {
	&GenUtil::writeToSocket($client,"Content-type: multipart/x-mixed-replace;boundary=---981251251---\r\n");
	&GenUtil::writeToSocket($client,"\r\n---981251251---\r\n");
	
	my $rec={};
	$rec->{client}=$client;
	$rec->{request}=$htmlfile;
	$rec->{alive}=1;
	push(@{$self->{html}},$rec);
	$self->sendHTML($rec);
	$closeclient=0;
      }
    } else {
      my $buffer="";
      $buffer.="Content-type: text/html\r\n\r\n";
      $buffer.="<HTML><HEAD>\r\n";
      $buffer.="<SCRIPT LANGUAGE=\"JavaScript\">\r\n";
      $buffer.=sprintf(" open(\"http://%s:%s/data?%s\",\"m%d\",\"width=450,height=%d,toolbar=no,location=no,directories=no,status=no,menubar=no,resizable=yes\");\r\n",
		       hostname,$self->{serverPort},$self->{serverID},int(rand(100000)),500);
      $buffer.="history.back();\r\n";
      $buffer.="</SCRIPT>\r\n";
      $buffer.="</HEAD><BODY></BODY></HTML>\r\n";
      &GenUtil::writeToSocket($client,$buffer);
    }
  } else {
    &GenUtil::writeToSocket($client,"ERROR: Command syntax error\n");
  }
  
  close $client if ($closeclient);
  
  return 0;
}

sub serviceClient {
  my $self=shift;
  my $clientid=shift;

  &GenUtil::log("Server::serviceClient","reading data from $clientid");

  my $line=&GenUtil::readFromSocket($self->{client}->{$clientid}->{socket},"[\r\n]\$");

  if (!defined $line || $line=~/^DONE/) {
    my $qt=$self->exitClient($clientid);
    if (defined $line) {
      &GenUtil::log("Server::serviceClient",
		    "client $clientid is done",0);
    } else {
      &GenUtil::log("Server::serviceClient",
		    "socket connection terminated at client $clientid",1);
    }

    vec($self->{readselect},fileno($self->{client}->{$clientid}->{socket}),1)=0;
    close $self->{client}->{$clientid}->{socket};
    $self->{client}->{$clientid}->{socket}=undef;
    return $qt;
  } elsif ($line=~/^DATA (.*)\n/) {
    &GenUtil::log("Server::serviceClient","$clientid sent data: $1");
    $self->processClientResponse($clientid,$1);
  } elsif ($line=~/^FILE (.*)\n/) {
    my $localfile=$1;
    my @farg=split(/\//,$localfile);
    if ($farg[$#farg] eq $MONSSTER::fileName{finalchain} ||
	$farg[$#farg] eq $MONSSTER::fileName{output} ||
	$farg[$#farg] eq $ReXServer::restartfile || 
	$farg[$#farg] eq $ReXServer::outpdb ||
	$farg[$#farg] eq $ReXServer::outcrd ||
	$farg[$#farg] eq $ReXServer::outphmd ||
	$farg[$#farg] eq $AmberReXServer::ambertrajout ||
	$farg[$#farg] eq $ReXServer::trajout ||
	$farg[$#farg] =~ /\.pdb$/ ||
	$farg[$#farg] =~ /\.crd$/ ||
	$farg[$#farg] =~ /\.lamb$/ ||
	$farg[$#farg] =~ /\.psf$/ ||
	$farg[$#farg] =~ /\.elog$/) {
      &GenUtil::writeToSocket($self->{client}->{$clientid}->{socket},"accepted\n");
      my $data=&GenUtil::recvFixedData($self->{client}->{$clientid}->{socket});
      my $out=&GenUtil::getOutputFile($localfile);
      print $out $data;
      close $out;
      &GenUtil::log("Server::serviceClient",
		    "$clientid sent file $localfile");
    } else {
      &GenUtil::log("Server::serviceClient",
		    "file $localfile not accepted by server",1);
      &GenUtil::writeToSocket($self->{client}->{$clientid}->{socket},"rejected\n");
    }
  } else {
    &GenUtil::log("Server::serviceClient",
		  "syntax error in response from client $clientid: $line",1);
    &GenUtil::writeToSocket($self->{client}->{$clientid}->{socket},"syntax error\n");
  }

  return 0;
}
  
sub serverTask {
  my $self=shift;
  return 0;
}

sub serverSetup {
  my $self=shift;
}

sub serverFinish {
  my $self=shift;

  foreach my $cid ( keys %{$self->{client}} ) {
    if (defined $self->{client}->{$cid}->{socket}) {
      &GenUtil::writeToSocket($self->{client}->{$cid}->{socket},"DONE\n");
      close $self->{client}->{$cid}->{socket};
    }
  }
  close $self->{serverSocket};
  &GenUtil::log("Server::finish","all done");
}

sub maxClients {
  my $self=shift;
  return 9999;
}

sub initClient {
  my $self=shift;
  my $clientid=shift;
  return "welcome";
}

sub exitClient {
  my $self=shift;
  my $clientid=shift;
  return 1;
}

sub getInfo {
  my $self=shift;
  my $name=shift;
  return "0\n";
}

sub processClientResponse {
  my $self=shift;
  my $clientid=shift;
  my $data=shift;
}

sub sendHTML {
  my $self=shift;
  my $rec=shift;

  if ($rec->{alive}) {
    my $client=$rec->{client};
    my $ret;
    if (($ret=&GenUtil::writeToSocket($client,"Content-type: text/html\r\n\r\n"))>0) {
      $ret=&GenUtil::writeToSocket($client,join("\r\n",@{$self->getHTMLPage($rec->{request})}).
				   "\r\n---981251251---\r\n");
    }
    $rec->{alive}=0 if ($ret<=0);
  }
}

sub getHTMLPage {
  my $self=shift;
  my $tag=shift;
}

sub pushHTML {
  my $self=shift;

  if (defined $self->{html} && $#{$self->{html}}>=0) {
    foreach my $hc ( @{$self->{html}} ) {
      $self->sendHTML($hc);
    }
  }
}

sub respondToClient {
  my $self=shift;
  my $clientid=shift;
  my $data=shift;

  &GenUtil::log("Server::respondToClient","data: $data");
  if (defined $self->{client}->{$clientid}->{socket}) {
    &GenUtil::writeToSocket($self->{client}->{$clientid}->{socket},$data."\n");
  } else {
    &GenUtil::log("Server::respondToClient",
		  "cannot respond to client $clientid",1);
  }
}
