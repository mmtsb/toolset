# JobServer package
# provide server for parallel ensemble jobs 
#
# http://mmtsb.scripps.edu/doc/JobServer.pm.html
# 2000, Michael Feig, Brooks group, TSRI
#
# derived from Server

package JobServer;

require 5.004;

use strict;

use IO::Socket;
use Net::hostent;

use GenUtil;
use Server;
use Ensemble;
use MONSSTER;

use vars qw ( @ISA );

@ISA = ( "Server" );

## constructor: new(joblist,ensemble[,updatefrq])
## creates a new JobServer object. The job list is given in 
## <mark>joblist</mark>.
## An Ensemble object <mark>ensemble</mark> is used to update
## properties for completed jobs. The update frequency can be given
## by <mark>updatefrq</mark>.

sub new {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $joblist=shift;
  my $ens=shift;
  my $updatefrq=shift;

  my $self=$class->SUPER::new();

  $self->{joblist}=();
  if (defined $joblist) {
    for (my $i=0; $i<=$#{$joblist}; $i++) {
      my $newrec={};
      $newrec->{job}=$joblist->[$i];
      $newrec->{status}="ready";
      $newrec->{client}=undef;
      push (@{$self->{joblist}},$newrec);
    }
  }

  $self->{status}->{activeclients}=0;
  $self->{status}->{finishedclients}=0;
  $self->{status}->{jobsready}=$#{$joblist}+1;
  $self->{status}->{jobstodo}=$#{$joblist}+1;
  $self->{status}->{jobsrunning}=0;
  $self->{status}->{jobsdone}=0;

  $self->{updatefrq}=(defined $updatefrq)?$updatefrq:20;
  $self->{ensemble}=$ens;
  
  return $self;
}

sub serverSetup {
  my $self=shift;

  &GenUtil::log("JobServer::serverSetup",
		"$self->{status}->{jobstodo} jobs to do");
}

sub serverTask {
  my $self=shift;
  return ($self->{status}->{jobsdone} == $self->{status}->{jobstodo} &&
          $self->{status}->{activeclients} == 0);
}

sub serverFinish {
  my $self=shift;
  
  if (defined $self->{ensemble}) {
    $self->{ensemble}->saveProp();
    $self->{ensemble}->save();
  }

  $self->SUPER::serverFinish();
}

sub processCommand {
  my $self=shift;
  my $client=shift;
  my $cmd=shift;
  my $htmlfile=shift;
  my @arg=@_;

  if ($cmd eq ".MKDIR") {
    if ($arg[0]!~/^\// && $arg[0]!~/\.\./) {
      &GenUtil::makeDir($arg[0]);
    } else {
      &GenUtil::log("JobServer::processCommand",
		    "directory name $arg[0] not accepted by server",1);
    }
  } else {
    return $self->SUPER::processCommand($client,$cmd,$htmlfile,@arg);
  }

  close $client;
  return 0;

}

sub initClient {
  my $self=shift;
  my $clientid=shift;

  my $sc=$self->{client}->{$clientid};

  $sc->{jobsreq}=0;
  $sc->{jobscompleted}=0;

  $sc->{lasttime}=undef;
  $sc->{inittime}=time;
  $sc->{exittime}=undef;
  $sc->{status}="active";

  $sc->{current}=undef;

  $self->{status}->{activeclients}++;

  return "";
}

sub exitClient {
  my $self=shift;
  my $clientid=shift;

  my $sc=$self->{client}->{$clientid};

  $sc->{exittime}=time;
  $sc->{status}="finished";
  
  $self->{status}->{activeclients}--;
  $self->{status}->{finishedclients}++;
      
  if (defined $sc->{current}) {
    $sc->{current}->{status}="ready";
    $self->{status}->{jobsrunning}--;
    $self->{status}->{jobsready}++;
    undef $sc->{current}->{job};
  }
}

sub processClientResponse {
  my $self=shift;
  my $clientid=shift;
  my $data=shift;

  my $sc=$self->{client}->{$clientid};

  if (defined $data) {
    my $job=$sc->{current}->{job};

    if (defined $job) {
      $sc->{lasttime}=time;    
      $sc->{jobscompleted}++;
      $sc->{current}->{status}="done";
      $sc->{current}=undef;

      $self->{status}->{jobsrunning}--;
      $self->{status}->{jobsdone}++;
      
      &GenUtil::log("JobServer::processClientResponse","job $job is done");


      if (defined $self->{ensemble}) {
	$self->{ensemble}->setPropsFromString($job,$data);
	$self->{ensemble}->cleanUp($job);
	$self->{ensemble}->saveProp() 
	  if ($self->{status}->{jobsdone}%$self->{updatefrq}==0);
      }
    }
  }

  my $found=undef;
  if ($self->{status}->{jobsready}>0) {
    for (my $i=0; !defined $found && $i<=$#{$self->{joblist}}; $i++) {
      $found=$i 
	if ($self->{joblist}->[$i]->{status} eq "ready");
    }

    if (defined $found) {
      $sc->{jobsreq}++;
      $sc->{lasttime}=time;
      $sc->{current}=$self->{joblist}->[$found];
      $sc->{current}->{status}="running";
      $sc->{current}->{client}=$sc;
      $self->{status}->{jobsready}--;
      $self->{status}->{jobsrunning}++;
	  
      my $jobnum=$sc->{current}->{job};

      &GenUtil::log("JobServer::processClientResponse","sending job $jobnum");

      $self->respondToClient($clientid,"$jobnum");
    } else {
      $self->{status}->{jobsready}=0;
    }
  }

  if (!defined $found) {
    &GenUtil::log("JobServer","No more jobs left to do");
    $self->respondToClient($clientid,"NO MORE JOBS");
  }
}

sub getInfo {
  my $self=shift;
  my $name=shift;
  my $arg=shift;
  
  if ($name eq "ACTIVE") {
    return sprintf("%d\n",$self->{status}->{activeclients});
  } elsif ($name eq "FINISHED") {
    return sprintf("%d\n",$self->{status}->{finishedclients});
  } elsif ($name eq "READY") {
    return sprintf("%d\n",$self->{status}->{jobsready});
  } elsif ($name eq "RUNNING") {
    return sprintf("%d\n",$self->{status}->{running});
  } elsif ($name eq "DONE") {
    return sprintf("%d\n",$self->{status}->{jobsdone});
  } elsif ($name eq "JOBINFO") {
    my $buffer="";
    foreach my $j ( @{$self->{joblist}} ) {
      if (!defined $arg ||  $j->{job} eq $arg) {
	$buffer.=sprintf("%s %s %d\n",
			 $j->{job},$j->{status},
			 $j->{client}->{id});
      }
    }
    return $buffer;
  } elsif ($name eq "CLIENTINFO") {
    my $buffer="";
    foreach my $c ( keys %{$self->{client}} ) {
      if (!defined $arg ||  $c->{id} eq $arg) {
	$buffer.=sprintf("%s %s %s %d %d %d %d %d\n",
	$c->{id},$c->{hostname},$c->{status},
	$c->{jobsreq},$c->{jobscompleted},
	(defined $c->{lasttime})?$c->{lasttime}:-1,
	(defined $c->{inittime})?$c->{inittime}:-1,
	(defined $c->{exittime})?$c->{exittime}:-1);
      }
    }
    return $buffer;
  } elsif ($arg eq "SUMMARY") {
    return sprintf("%d %d %d %d %d\n",
		   $self->{status}->{activeclients},
		   $self->{status}->{finishedclients},
		   $self->{status}->{jobsready}, 
		   $self->{status}->{jobsrunning},
		   $self->{status}->{jobsdone});
  } else {
    my $buffer="";
    $buffer.=sprintf("active clients  : %d\n",$self->{status}->{activeclients});
    $buffer.=sprintf("finished clients: %d\n",$self->{status}->{finishedclients});
    $buffer.=sprintf("jobs ready      : %d\n",$self->{status}->{jobsready});
    $buffer.=sprintf("jobs running    : %d\n",$self->{status}->{jobsrunning});
    $buffer.=sprintf("jobs done       : %d\n",$self->{status}->{jobsdone});
    return $buffer;
  }
}

1;
