#!/usr/bin/env perl
#
# Startup Job Server
#
# http://mmtsb.scripps.edu/doc/jobserver.pl.html
# 2000, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:    jobserver.pl [options] tag\n";
  printf STDERR "options:  [-dir workdir]\n";
  printf STDERR "          [-runs [from:]to]\n";
  printf STDERR "          [-use proptag]\n";
  printf STDERR "          [-log logfile]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use IO::Handle;
use IO::File;
use Sys::Hostname;

use JobServer;
use GenUtil;
use Ensemble;

my $from;
my $to;
my $dir=".";
my $proptag="etot";
my $tag;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-runs") {
    shift @ARGV;
    my ($n1,$n2)=split(/:/,shift @ARGV);
    if (defined $n2) {
      $from=$n1;
      $to=$n2;
    } else {
      $from=1;
      $to=$n1;
    }
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-use") {
    shift @ARGV;
    $proptag=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    &GenUtil::setLogFile(shift @ARGV);
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option $ARGV[0]\n";
    &usage();
  } else {
    $tag=shift @ARGV
      if (!defined $tag);
  }
}

&usage()
  if (!defined $tag);

my $ens=Ensemble->new($tag,$dir);

my $jobServer=JobServer->new($ens->jobList($from,$to,$proptag),$ens,10);
my ($port,$sid,$pid)=$jobServer->run(4000);

printf "%s:%d:%d started\n",hostname,$port,$sid;

waitpid($pid,0);

