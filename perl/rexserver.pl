#!/usr/bin/env perl
#
# All-atom replica exchange server
#
# http://mmtsb.scripps.edu/doc/rexserver.pl.html
# 2001, Michael Feig, Brooks group, TSRI
# 2001, John Karanicolas, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   rexserver.pl [options] [initfiles]\n";
  printf STDERR "options: [-n runs]\n";
  printf STDERR "         [-f listfile]\n";
  printf STDERR "         [-dir workdir]\n";
  printf STDERR "         [-monsster]\n";
  printf STDERR "         [-ens tag] [-ensdir dir]\n";
  printf STDERR "         [-par initruns=value,equilruns=value,\n";
  printf STDERR "               [no]save,savebestfreq=value,\n";
  printf STDERR "               natpdb=file,seq=file]\n";
  printf STDERR "         [-l refPDB min:max[=min:max ...]]\n";
  printf STDERR "         [-temp nwin:min:max]\n";
  printf STDERR "         [-condfile file]\n";
  printf STDERR "         [-log logfile]\n";
  printf STDERR "         [-serverid num]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use Sys::Hostname;

use LatSimReXServer;
use ReXServer;
use Ensemble;
use LatEnsemble;
use GenUtil;

my $presetid;

my $monsster=0;

my %par=();
my %defpar = (
   initruns      => 5,
   equilruns     => 20, 
   save          => 1,
   savebestfreq  => 1,
);

my $nwindows;
my $mintemp;
my $maxtemp;

my $nruns=100;

my $dir=".";

my @initfile;
my $listfile;

my $condfile;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-n") {
    shift @ARGV;
    $nruns=shift @ARGV;
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-f") {
    shift @ARGV;
    $listfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-monsster") {
    shift @ARGV;
    $monsster=1;
  } elsif ($ARGV[0] eq "-ens") {
    shift @ARGV;
    $par{enstag}=shift @ARGV;
  } elsif ($ARGV[0] eq "-ensdir") {
    shift @ARGV;
    $par{ensdir}=shift @ARGV;
  } elsif ($ARGV[0] eq "-serverid") {
    shift @ARGV;
    $presetid=shift @ARGV;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      if (defined $val) {
	$par{$key}=$val;
      } else {
	if ($key=~/^no(.+)$/) {
	  $par{$1}=0;
	} else {
	  $par{$key}=1;
	}
      }
    }
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $par{fragref}=shift @ARGV;
    $par{fraglist}=shift @ARGV;
  } elsif ($ARGV[0] eq "-temp") {
    shift @ARGV;
    ($nwindows,$mintemp,$maxtemp)=split(/:/,shift @ARGV);
  } elsif ($ARGV[0] eq "-condfile") {
    shift @ARGV;
    $condfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    my $logfile=shift @ARGV;
    &GenUtil::setLogFile("$logfile");
  } elsif ($ARGV[0]=~/^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    push(@initfile, shift @ARGV);
  }
}

&GenUtil::log("aarexserver","starting");

if (defined $listfile && &GenUtil::checkFile($listfile)) {
  my $lfile=&GenUtil::getInputFile($listfile);
  while (<$lfile>) {
    chomp;
    push(@initfile,$_);
  }
  close $lfile;
}

my $server=($monsster?
	    LatSimReXServer->new($nruns,$dir,\@initfile,%par):
	    ReXServer->new($nruns,$dir,\@initfile,%par));

if (defined $server->{par}->{enstag}) {
  my $ens=($monsster?
	   LatEnsemble->new($server->{par}->{enstag},$server->{par}->{ensdir}):
	   Ensemble->new($server->{par}->{enstag},$server->{par}->{ensdir}));
  $server->setEnsemble($ens);
  $server->{par}->{savebestfreq}=1 if ($server->{par}->{savebestfreq}<=0);
}

foreach my $p ( keys %defpar ) {
  if (!defined $server->{par} || !defined $server->{par}->{$p}) {
    $server->set($p=>$defpar{$p});
  }
}

if (defined $nwindows && !defined $condfile) {
  my $tlist=&ReXServer::expTempList($nwindows,$mintemp,$maxtemp);
  my $out=&GenUtil::getOutputFile("tmp-cond$$");
  foreach my $t ( @{$tlist} ) {
    printf $out "%f\n",$t;
  }
  close $out;
  $server->setup("tmp-cond$$");
  &GenUtil::remove("tmp-cond$$");
} else {
  $server->setup($condfile);
}

$server->save();

my ($port,$sid,$pid)=
  $server->run(4100,$presetid);

printf "%s:%d:%d started\n",hostname,$port,$sid;
close STDOUT;

waitpid($pid,0);
