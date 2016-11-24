#!/usr/bin/env perl
#
# clusters a set of structures based on RMSD or residue contact map
#
# http://mmtsb.scripps.edu/doc/cluster.pl.html
# 2001, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   cluster.pl [options] [files]\n";
  printf STDERR "options: [-jclust] [-kclust]\n";
  printf STDERR "         [-maxnum value] [-minsize value] [-maxlevel value]\n";
  printf STDERR "         [-radius value] [-[no]iterate]\n";
  printf STDERR "         [-mode rmsd|contact|phi|psi|phipsi|mix]\n";
  printf STDERR "         [-contmaxdist value] [-mixfactor value]\n";
  printf STDERR "         [-pdb | -sicho]\n";
  printf STDERR "         [-selmode ca|cb|cab|heavy|all]\n";
  printf STDERR "         [-l min:max[=min:max ...]] [-fitxl]\n";
  printf STDERR "         [-[no]lsqfit]\n";
  printf STDERR "         [-centroid] [-centout template]\n";
  printf STDERR "         [-log file]\n";
  exit 1;
}

require 5.004;

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use Cluster;
use GenUtil;

my $maxnum=4;
my $minsize=undef;
my $maxlevel=999;
my $clmode="rmsd";
my $inpmode="pdb";
my $method="jclust";
my $radius=2.5;
my $mixfactor=0.3;
my $selmode="cab";
my $iterate=1;
my $lsqfit=1;

my $contmaxdist;

my $filelist=();

my $fraglist;
my $fitfraglist=undef;

my $centroid=0;
my $centout="cent";
my $logfile;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-maxnum") {
    shift @ARGV;
    $maxnum=shift @ARGV;
  } elsif ($ARGV[0] eq "-jclust") {
    shift @ARGV;
    $method="jclust";
  } elsif ($ARGV[0] eq "-kclust") {
    shift @ARGV;
    $method="kclust";
  } elsif ($ARGV[0] eq "-minsize") {
    shift @ARGV;
    $minsize=shift @ARGV;
  } elsif ($ARGV[0] eq "-maxlevel") {
    shift @ARGV;
    $maxlevel=shift @ARGV;
  } elsif ($ARGV[0] eq "-mode") {
    shift @ARGV;
    $clmode=shift @ARGV;
  } elsif ($ARGV[0] eq "-selmode") {
    shift @ARGV;
    $selmode=shift @ARGV;
  } elsif ($ARGV[0] eq "-lsqfit") {
    shift @ARGV;
    $lsqfit=1;
  } elsif ($ARGV[0] eq "-nolsqfit") {
    shift @ARGV;
    $lsqfit=0;
  } elsif ($ARGV[0] eq "-radius") {
    shift @ARGV;
    $radius=shift @ARGV;
  } elsif ($ARGV[0] eq "-iterate") {
    shift @ARGV;
    $iterate=1;
  } elsif ($ARGV[0] eq "-noiterate") {
    shift @ARGV;
    $iterate=0;
  } elsif ($ARGV[0] eq "-mixfactor") {
    shift @ARGV;
    $mixfactor=shift @ARGV;
  } elsif ($ARGV[0] eq "-contmaxdist") {
    shift @ARGV;
    $contmaxdist=shift @ARGV;
  } elsif ($ARGV[0] eq "-pdb") {
    shift @ARGV;
    $inpmode="pdb";
  } elsif ($ARGV[0] eq "-sicho") {
    shift @ARGV;
    $inpmode="sicho";
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fraglist=shift @ARGV;
    $fitfraglist=$fraglist unless (defined $fitfraglist);
  } elsif ($ARGV[0] eq "-fitxl") {
    shift @ARGV;
    $fitfraglist="fitxl";
  } elsif ($ARGV[0] eq "-centroid") {
    shift @ARGV;
    $centroid=1;
  } elsif ($ARGV[0] eq "-centout") {
    shift @ARGV;
    $centout=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] =~ /^-/) {
    die "unknown option $ARGV[0]";
  } else {
    my $fname=shift @ARGV;
    push(@{$filelist},$fname) 
      if (-r $fname && !-z $fname);
  }
}

$fitfraglist=undef if ($fitfraglist eq "fitxl");

die "contact clustering not available with kclust"
  if ($method eq "kclust" && $clmode eq "contact");

die "$clmode clustering not available with jclust"
  if ($method eq "jclust" && $clmode eq "mix");

&GenUtil::setLogFile($logfile) if (defined $logfile);

if (!defined $filelist || $#{$filelist}<0) {
  while (<>) {
    chomp;
    push(@{$filelist},$_) 
      if (-r $_ && !-z $_);
  }
}

my $numfiles=$#{$filelist}+1;
$minsize=$numfiles if (!defined $minsize);

&GenUtil::log("cluster.pl","inpmode: $inpmode, clmode: $clmode, maxnum: $maxnum, minsize: $minsize, files: $numfiles");

my $cluster=Cluster::new(filetype=>$inpmode, clustermode=>$clmode,
			 contmaxdist=>$contmaxdist, fraglist=>$fraglist,
                         fitfraglist=>$fitfraglist, 
                         maxnum=>$maxnum, minsize=>$minsize, radius=>$radius,
			 mixfactor=>$mixfactor,iterate=>$iterate,lsqfit=>$lsqfit,
			 method=>$method,selmode=>$selmode);

$cluster->setFileList($filelist);
$cluster->clusterHierarchy($maxlevel);

$cluster->writeFile("-",$centroid,$centout);

exit 0;

