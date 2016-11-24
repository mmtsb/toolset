#!/usr/bin/env perl

# calculates RMSD between two protein structures
#
# http://mmtsb.scripps.edu/doc/rms.pl.html
# 2000, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   rms.pl [options] refPDB cmpPDB\n";
  printf STDERR "options: [-out ca|cab|c|o|n|side|back|all...]\n";
  printf STDERR "         [-detailed] [-chains]\n";
  printf STDERR "         [-l min:max[...]]\n";
  printf STDERR "         [-fit] [-fitxl]\n";
  printf STDERR "         [-fitl min:max[...]] [-fitx min:max[...]]\n";
  printf STDERR "         [-fitsel cab|ca|cb|heavy]\n";
  printf STDERR "         [-nowarn]\n";
  printf STDERR "         [-resnumonly] [-useseg]\n";
  printf STDERR "         [-align fasta]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use GenUtil;
use Molecule;
use Analyze;

my @outmode;

my $fitflag=0;
my $fitxlflag=0;

my $fraglist;

my $fitfraglist;
my $fitexclmode;

my $fitselmode="cab";

my ($refpdb,$cmppdb);

my $warn=1;
my $chains=0;
my $detailed=0;

my $resnumonly;

my $alignfile;

my $useseg=0;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-out") {
    shift @ARGV;
    foreach my $om ( split(/,/,shift @ARGV) ) {
      push(@outmode,$om);
    }
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-fit") {
    shift @ARGV;
    $fitflag=1;
    $fitexclmode=0;
  } elsif ($ARGV[0] eq "-fitxl") {
    shift @ARGV;
    $fitflag=1;
    $fitexclmode=1;
  } elsif ($ARGV[0] eq "-fitl") {
    shift @ARGV;
    $fitflag=1;
    $fitexclmode=0;
    $fitfraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-fitx") {
    shift @ARGV;
    $fitflag=1;
    $fitexclmode=1;
    $fitfraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-fitsel") {
    shift @ARGV;
    $fitselmode=shift @ARGV;
  } elsif ($ARGV[0] eq "-nowarn") {
    shift @ARGV;
    $warn=0;
  } elsif ($ARGV[0] eq "-useseg") {
    shift @ARGV;
    $useseg=1;
  } elsif ($ARGV[0] eq "-resnumonly") {
    shift @ARGV;
    $resnumonly=1;
    $warn=0;
  } elsif ($ARGV[0] eq "-align") {
    shift @ARGV;
    $resnumonly=1;
    $warn=0;
    $alignfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-chains") {
    shift @ARGV;
    $chains=1;
  } elsif ($ARGV[0] eq "-detailed") {
    shift @ARGV;
    $detailed=1;
  } elsif ($ARGV[0] =~ /^-.+/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $refpdb = shift @ARGV;
    $cmppdb = shift @ARGV;
    $done=1;
  }
}

my $refmol=Molecule::new($refpdb);
my $cmpmol=Molecule::new();
$cmpmol->readPDB($cmppdb);

my $analyze=Analyze::new($refmol);

if ($fitflag) {
  $fitfraglist=$fraglist
    if (defined $fraglist && !defined $fitfraglist);

  $cmpmol->setValidResidues($fitfraglist,$fitexclmode)
    if (defined $fitfraglist);

  $analyze->lsqfit($cmpmol,$fitselmode,$warn,$resnumonly,$alignfile,$useseg);

  $cmpmol->resetValidResidues()
    if (!defined $fraglist && defined $fitfraglist);
}

$cmpmol->setValidResidues($fraglist,0)
  if (defined $fraglist);

my $rmsd=$analyze->rmsd($cmpmol,$warn,undef,$resnumonly,$alignfile,$useseg);

my @keys=sort keys %{$rmsd};

my @aksel=sort map (/all(_.*)/,@keys);
my @allkeys=();

if ($chains) {
  @allkeys=@aksel;
} else {
  push(@allkeys,"");
}

push(@outmode,"all") 
  if ($#outmode<0);

my $lastchain;
foreach my $ak ( @allkeys ) {
  my $ctag=$ak;
  my $chain=$ctag;
  $chain=~s/_//;

  if ($detailed) {
    printf STDOUT "\n" if (defined $lastchain);
    printf STDOUT "RMSD values for chain $chain:\n"
      if ($ctag ne "");

    printf STDOUT "all      %9.4f\n",$rmsd->{"all$ctag"};
    printf STDOUT " back    %9.4f\n",$rmsd->{"back$ctag"};
    printf STDOUT "  CA     %9.4f\n",$rmsd->{"CA$ctag"};
    printf STDOUT "  N      %9.4f\n",$rmsd->{"N$ctag"};
    printf STDOUT "  C      %9.4f\n",$rmsd->{"C$ctag"};
    printf STDOUT "  O      %9.4f\n",$rmsd->{"O$ctag"};
    printf STDOUT " side    %9.4f\n",$rmsd->{"side$ctag"};
    printf STDOUT "  CAB    %9.4f\n",$rmsd->{"CAB$ctag"};
    printf STDOUT "  CB     %9.4f\n",$rmsd->{"CB$ctag"};

    foreach my $r ( grep (/^side:[A-Z]+$ctag$/,@keys) ) {
      if ($r !~ /side:GLY/) {
	$r=~s/_(.*)$//;
	my @res=split(/:/,$r);
	printf STDOUT "  %-7s%9.4f\n",$res[1],$rmsd->{"$r$ctag"};
      }
    }


    if ($chains || $#aksel<=0) {
      $ctag="_" if ($ctag eq "");
      foreach my $n ( sort { $a <=> $b } map(/all:([0-9]+)$ctag/,@keys) ) {
	my $r=$cmpmol->getResidue($n,$chain);
	my $resname=$r->{name};
	printf STDOUT " %-8s%9.4f%9.4f%9.4f\n",
	"$resname:$n",$rmsd->{"all:$n$ctag"},$rmsd->{"back:$n$ctag"},$rmsd->{"side:$n$ctag"};
      }
      $lastchain=$chain;
    }
  } else {
    foreach my $om (@outmode) {
      $om=uc $om if ($om =~ /^(ca|c|cb|cab|o|n)$/);
      die "no RMSD value for $om available" 
	unless (defined $rmsd->{"$om$ctag"});
      printf STDOUT "%9.4f $om $chain\n",$rmsd->{"$om$ctag"};
    }
  }
}

