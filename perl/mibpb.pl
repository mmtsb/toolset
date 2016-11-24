#!/usr/bin/env perl

# get scalar quantity from CHARMM
#
# Michael Feig/Duan Chen, 2007, MSU

sub usage {
  printf STDERR "usage:   mibpb.pl [options] [PDBfile]\n";
  printf STDERR "options: [-par CHARMMoptions]\n";
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
use CHARMM;

my %par;

my $inpfile="-";

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    foreach my $p ( split(/,/,shift @ARGV) ) {
      my ($key,$val)=split(/=/,$p);
      $par{$key}=(defined $val)?$val:1;
    }
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    $inpfile=(shift @ARGV);
    $done=1;
  }
}

my $charmm=&CHARMM::new();
$charmm->loadParameters(%par);
$charmm->setupFromPDB($inpfile);

my $charges=$charmm->getScalar("charge");
my $radii=$charmm->getScalar("radius");

open XYZR,">$$.xyzr";
open PQR,">$$.pqr";

my $mol=$charmm->{molecule};

my $i=0;
foreach my $c ( @{$mol->{chain}} ) {
  foreach my $a ( @{$c->{atom}} ) {
    printf XYZR "%11.5f %11.5f %11.5f %9.5f\n",$a->{xcoor},$a->{ycoor},$a->{zcoor},$radii->[$i];
    printf PQR "%9.5f %9.5f %9.5f %8.5f\n",$a->{xcoor},$a->{ycoor},$a->{zcoor},$charges->[$i];
    $i++;
  }
}

close PQR;
close XYZR;

$charmm->finish();

exit 0;

