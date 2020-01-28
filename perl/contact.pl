#!/usr/bin/env perl

# calculates number of native contacts for a given 
# protein structure compared to its native structure
#
# http://mmtsb.scripps.edu/doc/contact.pl.html
# 2000, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   contact.pl [options] [refPDB] cmpPDB\n";
  printf STDERR "options: [-l min:max[...]]\n";
  printf STDERR "         [-read file]\n";
  printf STDERR "         [-list]\n";
  printf STDERR "         [-mindist value]\n";
  printf STDERR "         [-betweenchains]\n";
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

my $list=0;

my $fraglist;
my ($pdb1,$pdb2);
my $readfile;

my $mindist=4.2;

my $betweenchains=0;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-list") {
    shift @ARGV;
    $list=1;
  } elsif ($ARGV[0] eq "-read") {
    shift @ARGV;
    $readfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-mindist") {
    shift @ARGV;
    $mindist=shift @ARGV;
  } elsif ($ARGV[0] eq "-betweenchains") {
    shift @ARGV;
    $betweenchains=1;
  } elsif ($ARGV[0] =~ /^-.+/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $pdb1 = shift @ARGV;
    $pdb2 = shift @ARGV;
    $done=1;
  }
}

my $analyze;
my $mol;

if (defined $pdb1 && defined $pdb2) {
  $mol=Molecule::new($pdb2);
  my $refmol=Molecule::new($pdb1);
  $analyze=Analyze::new($refmol);
} elsif (defined $readfile) {
  $mol=Molecule::new();
  $mol->readPDB($pdb1);
  $analyze=Analyze::new();
  $analyze->readContacts($readfile);
} else {
  $mol=Molecule::new();
  $mol->readPDB($pdb1,splitseg=>1);
  $mol->setValidResidues($fraglist) if (defined $fraglist);
  $analyze=Analyze::new();
  $analyze->{contactReferenceList}=$analyze->contactList($mol,$mindist);
  $analyze->writeContacts("-",$betweenchains);
  exit 0;
} 

$mol->setValidResidues($fraglist) if (defined $fraglist);
my ($natcont, $natcontfrac, $rnatcont, $contlist)=$analyze->contacts($mol);

printf STDOUT "%d contacts ( fraction: %f ), rho: %f\n",
  $natcont,$natcontfrac,$rnatcont;
  
if ($list) {
  foreach my $c ( @{$contlist}) {
    my $r1=$mol->getResidue($c->{res1},$c->{chain1});
    my $r2=$mol->getResidue($c->{res2},$c->{chain2});
    if ($betweenchains == 0 || $c->{chain1} ne $c->{chain2}) {
     printf STDOUT " %1s %-9s - %-9s   %f %f\n",
     (($c->{compd}<4.2)?"*":" "),
     "$r1->{name}:$r1->{chain}$r1->{num}","$r2->{name}:$r2->{chain}$r2->{num}",
     $c->{d},$c->{compd};
    }
  }
}


