#!/usr/bin/env perl

# compares dihedral distribution of a 
# protein structure compared to a reference structure
#
# http://mmtsb.scripps.edu/doc/dihed.pl.html
# 2001, Michael Feig, Brooks group, TSRI
#
## 2009-06-13: Modified by Lev Gorenstein <lev@purdue.edu>
## 	Added capability to report Chi2, Chi3, Chi4 and Chi5 backbone
## 	torsions.  This also required adding corresponding helper methods
## 	to Analyze.pm module.  Corresponding dihedral definitions are from
## 	http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
## 	Note that now Ensemble.pm could (and probably should) use these
## 	methods, but I haven't modified it.

sub usage {
  printf STDERR "usage:   dihed.pl [options] [refPDB [cmpPDB]]\n";
  printf STDERR "options: [-l min:max[...]]\n";
  printf STDERR "         [-list phi|psi|chi1|chi2|chi3|chi4|chi5|omega]\n";
  printf STDERR "         [-atoms [chain]res:name [chain]res:name [chain]res:name [chain]res:name\n";
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

my %list = { phi => 0,
	     psi => 0,
	     chi1 => 0,
	     chi2 => 0,
	     chi3 => 0,
	     chi4 => 0,
	     chi5 => 0,
             omega =>0 };

my $fraglist;
my $refpdb;

my $aname1;
my $aname2;
my $aname3;
my $aname4;


my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-list") {
    shift @ARGV;
    foreach my $tl ( split(/,/,shift @ARGV) ) {
      $list{$tl}=1;
    }
  } elsif ($ARGV[0] eq "-atoms") {
    shift @ARGV;
    $aname1=shift @ARGV;
    $aname2=shift @ARGV;
    $aname3=shift @ARGV;
    $aname4=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $refpdb = shift @ARGV;
    $done=1;
  }
}

my $analyze;

if (defined $aname1) {
  my $mol=Molecule::new();
  $mol->readPDB($refpdb);
  
  my $a1=&getAtom($mol,$aname1);
  my $a2=&getAtom($mol,$aname2);
  my $a3=&getAtom($mol,$aname3);
  my $a4=&getAtom($mol,$aname4);

  printf "%f\n",&GenUtil::dihedral($a1,$a2,$a3,$a4);

} elsif ($list{phi} || $list{psi} || $list{chi1} || $list{chi2} || $list{chi3} || $list{chi4} || $list{chi5} || $list{omega}) {
  my $mol=Molecule::new();
  $mol->readPDB($refpdb);
  $mol->setValidResidues($fraglist,0)
    if (defined $fraglist);

  &Analyze::phipsi($mol) if ($list{phi} || $list{psi} || $list{omega});
  &Analyze::chi1($mol) if ($list{chi1});
  &Analyze::chi2($mol) if ($list{chi2});
  &Analyze::chi3($mol) if ($list{chi3});
  &Analyze::chi4($mol) if ($list{chi4});
  &Analyze::chi5($mol) if ($list{chi5});

  foreach my $c ( @{$mol->activeChains()} ) {
    foreach my $r ( @{$c->{res}} ) {
      if ($r->{valid}) {
	printf "%s%d:%s",$r->{name},$r->{num},$r->{chain};
	printf " %8.3f",$r->{phi}  if ($list{phi});
	printf " %8.3f",$r->{psi}  if ($list{psi});
	printf " %8.3f",$r->{omega}  if ($list{omega});
	printf " %8.3f",$r->{chi1} if ($list{chi1});
	printf " %8.3f",$r->{chi2} if ($list{chi2});
	printf " %8.3f",$r->{chi3} if ($list{chi3});
	printf " %8.3f",$r->{chi4} if ($list{chi4});
	printf " %8.3f",$r->{chi5} if ($list{chi5});
	printf "\n";
      }
    }
  }
} else {
  my $refmol=Molecule::new();
  $refmol->readPDB($refpdb);
  $analyze=Analyze::new($refmol);

  my $mol=Molecule::new();
  $mol->readPDB(shift @ARGV);
  $mol->setValidResidues($fraglist,0)
    if (defined $fraglist);

  my ($phi,$psi,$cphi,$cpsi)=$analyze->phipsiRMSD($mol);
  my ($chi1,$cchi1)=$analyze->chi1RMSD($mol);
  my ($chi2,$cchi2)=$analyze->chi2RMSD($mol);
  my ($chi3,$cchi3)=$analyze->chi3RMSD($mol);
  my ($chi4,$cchi4)=$analyze->chi4RMSD($mol);
  my ($chi5,$cchi5)=$analyze->chi5RMSD($mol);

  printf STDOUT "phi: %1.3f ( %1.2f \% ), psi: %1.3f ( %1.2f \% ), chi1: %1.3f ( %1.2f \% ), chi2: %1.3f ( %1.2f \% ), chi3: %1.3f ( %1.2f \% ), chi4: %1.3f ( %1.2f \% ), chi5: %1.3f ( %1.2f \% )\n",
  $phi,$cphi,$psi,$cpsi,$chi1,$cchi1,$chi2,$cchi2,$chi3,$cchi3,$chi4,$cchi4,$chi5,$cchi5;
}  

1;

sub getAtom {
  my $m=shift;
  my $s=shift;

  my $chain;
  my $resid;
  my $aname="";

  my @f=split(/:/,$s);
  if ($#f==0) {
    $aname=$s;
  } else {
    $aname=$f[1];
    if ($f[0]=~/([A-Z\+])([0-9]+)/) {
      $chain=$1;
      $resid=$2;
    } else {
      $resid=$f[0];
    }
  }

  my $c=$m->getChain($chain);
  my $r;
  if (defined $resid) {
    $r=$m->getResidueInChain($resid,$c);
  } else {
    $r=$c->{res}->[0];
  }
  my $a=$m->getAtomInResidue($r,$aname,$c);
  return $a;
}
