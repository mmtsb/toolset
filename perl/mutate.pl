#!/usr/bin/env perl

# mutates residues in a given PDB files
#
# http://mmtsb.scripps.edu/doc/mutate.pl.html
# 2003, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   mutate.pl [options] [PDBfile]\n";
  printf STDERR "options: [-seq index:sequence[=index:sequence]]\n";
  printf STDERR "         [-minimize]\n";
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
use Sequence;
use SICHO;
use CHARMM;

my $mutlist;
my $template;
my $fname="-";
my $minimize;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-seq") {
    shift @ARGV;
    foreach my $t ( split(/=/,shift @ARGV) ) {
      my ($tinx,$tseq)=split(/:/,$t);
      my $srec={};
      $srec->{inx}=$tinx;
      $srec->{seq}=$tseq;
      &usage() if (!defined $tseq || $tseq eq "");
      push(@{$mutlist},$srec);    
    }
  } elsif ($ARGV[0] eq "-minimize") {
    shift @ARGV;
    $minimize=1;
  } else {
    $fname = shift @ARGV;
  }
}

my $mol=Molecule::new();
$mol->readPDB($fname);

my $ssbonds=$mol->getSSBonds();

if ($#{$mol->{chain}}>0) {
  printf STDERR "Warning! Only first chain will be used.\n";
}

my $savechainid=$mol->{chain}->[0]->{id};
#printf STDERR "chainid: >%s<\n",$savechainid;

$mol->selectChain("");

$mol->setChain(" ");

if (!defined $mutlist) {
  $mol->writePDB(\*STDOUT,translate=>"GENERIC");
  exit(0);
}

my $seq=Sequence::new($mol,$mutlist);

my $keeparr=();
foreach my $m ( @{$mutlist} ) {
  for (my $i=$m->{inx}; $i<$m->{inx}+length($m->{seq}); $i++) {
    push(@{$keeparr},$i);
  }
}

my $fraglist=&GenUtil::fragListFromArray($keeparr);
my $fragoption=&GenUtil::fragOptionFromList($fraglist);
#printf STDERR "new sequence: %s, fraglist: %s\n",$seq->abbrevSeq(),$fragoption;

my $sicho=SICHO::new(offsetx=>0,offsety=>0,offsetz=>0,resolution=>1.0,intflag=>0);
$sicho->genSimpleFromAllAtom($mol,ca=>1);

my $refpdb="tmp-$$.pdb";
$mol->writePDB($refpdb,ssbond=>0);

my $rmol=Molecule::new();
$rmol->rebuildFromSICHO($seq,$sicho,$fragoption,$refpdb,1);

$rmol->setValidResidues($fraglist);
$mol->merge($rmol->clone(1));

&GenUtil::remove($refpdb);

if (defined $minimize && $minimize) {
  my $charmm=&CHARMM::new();
  $charmm->loadParameters(dielec=>"RDIE",epsilon=>4.0,cuton=>10.0,cutoff=>12.0,cutnb=>15.0,
			  minsteps=>100,sdsteps=>50);
  $mol->generateSegNames();
  $charmm->setupFromMolecule($mol);

  $charmm->setupEnergy();
  $charmm->shake();

  my $cons=();
  my $c={};
  $c->{sel}="heavy";
  $c->{type}="self";
  $c->{force}=10.0;
  $c->{list}=$fraglist;
  $c->{exclmode}=1;
  push (@{$cons},$c);
  $charmm->setupRestraints(1.0,$cons);

  if ($charmm->{par}->{sdsteps}>0) {
    $charmm->minimizeSD();
  }

  $charmm->minimize();
  
  my $chmoutpdb=lc "pdb$$-out";
  $charmm->writePDB($chmoutpdb);

  my $outmol=Molecule::new();
  $outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),chainfromseg=>1);
  $outmol->setChain($savechainid);
  $outmol->setSSBonds($ssbonds);
  foreach my $r ( @{$mol->{chain}->[0]->{res}} ) {
    $r->{seg}=$mol->{chain}->[0]->{res}->[0]->{seg};
  }
  foreach my $a ( @{$mol->{chain}->[0]->{atom}} ) {
    $a->{seg}=$mol->{chain}->[0]->{res}->[0]->{seg};
  }
  $outmol->writePDB(\*STDOUT,translate=>"GENERIC");
  &GenUtil::remove($chmoutpdb);
  $charmm->finish();
} else {
  $mol->setChain($savechainid);
  $mol->setSSBonds($ssbonds);
  foreach my $r ( @{$mol->{chain}->[0]->{res}} ) {
    $r->{seg}=$mol->{chain}->[0]->{res}->[0]->{seg};
  }
  foreach my $a ( @{$mol->{chain}->[0]->{atom}} ) {
    $a->{seg}=$mol->{chain}->[0]->{res}->[0]->{seg};
  }
  $mol->writePDB(\*STDOUT,translate=>"GENERIC");
}
