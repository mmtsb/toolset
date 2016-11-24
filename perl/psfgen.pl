#!/usr/bin/env perl

# generates a PSF file using NAMD's psfgen
#

sub usage {
  printf STDERR "usage: psfgen.pl -par CHARMMparams -psfout file -pdbout file [pdbFile]\n";
  exit 1;
}

use vars qw ( $perllibdir $exec );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use GenUtil;
use Molecule;
use NAMD;

if ($ENV{'PSFGENEXEC'} ne "") {
  $exec=$ENV{'PSFGENEXEC'};
} else {
  $exec=&GenUtil::findExecutable("psfgen");
}
my %par;

my $fname;
my $psfout;
my $pdbout;
my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-psfout") {
    shift @ARGV;
    $psfout=shift @ARGV;
  } elsif ($ARGV[0] eq "-pdbout") {
    shift @ARGV;
    $pdbout=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $fname = shift @ARGV;
    $done=1;
  }
}

my $namd=&NAMD::new();
$namd->setParameter(%par);

my $mol=Molecule::new();
$mol->readPDB($fname);

#open OUT,"| $exec";
open OUT,">psfinp";

foreach my $n ( split(/:/,$namd->{par}->{xtop}) ) {
  printf OUT "topology %s\n",$n;
}

printf OUT "pdbalias atom ILE CD1 CD\n";
printf OUT "pdbalias residue HOH TIP3\n";
printf OUT "pdbalias atom HOH O OH2\n";

my @pdbs=();
my $inx=1;
my $segnames=$mol->getSegNames();
for my $s ( @{$segnames} ) {
  $mol->setValidSelection($s->{name}.":");
  my $smol=$mol->clone(1);
  my $pdbfile=sprintf("t%d%d",$$,$inx);
  push(@pdbs,$pdbfile);
  if ($smol->{chain}->[0]->{res}->[0]->{name} =~ /HOH|TIP3/) {
    $smol->renumber(0);
  } 
  my $ainx=1;
  foreach my $a ( @{$smol->{chain}->[0]->{atom}} ){
    $a->{atominx}=$ainx++;
  }
  $smol->writePDB($pdbfile,translate=>"CHARMM22",genresno=>1);
  printf OUT "segment %s {\n",$s->{name};
  printf OUT " auto none\n" if ($smol->{chain}->[0]->{res}->[0]->{name} =~ /HOH|TIP3/);
  printf OUT " pdb %s\n",$pdbfile;
  if ($s->{name}=~/^N/) {
     printf OUT  " first 5TER\n";
     printf OUT  " last 3TER\n";
  }
  printf OUT "}\n";

  if (defined $namd->{par}->{patch} && $namd->{par}->{patch} ne "") {
#  printf STDERR "%s\n",$namd->{par}->{patch};
   my @ppl=split(/_/,$namd->{par}->{patch});
   foreach my $pp ( @ppl) {
    if ($pp=~/$s->{name}/) {
     my @pl=split(/:/,$pp);
#     printf STDERR "pp: %s\n",$pp;
     printf OUT "patch %s ",$pl[0];
     if ($#pl>=1) {
      $pl[1]=~s/\./:/; 
      $pl[2]=~s/\./:/; 
      printf OUT "%s %s\n",$pl[1],$pl[2];
     } else {
      $pl[1]=~s/\./:/; 
      printf OUT "%s\n",$pl[1];
     }
    }  
   }
  }
  printf OUT "coordpdb %s %s\n",$pdbfile,$s->{name};
  printf OUT "\n";
  $inx++;
}

printf OUT "writepsf $psfout\n" if (defined $psfout);
printf OUT "writepdb $pdbout\n" if (defined $pdbout);


close OUT;

foreach my $f ( @pdbs ) {
#  system "rm -f $f\n";
}

1;
