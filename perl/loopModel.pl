#!/usr/bin/env perl

# build loops with Modeller
#
# 2009, Michael Feig, Michigan State University

sub usage {
  printf STDERR "usage:   loopModel.pl [options] [PDBfile]\n";
  printf STDERR "options:              [-loop inx:sequence[:secstr][=inx:sequence[:secstr]]]\n";
  printf STDERR "                      [-models number]\n";
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

my $models=5; 

my $pdbname;
my $loopinp;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-loop") {
    shift @ARGV;
    $loopinp=shift @ARGV;
  } elsif ($ARGV[0] eq "-models") {
    shift @ARGV;
    $models=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option $ARGV[0]\n";
    &usage();
  } else {
    $pdbname = shift @ARGV;
    $done=1;
  }
}

die "no loops specified" if (!defined $loopinp);

my $tag=$$;

my $tmpdir=sprintf("$$-loops");
&GenUtil::makeDir($tmpdir);

my $mol=Molecule::new();
$mol->readPDB($pdbname);

my $loops=();

foreach my $f ( split(/=/,$loopinp) ) {
  my @ff=split(/:/,$f);
  my $sinx=$ff[0];
  my $trec={};
  if ($sinx=~/([A-Z])([0-9]+)/) {
     $trec->{inx}=$2;
     $trec->{chain}=$1;
  } else {
     $trec->{inx}=$sinx;
     $trec->{chain}=$mol->{chain}->[0]->{id};
  } 
  $trec->{seq}=$ff[1];
  $trec->{from}=$trec->{inx};
  $trec->{to}=$trec->{inx}+length($trec->{seq})-1;
  $trec->{secstr}=$ff[2] if ($#ff>1);
  push(@{$loops},$trec);
}

my $modinp=sprintf("mod%d.py",$tag);

my @loopstr=();
my $looprangestr="";
my $secstrrestraints="";

foreach my $l ( @{$loops} ) {
  my $nmol=Molecule::new();
  $nmol->fromSequence(sprintf("%s%d",$l->{chain},$l->{inx}),"CA",$l->{seq});
  $nmol->merge($mol);
  $mol=$nmol->clone(1);

  if ($l->{chain} ne "") {
    push(@loopstr,
         sprintf("                         self.residue_range('%d:%s', '%d:%s')",$l->{from},$l->{chain},$l->{to},$l->{chain}));
  } else {
    push(@loopstr,
         sprintf("                         self.residue_range('%d', '%d')",$l->{from},$l->{to}));
  }
  if (exists $l->{secstr}) {
    for (my $is=0; $is<length($l->{secstr}); $is++) {
      my $s=substr($l->{secstr},$is,1);
      if ($s eq "H") {
        if ($l->{chain} ne "") {
	  $secstrrestraints.=
             sprintf("       rsr.add(secondary_structure.alpha(self.residue_range('%d:%s', '%d:%s')))\n",
		  $l->{from}+$is,$l->{chain},$l->{from}+$is,$l->{chain});
        } else {
  	  $secstrrestraints.=
	     sprintf("       rsr.add(secondary_structure.alpha(self.residue_range('%d', '%d')))\n",
	 	  $l->{from}+$is,$l->{from}+$is);
        }
      } elsif ($s eq "E") {
        if ($l->{chain} ne "") {
	  $secstrrestraints.=
             sprintf("       rsr.add(secondary_structure.strand(self.residue_range('%d:%s', '%d:%s')))\n",
		  $l->{from}+$is,$l->{chain},$l->{from}+$is,$l->{chain});
        } else {
  	  $secstrrestraints.=
	     sprintf("       rsr.add(secondary_structure.strand(self.residue_range('%d', '%d')))\n",
	 	  $l->{from}+$is,$l->{from}+$is);
        }
      }
    }
  }
}

if ($#loopstr>=0) {
  $looprangestr.=join(",\n",@loopstr);
}

if (length($secstrrestraints)>0) {
   $secstrrestraints="    def special_restraints(self,aln):\n       rsr = self.restraints\n       at = self.restraints\n".$secstrrestraints."\n";
}

$mol->writePDB("$tmpdir/$tag.pdb",translate=>"GENERIC");

open OUT,">$tmpdir/$modinp";
print OUT <<ENDinp;
from modeller import *
from modeller.automodel import *
env = environ()
class MyLoop(loopmodel):
    def select_loop_atoms(self):
       return selection(
$looprangestr
                       )
$secstrrestraints
m = MyLoop(env,
           inimodel='$tag.pdb', 
           sequence='$tag',
           loop_assess_methods=(assess.DOPE)) 
m.loop.starting_model= 1
m.loop.ending_model  = $models
m.loop.md_level = refine.fast 
m.make()
ENDinp
close OUT;

my $modexec;
if ($ENV{'MODELLEREXEC'} ne "") {
   $modexec=$ENV{'MODELLEREXEC'};
} else {
   $modexec=&GenUtil::findExecutable("mod9v8");
   $modexec=&GenUtil::findExecutable("mod9v7") if (!defined $modexec);
   $modexec=&GenUtil::findExecutable("mod9v6") if (!defined $modexec);
   $modexec=&GenUtil::findExecutable("mod9v5") if (!defined $modexec);
   $modexec=&GenUtil::findExecutable("mod9v4") if (!defined $modexec);
   $modexec=&GenUtil::findExecutable("mod9v3") if (!defined $modexec);
   $modexec=&GenUtil::findExecutable("mod9v2") if (!defined $modexec);
   $modexec=&GenUtil::findExecutable("mod9v1") if (!defined $modexec);
}

die "cannot find MODELLER exectuable" if (!defined $modexec);
 
system "cd $tmpdir; $modexec $modinp";

my @model=();
for (my $i=1; $i<=$models; $i++) {
  system sprintf("cp $tmpdir/$$.BL%04d0001.pdb model.%d.pdb\n",$i,$i);
  push(@model,0);
}
  
my $mi=0;
open INP,"$tmpdir/mod$tag.log";
while (<INP>) {
  chomp;
  if (/Filename.+DOPE score/) {
    $mi++;
  } elsif ($mi>0 && /BL.+pdb\s+\S+\s+(\S+)/) {
    $model[$mi-1]=$1;
    $mi++;
  }
}
close INP;

my $minscore=999999999;
my $mininx=-1;
for (my $i=1; $i<=$models; $i++) {
  printf STDOUT "model %d: %s\n",$i,$model[$i-1];
  if ($model[$i-1]<$minscore && ($model[$i-1]>0.0001 || $model[$i-1]<-0.0001) &&
     -r "model.$i.pdb") {
    $minscore=$model[$i-1];
    $mininx=$i;
  }
}

system "/bin/rm -rf $tmpdir" if (-d $tmpdir);

1;

