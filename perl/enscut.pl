#!/usr/bin/env perl

# generate structure ensemble with cut out protein fragments 
# for loop modeling
#
# http://mmtsb.scripps.edu/doc/enscut.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   enscut.pl [options] intag outtag\n";
  printf STDERR "options: [-dir workdir]\n";
  printf STDERR "         [-run [from:]to]\n";
  printf STDERR "         [-opt file[:file]]\n";
  printf STDERR "         [-l min:max[=min:max=...]]\n";
  printf STDERR "         [-hard value[:force]]\n";
  printf STDERR "         [-soft value[:force]]\n";
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

use Sys::Hostname;

use GenUtil;
use Molecule;
use Ensemble;

my $intag;
my $outtag;
my $dir=".";
my ($from,$to);

my $fraglist;

my %par;

my %defpar = ( 
  hardcutoff      => 16.0,
  softcutoff      => 12.0,
  hardcutforce    => 10.0,
  softcutforce    => 0.5
	    );

my $optfile;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-run") {
    shift @ARGV;
    ($from,$to)=split(/:/,shift @ARGV);
    if (!defined $to) {
      $to=$from;
      $from=1;
    }
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fraglist=shift @ARGV;
  } elsif ($ARGV[0] eq "-hard") {
    shift @ARGV;
    ($par{hardcutoff},$par{hardcutforce})=split(/:/,shift @ARGV);
  } elsif ($ARGV[0] eq "-soft") {
    shift @ARGV;
    ($par{softcutoff},$par{softcutforce})=split(/:/,shift @ARGV);
  } elsif ($ARGV[0] eq "-opt") {
    shift @ARGV;
    $optfile=shift @ARGV;
  } else {
    if (!defined $intag) {
      $intag=shift @ARGV;
    } elsif (!defined $outtag) {
      $outtag=shift @ARGV;
    }
  }
}

&usage() if (!defined $intag || !defined $outtag);

my $ens=Ensemble->new($outtag,$dir);

if (defined $optfile) {
  foreach my $o ( split(/:/,$optfile) ) {
    $ens->readOptions($o);
  }
}

$ens->set(fraglist=>$fraglist);
$ens->setOption(%par);

foreach my $p ( keys %defpar ) {
  if (!defined $ens->{opt} || !defined $ens->{opt}->{$p}) {
    $ens->setOption($p=>$defpar{$p});
  }
}

$ens->save();

my $farcutoff=24.0;
my $filelist=();

my $jlist=$ens->jobList($from,$to,"none");

my $flist;
my $clist;

foreach my $j (@{$jlist}) {
  my $datadir=$ens->{dir}."/".&GenUtil::dataDir($j);
  my $f=$datadir."/".$intag.".pdb";

  my $mol=Molecule::new($f);

  if (!defined $flist) {
    $flist=();
    $clist=();
    $mol->setValidResidues($ens->getFragList());
    
    for (my $ic=0; $ic<=$#{$mol->{chain}}; $ic++) {
      my $c=$mol->{chain}->[$ic];
      my $newfrec={};
      $newfrec->{inx}=$ic;
      $newfrec->{res}=();

      my $newcrec={};
      $newcrec->{nres}=$#{$c->{res}}+1;
      $newcrec->{id}=$c->{id};
      $newcrec->{stat}=();
      $newcrec->{num}=();
      
      for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
	my $r=$c->{res}->[$ir];
	if ($r->{valid}) {
	  my $fresrec={};
	  $fresrec->{inx}=$ir;
	  push(@{$newfrec->{res}},$fresrec);
	  $newcrec->{stat}->[$ir]=5;
	} else {
	  $newcrec->{stat}->[$ir]=0;
	} 
	$newcrec->{num}->[$ir]=$r->{num};
      }

      push(@{$flist},$newfrec) if ($#{$newfrec->{res}}>=0);
      push(@{$clist},$newcrec);
    }
  }

  for (my $icc=0; $icc<=$#{$clist}; $icc++) {
    my $cchain=$mol->{chain}->[$icc];
    my $cstat=$clist->[$icc]->{stat};
    for (my $icr=0; $icr<$clist->[$icc]->{nres}; $icr++) {
      my $cres=$cchain->{res}->[$icr];
      if ($cstat->[$icr]<3) {
      DONE:
	foreach my $fc ( @{$flist} ) {
	  my $fchain=$mol->{chain}->[$fc->{inx}];
	  foreach my $fr ( @{$fc->{res}} ) {
	    my $fres=$fchain->{res}->[$fr->{inx}];
	    my $mind=$mol->minDistance($cres,$fres);
	    if ($mind<$ens->{opt}->{softcutoff}) {
	      $cstat->[$icr]=3;
	      last DONE;
	    } elsif ($mind<$ens->{opt}->{hardcutoff}) {
	      $cstat->[$icr]=1;
	    }
	  }
	}
      }
    }
  }
}

my $havereslist=();
my $havecons3list=();
my $havecons1list=();

for (my $icc=0; $icc<=$#{$clist}; $icc++) {
  my $cstat=$clist->[$icc]->{stat};
  my $cnum=$clist->[$icc]->{num};

  my $haveres=();
  my $havecons1=();
  my $havecons3=();

  for (my $icr=0; $icr<$clist->[$icc]->{nres}; $icr++) {
    $cstat->[$icr]=1 if ($cstat->[$icr]==0 &&
			 (($icr>0 && $cnum->[$icr-1]==$cnum->[$icr]-1 && 
			   $cstat->[$icr-1]==3) ||
			  ($icr<$clist->[$icc]->{nres}-1 && 
			   $cnum->[$icr+1]==$cnum->[$icr]+1 &&
			   $cstat->[$icr+1]==3)));

    push(@{$haveres},$cnum->[$icr]) if ($cstat->[$icr]>0);
    push(@{$havecons1},$cnum->[$icr]) if ($cstat->[$icr]==1);
    push(@{$havecons3},$cnum->[$icr]) if ($cstat->[$icr]==3);
  }

  if ($#{$haveres}>=0) {
    foreach my $l ( @{&GenUtil::fragListFromArray($haveres,$clist->[$icc]->{id})} ) {
      push(@{$havereslist},$l);
    }
  }
  
  if ($#{$havecons3}>=0) {
    foreach my $l ( @{&GenUtil::fragListFromArray($havecons3,$clist->[$icc]->{id})} ) {
      push(@{$havecons3list},$l);
    }
  }

  if ($#{$havecons1}>=0) {
    foreach my $l ( @{&GenUtil::fragListFromArray($havecons1,$clist->[$icc]->{id})} ) {
      push(@{$havecons1list},$l);
    }
  }
}

foreach my $j (@{$jlist}) {
  my $datadir=$ens->{dir}."/".&GenUtil::dataDir($j);
  my $inf=$datadir."/".$intag.".pdb";
  my $outf=$datadir."/".$ens->{tag}.".pdb";

  my $mol=Molecule::new($inf);
  $mol->setValidResidues($havereslist);
  $mol=$mol->clone(1);

  $mol->writePDB($outf,translate=>"CHARMM22");

  $ens->cleanUp($j);
}

my $consstr;
if ($#{$havecons1list}>=0) {
  $consstr="heavy,self,".
    &GenUtil::fragOptionFromList($havecons1list,$ens->{opt}->{hardcutforce})."=".
    &GenUtil::fragOptionFromList($havecons3list,$ens->{opt}->{softcutforce}).",";
} else {
  $consstr="heavy,self,".
    &GenUtil::fragOptionFromList($havecons3list,$ens->{opt}->{softcutforce}).",";
}

$ens->setOption(cons=>$consstr);
$ens->save();

exit 0;
