#!/usr/bin/env perl

# determines residues in vicinity of given fragment
#
# http://mmtsb.scripps.edu/doc/vicinity.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   vicinity.pl [options] [file(s)]\n";
  printf STDERR "options: [-l min:max[=min:max=...]]\n";
  printf STDERR "         [-hard value[:force]]\n";
  printf STDERR "         [-soft value[:force]]\n";
  printf STDERR "         [-f filelist]\n";
  printf STDERR "         [-[no]margin]\n";
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

my $fraglist;
my $cutoff=16.0;
my $soft=12.0;
my $filelist=();
my $forcesoft=0.5;
my $forcehard=10.0;
my $listfile;
my $margin=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-hard") {
    shift @ARGV;
    my @f=split(/:/,shift @ARGV);
    $cutoff=$f[0];
    $forcehard=$f[1] if (defined $f[1]);
  } elsif ($ARGV[0] eq "-soft") {
    shift @ARGV;
    my @f=split(/:/,shift @ARGV);
    $soft=$f[0];
    $forcesoft=$f[1] if (defined $f[1]);
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-f") {
    shift @ARGV;
    $listfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-margin") {
    shift @ARGV;
    $margin=1;
  } elsif ($ARGV[0] eq "-nomargin") {
    shift @ARGV;
    $margin=0;
  } elsif ($ARGV[0] =~/^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    push(@{$filelist},shift @ARGV);
  }    
}

if (defined $listfile || -r $listfile) {
  my $lfile=&GenUtil::getInputFile($listfile);
  while (<$lfile>) {
    chomp;
    push(@{$filelist},$_);
  }
  close $lfile;
}

&usage() if ($#{$filelist}<0 || !defined $fraglist);

my $flist;
my $clist;

foreach my $f (@{$filelist}) {
  my $mol=Molecule::new($f);

  if (!defined $flist) {
    $flist=();
    $clist=();
    $mol->setValidResidues($fraglist);
    
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
	    if ($mind<$soft) {
	      $cstat->[$icr]=3;
	      last DONE;
	    } elsif ($mind<$cutoff) {
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
    $cstat->[$icr]=1 if ($cstat->[$icr]==0 && $margin &&
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

printf "%s\n",
  &GenUtil::fragOptionFromList($havereslist);

if ($#{$havecons1list}>=0) {
  printf "%s\n",
  &GenUtil::fragOptionFromList($havecons1list,$forcehard)."=".
  &GenUtil::fragOptionFromList($havecons3list,$forcesoft);
} else {
  printf "%s\n",
  &GenUtil::fragOptionFromList($havecons3list,$forcesoft);
}

  
