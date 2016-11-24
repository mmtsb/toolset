#!/usr/bin/env perl

# read file from archive
#
# http://mmtsb.scripps.edu/doc/readArchive.pl.html
# 2002, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:    readArchive.pl archivefile\n";
  printf STDERR "options:  [-showheader]\n";
  printf STDERR "          [-inx index[:to]] [-step n]\n";
  printf STDERR "          [-pdb template]\n";
  printf STDERR "          [-xtract] [-prefix name] [-listfiles]\n";
  printf STDERR "          [-center] [-fit ref]\n";
  printf STDERR "          [-ensdir dir] [-ens tag] [-at inx]\n";
  printf STDERR "          [-apply cmd]\n";
  printf STDERR "          [-rms CA|CAB|C|O|N|side|back|all ... ref]\n";
  printf STDERR "          [-qscore ref]\n";
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
use Ensemble;
use Analyze;

my $from=1;
my $to=1;
my $step=1;
my $pdbtemplate;
my $showheader=0;
my $arfile;
my $xtract=0;
my $center=0;
my $fitref;
my $listfiles=0;
my $prefix="set";
my $apply;
my $rmsmode;
my $ref;
my $ensdir=".";
my $enstag;
my $ensat;
my $qscore;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-inx") {
    shift @ARGV;
    my @f=split(/:/,shift @ARGV);
    $from=$f[0];
    $to=(defined $f[1])?$f[1]:$from;
  } elsif ($ARGV[0] eq "-step") {
    shift @ARGV;
    $step=shift @ARGV;
  } elsif ($ARGV[0] eq "-pdb") {
    shift @ARGV;
    $pdbtemplate=shift @ARGV;
  } elsif ($ARGV[0] eq "-showheader") {
    shift @ARGV;
    $showheader=1;
  } elsif ($ARGV[0] eq "-xtract") {
    shift @ARGV;
    $xtract=1;
  } elsif ($ARGV[0] eq "-center") {
    shift @ARGV;
    $center=1;
  } elsif ($ARGV[0] eq "-fit") {
    shift @ARGV;
    $fitref=shift @ARGV;
  } elsif ($ARGV[0] eq "-prefix") {
    shift @ARGV;
    $prefix=shift @ARGV;
  } elsif ($ARGV[0] eq "-listfiles") {
    shift @ARGV;
    $listfiles=1;
  } elsif ($ARGV[0] eq "-ens") {
    shift @ARGV;
    $enstag=shift @ARGV;
  } elsif ($ARGV[0] eq "-at") {
    shift @ARGV;
    $ensat=shift @ARGV;
  } elsif ($ARGV[0] eq "-ensdir") {
    shift @ARGV;
    $ensdir=shift @ARGV;
  } elsif ($ARGV[0] eq "-apply") {
    shift @ARGV;
    $apply=shift @ARGV;
  } elsif ($ARGV[0] eq "-rms") {
    shift @ARGV;
    $rmsmode=shift @ARGV;
    $ref=shift @ARGV;
  } elsif ($ARGV[0] eq "-qscore") {
    shift @ARGV;
    $qscore=1;
    $ref=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    $arfile=(shift @ARGV);
  }
}

my ($len,$count,$ft,$aux)=&GenUtil::readArchiveHeader($arfile);

if ($showheader) {
  printf "data length: $len\n";
  printf "data items : $count\n";
  printf "file type  : $ft\n";
  printf "aux        : $aux\n";
} elsif (defined $enstag) {
  die "can generate ensembles only for PDB files" 
    if (!defined $pdbtemplate);

  my $tmol=Molecule::new();
  $tmol->readPDB($pdbtemplate);

  my $ens=Ensemble->new($enstag,$ensdir);
  $ens->readFileList();
  my $at=(defined $ensat)?$ensat:$ens->{par}->{runs}+1;
  
  for (my $i=$from; $i<=$to && $i<=$count; $i+=$step) {
    my $data=&GenUtil::readArchiveFile($arfile,$i);
    $ens->setFileList($at,"$arfile/$i");

    my $start=0;
    foreach my $c ( @{$tmol->{chain}} ) {
      foreach my $a ( @{$c->{atom}} ) {
	$a->{xcoor}=substr($data,$start,8)+0.0;
	$a->{ycoor}=substr($data,$start+8,8)+0.0;
	$a->{zcoor}=substr($data,$start+16,8)+0.0;
	$start+=24;
      }
    }

    $ens->checkinPDB($at++,$tmol,undef,"");
  }
  $ens->save();
} elsif ((defined $rmsmode || defined $qscore) && defined $ref) {
  die "can apply rms/qscore command only for PDB files" 
    if (!defined $pdbtemplate);
  
  my $refmol=Molecule::new($ref);
  my $analyze=Analyze::new($refmol);
  my $cmpmol=Molecule::new();
  $cmpmol->readPDB($pdbtemplate);

  for (my $i=$from; $i<=$to; $i+=$step) {
    my $data=&GenUtil::readArchiveFile($arfile,$i);

    my $start=0;
    foreach my $c ( @{$cmpmol->{chain}} ) {
      foreach my $a ( @{$c->{atom}} ) {
	$a->{xcoor}=substr($data,$start,8)+0.0;
	$a->{ycoor}=substr($data,$start+8,8)+0.0;
	$a->{zcoor}=substr($data,$start+16,8)+0.0;
	$start+=24;
      }
    }
    
    if (defined $rmsmode) {
      $analyze->lsqfit($cmpmol,"cab",0,1);
      my $rmsd=$analyze->rmsd($cmpmol,0,undef,1);
      printf "%i %f\n",$i,$rmsd->{$rmsmode};
    } elsif (defined $qscore) {
      my $qsc=$analyze->qscore($cmpmol,1);
      printf "%i %f %f %f %f\n",$i,$qsc->{all},$qsc->{short},$qsc->{medium},$qsc->{long};
    }
  }
} else {
  my @pdbtemp;
  if (defined $pdbtemplate) {
    my $templ=&GenUtil::getInputFile($pdbtemplate);
    @pdbtemp=<$templ>;
    close $templ;
  }

  if (defined $apply) {
    die "can apply commands only for PDB files" 
      if (!defined $pdbtemplate);

    for (my $i=$from; $i<=$to; $i+=$step) {
      my $data=&GenUtil::readArchiveFile($arfile,$i);
      open OUT,"|$apply";
      my $start=0;
      foreach my $line (@pdbtemp) {
	if ($line=~/ATOM/) {
	  substr($line,30,24)=substr($data,$start,24);
	  $start+=24;
	}
	print OUT $line;
      }
      close OUT;
    }
  } elsif ($xtract) {
    my $refmol;
    my $analyze;
    my $cmpmol;
    if ($center || (defined $fitref && -r $fitref)) {
      die "can center/fit command only for PDB files" 
	if (!defined $pdbtemplate);
      if (defined $fitref && -r $fitref) {
	$refmol=Molecule::new($fitref);
	$analyze=Analyze::new($refmol);
      } else {
        undef $fitref;
      }
      $cmpmol=Molecule::new();
      $cmpmol->readPDB($pdbtemplate);
    }

    for (my $i=$from; $i<=$to && $i<=$count; $i+=$step) {
      my $data=&GenUtil::readArchiveFile($arfile,$i);
      my $fname=sprintf("%s%d.%s",$prefix,$i,(defined $pdbtemplate)?"pdb":"dat");

      if ($center || defined $fitref) {
	my $start=0;
	foreach my $c ( @{$cmpmol->{chain}} ) {
	  foreach my $a ( @{$c->{atom}} ) {
	    $a->{xcoor}=substr($data,$start,8)+0.0;
	    $a->{ycoor}=substr($data,$start+8,8)+0.0;
	    $a->{zcoor}=substr($data,$start+16,8)+0.0;
	    $start+=24;
	  }
	}
    
	$cmpmol->center() if ($center);
        $analyze->lsqfit($cmpmol,"cab",0,1) if (defined $fitref);
	$cmpmol->writePDB($fname);
      } else {
	open OUT,">$fname";
	if (defined $pdbtemplate) {
	  my $start=0;
	  foreach my $line (@pdbtemp) {
	    if ($line=~/ATOM/) {
	      substr($line,30,24)=substr($data,$start,24);
	      $start+=24;
	    }
	    print OUT $line;
	  }
	} else {
	  syswrite(OUT,$data,length($data));
	}
	close OUT;
      }
      printf "%s\n",$fname if ($listfiles);
    }
  } else {
    my $data=&GenUtil::readArchiveFile($arfile,$from);
    if (defined $pdbtemplate) {
      my $start=0;
      foreach my $line (@pdbtemp) {
	if ($line=~/ATOM/) {
	  substr($line,30,24)=substr($data,$start,24);
	  $start+=24;
	}
	print $line;
      }
    } else {
      syswrite(STDOUT,$data,length($data));
    }
  }
}
