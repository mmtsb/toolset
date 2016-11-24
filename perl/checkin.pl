#!/usr/bin/env perl
#
# checks chain files or PDB files into Ensemble directory structure
#
# http://mmtsb.scripps.edu/doc/checkin.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   checkin.pl [options] tag [file(s)]\n";
  printf STDERR "options: [-sicho] [-seq file]\n";
  printf STDERR "         [-at index]\n";
  printf STDERR "         [-dir datadir]\n";
  printf STDERR "         [-f filelist]\n";
  printf STDERR "         [-links]\n";
  printf STDERR "         [-verbose]\n";
  printf STDERR "         [-[no]compress]\n";
  printf STDERR "         [-dcd PDBtemplate CHARMMtrajectory]\n";
  printf STDERR "         [-dcdskip value]\n";
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
use Ensemble;
use LatEnsemble;

my $tag;
my $at=-1;
my $filelist=();
my $dir=".";

my $sicho=0;

my $listfile;

my $dcd;
my $pdbtemplate;
my $dcdskip=1;

my $verbose=0;

my $seqfile;

my $links=0;

my $compress;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-sicho") {
    shift @ARGV;
    $sicho=1;
  } elsif ($ARGV[0] eq "-seq") {
    shift @ARGV;
    $seqfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-at") {
    shift @ARGV;
    $at=shift @ARGV;
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-verbose") {
    shift @ARGV;
    $verbose=1;
  } elsif ($ARGV[0] eq "-f") {
    shift @ARGV;
    $listfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-links") {
    shift @ARGV;
    $links=1;
 } elsif ($ARGV[0] eq "-nocompress") {
    shift @ARGV;
    $compress=0;
  } elsif ($ARGV[0] eq "-compress") {
    shift @ARGV;
    $compress=1;
  } elsif ($ARGV[0] eq "-dcd") {
    shift @ARGV;
    $pdbtemplate=shift @ARGV;
    $dcd=shift @ARGV;
  } elsif ($ARGV[0] eq "-dcdskip") {
    shift @ARGV;
    $dcdskip=shift @ARGV;
  } elsif ($ARGV[0] =~/^-.+/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    if (!defined $tag) {
      $tag=shift @ARGV;
    } else {
      push(@{$filelist},shift @ARGV);
    }
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

&usage() if (($#{$filelist}<0 && !defined $dcd) || !defined $tag);

my $ens=($sicho)?LatEnsemble->new($tag,$dir):Ensemble->new($tag,$dir);

$ens->set(compress=>$compress) if (defined $compress);
$ens->readFileList();

$at=$ens->{par}->{runs}+1 if ($at<0);

if (defined $dcd && &GenUtil::checkFile($dcd) &&
    defined $pdbtemplate && &GenUtil::checkFile($pdbtemplate)) {
  my $tmol=Molecule::new();
  $tmol->readPDB($pdbtemplate);

  my $dcdfile=&GenUtil::getInputFile($dcd);
  
  my $buffer;
  my $len;
  ($buffer,$len)=&GenUtil::readFortran($dcdfile);
  my ($tag,@icontrol)=unpack("A4L*",$buffer);

  ($buffer,$len)=&GenUtil::readFortran($dcdfile);
  ($buffer,$len)=&GenUtil::readFortran($dcdfile);
  my $natom=unpack("L",$buffer);
  
  my $nfiles=$icontrol[0];
  my $crystal=$icontrol[10];

#  printf STDERR "natom: $natom, nfiles: $nfiles\n";

  my ($xbuf,$ybuf,$zbuf);
  for (my $i=1; $i<=$nfiles; $i+=$dcdskip) {
    if ($crystal) {
	my ($tbuf,$tlen)=&GenUtil::readFortran($dcdfile);
    }
    ($xbuf,$len)=&GenUtil::readFortran($dcdfile);
    ($ybuf,$len)=&GenUtil::readFortran($dcdfile);
    ($zbuf,$len)=&GenUtil::readFortran($dcdfile);

    my @xcoor=unpack("f*",$xbuf);
    my @ycoor=unpack("f*",$ybuf);
    my @zcoor=unpack("f*",$zbuf);

    my $n=0;
    foreach my $c ( @{$tmol->{chain}} ) {
      foreach my $a ( @{$c->{atom}} ) {
	$a->{xcoor}=$xcoor[$n];
	$a->{ycoor}=$ycoor[$n];
	$a->{zcoor}=$zcoor[$n];
	$n++;
      }
    }
    $ens->setFileList($at,"$dcd/$i");
    $ens->checkinPDB($at++,$tmol,undef,"");
  }
} else {
  foreach my $f ( @{$filelist} ) {
    if ($verbose) {
      printf STDERR "checking in %s at %d\n",$f,$at;
    }
     
    $ens->setFileList($at,$f);
    if ($sicho) {
      $ens->set(seq=>$seqfile) if (defined $seqfile);
      my $chain=SICHO::new();
      $chain->readChain($f);
      $ens->checkinSICHO($at++,$chain);
    } else {
      my $mol;

      if ($links) {
        if ($f=~/^\//) {
	  $mol=$f;
	} else {
	  my $pwd=`pwd`;
          chomp $pwd;
	  $mol=$pwd."/".$f;
	}
      } else {
	$mol=Molecule::new($f);
      }
      if ($links || !$mol->empty()) {
        $ens->checkinPDB($at++,$mol,undef,"",$links);
      } else {
        printf STDERR "molecule from $f is empty and will be ignored\n";
      }
    }
  }
}

$ens->save();


