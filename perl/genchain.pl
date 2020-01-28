#!/usr/bin/env perl

# generates SICHO chain
#
# http://mmtsb.scripps.edu/doc/genchain.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   genchain.pl [options] [[-m | -s | -pdb | -primo | -primo2 | -htrna] file] | [-rnd num]\n";
  printf STDERR "options: [-r resolution] [-g gridsize]\n";
  printf STDERR "         [-float] [-center] [-ca]\n";
  printf STDERR "         [-o offsetx offsety offsetz]\n";
  printf STDERR "         [-l min:max[=min:max=...]]\n";
  printf STDERR "         [-setres num:name[,num:name]]\n";
  printf STDERR "         [-seq seqfile]\n";
  printf STDERR "         [-dcd dcdinp dcdout]\n";
  printf STDERR "         [-nsel selection]\n";
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

my $resolution=1.45;
my $gridsize=100;
my $offsetx=50.0;
my $offsety=50.0;
my $offsetz=50.0;
my $fraglist;
my $mode="monsster";
my $filename="";
my $randomnum=-1;
my $center;
my $wantca=0;
my $intflag=1;
my $newresnames;
my $seqfile;
my $primo=0;
my $primo2=0;
my $dcdinp=undef;
my $dcdout=undef;
my $sel=undef;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-m") {
    shift @ARGV;
    $center=1 if (!defined $center);
    $mode="monsster";
  } elsif ($ARGV[0] eq "-rnd") {
    shift @ARGV;
    $center=1 if (!defined $center);
    $mode="random";
    $randomnum=shift @ARGV;
  } elsif ($ARGV[0] eq "-s") {
    shift @ARGV;
    $center=0 if (!defined $center);
    $mode="simple";
  } elsif ($ARGV[0] eq "-pdb") {
    shift @ARGV;
    $center=0 if (!defined $center);
    $mode="pdb";
    $resolution=0;
    $offsetx=$offsety=$offsetz=0;
  } elsif ($ARGV[0] eq "-primo") {
    shift @ARGV;
    $mode="primo";
  } elsif ($ARGV[0] eq "-primo2") {
    shift @ARGV;
    $mode="primo2";
  } elsif ($ARGV[0] eq "-htrna") {
    shift @ARGV;
    $mode="htrna";
  } elsif ($ARGV[0] eq "-r") {
    shift @ARGV;
    $resolution=shift @ARGV;
  } elsif ($ARGV[0] eq "-g") {
    shift @ARGV;
    $gridsize=shift @ARGV;
    $offsetx=$offsety=$offsetz=int($gridsize/2);
  } elsif ($ARGV[0] eq "-center") {
    shift @ARGV;
    $center=1;
  } elsif ($ARGV[0] eq "-nocenter") {
    shift @ARGV;
    $center=0;
  } elsif ($ARGV[0] eq "-float") {
    shift @ARGV;
    $intflag=0;
  } elsif ($ARGV[0] eq "-setres") {
    shift @ARGV;
    $newresnames=shift @ARGV;
  } elsif ($ARGV[0] eq "-seq") {
    shift @ARGV;
    $seqfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-o") {
    shift @ARGV;
    $offsetx=shift @ARGV;
    $offsety=shift @ARGV;
    $offsetz=shift @ARGV;
  } elsif ($ARGV[0] eq "-l") {
    shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
    $center=0;
  } elsif ($ARGV[0] eq "-ca") {
    shift @ARGV;
    $wantca=1;
  } elsif ($ARGV[0] eq "-dcd") {
    shift @ARGV;
    $dcdinp=shift @ARGV;
    $dcdout=shift @ARGV;
  } elsif ($ARGV[0] eq "-nsel"){
    shift @ARGV;
    $sel=shift @ARGV;
  } elsif ($ARGV[0] =~/^-/) {
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    $filename=shift @ARGV;
  }    
}

if ($mode eq "none") {
  printf STDERR "Unknown mode: don't know what to do!\n";
  &usage();
}


my $mol;
my $atoms;
$atoms->{max}=0;
my $a;
my $c;
if ($mode eq "primo" || $mode eq "primo2" || $mode eq "htrna") {
  $mol=Molecule::new();
  $mol->readPDB($filename);
  if (defined $sel){
    $mol->setValidSelection($sel);
    $mol=$mol->clone(1);
  }
  for $c (@{$mol->activeChains()}){
      foreach $a (@{$c->{atom}}){
	  $atoms->{$a->{atominx}}=1;
	  $atoms->{max}=$a->{atominx} if ($a->{atominx} > $atoms->{max});
      }
  }
  my $cg;
  if ($mode eq "primo") {
    $cg=$mol->genPRIMO();
  } elsif ($mode eq "primo2") {
    $cg=$mol->genPRIMO2();
  } elsif ($mode eq "htrna") {
    $cg=$mol->genHTRNA();
  } else {
    die "unknown mode\n";
  }
  my $natom=0;
  for $c (@{$cg->activeChains()}){
    foreach $a (@{$c->{atom}}){
      $natom++;
    }
  }
  if (!defined $dcdinp && !defined $dcdout){
    $cg->writePDB(\*STDOUT);
  }
  else{
    #Read All-Atom DCD Header
    $dcdinp=&GenUtil::getInputFile($dcdinp);
    binmode $dcdinp;
    my $buffer;
    my $len;
    my ($xbuf, $ybuf, $zbuf);
    ($buffer,$len)=&GenUtil::readFortran($dcdinp);
    my ($tag,@icontrol)=unpack("A4L*",$buffer);

    ($buffer,$len)=&GenUtil::readFortran($dcdinp);
    ($buffer,$len)=&GenUtil::readFortran($dcdinp);
    my $NATOM=unpack("L",$buffer);
    
    my $tstep=unpack("f",pack("L",$icontrol[9]))*4.88882129E-02;
    my $nfiles=$icontrol[0];
    my $first=$icontrol[1];
    my $delta=$icontrol[2];
    my $deltat=$icontrol[2]*$tstep;
    my $crystal=$icontrol[10];
    my $fixed=$icontrol[8];

    $icontrol[7]=0; #DEGREES OF FREEDOM
    $icontrol[10]=0; #Crystal, Periodic Boundaries
    $icontrol[11]=0; #4-Dimensions
    $icontrol[12]=0; #Fluctuation Charges

    my $firstframe=$first/$delta;
    
    #Write PRIMO DCD Header
    $dcdout=&GenUtil::getOutputFile($dcdout);
    binmode $dcdout;

    &writeFortran($dcdout,pack("A4L*","CORD",@icontrol));
    &writeFortran($dcdout,pack("LA80A80",2,"* TITLE","* CREATED by genchain.pl"));
    &writeFortran($dcdout,pack("L",$natom));
    
    for (my $i=1; $i<=$nfiles; $i++){
      #Read All-Atom DCD Coordinates
      if ($crystal){
	my ($tbuf,$tlen)=&GenUtil::readFortran($dcdinp);
      }
      if ($atoms->{max} == $NATOM){
	  ($xbuf,$len)=&GenUtil::readFortran($dcdinp);
	  ($ybuf,$len)=&GenUtil::readFortran($dcdinp);
	  ($zbuf,$len)=&GenUtil::readFortran($dcdinp); 
      }
      else{
	  ($xbuf,$len)=&readBinary($dcdinp,$atoms->{max},$NATOM-$atoms->{max});
	  ($ybuf,$len)=&readBinary($dcdinp,$atoms->{max},$NATOM-$atoms->{max});
	  ($zbuf,$len)=&readBinary($dcdinp,$atoms->{max},$NATOM-$atoms->{max});
      }
      my @xcoor=unpack("f*",$xbuf);
      my @ycoor=unpack("f*",$ybuf);
      my @zcoor=unpack("f*",$zbuf);
      my $start=0;
      for $c (@{$mol->{chain}}){
	  foreach $a (@{$c->{atom}}){
	      if ($atoms->{$a->{atominx}}){
		  $a->{xcoor}=$xcoor[$start];
		  $a->{ycoor}=$ycoor[$start];
		  $a->{zcoor}=$zcoor[$start];
	      }
	      elsif ($a->{atominx} > $atoms->{max}){
		  last;
	      }
	      $start++;
	  }
      }

      if ($mode eq "primo") {
        $cg=$mol->genPRIMO();
      } elsif ($mode eq "primo2") {
        $cg=$mol->genPRIMO2();
      } elsif ($mode eq "htrna") {
        $cg=$mol->genHTRNA();
      } else {
        die "unknown mode\n";
      }
      
      #Write PRIMO DCD Coordinates
      $xbuf=undef;
      $ybuf=undef;
      $zbuf=undef;
      foreach $c (@{$cg->{chain}}){
	  $xbuf.=pack("f*",@{$c->{xcoor}});
	  $ybuf.=pack("f*",@{$c->{ycoor}});
	  $zbuf.=pack("f*",@{$c->{zcoor}});
      }
      &writeFortran($dcdout,$xbuf);
      &writeFortran($dcdout,$ybuf);
      &writeFortran($dcdout,$zbuf);
    }
  }
} else {
  my $sicho=SICHO::new(gridsize => $gridsize,
		       offsetx  => $offsetx,
		       offsety  => $offsety,
		       offsetz  => $offsetz,
		       resolution => $resolution,
		       intflag    => $intflag);
  

  if ($mode eq "monsster") {
    $mol=Molecule::new();
    $mol->readPDB($filename);
    $mol->selectChain("");
    $mol->center() if ($center);
    $sicho->genMONSSTERFromAllAtom($mol, fraglist => $fraglist);
  } elsif ($mode eq "random") {
    $sicho->genRandomMONSSTER($randomnum);  
  } elsif ($mode eq "simple" || $mode eq "pdb") {
    $mol=Molecule::new();
    $mol->readPDB($filename);
    $mol->selectChain("");
    $mol->center() if ($center);
    $sicho->genSimpleFromAllAtom($mol,ca => $wantca);
  }

  if ($#{$sicho->{sidechain}}<0) {
    printf STDERR "Could not generate chain\n";
  } else {
    if ($mode eq "pdb") {
      my $outmol=Molecule::new();
      my $seq;
      if (defined $seqfile) {
	$seq=Sequence::new();
	$seq->readMONSSTER($seqfile);
      } else {
	$seq=Sequence::new($mol);
	$seq->modifyResidueName($newresnames)
	  if (defined $newresnames);
      }
      $outmol->fromSICHO($seq,$sicho);
      $outmol->writePDB(\*STDOUT,ssbond=>0);
    } else {
      $sicho->writeChain(\*STDOUT);
    }
  }
}

sub writeFortran {
  my $handle=shift;
  my $buffer=shift;
  my $len=length($buffer);

  syswrite($handle,pack("L",$len),4);
  syswrite($handle,$buffer,$len);
  syswrite($handle,pack("L",$len),4);
}
## function: readFortran(handle)

sub readBinary {
  my $handle=shift;
  my $maxatoms=shift;
  my $seekatoms=shift;

  $maxatoms=$maxatoms*4;
  $seekatoms=$seekatoms*4;

  my $dat;
  my $tdat;

  my $len;
  if ($GenUtil::fortranheadersize<0) {
    read($handle,$tdat,8) || die "cannot read data";
    if (substr($tdat,4,1) eq "C") {
       $GenUtil::fortranheadersize=4;
       $len=unpack("L",substr($tdat,0,4));
       read($handle,$dat,$len-4) || die "cannot read data";
       $dat=substr($tdat,4,4).$dat;
    } else {
       $GenUtil::fortranheadersize=8;
       $len=unpack("L",$tdat);
       read($handle,$dat,$len) || die "cannot read data";
    } 
  } else {
    read($handle,$tdat,$GenUtil::fortranheadersize) || die "cannot read data";
    $len=unpack("L",$tdat);
    read($handle,$dat,$maxatoms) || die "cannot read data";
    seek($handle,$seekatoms,1) || die "cannot read data";
  }
  read($handle,$tdat,$GenUtil::fortranheadersize) || die "cannot read data";
 
  return ($dat,$len);
}
