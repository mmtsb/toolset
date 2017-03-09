#!/usr/bin/env perl

# generates trajectory with random vectors to 
# be rotated with frame for rotational correlation 
# time analysis
#
# 2017, Michael Feig

sub usage {
  printf STDERR "usage: vecdcd.pl nvec nframes dcdfile\n";
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

my $nvec=shift @ARGV;
my $nframes=shift @ARGV;
my $dcdname=shift @ARGV;

my @xyz=();
my $i=0;
do {
  my $x1=rand(2.0)-1.0;
  my $x2=rand(2.0)-1.0;
  my $x12=$x1*$x1;
  my $x22=$x2*$x2;

  if ($x12+$x22<1) {
    my $trec={};
    $trec->{x}=2.0*$x1*sqrt(1.0-$x12-$x22);
    $trec->{y}=2.0*$x2*sqrt(1.0-$x12-$x22);
    $trec->{z}=1.0-2.0*($x12+$x22);
    push(@xyz,$trec);
    $i++;
  }
} while ($i<$nvec);

my $mol=Molecule::new();
$mol->fromXYZ(\@xyz,"CES","CES");

my ($xbuf,$ybuf,$zbuf);
foreach my $c ( @{$mol->{chain}} ) {
  $xbuf.=pack("f*",@{$c->{xcoor}});
  $ybuf.=pack("f*",@{$c->{ycoor}});
  $zbuf.=pack("f*",@{$c->{zcoor}});
}

my $out=&GenUtil::getOutputFile($dcdname);
binmode $out;

&writeFortran($out,pack("A4L*","CORD",@{&getIControl(nfiles=>$nframes)}));
&writeFortran($out,pack("LA80A80",2,"* TITLE","* CREATED by vecdcd.pl"));
&writeFortran($out,pack("L",$nvec));

for (my $i=0; $i<$nframes; $i++) {
  &writeFortran($out,$xbuf);
  &writeFortran($out,$ybuf);
  &writeFortran($out,$zbuf);
}

undef $out;

$mol->writePDB("-");

exit 0;

sub getIControl {
  my %par=@_;
  
  $par{nfiles}=1 if (!defined $par{nfiles});
  $par{istep}=1 if (!defined $par{istep});
  $par{interval}=1 if (!defined $par{interval});
  $par{nstep}=$par{nfiles}*$par{interval} if (!defined $par{nstep});
  $par{nsavv}=0 if (!defined $par{nsavv});
  $par{ndegf}=0 if (!defined $par{ndegf});
  $par{delta}=1.0 if (!defined $par{delta});
  $par{delta}*=20.45482667;
  $par{version}=30 if (!defined $par{version});

  my @icontrol;
  $icontrol[0]=$par{nfiles};
  $icontrol[1]=$par{istep};
  $icontrol[2]=$par{interval};
  $icontrol[3]=$par{nstep};
  $icontrol[4]=$par{nsavv};
  $icontrol[5]=0;
  $icontrol[6]=0;
  $icontrol[7]=$par{ndegf};
  $icontrol[8]=0;
  $icontrol[9]=unpack("L",pack("f",$par{delta}));
  $icontrol[10]=0;
  $icontrol[11]=0;
  $icontrol[12]=0;
  $icontrol[13]=0;
  $icontrol[14]=0;
  $icontrol[15]=0;
  $icontrol[16]=0;
  $icontrol[17]=0;
  $icontrol[18]=0;
  $icontrol[19]=$par{version};
  return \@icontrol;
}

sub writeFortran {
  my $handle=shift;
  my $buffer=shift;
  my $len=length($buffer);

  syswrite($handle,pack("L",$len),4);
  syswrite($handle,$buffer,$len);
  syswrite($handle,pack("L",$len),4);
}

