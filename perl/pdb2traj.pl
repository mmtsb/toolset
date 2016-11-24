#!/usr/bin/env perl

# MMTSB Tool Set
#
# Michael Feig, 2005, MSU
#

sub usage {
  printf STDERR "usage:   pdb2traj.pl [options] [PDBfiles]\n";
  printf STDERR "options: [-f listfile]\n";
  printf STDERR "         [-out file]\n";
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
use Molecule;

my $listfile;
my @filelist;
my $outfile="traj.dcd";
my $sel=undef;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-f") {
    shift @ARGV;
    $listfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-out") {
    shift @ARGV;
    $outfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-nsel"){
    shift @ARGV;
    $sel=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option %s\n",$ARGV[0];
    &usage();
  } else {
    push(@filelist,shift @ARGV);
  }
}

if (defined $listfile) {
  my $linp=&GenUtil::getInputFile($listfile);
  while (<$linp>) {
    chomp;
    push(@filelist,$_);
  }
  undef $linp;
}

my $first=Molecule::new($filelist[0]);
if (defined $sel){
  $first->setValidSelection($sel);
  $first=$first->clone(1);
}
		  
my $natom=0;
foreach my $c ( @{$first->{chain}} ) {
  $natom+=$#{$c->{atom}}+1;
}

my $out=&GenUtil::getOutputFile($outfile);

binmode $out;

&writeFortran($out,pack("A4L*","CORD",@{&getIControl(nfiles=>$#filelist+1)}));
&writeFortran($out,pack("LA80A80",2,"* TITLE","* CREATED by pdb2traj.pl"));
&writeFortran($out,pack("L",$natom));

foreach my $f ( @filelist ) {
#  printf STDERR "$f\n";
  my $mol=Molecule::new($f);
  if (defined $sel){
    $mol->setValidSelection($sel);
    $mol=$mol->clone(1);
  }
		    
  my ($xbuf,$ybuf,$zbuf);
  foreach my $c ( @{$mol->{chain}} ) {
    $xbuf.=pack("f*",@{$c->{xcoor}});
    $ybuf.=pack("f*",@{$c->{ycoor}});
    $zbuf.=pack("f*",@{$c->{zcoor}});
  }
  &writeFortran($out,$xbuf);
  &writeFortran($out,$ybuf);
  &writeFortran($out,$zbuf);
}

undef $out;

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


