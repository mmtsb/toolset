#!/usr/bin/env perl

# read DCD file
#
#

sub usage {
  printf STDERR "usage:    processDCD.pl [template] [dcdfile(s)]\n";
  printf STDERR "options:  [-inx index[:to]] [-step n]\n";
  printf STDERR "          [-multi from:to]\n";
  printf STDERR "          [-apply cmd]\n";
  printf STDERR "          [-function file]\n";
  printf STDERR "          [-extract name]\n";
  printf STDERR "          [-ensdir dir] [-ens tag]\n";
  printf STDERR "          [-rms CA|CAB|C|O|N|side|back|all ref] [-useseg]\n";
  printf STDERR "          [-rmsl min:max[...]]\n";
  printf STDERR "          [-qscore ref] [-boxsize]\n";
  printf STDERR "          [-wrapseg]\n";
  printf STDERR "          [-average] [-fit] [-fitsel cab|ca|cb|heavy] [-fitresnumonly]\n";
  printf STDERR "          [-fitl min:max[...]]\n";
  printf STDERR "          [-ref ref]\n";
  printf STDERR "          [-psf file]\n";
  printf STDERR "          [-tag value]\n";
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
my $to=99999999;
my $step=1;
my $pdbtemplate;
my @dcdfiles=();
my $apply;
my $ffile;
my $rmsmode;
my $ref;
my $extract;
my $qscore;
my $avg=0;
my $psffile;
my $enstag;
my $ensdir=".";
my $multi;
my $boxa;
my $boxb;
my $boxc;
my $starttag;

my $frames;

my $warn=1;
my $resnumonly=undef;

my $selmode="cab";

my $boxsize=0;
my $wrapseg=0;

my $useseg=0;
my $lsqfit=0;
my $fraglist;
my $fitfraglist;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-inx") {
    shift @ARGV;
    my @f=split(/:/,shift @ARGV);
    $from=$f[0];
    $to=(defined $f[1])?$f[1]:$from;
    $frames=sprintf("%d-%d",$from,$to);
  } elsif ($ARGV[0] eq "-step") {
    shift @ARGV;
    $step=shift @ARGV;
  } elsif ($ARGV[0] eq "-apply") {
    shift @ARGV;
    $apply=shift @ARGV;
  } elsif ($ARGV[0] eq "-tag") {
    shift @ARGV;
    $starttag=shift @ARGV;
  } elsif ($ARGV[0] eq "-function") {
    shift @ARGV;
    $ffile=shift @ARGV;
  } elsif ($ARGV[0] eq "-rms") {
    shift @ARGV;
    $rmsmode=shift @ARGV;
    $ref=shift @ARGV;
  } elsif ($ARGV[0] eq "-useseg") {
    shift @ARGV;
    $useseg=1;
  } elsif ($ARGV[0] eq "-fit") {
    shift @ARGV;
    $lsqfit=1;
  } elsif ($ARGV[0] eq "-rmsl") {
    shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-fitl") {
    shift @ARGV;
    $fitfraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-multi") {
    shift @ARGV;
    $multi=shift @ARGV;
  } elsif ($ARGV[0] eq "-ensdir") {
    shift @ARGV;
    $ensdir=shift @ARGV;
  } elsif ($ARGV[0] eq "-ens") {
    shift @ARGV;
    $enstag=shift @ARGV;
  } elsif ($ARGV[0] eq "-qscore") {
    shift @ARGV;
    $qscore=1;
    $ref=shift @ARGV;
  } elsif ($ARGV[0] eq "-psf") {
    shift @ARGV;
    $psffile=shift @ARGV;
  } elsif ($ARGV[0] eq "-average") {
    shift @ARGV;
    $avg=1;
  } elsif ($ARGV[0] eq "-boxsize") {
    shift @ARGV;
    $boxsize=1;
  } elsif ($ARGV[0] eq "-wrapseg") {
    shift @ARGV;
    $wrapseg=1;
  } elsif ($ARGV[0] eq "-ref") {
    shift @ARGV;
    $ref=shift @ARGV;
  } elsif ($ARGV[0] eq "-fitsel") {
    shift @ARGV;
    $selmode=shift @ARGV;
  } elsif ($ARGV[0] eq "-fitresnumonly") {
    shift @ARGV;
    $resnumonly=1;
  } elsif ($ARGV[0] eq "-extract") {
    shift @ARGV;
    $extract=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    if (!defined $psffile && !defined $pdbtemplate && !$boxsize) {
      $pdbtemplate=shift @ARGV;
    } else {
      push(@dcdfiles,(shift @ARGV));
    }
  }
}

require "$ffile" if (defined $ffile);

my $nn=0;
my $mfrom=1;
my $mto=1;
if (defined $multi) {
  ($mfrom,$mto)=split(/:/,$multi);
}

my $fromframe=0;

my $ensat=-1;

my $ens;
    
if (defined $enstag) {
  $ens=Ensemble->new($enstag,$ensdir);
  $ensat=$ens->{par}->{runs}+1 if ($ensat<0);
} 

my $navg=0;
my $cmpmol;
my $avgmol;
my $refmol;
my $analyze;

if (defined $extract || $avg || defined $enstag) {
  $cmpmol=Molecule::new();
  if (defined $psffile) {
    $cmpmol->readPSF($psffile);
  } else {
    $cmpmol->readPDB($pdbtemplate);
  }

  if ($avg) {
    $avgmol=Molecule::new();
    if (defined $psffile) {
      $avgmol->readPSF($psffile);
    } else {
      $avgmol->readPDB($pdbtemplate);
    }
  }

  if (defined $ref) {
    $refmol=Molecule::new($ref);
    $analyze=Analyze::new($refmol);
  }

  foreach my $c ( @{$avgmol->{chain}} ) {
    foreach my $a ( @{$c->{atom}} ) {
      $a->{xcoor}=0.0;
      $a->{ycoor}=0.0;
      $a->{zcoor}=0.0;
    }
  }
}

my $i;

$from=1 if ($from<0);  #$nfiles+$from+1 if ($from<0);
$to=9999999 if ($to<0); # $nfiles+$to+1 if ($to<0);

my $itot=0;

&start($starttag) if (defined &start && defined $ffile);

foreach my $dcd ( @dcdfiles ) {
 for (my $im=$mfrom; $im<=$mto; $im++) {
  my $dcdfile;
  
  my $fdcd;
  if (defined $multi) {
    $fdcd=sprintf("$dcd",$im);
#    printf STDERR "reading from %s\n",$fdcd;
  } else {
    $fdcd=$dcd;
#    printf STDERR "reading from %s\n",$fdcd;
  }

  $dcdfile=&GenUtil::getInputFile($fdcd);

  binmode $dcdfile;

  my $buffer;
  my $len;
  ($buffer,$len)=&GenUtil::readFortran($dcdfile);
  my ($tag,@icontrol)=unpack("A4L*",$buffer);
  
  ($buffer,$len)=&GenUtil::readFortran($dcdfile);
  ($buffer,$len)=&GenUtil::readFortran($dcdfile);
  my $natom=unpack("L",$buffer);
  
  my $tstep=unpack("f",pack("L",$icontrol[9]))*4.88882129E-02;
  my $nfiles=$icontrol[0];
  my $first=$icontrol[1];
  my $delta=$icontrol[2];
  my $deltat=$icontrol[2]*$tstep;
  my $crystal=$icontrol[10];
  my $fixed=$icontrol[8];

  my $firstframe=$first/$delta;

  my @fixedinx=();
  if ($fixed>0) {
    my ($ibuf,$ilen)=&GenUtil::readFortran($dcdfile);
    @fixedinx=unpack("i*",$ibuf);
  }

  my ($xbuf,$ybuf,$zbuf);

  my @xcoor;
  my @ycoor;
  my @zcoor;

  my ($tbuf,$tlen);
  my @txcoor;
  my @tycoor;
  my @tzcoor;
  my $aa;
  my $ch;
  my $at;

  if (defined $extract || $avg || defined $enstag || ($boxsize && $crystal)) {
    
    foreach $ch ( @{$cmpmol->{chain}} ) {
      foreach $at ( @{$ch->{atom}} ) {
	$at->{xcoor}=0.0;
	$at->{ycoor}=0.0;
	$at->{zcoor}=0.0;
      }
    }

    for ($i=1; $itot<=$to && $i<=$nfiles; $i++) {
      $itot++;
    
      my $a=0;
      my $b=0;
      my $c=0; 
       
      if ($crystal) {
	($tbuf,$tlen)=&GenUtil::readFortran($dcdfile);
        if ($boxsize || $wrapseg) {
          my @cdat=unpack("d*",$tbuf);

	  $a=sqrt($cdat[0]*$cdat[0]+$cdat[1]*$cdat[1]+$cdat[3]*$cdat[3]);
	  $b=sqrt($cdat[1]*$cdat[1]+$cdat[2]*$cdat[2]+$cdat[4]*$cdat[4]);
	  $c=sqrt($cdat[3]*$cdat[3]+$cdat[4]*$cdat[4]+$cdat[5]*$cdat[5]);

	  my $ab=$cdat[1]*($cdat[0]+$cdat[2])+$cdat[3]*$cdat[4];
	  my $bc=$cdat[4]*($cdat[2]+$cdat[5])+$cdat[1]*$cdat[3];
	  my $ca=$cdat[3]*($cdat[0]+$cdat[5])+$cdat[1]*$cdat[4];

	  my $a1=&GenUtil::acos($bc/($b*$c))*180.0/$GenUtil::pi;
	  my $a2=&GenUtil::acos($ca/($c*$a))*180.0/$GenUtil::pi;
	  my $a3=&GenUtil::acos($ab/($a*$b))*180.0/$GenUtil::pi;
  
          if ($itot>=$from && $boxsize) { 
            printf "%f %f %f - %f %f %f\n",$a,$b,$c,$a1,$a2,$a3;
          }
        }
      }
      
      ($xbuf,$len)=&GenUtil::readFortran($dcdfile); # printf STDERR "%d ",$len;
      ($ybuf,$len)=&GenUtil::readFortran($dcdfile); # printf STDERR "%d ",$len;
      ($zbuf,$len)=&GenUtil::readFortran($dcdfile); # printf STDERR "%d \n",$len;

      if ($fixed>0 && $i==1) {
	@xcoor=unpack("f*",$xbuf);
	@ycoor=unpack("f*",$ybuf);
	@zcoor=unpack("f*",$zbuf);
      }
      
      if ($itot>=$from && $itot<=$to && ($itot%$step)==0 && !$boxsize) {
	if ($fixed>0 && $i>1) {
	  @txcoor=unpack("f*",$xbuf);
	  @tycoor=unpack("f*",$ybuf);
	  @tzcoor=unpack("f*",$zbuf);
	  
	  for (my $in=0; $in<=$#fixedinx; $in++) {
	    $xcoor[$fixedinx[$in]-1]=$txcoor[$in];
	    $ycoor[$fixedinx[$in]-1]=$tycoor[$in];
	    $zcoor[$fixedinx[$in]-1]=$tzcoor[$in];
	  }
	} elsif ($fixed==0) {
	  @xcoor=unpack("f*",$xbuf);
	  @ycoor=unpack("f*",$ybuf);
	  @zcoor=unpack("f*",$zbuf);
	}

	my $start=0;
	foreach $ch ( @{$cmpmol->{chain}} ) {
          $ch->{avgx}=0.0;
          $ch->{avgy}=0.0;
          $ch->{avgz}=0.0;
          $ch->{navg}=0;
	  foreach $at ( @{$ch->{atom}} ) {
	    $at->{xcoor}=$xcoor[$start];
	    $at->{ycoor}=$ycoor[$start];
	    $at->{zcoor}=$zcoor[$start];
	    $start++;
        
            $ch->{avgx}+=$at->{xcoor};
            $ch->{avgy}+=$at->{ycoor};
            $ch->{avgz}+=$at->{zcoor};
            $ch->{navg}++;
	  }
          $ch->{avgx}/=$ch->{navg};
          $ch->{avgy}/=$ch->{navg};
          $ch->{avgz}/=$ch->{navg};
	}

        if ($wrapseg) {
          my $chx=$cmpmol->{chain}->[0]->{avgx};
          my $chy=$cmpmol->{chain}->[0]->{avgy};
          my $chz=$cmpmol->{chain}->[0]->{avgz};
          for (my $ich=1; $ich<=$#{$cmpmol->{chain}}; $ich++) {
            my $ch=$cmpmol->{chain}->[$ich];
            my $dx=$ch->{avgx}-$chx;
            my $dy=$ch->{avgy}-$chy;
            my $dz=$ch->{avgz}-$chz;
            my $delx=(($dx>0)?int($dx/$a+0.5):int($dx/$a-0.5))*$a; 
            my $dely=(($dy>0)?int($dy/$b+0.5):int($dy/$b-0.5))*$b; 
            my $delz=(($dz>0)?int($dz/$c+0.5):int($dz/$c-0.5))*$c; 
	    foreach $at ( @{$ch->{atom}} ) {
	      $at->{xcoor}-=$delx;
	      $at->{ycoor}-=$dely;
	      $at->{zcoor}-=$delz;
            }
          }
        }

	if (defined $refmol) {
	  $analyze->lsqfit($cmpmol,$selmode,$warn,$resnumonly);
	}
	
	if (!$avg) {
	  if (defined $enstag) {
	    $ens->checkinPDB($ensat++,$cmpmol,undef,"");
	  } else {
	    my $fname=sprintf("%s.%d.pdb",$extract,++$nn);
	    $cmpmol->writePDB($fname);
	  }
	} else {
	  for (my $ic=0; $ic<=$#{$cmpmol->{chain}}; $ic++) {
	    my $c=$cmpmol->{chain}->[$ic];
	    my $ac=$avgmol->{chain}->[$ic];
	    for (my $ia=0; $ia<=$#{$c->{atom}}; $ia++) {
	      $a=$c->{atom}->[$ia];
	      $aa=$ac->{atom}->[$ia];
	      
	      $aa->{xcoor}+=$a->{xcoor};
	      $aa->{ycoor}+=$a->{ycoor};
	      $aa->{zcoor}+=$a->{zcoor};
	    }
	  }
	  
	  $navg++;
	}
      }
    }
  } elsif ((defined $rmsmode || defined $qscore) && defined $ref) {
    die "can apply rms/qscore command only for PDB files" 
      if (!defined $pdbtemplate && !defined $psffile);
    
    my $refmol=Molecule::new($ref);
    my $analyze=Analyze::new($refmol);
    my $cmpmol=Molecule::new();
    if (defined $psffile) {
      $cmpmol->readPSF($psffile);
    } else {
      $cmpmol->readPDB($pdbtemplate);
    }


    for ($i=1; $itot<=$to && $i<=$nfiles; $i++) {
      $itot++;

      my $a=0;
      my $b=0;
      my $c=0; 
      if ($crystal) {
	my ($tbuf,$tlen)=&GenUtil::readFortran($dcdfile);
        if ($wrapseg) {
          my @cdat=unpack("d*",$tbuf);
	  $a=sqrt($cdat[0]*$cdat[0]+$cdat[1]*$cdat[1]+$cdat[3]*$cdat[3]);
	  $b=sqrt($cdat[1]*$cdat[1]+$cdat[2]*$cdat[2]+$cdat[4]*$cdat[4]);
	  $c=sqrt($cdat[3]*$cdat[3]+$cdat[4]*$cdat[4]+$cdat[5]*$cdat[5]);
        }
      }
      
      ($xbuf,$len)=&GenUtil::readFortran($dcdfile); #printf STDERR "%d ",$len;
      ($ybuf,$len)=&GenUtil::readFortran($dcdfile); #printf STDERR "%d ",$len;
      ($zbuf,$len)=&GenUtil::readFortran($dcdfile); #printf STDERR "%d \n",$len;

      if ($fixed>0 && $i==1) {
	@xcoor=unpack("f*",$xbuf);
	@ycoor=unpack("f*",$ybuf);
	@zcoor=unpack("f*",$zbuf);
      }

      if ($itot>=$from && $itot<=$to && ($itot%$step)==0) {
	if ($fixed>0 && $i>1) {
	  my @txcoor=unpack("f*",$xbuf);
	  my @tycoor=unpack("f*",$ybuf);
	  my @tzcoor=unpack("f*",$zbuf);
	  
	  for (my $in=0; $in<=$#fixedinx; $in++) {
	    $xcoor[$fixedinx[$in]-1]=$txcoor[$in];
	    $ycoor[$fixedinx[$in]-1]=$tycoor[$in];
	    $zcoor[$fixedinx[$in]-1]=$tzcoor[$in];
	  }
	} elsif ($fixed==0) {
	  @xcoor=unpack("f*",$xbuf);
	  @ycoor=unpack("f*",$ybuf);
	  @zcoor=unpack("f*",$zbuf);
	}

	my $start=0;
	foreach $ch ( @{$cmpmol->{chain}} ) {
          $ch->{avgx}=0.0;
          $ch->{avgy}=0.0;
          $ch->{avgz}=0.0;
          $ch->{navg}=0;
	  foreach $at ( @{$ch->{atom}} ) {
	    $at->{xcoor}=$xcoor[$start];
	    $at->{ycoor}=$ycoor[$start];
	    $at->{zcoor}=$zcoor[$start];
	    $start++;
        
            $ch->{avgx}+=$at->{xcoor};
            $ch->{avgy}+=$at->{ycoor};
            $ch->{avgz}+=$at->{zcoor};
            $ch->{navg}++;
	  }
          $ch->{avgx}/=$ch->{navg};
          $ch->{avgy}/=$ch->{navg};
          $ch->{avgz}/=$ch->{navg};
	}

        if ($wrapseg) {
          my $chx=$cmpmol->{chain}->[0]->{avgx};
          my $chy=$cmpmol->{chain}->[0]->{avgy};
          my $chz=$cmpmol->{chain}->[0]->{avgz};
          for (my $ich=1; $ich<=$#{$cmpmol->{chain}}; $ich++) {
            my $ch=$cmpmol->{chain}->[$ich];
            my $dx=$ch->{avgx}-$chx;
            my $dy=$ch->{avgy}-$chy;
            my $dz=$ch->{avgz}-$chz;
            my $delx=(($dx>0)?int($dx/$a+0.5):int($dx/$a-0.5))*$a; 
            my $dely=(($dy>0)?int($dy/$b+0.5):int($dy/$b-0.5))*$b; 
            my $delz=(($dz>0)?int($dz/$c+0.5):int($dz/$c-0.5))*$c; 
	    foreach $at ( @{$ch->{atom}} ) {
	      $at->{xcoor}-=$delx;
	      $at->{ycoor}-=$dely;
	      $at->{zcoor}-=$delz;
            }
          }
        }

	if (defined $rmsmode) {
          $fitfraglist=$fraglist
            if (defined $fraglist && !defined $fitfraglist);

          $cmpmol->setValidResidues($fitfraglist,0) if (defined $fitfraglist);
          $analyze->lsqfit($cmpmol,$selmode,0,0,undef,$useseg) if ($lsqfit);
          $cmpmol->resetValidResidues() if (!defined $fraglist && defined $fitfraglist);
          $cmpmol->setValidResidues($fraglist,0) if (defined $fraglist);

          my $rmsd=$analyze->rmsd($cmpmol,0,undef,0,undef,$useseg);
          printf "%d %f %f\n",($itot),($itot)*$deltat,$rmsd->{$rmsmode};
	} elsif (defined $qscore) {
	  my $qsc=$analyze->qscore($cmpmol,1);
	  printf "%d %f %f %f %f %f\n",($itot),($itot)*$deltat,
	    $qsc->{all},$qsc->{short},$qsc->{medium},$qsc->{long};
	}
      }
    }
  } else {
    my $cmpmol=Molecule::new();
    if (defined $psffile) {
      $cmpmol->readPSF($psffile);
    } else {
      $cmpmol->readPDB($pdbtemplate);
    }
    
    if (defined $apply || defined $ffile) {
      die "can apply commands only for PDB files" 
	if (!defined $pdbtemplate && !defined $psffile);
      
      for ($i=1; $itot<=$to && $i<=$nfiles; $i++) {
	$itot++;

        my $a=0;
        my $b=0;
        my $c=0; 

	if ($crystal) {
	  my ($tbuf,$tlen)=&GenUtil::readFortran($dcdfile);
          my @cdat=unpack("d*",$tbuf);
          $a=sqrt($cdat[0]*$cdat[0]+$cdat[1]*$cdat[1]+$cdat[3]*$cdat[3]);
	  $b=sqrt($cdat[1]*$cdat[1]+$cdat[2]*$cdat[2]+$cdat[4]*$cdat[4]);
	  $c=sqrt($cdat[3]*$cdat[3]+$cdat[4]*$cdat[4]+$cdat[5]*$cdat[5]);
          $boxa=$a;
          $boxb=$b;
          $boxc=$c;
	}
	
	($xbuf,$len)=&GenUtil::readFortran($dcdfile); #printf STDERR "%d ",$len;
	($ybuf,$len)=&GenUtil::readFortran($dcdfile); #printf STDERR "%d ",$len;
	($zbuf,$len)=&GenUtil::readFortran($dcdfile); #printf STDERR "%d \n",$len;

	if ($fixed>0 && $i==1) {
	  @xcoor=unpack("f*",$xbuf);
	  @ycoor=unpack("f*",$ybuf);
	  @zcoor=unpack("f*",$zbuf);
	}

	if ($itot>=$from && $itot<=$to && ($itot%$step)==0) {	
	  if ($fixed>0 && $i>1) {
	    my @txcoor=unpack("f*",$xbuf);
	    my @tycoor=unpack("f*",$ybuf);
	    my @tzcoor=unpack("f*",$zbuf);
	    
	    for (my $in=0; $in<=$#fixedinx; $in++) {
	      $xcoor[$fixedinx[$in]-1]=$txcoor[$in];
	      $ycoor[$fixedinx[$in]-1]=$tycoor[$in];
	      $zcoor[$fixedinx[$in]-1]=$tzcoor[$in];
	    }
	  } elsif ($fixed==0) {
	    @xcoor=unpack("f*",$xbuf);
	    @ycoor=unpack("f*",$ybuf);
	    @zcoor=unpack("f*",$zbuf);
	  }
	
	my $start=0;
	foreach $ch ( @{$cmpmol->{chain}} ) {
          $ch->{avgx}=0.0;
          $ch->{avgy}=0.0;
          $ch->{avgz}=0.0;
          $ch->{navg}=0;
	  foreach $at ( @{$ch->{atom}} ) {
	    $at->{xcoor}=$xcoor[$start];
	    $at->{ycoor}=$ycoor[$start];
	    $at->{zcoor}=$zcoor[$start];
	    $start++;
        
            $ch->{avgx}+=$at->{xcoor};
            $ch->{avgy}+=$at->{ycoor};
            $ch->{avgz}+=$at->{zcoor};
            $ch->{navg}++;
	  }
          $ch->{avgx}/=$ch->{navg};
          $ch->{avgy}/=$ch->{navg};
          $ch->{avgz}/=$ch->{navg};
	}

        if ($wrapseg) {
          my $chx=$cmpmol->{chain}->[0]->{avgx};
          my $chy=$cmpmol->{chain}->[0]->{avgy};
          my $chz=$cmpmol->{chain}->[0]->{avgz};
          for (my $ich=1; $ich<=$#{$cmpmol->{chain}}; $ich++) {
            my $ch=$cmpmol->{chain}->[$ich];
            my $dx=$ch->{avgx}-$chx;
            my $dy=$ch->{avgy}-$chy;
            my $dz=$ch->{avgz}-$chz;
            my $delx=(($dx>0)?int($dx/$a+0.5):int($dx/$a-0.5))*$a; 
            my $dely=(($dy>0)?int($dy/$b+0.5):int($dy/$b-0.5))*$b; 
            my $delz=(($dz>0)?int($dz/$c+0.5):int($dz/$c-0.5))*$c; 
	    foreach $at ( @{$ch->{atom}} ) {
	      $at->{xcoor}-=$delx;
	      $at->{ycoor}-=$dely;
	      $at->{zcoor}-=$delz;
            }
          }
        }

          if (defined $apply) {	  
	    printf "%d %1.4f ",($i+$firstframe-1),($i+$firstframe-1)*$deltat;
	    $cmpmol->writePDB("tmp-$$");
	    system "cat tmp-$$ | $apply";
	    &GenUtil::remove("tmp-$$");
	  } elsif (defined $ffile) {
	    my @res=&analyze($cmpmol,$boxa,$boxb,$boxc);
            if (defined $res[0]) {
  	      printf "%d %1.4f ",($i+$firstframe-1),($i+$firstframe-1)*$deltat;
	      printf "%s\n",join(" ",@res);
            }
	  }
	}
      }
    }
  }
 }
}

$ens->save() if (defined $enstag);

if ($avg && $navg>0) {
  printf STDERR "navg: %d\n",$navg;
  foreach my $c ( @{$avgmol->{chain}} ) {
    foreach my $a ( @{$c->{atom}} ) {
      $a->{xcoor}/=$navg;
      $a->{ycoor}/=$navg;
      $a->{zcoor}/=$navg;
    }
  }
  $avgmol->writePDB();
}

&end() if (defined &end && defined $ffile);
  
1;
