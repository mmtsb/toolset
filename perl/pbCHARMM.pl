#!/usr/bin/env perl

# get Poisson-Boltzmann energy from PDB file 
#
# http://mmtsb.scripps.edu/doc/pbCHARMM.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   pbCHARMM.pl [options] [PDBfile]\n";
  printf STDERR "options: [-par param=19|22,hsd=list,hse=list,scalerad,\n";
  printf STDERR "               smooth,dcel=value,epsp=value,epsw=value,epsr=value,\n";
  printf STDERR "               pbionconc=value,pbtemp=value,pbionr=value]\n";
  printf STDERR "         [-psf PSFfile CRDfile]\n";
  printf STDERR "         [-mol2 MOL2file]\n";
  printf STDERR "         [-radii file]\n";
  printf STDERR "         [-partial file] [-threads n]\n";
  printf STDERR "         [-atomic] [-pairs] [-keepcharge]\n";
  printf STDERR "         [-log logFile] [-cmd logFile]\n";
  printf STDERR "         [-custom file]\n";
  printf STDERR "         [-epsgrid file] [-epssize num]\n";
  printf STDERR "         [-grid phi|phix|chrg|epsx|epsy|epsz file] [-dx]\n";
  printf STDERR "         [-emap file]\n"; 
  printf STDERR "         [-nocenter]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use IO::File;
use IO::Handle;

use GenUtil;
use Molecule;
use CHARMM;

my %par;

my $logFile;
my $cmdlog;

my $inpfile="-";
my $base="";

my $single;
my $double;
my $keepcharge;
my $threads=1;

my $partial;

my $psffile;
my $crdfile;
my $mol2file;

my $customfile;

my $epsgrid;
my $epssize=12;

my $dx=0;
my $gridname;

my $emap;

my $center=1;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=(shift @ARGV);
  } elsif ($ARGV[0] eq "-atomic") {
    shift @ARGV;
    $single=1;
  } elsif ($ARGV[0] eq "-pairs") {
    shift @ARGV;
    $double=1;
  } elsif ($ARGV[0] eq "-nocenter") {
    shift @ARGV;
    $center=0;
  } elsif ($ARGV[0] eq "-radii") {
    shift @ARGV;
    $par{scalerad}="";
    open INP,shift @ARGV;
    while (<INP>) {
      chomp;
      $par{scalerad}.=":" if ($par{scalerad} ne "");
      $par{scalerad}.=$_;
    }
    close INP;
  } elsif ($ARGV[0] eq "-threads") {
    shift @ARGV;
    $threads=shift @ARGV;
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-partial") {
    shift @ARGV;
    $partial=shift @ARGV;
  } elsif ($ARGV[0] eq "-keepcharge") {
    shift @ARGV;
    $keepcharge=1;
  } elsif ($ARGV[0] eq "-psf") {
    shift @ARGV;
    $psffile=shift @ARGV;
    $crdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-mol2") {
    shift @ARGV;
    $mol2file=shift @ARGV;
  } elsif ($ARGV[0] eq "-epsgrid") {
    shift @ARGV;
    $epsgrid=shift @ARGV;
  } elsif ($ARGV[0] eq "-emap") {
    shift @ARGV;
    $emap=shift @ARGV;
  } elsif ($ARGV[0] eq "-epssize") {
    shift @ARGV;
    $epssize=shift @ARGV;
  } elsif ($ARGV[0] eq "-grid") {
    shift @ARGV;
    $gridname=shift @ARGV;
    $epsgrid=shift @ARGV;
  } elsif ($ARGV[0] eq "-dx") {
    shift @ARGV;
    $dx=1;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    $inpfile=(shift @ARGV);
    $done=1;
  }
}

if (defined $single || defined $double) {
  my $tag=$$;
  my $mol=Molecule::new();

  if (defined $psffile) {
    $mol->readCRD($crdfile);
  } elsif (defined $mol2file) {
    $mol->readMol2($mol2file);
    $mol->generateSegNames();
  } else {
    $mol->readPDB($inpfile);
    $mol->fixHistidine($par{hsd},$par{hse});
    $mol->translate(&CHARMM::getConvType($par{param}));
    $mol->generateSegNames();
  }

  my %pb;

  if (defined $partial && -r $partial) {
    open INP,"$partial";
    while (<INP>) {
      chomp;
      my ($inx,$pbe)=split(/ +/);
      $pb{$inx}=$pbe;
    }
    close INP;
  }

  my $natom=$#{$mol->{atom}}+1;

  my @todolist=();
  my @allist=();
  my $inx=0;

  if (defined $double) {
    my @atomlist=();
    foreach my $c ( @{$mol->{chain}} ) {
      foreach my $a ( @{$c->{atom}} ) { 
       my $jnx=$c->{id}.":".$a->{resnum}.":".$a->{atomname}.":".$inx;
       $inx++;
       push(@atomlist,$jnx);
      }
    }

    for (my $i=0; $i<$#atomlist; $i++) {
     for (my $j=$i; $j<=$#atomlist; $j++) {
      my $jnx=$atomlist[$i]."=".$atomlist[$j];  
      if (!exists $pb{$jnx}) {
	push(@todolist,$jnx);
      }
      push(@allist,$jnx);
     }
    }
  } elsif (defined $single) { 
   foreach my $c ( @{$mol->{chain}} ) {
    foreach my $a ( @{$c->{atom}} ) { 
      my $j=$c->{id}.":".$a->{resnum}.":".$a->{atomname}.":".$inx;
      $inx++;
      if (!exists $pb{$j}) {
	push(@todolist,$j);
      }
      push(@allist,$j);
    }
   }
  }
  
  my @pidlist;
  for (my $i=1; $i<=$threads; $i++) {
    my $pid=fork();
    if (!$pid) {
      $logFile="$i-$logFile" if (defined $logFile);
      my $outfile="$tag-$i-pbout";
      my $handle=&GenUtil::getOutputFile($outfile);
      $handle->autoflush(1);
      my $charmm=&CHARMM::new($logFile,$cmdlog);

      $charmm->loadParameters(%par);

      if (defined $psffile) {
	$charmm->setupFromPSF($psffile,$crdfile);
      } elsif (defined $mol2file) {
	$charmm->setupFromMol2($mol2file);
      } else {
	$charmm->setupFromMolecule($mol);
      }

      if (defined $customfile && &GenUtil::checkFile($customfile)) {
	my $custom=&GenUtil::readData(&GenUtil::getInputFile($customfile));
	$charmm->stream($custom);
      }

      $charmm->orient() unless (!$center);
      for (my $ia=$i-1; $ia<=$#todolist; $ia+=$threads) {
	my $inx=$todolist[$ia];
	my $pbener=$charmm->atomPoissonBoltzmann($inx,$keepcharge);
	printf $handle "%s %f\n",$inx,$pbener;
	printf STDERR "%s %f\n",$inx,$pbener;
      }
      $charmm->finish();
      close $handle;
      exit 0;
    } else {
      push (@pidlist,$pid);
    }
  } 
  foreach my $p (@pidlist) {
    waitpid($p,0);
  }

  for (my $i=1; $i<=$threads; $i++) {
    open INP,"$tag-$i-pbout";
    while (<INP>) {
      chomp;
      my ($inx,$pbe)=split(/ +/);
      $pb{$inx}=$pbe;
    }
    close INP;
    unlink "$tag-$i-pbout";
  }

  foreach my $pkey (@allist) { 
    printf "%s %f %f\n",$pkey,$pb{$pkey},-163.9603525/$pb{$pkey}
      if ($pb{$pkey} != 0.0);
  }
} else {
  my $charmm=&CHARMM::new($logFile,$cmdlog);
  $charmm->loadParameters(%par);

  if (defined $psffile) {
    $charmm->setupFromPSF($psffile,$crdfile);
  } elsif (defined $mol2file) {
    $charmm->setupFromMol2($mol2file);
  } else {
    $charmm->setupFromPDB($inpfile);
  }

  if (defined $customfile && &GenUtil::checkFile($customfile)) {
    my $custom=&GenUtil::readData(&GenUtil::getInputFile($customfile));
    $charmm->stream($custom);
  }

  if (defined $emap) {
    my $saveeps=$charmm->{par}->{epsw};
    $charmm->{par}->{epsw}=1;
    my ($tmaxx,$tmaxy,$tmaxz,$nxcel,$nycel,$nzcel,$xcent,$ycent,$zcent)=$charmm->pbgrid("$$.eps1","phi"); 
    $charmm->{par}->{epsw}=$saveeps;
    ($tmaxx,$tmaxy,$tmaxz,$nxcel,$nycel,$nzcel,$xcent,$ycent,$zcent)=$charmm->pbgrid("$$.epsw","phi"); 
    $charmm->finish();

    &convert2DX($tmaxx,$tmaxy,$tmaxz,$nxcel,$nycel,$nzcel,$xcent,$ycent,$zcent,$charmm->{par}->{dcel},"$$.eps1");
    &convert2DX($tmaxx,$tmaxy,$tmaxz,$nxcel,$nycel,$nzcel,$xcent,$ycent,$zcent,$charmm->{par}->{dcel},"$$.epsw");
    &diffDX("$$.eps1","$$.epsw",$emap);
    &GenUtil::remove("$$.eps1");
    &GenUtil::remove("$$.epsw");
  } elsif (defined $epsgrid) {
    if (defined $gridname) {
      my ($tmaxx,$tmaxy,$tmaxz,$nxcel,$nycel,$nzcel,$xcent,$ycent,$zcent)=$charmm->pbgrid($epsgrid,($gridname eq "eps")?"epsx":$gridname); 
      $charmm->finish();

      my $offsetx=0;
      my $offsety=0;
      my $offsetz=0;
      
      $offsetx=$charmm->{par}->{dcel}/2.0 if ($gridname eq "eps");

      &convert2DX($tmaxx,$tmaxy,$tmaxz,$nxcel,$nycel,$nzcel,$xcent,$ycent,$zcent,$charmm->{par}->{dcel},$epsgrid,$offsetx,$offsety,$offsetz) if ($dx);
    } else {
      $charmm->epsgrid($epsgrid,$epssize);
      $charmm->finish();
    }

  } else {
    $charmm->orient() unless (!$center);
    my $pbener=$charmm->poissonBoltzmann();
    $charmm->finish();
    printf "%f\n",$pbener;
  }
}

exit 0;

sub convert2DX {
  my $tx=shift;
  my $ty=shift;
  my $tz=shift;
  my $nx=shift;
  my $ny=shift;
  my $nz=shift;
  my $cx=shift;
  my $cy=shift;
  my $cz=shift;
  my $dcel=shift;
  my $f=shift;
  my $offx=shift;
  my $offy=shift;
  my $offz=shift;

  $offx=0 if (!defined $offx);
  $offy=0 if (!defined $offy);
  $offz=0 if (!defined $offz);


  my $ftmp="tmp$$.pot";

  open OUT,">$ftmp";

  printf OUT "# Data from CHARMM, PBEQ module\n";
  printf OUT "#\n";
  printf OUT "object 1 class gridpositions counts $nx $ny $nz\n";
  printf OUT "origin %e %e %e\n",-($nx/2.0-0.5)*$dcel+$offx+$cx,-($ny/2.0-0.5)*$dcel+$offy+$cy,-($nz/2.0-0.5)*$dcel+$offz+$cz;
  printf OUT "delta %f 0.0 0.0\n",$dcel;
  printf OUT "delta 0.0 %f 0.0\n",$dcel;
  printf OUT "delta 0.0 0.0 %f\n",$dcel;
  printf OUT "object 2 class gridconnections counts $nx $ny $nz\n";
  printf OUT "object 3 class array type double rank 0 items %d data follows\n",$nx*$ny*$nz;
  
  my $cnt=0;
  open INP,$f;
  while (<INP>) {
    my @f=split(/ +/);
    printf OUT "%e ",$f[4];
    printf OUT "\n" if (++$cnt%3==0);
  }
  close INP;

  close OUT;
  system "mv $ftmp $f";
}

sub diffDX {
  my $file1=shift;
  my $file2=shift;
  my $file3=shift;

  open OUT,">$file3";
  
  my $flag=0;
  open INP1,$file1;
  open INP2,$file2;
  while (<INP1>) {
    if (!$flag) {
      $flag=1 if (/object 3/);
      print OUT;
      <INP2>;
    } else {
      s/^\s+//g;
      s/\s+$//g;
      my @f=split(/\s+/);
      my $line=<INP2>;
      $line=~s/^\s+//g;
      $line=~s/\s+$//g;
      my @g=split(/\s+/,$line);
      for (my $i=0; $i<=$#f; $i++) {
	printf OUT "%e ",$f[$i]-$g[$i];
      }
      printf OUT "\n";
    }
  }
  close INP1;
  close INP2;
}
