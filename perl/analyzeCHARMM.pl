#!/usr/bin/env perl

# CHARMM trajectory analysis 
#
# 2005, Michael Feig, Michigan State University
#
## 2008-12-04: Modified by Lev Gorenstein <lev@purdue.edu>
## 	Fixed a bug in handling '-inx from:to' and '-time from:to' limits
## 	(stop looping over frames once we got outside of the upper limit).
##      The original file is saved as analyzeCHARMM.pl.orig.
##
## 2008-12-06: Modified by Lev Gorenstein <lev@purdue.edu>
## 	Added a '-custom1 file' command line option.  It's similar to
## 	'-custom', but the specified file is only streamed *once* (after
## 	molecule setup and before any trajectory reads).  This allows 
## 	to set up whatever parameters you want for your '-custom' analysis
## 	without them being reread every time.
##
## 2009-04-01: Modified by Lev Gorenstein <lev@purdue.edu>
## 	Fixed a bug in filelist handling (line 207 in original file) that
## 	was hinderng usage of '-f filelist' option. 
##
## 2009-06-13: Modified by Lev Gorenstein <lev@purdue.edu>
## 	Added '-chi2|-chi3|-chi4|-chi5' command line options to complement
## 	other backbone torsions.  Corresponding dihedral definitions are from
## 	http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html


sub usage {
  printf STDERR "usage:   analyzeCHARMM.pl [options] [DCDfile(s)]\n";
  printf STDERR "options: [-inx from:to] [-time from:to] [-step size]\n";
  printf STDERR "         [-psf PSFfile] [-pdb PDBfile] [-crd CRDfile]\n";
  printf STDERR "         [-comp PDBfile] [-compcrd CRDfile]\n";
  printf STDERR "         [-timeonly]\n";
  printf STDERR "         [-out outputfile] [-cout outputfile]\n";
  printf STDERR "         [-log logFile] [-cmd logFile]\n";
  printf STDERR "         [-f listfile]\n";
  printf STDERR "         [-par CHARMMparams]\n";
  printf STDERR "         [-custom file]\n";
  printf STDERR "         [-custom1 file]\n"; 		# Added by LG
  printf STDERR "         [-sel Selection] [-fitsel Selection] [-selout Selection]\n";
  printf STDERR "         [-dsel Selection Selection]\n";
  printf STDERR "         [-tsel Selection Selection Selection]\n";
  printf STDERR "         [-qsel Selection Selection Selection Selection]\n";
  printf STDERR "         [-mass] [-[no]fit] [-ref x y z] [-cofm] [-ndens]\n";
  printf STDERR "         [-rgyr] [-center] [-rms]\n";
  printf STDERR "         [-mindist] [-maxdist] [-alldist] [-cutoff value]\n";
  printf STDERR "         [-dist] [-angle] [-dihedral] [-pucker] [-glycosidic]\n";
  printf STDERR "         [-phi] [-psi] [-omega] [-chi1] [-chi2] [-chi3] [-chi4] [-chi5] \n"; 	# By LG
  printf STDERR "         [-epsilon] [-zeta] [-delta] [-alpha] [-beta] [-gamma]\n";
  printf STDERR "         [-volume] [-ientropy] [-surface]\n";
  printf STDERR "         [-rdist] [-diffusion]\n";
  printf STDERR "         [-rmsf] [-s2] [-aniso] [-rmsdyn] [-rotcor]\n";
  printf STDERR "         [-hbond] [-hbdist] [-hbangle] [-hbtime] [-hbunit]\n";
  printf STDERR "         [-avg] [-orient] [-recenter] [-norotate]\n";
  printf STDERR "         [-cpath] [-periodic] [-energy]\n";
  printf STDERR "         [-crot]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use POSIX ":sys_wait_h";

use Fcntl;

use GenUtil;
use Molecule;
use CHARMM;

my %par;

my $logFile;
my $cmdlog;

my $pdbfile;
my $psffile;
my $crdfile;

my $listfile;

my $timeonly=0;

my $ifrom;
my $ito;
my $tfrom;
my $tto;

my $periodic=0;

my $stepsize=1;

my $pdbcomp;
my $crdcomp;

my $customfile;
my @customfile1; 	# Added by LG. Note: array to allow multiple '-custom1'

my @dcdfiles;

my $selection1;
my $selection2;
my $selection3;
my $selection4;

my $mass=0;
my $cofm=0;
my $refx;
my $refy;
my $refz;

my $ndens=0;

my $fit=1;
my $fitsel;

my $outfile;
my $coutfile;
my $selout;

my $rmsdyn=0;

my $orient=0;
my $recenter=0;
my $rotate=1;
my $hbdist=0; #in anstroms
my $hbtime=0; #in picoseconds
my $hbang=999; #Angle
my $hbunit=51; #Starting unit to read from
my @hbdcds=(); #DCD files

my $cutoff=9999;

my $setupenergy=0;

my @analyze=();

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logFile=(shift @ARGV);
  } elsif ($ARGV[0] eq "-cmd") {
    shift @ARGV;
    $cmdlog=(shift @ARGV);
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-custom1") {
    # Added by LG.
    # Note: unlike '-custom' (which can be specified only once), @customfile1
    # is an array, and thus multiple '-custom1' options are allowed
    # (with corresponding files streamed sequentially).
    shift @ARGV;
    push(@customfile1, shift @ARGV);
  } elsif ($ARGV[0] eq "-f") {
    shift @ARGV;
    $listfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-pdb") {
    shift @ARGV;
    $pdbfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-crd") {
    shift @ARGV;
    $crdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-comp") {
    shift @ARGV;
    $pdbcomp=shift @ARGV;
  } elsif ($ARGV[0] eq "-compcrd") {
    shift @ARGV;
    $crdcomp=shift @ARGV;
  } elsif ($ARGV[0] eq "-psf") {
    shift @ARGV;
    $psffile=shift @ARGV;
  } elsif ($ARGV[0] eq "-out") {
    shift @ARGV;
    $outfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-cout") {
    shift @ARGV;
    $coutfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-timeonly") {
    shift @ARGV;
    $timeonly=1;
  } elsif ($ARGV[0] eq "-inx") {
    shift @ARGV;
    my @f=split(/:/,shift @ARGV);
    $ifrom=$f[0];
    $ito=($#f>0)?$f[1]:$ifrom;
  } elsif ($ARGV[0] eq "-time") {
    shift @ARGV;
    my @f=split(/:/,shift @ARGV);
    $tfrom=$f[0];
    $tto=($#f>0)?$f[1]:$tfrom;
  } elsif ($ARGV[0] eq "-step") {
    shift @ARGV;
    $stepsize=shift @ARGV;
  } elsif ($ARGV[0] eq "-sel") {
    shift @ARGV;
    $selection1=shift @ARGV;
  } elsif ($ARGV[0] eq "-sel2") {
    shift @ARGV;
    $selection2=shift @ARGV;
  } elsif ($ARGV[0] eq "-sel3") {
    shift @ARGV;
    $selection3=shift @ARGV;
  } elsif ($ARGV[0] eq "-sel4") {
    shift @ARGV;
    $selection4=shift @ARGV;
  } elsif ($ARGV[0] eq "-dsel") {
    shift @ARGV;
    $selection1=shift @ARGV;
    $selection2=shift @ARGV;
  } elsif ($ARGV[0] eq "-tsel") {
    shift @ARGV;
    $selection1=shift @ARGV;
    $selection2=shift @ARGV;
    $selection3=shift @ARGV;
  } elsif ($ARGV[0] eq "-qsel") {
    shift @ARGV;
    $selection1=shift @ARGV;
    $selection2=shift @ARGV;
    $selection3=shift @ARGV;
    $selection4=shift @ARGV;
  } elsif ($ARGV[0] eq "-mass") {
    shift @ARGV;
    $mass=1;
  } elsif ($ARGV[0] eq "-fit") {
    shift @ARGV;
    $fit=1;
  } elsif ($ARGV[0] eq "-nofit") {
    shift @ARGV;
    $fit=0;
  } elsif ($ARGV[0] eq "-cutoff") {
    shift @ARGV;
    $cutoff=shift @ARGV;
  } elsif ($ARGV[0] eq "-ref") {
    shift @ARGV;
    $refx=shift @ARGV;
    $refy=shift @ARGV;
    $refz=shift @ARGV;
  } elsif ($ARGV[0] eq "-cofm") {
    shift @ARGV;
    $cofm=1;
  } elsif ($ARGV[0] eq "-ndens") {
    shift @ARGV;
    $ndens=1;
  } elsif ($ARGV[0] eq "-fitsel") {
    shift @ARGV;
    $fitsel=shift @ARGV;
  } elsif ($ARGV[0] eq "-selout") {
    shift @ARGV;
    $selout=shift @ARGV;
  } elsif ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-recenter") {
    shift @ARGV;
    $recenter=1; 
  } elsif ($ARGV[0] eq "-norotate") {
    shift @ARGV;
    $rotate=0;
  } elsif ($ARGV[0] eq "-orient") {
    shift @ARGV;
    $orient=1;
  } elsif ($ARGV[0] eq "-rmsdyn"){
    shift @ARGV;
    $rmsdyn=1;
  } elsif ($ARGV[0] eq "-hbdist"){
    shift @ARGV;
    $hbdist=shift @ARGV;
  } elsif ($ARGV[0] eq "-hbangle"){
    shift @ARGV;
    $hbang=shift @ARGV;
  } elsif ($ARGV[0] eq "-hbtime"){
    shift @ARGV;
    $hbtime=shift @ARGV;
  } elsif ($ARGV[0] eq "-hbunit"){
    shift @ARGV;
    $hbunit=shift @ARGV;
  } elsif ($ARGV[0] eq "-periodic") {
    shift @ARGV;
    $periodic=1;
  } elsif ($ARGV[0] eq "-energy") {
    shift @ARGV;
    push(@analyze,"energy");
    $setupenergy=1;
  } elsif ($ARGV[0] eq "-setupenergy") {
    shift @ARGV;
    $setupenergy=1;
  } elsif ($ARGV[0] =~ /^-([a-zA-Z0-9]+)$/) {
    shift @ARGV;
    if ($1 eq "center") { $fit=0; }
    push(@analyze,$1);
  } else {
    die "Unknown option $ARGV[0]" if ($ARGV[0]=~/^-/);
    my $fname=shift @ARGV;
    push(@dcdfiles,$fname) if (-r $fname);
  }
}


push(@analyze,"transform") if ($recenter || $orient || $selout || $rmsdyn);

if (-r $listfile) {
  open INP,"$listfile";
  while (<INP>) {
    chomp;
    push(@dcdfiles,$_) if (-r $_);
  }
  close INP;
}

if ($#dcdfiles<0) {
  printf STDERR "no DCD files given\n";
  usage();
}

my $charmm=CHARMM::new($logFile,$cmdlog);

$charmm->loadParameters(%par);

if (defined $psffile) {
  $charmm->setupFromPSF($psffile,$crdfile);
  $charmm->readFromPDB($pdbfile) if (defined $pdbfile && !defined $crdfile);
} else {
  $charmm->setupFromPDB($pdbfile);
}

# Added by LG.
# Customfile1 is almost the same as '-custom', but is read only once.
# Note: if your initial setup includes some energy settings 
# 	(i.e. UPDATE ;-), then you'd better be sure to use 
# 	'-crd' or '-pdb' option on the command line!
# Note: unlike '-custom' (which can be specified only once), @customfile1
# 	is an array, and thus multiple '-custom1' options are allowed
# 	(with corresponding files streamed sequentially).
foreach my $file1 ( @customfile1 ) {
   if (defined $file1 && &GenUtil::checkFile($file1)) {
      $charmm->_sendCommand();
      my $custom1=&GenUtil::readData(&GenUtil::getInputFile($file1));
      $charmm->stream($custom1);
   }
}

if (($periodic || grep(/rdist|diffusion|crot|energy/,@analyze)) && !defined $charmm->{par}->{boxsize} && 
    !defined $charmm->{par}->{boxx}) {
  my ($bx,$by,$bz,$shape)=&_getBoxSizeFromDCD($dcdfiles[0]);
  if (defined $bx) {
    $charmm->setParameter(boxx=>$bx, boxy=>$by, boxz=>$bz, boxshape=>$shape);
  }
}

#$charmm->setupEnergy() if ((defined $pdbfile || defined $crdfile) && grep(/rdist/,@analyze));

$charmm->setupEnergy() if ($setupenergy);

$charmm->periodicBoundaries() if ($periodic);

my $comp=0;
if (defined $pdbcomp) {
  $charmm->loadReference($pdbcomp,9001.0);
  $comp=1;
}

if (defined $crdcomp) {
  $charmm->loadReferenceCRD($crdcomp,9001.0);
  $comp=1;
}

my $sel1;
my $sel2;
my $sel3;
my $sel4;
my $selft;
my $selo;

my @selections=();
my @puckerlist=();

if (grep(/phi|psi|omega|chi|alpha|beta|gamma|delta|epsilon|zeta|pucker|glycosidic/,@analyze)) {
  push(@analyze,"dihedral") if (grep(/phi|psi|omega|chi|alpha|beta|gamma|delta|epsilon|zeta|glycosidic/,@analyze)); 
  my $m=$charmm->{molecule};
  $m->setValidSelection($selection1);
  foreach my $cm ( @{$m->{chain}} ) {
    for (my $ir=0; $ir<=$#{$cm->{res}}; $ir++) {
      my $rm=$cm->{res}->[$ir];
      if ($rm->{valid}) {
       if ($rm->{name}=~/ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HSP|HID|HIE|HIP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX/) {
        foreach my $a ( @analyze ) {
         if ($a eq "phi") {
	   if ($ir>0 && $rm->{num}-1==$cm->{res}->[$ir-1]->{num}) {
	     my $srec={};
	     $srec->{sel1}=sprintf("ATOM %s %d C ",$rm->{seg},$rm->{num}-1);
	     $srec->{sel2}=sprintf("ATOM %s %d N ",$rm->{seg},$rm->{num});
	     $srec->{sel3}=sprintf("ATOM %s %d CA",$rm->{seg},$rm->{num});
	     $srec->{sel4}=sprintf("ATOM %s %d C ",$rm->{seg},$rm->{num});
	     push(@selections,$srec);
	   } else {
	     my $srec={};
	     $srec->{sel1}=sprintf("ATOM %s %d CY",$rm->{seg},$rm->{num});
	     $srec->{sel2}=sprintf("ATOM %s %d N ",$rm->{seg},$rm->{num});
	     $srec->{sel3}=sprintf("ATOM %s %d CA",$rm->{seg},$rm->{num});
	     $srec->{sel4}=sprintf("ATOM %s %d C ",$rm->{seg},$rm->{num});
	     push(@selections,$srec);
	   }
	 }
      
         if ($a eq "psi") {
	   if ($ir<$#{$cm->{res}} && $rm->{num}+1==$cm->{res}->[$ir+1]->{num}) {
	     my $srec={};
	     $srec->{sel1}=sprintf("ATOM %s %d N ",$rm->{seg},$rm->{num});
	     $srec->{sel2}=sprintf("ATOM %s %d CA ",$rm->{seg},$rm->{num});
	     $srec->{sel3}=sprintf("ATOM %s %d C ",$rm->{seg},$rm->{num});
	     $srec->{sel4}=sprintf("ATOM %s %d N ",$rm->{seg},$rm->{num}+1);
	     push(@selections,$srec);
	   } else {
	     my $srec={};
	     $srec->{sel1}=sprintf("ATOM %s %d N ",$rm->{seg},$rm->{num});
	     $srec->{sel2}=sprintf("ATOM %s %d CA ",$rm->{seg},$rm->{num});
	     $srec->{sel3}=sprintf("ATOM %s %d C ",$rm->{seg},$rm->{num});
	     $srec->{sel4}=sprintf("ATOM %s %d NT",$rm->{seg},$rm->{num});
	     push(@selections,$srec);
	   }
	 }
         if ($a eq "omega" && $ir<$#{$cm->{res}} && $rm->{num}+1==$cm->{res}->[$ir+1]->{num}) {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d CA ",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d C ",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d N ",$rm->{seg},$rm->{num}+1);
          $srec->{sel4}=sprintf("ATOM %s %d CA ",$rm->{seg},$rm->{num}+1);
          push(@selections,$srec);
         } 
         if ($a eq "chi1") {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d N ",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d CA ",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d CB ",$rm->{seg},$rm->{num});
          $srec->{sel4}=sprintf("ATOM %s %d CG ",$rm->{seg},$rm->{num}) 
            if ($rm->{name}=~/GLU|PRO|LYS|GLN|ARG|LEU|PHE|TYR|TRP|ASN|ASP|MET|HIS|HSD|HSE|HSP/);
          $srec->{sel4}=sprintf("ATOM %s %d OG ",$rm->{seg},$rm->{num}) 
            if ($rm->{name} eq "SER");
          $srec->{sel4}=sprintf("ATOM %s %d OG1 ",$rm->{seg},$rm->{num}) 
            if ($rm->{name} eq "THR");
          $srec->{sel4}=sprintf("ATOM %s %d CG1 ",$rm->{seg},$rm->{num}) 
            if ($rm->{name}=~/VAL|ILE/);
          $srec->{sel4}=sprintf("ATOM %s %d SG ",$rm->{seg},$rm->{num}) 
            if ($rm->{name} eq "CYS");
          push(@selections,$srec) if (defined $srec->{sel4});
         } 

         # Added by LG.
         if ($a eq "chi2") {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d CA ",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d CB ",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d CG ",$rm->{seg},$rm->{num})
            if ($rm->{name}=~/ARG|ASN|ASP|GLN|GLU|HIS|HSD|HSE|HSP|LEU|LYS|MET|PHE|PRO|TRP|TYR/);
          $srec->{sel3}=sprintf("ATOM %s %d CG1 ",$rm->{seg},$rm->{num})
            if ($rm->{name}=~/ILE/);
          $srec->{sel4}=sprintf("ATOM %s %d CD ",$rm->{seg},$rm->{num}) 
            if ($rm->{name}=~/ARG|GLN|GLU|ILE|LYS|PRO/);
          $srec->{sel4}=sprintf("ATOM %s %d CD1 ",$rm->{seg},$rm->{num}) 
            if ($rm->{name}=~/LEU|PHE|TRP|TYR/);
          $srec->{sel4}=sprintf("ATOM %s %d OD1 ",$rm->{seg},$rm->{num}) 
            if ($rm->{name}=~/ASN|ASP/);
          $srec->{sel4}=sprintf("ATOM %s %d ND1 ",$rm->{seg},$rm->{num}) 
            if ($rm->{name}=~/HIS|HSD|HSE|HSP/);
          $srec->{sel4}=sprintf("ATOM %s %d SD ",$rm->{seg},$rm->{num}) 
            if ($rm->{name} eq "MET");
          push(@selections,$srec) 
            if (defined $srec->{sel3} && defined $srec->{sel4});
         } 

         # Added by LG.
         if ($a eq "chi3") {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d CB ",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d CG ",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d CD ",$rm->{seg},$rm->{num})
            if ($rm->{name}=~/ARG|GLN|GLU|LYS/);
          $srec->{sel3}=sprintf("ATOM %s %d SD ",$rm->{seg},$rm->{num})
            if ($rm->{name} eq "MET");
          $srec->{sel4}=sprintf("ATOM %s %d NE ",$rm->{seg},$rm->{num}) 
            if ($rm->{name} eq "ARG");
          $srec->{sel4}=sprintf("ATOM %s %d OE1 ",$rm->{seg},$rm->{num}) 
            if ($rm->{name}=~/GLN|GLU/);
          $srec->{sel4}=sprintf("ATOM %s %d CE ",$rm->{seg},$rm->{num}) 
            if ($rm->{name}=~/LYS|MET/);
          push(@selections,$srec) 
            if (defined $srec->{sel3} && defined $srec->{sel4});
         } 

         # Added by LG.
         if ($a eq "chi4") {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d CG ",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d CD ",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d NE ",$rm->{seg},$rm->{num})
            if ($rm->{name} eq "ARG");
          $srec->{sel3}=sprintf("ATOM %s %d CE ",$rm->{seg},$rm->{num})
            if ($rm->{name} eq "LYS");
          $srec->{sel4}=sprintf("ATOM %s %d CZ ",$rm->{seg},$rm->{num}) 
            if ($rm->{name} eq "ARG");
          $srec->{sel4}=sprintf("ATOM %s %d NZ ",$rm->{seg},$rm->{num}) 
            if ($rm->{name} eq "LYS");
          push(@selections,$srec) 
            if (defined $srec->{sel3} && defined $srec->{sel4});
         } 

         # Added by LG.
         if ($a eq "chi5") {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d CD ",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d NE ",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d CZ ",$rm->{seg},$rm->{num}) 
            if ($rm->{name} eq "ARG");
          $srec->{sel4}=sprintf("ATOM %s %d NH1 ",$rm->{seg},$rm->{num}) 
            if ($rm->{name} eq "ARG");
          push(@selections,$srec) 
            if (defined $srec->{sel3} && defined $srec->{sel4});
         } 
        }
       } elsif ($rm->{name}=~/ADE|GUA|CYT|THY|URA/) {
        foreach my $a ( @analyze ) {
         if ($a eq "pucker") {
           my $srec={};
           $srec->{segid}=$rm->{seg};
           $srec->{resid}=$rm->{num};
           push(@puckerlist,$srec);
         } 
         if ($a eq "zeta" && $ir<$#{$cm->{res}} && $rm->{num}+1==$cm->{res}->[$ir+1]->{num}) {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d C3' ",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d O3' ",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d P ",$rm->{seg},$rm->{num}+1);
          $srec->{sel4}=sprintf("ATOM %s %d O5' ",$rm->{seg},$rm->{num}+1);
          push(@selections,$srec);
         } 
         if ($a eq "epsilon" && $ir<$#{$cm->{res}} && $rm->{num}+1==$cm->{res}->[$ir+1]->{num}) {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d C4' ",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d C3' ",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d O3' ",$rm->{seg},$rm->{num});
          $srec->{sel4}=sprintf("ATOM %s %d P ",$rm->{seg},$rm->{num}+1);
          push(@selections,$srec);
         } 
         if ($a eq "delta") {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d C5' ",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d C4' ",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d C3' ",$rm->{seg},$rm->{num});
          $srec->{sel4}=sprintf("ATOM %s %d O3' ",$rm->{seg},$rm->{num});
          push(@selections,$srec);
         } 
         if ($a eq "gamma") {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d O5'",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d C5'",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d C4'",$rm->{seg},$rm->{num});
          $srec->{sel4}=sprintf("ATOM %s %d C3'",$rm->{seg},$rm->{num});
          push(@selections,$srec);
         } 
         if ($a eq "beta") {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d P",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d O5'",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d C5'",$rm->{seg},$rm->{num});
          $srec->{sel4}=sprintf("ATOM %s %d C4'",$rm->{seg},$rm->{num});
          push(@selections,$srec);
         } 
         if ($a eq "alpha" && $ir>0 && $rm->{num}-1==$cm->{res}->[$ir-1]->{num}) {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d O3'",$rm->{seg},$rm->{num}-1);
          $srec->{sel2}=sprintf("ATOM %s %d P",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d O5'",$rm->{seg},$rm->{num});
          $srec->{sel4}=sprintf("ATOM %s %d C5'",$rm->{seg},$rm->{num});
          push(@selections,$srec);
	 }
         if ($a eq "glycosidic" && $rm->{name}=~/ADE|GUA/ ) {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d O4'",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d C1'",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d N9",$rm->{seg},$rm->{num});
          $srec->{sel4}=sprintf("ATOM %s %d C4",$rm->{seg},$rm->{num});
          push(@selections,$srec);
         } 
         if ($a eq "glycosidic" && $rm->{name}=~/CYT|THY|URA/ ) {
          my $srec={};
          $srec->{sel1}=sprintf("ATOM %s %d O4'",$rm->{seg},$rm->{num});
          $srec->{sel2}=sprintf("ATOM %s %d C1'",$rm->{seg},$rm->{num});
          $srec->{sel3}=sprintf("ATOM %s %d N1",$rm->{seg},$rm->{num});
          $srec->{sel4}=sprintf("ATOM %s %d C2",$rm->{seg},$rm->{num});
          push(@selections,$srec);
         } 

        }
       }
      }
    }  
  }
} else {
  if (defined $selection1) {
    $sel1="SL1";
    $charmm->defineSelection($selection1,$sel1);
  }

  if (defined $selection2) {
    $sel2="SL2";
    $charmm->defineSelection($selection2,$sel2);
  }

  if (defined $selection3) {
    $sel3="SL3";
    $charmm->defineSelection($selection3,$sel3);
  }

  if (defined $selection4) {
    $sel4="SL4";
    $charmm->defineSelection($selection4,$sel4);
  }

  if (defined $fitsel) {
    $selft="SLF";
    $charmm->defineSelection($fitsel,$selft," .and. property xcomp .lt. 9000");
  }

  if (defined $selout){
    $selo="SLO";
    $charmm->defineSelection($selout,$selo);
  }
}

my $accu;
my $naccu=0;
foreach my $f ( @dcdfiles ) {
  foreach my $a ( @analyze ) {
    if ($a eq "rdist") {
       my ($n,$data)=$charmm->analyzeRadialDistribution($f,
         selection=>$sel1,xselection=>$sel2,refx=>$refx,refy=>$refy,refz=>$refz,
         cofm=>$cofm, ndens=>$ndens );

       if (!defined $accu) {
	 $accu=();
         foreach my $d (@{$data} ) {
           $d->{val}*=$n;
           push(@{$accu},$d);
         }
       } else {
         for (my $id=0; $id<=$#{$data}; $id++) {
           $accu->[$id]->{val}+=$n*$data->[$id]->{val};
         }
       }
       $naccu+=$n;
     } elsif ($a eq "crot") {
       my ($n,$data)=$charmm->analyzeRotationalCorrelation($f,selection=>$sel1);
       if (!defined $accu) {
	 $accu=();
         foreach my $d (@{$data} ) {
           $d->{val}*=$n;
           push(@{$accu},$d);
         }
       } else {
         for (my $id=0; $id<=$#{$data}; $id++) {
           $accu->[$id]->{val}+=$n*$data->[$id]->{val};
         }
       }
       $naccu+=$n;
     } elsif ($a eq "diffusion") {
       my ($n,$diff)=$charmm->analyzeDiffusion($f,
	selection=>$sel1,xselection=>$sel2,refx=>$refx,refy=>$refy,refz=>$refz,cofm=>$cofm);
       if (!defined $accu) {
	 $accu=();
	 my $d={};
         $d->{desc}="diffusion [A*A/ps]";
	 $d->{val}=$diff*$n;
	 push(@{$accu},$d);
       } else {
	 $accu->[0]->{val}+=$diff*$n;
       }
       $naccu+=$n;
     } elsif ($a eq "rmsf") {
       if ($naccu>0) {
         printf STDERR "only first trajectory is processed for RMSF calculation\n";
       } else {
        if (!defined $sel1) {
         $sel1="SL1";
         $selection1="CA";
         $charmm->defineSelection($selection1,$sel1);
        }
        if ($fit && !defined $selft) {
         $selft="SLF";
         $charmm->defineSelection("solute.heavy",$selft," .and. property xcomp .lt. 9000");
        }
        my ($alldata)=$charmm->analyzeRMSFluctuations($f,mass=>$mass,fit=>$fit,
	                                             selection=>$sel1,fitsel=>$selft);
        my $cmol=$charmm->{molecule}->clone();
        $cmol->setAtomPropFromList("aux1",$alldata);
        $cmol->setValidSelection($selection1);

        foreach my $c ( @{$cmol->{chain}} ) {
         foreach my $a ( @{$c->{atom}} ) {
           if ($a->{valid}) {
             my $trec={};
             $trec->{desc}=sprintf("%s:%s:%s:%s",$a->{seg},$a->{resname},$a->{resnum},$a->{atomname});
             $trec->{val}=$a->{aux1};
             push(@{$accu},$trec);
           }
         }
        }
        $naccu=1;
       }
     } elsif ($a eq "avg") {
       if ($naccu>0) {
         printf STDERR "only first trajectory is processed for average calculation\n";
       } else {
        if (!defined $sel1) {
         $sel1="SL1";
         $selection1="solute";
         $charmm->defineSelection($selection1,$sel1);
        }
        if ($fit && !defined $selft) {
         $selft="SLF";
         $charmm->defineSelection("solute.heavy",$selft," .and. property xcomp .lt. 9000");
        }
        $charmm->analyzeTrajectoryAverage($f,mass=>$mass,fit=>$fit,
	                                  selection=>$sel1,fitsel=>$selft);
        my $chmoutpdb=lc "pdb$$-out";
        $charmm->writePDB($chmoutpdb);
	my $tmol=Molecule::new();
	$tmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),
		       chainfromseg=>1);
	$tmol->setValidSelection($selection1);
	my $outmol=$tmol->clone(1);
	$outmol->setSSBonds($charmm->{molecule}->getSSBonds());
	$outmol->writePDB("-",translate=>"CHARMM22");
	&GenUtil::remove($chmoutpdb);

        $naccu=1;
       }
     } elsif ($a eq "transform") {
       if (!defined $outfile ) {
	 printf STDERR "need to specify output file with -out\n";
	 exit(1);
       }
       if ((!defined $pdbcomp || !-r $pdbcomp) && !defined $selout && !$rmsdyn) {
	 printf STDERR "need to specify reference file with -comp\n";
	 exit(1);
       }
       if ($naccu>0) {
         printf STDERR "only first trajectory is transformed\n";
       } else {
        if ($fit && !defined $selft) {
          $selft="SLF";
         $charmm->defineSelection("solute.heavy",$selft," .and. property xcomp .lt. 9000");
        }
	if (defined $selout){
	    $charmm->selectTrajectoryOutput($f,selout=>$selo,output=>$outfile,mass=>$mass,fit=>$fit,fitsel=>$selft,
					    recenter=>$recenter,rotate=>$rotate);
	}
  elsif ($rmsdyn){
      if (!defined $sel1){
        $sel1="SL1";
        $charmm->defineSelection("solute.heavy",$sel1);
      }
      $charmm->analyzeTrajectoryRmsdyn($f,sel=>$sel1,step=>$stepsize,output=>$outfile,mass=>$mass);
  }
	else{
	    $charmm->analyzeTrajectoryTransform($f,fitsel=>$selft,output=>$outfile,mass=>$mass,
                                            orient=>$orient,recenter=>$recenter,rotate=>$rotate);
	}
        $naccu=1;
       }
     } elsif ($a eq "s2") {
       if ($naccu>0) {
         printf STDERR "only first trajectory is processed for order parameter calculation\n";
       } else {
	 if (!defined $sel1) {
	   $sel1="SL1";
	   $selection1="solute.heavy";
	   $charmm->defineSelection($selection1,$sel1);
	 }
	 if ($fit && !defined $selft) {
	   $selft="SLF";
	   $charmm->defineSelection("solute.heavy",$selft," .and. property xcomp .lt. 9000");
	 }
	
	 if ($fit && defined $pdbcomp && -r $pdbcomp) {
	   $outfile="$$-traj.dcd";
	   $charmm->analyzeTrajectoryTransform($f,fitsel=>$selft,output=>$outfile,orient=>1,recenter=>$recenter,rotate=>$rotate);
	   $f=$outfile;
	 }
	 
	 my ($alldata)=$charmm->analyzeTrajectoryOrderParameters($f,selection=>$selection1);
	 
	 if ($fit && defined $pdbcomp && -r $pdbcomp) {
	   &GenUtil::remove($outfile);
	 }

	 foreach my $a (@{$alldata}) {
	   my $trec={};
	   $trec->{desc}=sprintf("%s:%s:%s",$a->{seg},$a->{resname},$a->{resnum});
	   $trec->{val}=$a->{s2};
	   push(@{$accu},$trec);
	 }
	 $naccu=1;
       }
     }
    elsif ($a eq "aniso" || $a eq "anisotropy"){
	if ($naccu>0) {
	    printf STDERR "only first trajectory is processed for anisotropy calculation\n";
	} else {
	    if (!defined $sel1) {
		printf STDERR "Error: selection not defined\n";
		&usage();
	    }
	}

	if ($fit && !defined $selft) {
	   $selft="SLF";
	   $charmm->defineSelection("solute.heavy",$selft," .and. property xcomp .lt. 9000");
	 }
	 if ($fit && defined $pdbcomp && -r $pdbcomp) {
	   $outfile="$$-traj.dcd";
	   $charmm->analyzeTrajectoryTransform($f,fitsel=>$selft,output=>$outfile,orient=>1,recenter=>$recenter,rotate=>$rotate);
	   $f=$outfile;
	 }
	($accu)=$charmm->analyzeTrajectoryAnisotropy($f,selection=>$selection1);
	if ($fit && defined $pdbcomp && -r $pdbcomp) {
	    &GenUtil::remove($outfile);
	}
	
#	foreach my $a (@{$accu}) {
#	    printf "%f %f\n", $a->{time}, $a->{P2};
#	}
	 
	$naccu=1;
    }
    elsif ($a eq "rotcor") {
	if ($naccu>0) {
	    printf STDERR "only first trajectory is processed for anisotropy calculation\n";
	} else {
	    if (!defined $sel1) {
		printf STDERR "Error: selection not defined\n";
		&usage();
	    }
	}
        my ($alldata)=$charmm->analyzeTrajectoryRotationalCorrelation($f,selection=>$selection1);
	foreach my $a (@{$alldata}) {
	   $a->{desc}=sprintf("%lf",$a->{t});
	   push(@{$accu},$a);
	 }
        $naccu=1;
    }
    elsif ($a eq "hbond"){
	if (!defined $sel1){
	    printf STDERR "Error: selection 1 not defined\n";
	    &usage();
	}
	elsif (!defined $sel2){
	    printf STDERR "Error: selection 2 not defined\n";
	}
	elsif ($#hbdcds == $#dcdfiles-2){
	    push(@hbdcds,$f);
	}
	else{
	    push(@hbdcds,$f);
	    $charmm->analyzeTrajectoryHbond(\@hbdcds,$hbunit,sel1=>$sel1,sel2=>$sel2,dist=>$hbdist,time=>$hbtime,angle=>$hbang);
	}
	$naccu=1;
    }
  }
}

if ($naccu>0) {
  my $oneline=0;
  foreach my $d ( @{$accu} ) {
    if (exists $d->{desc}) {
       printf "%s %1.10f\n",$d->{desc},$d->{val}/$naccu;
    } elsif (exists $d->{x}) {
      printf "%f %1.10f\n",$d->{x},$d->{val}/$naccu;
    } elsif (exists $d->{P2}){
	printf "%f %1.10f %1.10f\n", $d->{time}, $d->{P2}, 0.4*($d->{P2});
  #P2 is the second order Legendre Polynomial
  #0.4*P2 is the anisotropy value commonly found - r(t)
    } else {
      printf "%1.10f ",$d->{val}/$naccu;
      $oneline=1;
    }
  }
  printf "\n" if ($oneline);
} else {
 my $iframe=0;
 foreach my $f ( @dcdfiles ) {
  $charmm->initTrajectory($f,$coutfile);

  while ($charmm->nextFrame()>0) {
    my $timestep=($charmm->{_firstframe}+$charmm->{_readframes}-1)*
      $charmm->{_trajfreq}*$charmm->{_trajdelta};
    
    $iframe++;

    # Added by LG - no point of looping past requested limits (if any).
    last if (defined $ito && $iframe>$ito);
    last if (defined $tto && $timestep>$tto+0.0000001);

    if ((!defined $ifrom || ($iframe>=$ifrom && $iframe<=$ito)) && 
	(!defined $tfrom || ($timestep>=$tfrom-0.0000001 && $timestep<=$tto+0.0000001)) &&
	(!defined $stepsize || ($iframe%$stepsize==0))) {
      if (defined $customfile && &GenUtil::checkFile($customfile)) {
	my $custom=&GenUtil::readData(&GenUtil::getInputFile($customfile));
	$charmm->stream($custom);
      } elsif ($#analyze>=0) {
	my @results=();
	foreach my $a ( @analyze ) {
	  if ($a eq "energy") {
	    my $ener=$charmm->getEnergy();
	    push(@results,$ener->{total},$ener->{elec},$ener->{vdwaals},$ener->{bonds},$ener->{angles},$ener->{ureyb},$ener->{dihedrals},$ener->{impropers},$ener->{cmap},$ener->{asp},$ener->{gb},$ener->{pmf1d},$ener->{pmf2d},$ener->{primo},$ener->{user});
	  } elsif ($a eq "rgyr") {
	    push(@results,$charmm->analyzeRadiusOfGyration(selection=>$sel1, mass=>$mass));
	  } elsif ($a eq "center") {
	    push(@results,$charmm->analyzeCenterOfMass(selection=>$sel1, mass=>$mass,fitsel=>$selft,fit=>$fit));
	  } elsif ($a eq "cpath") {
	    push(@results,$charmm->analyzeCpath());
	  } elsif ($a eq "rms") {
	    push(@results,$charmm->analyzeRMS(selection=>$sel1, mass=>$mass, fitsel=>$selft, fit=>$fit));
	  } elsif ($a eq "alldist") {
	    push(@results,$charmm->analyzeAllDistances(selection=>$sel1, xselection=>$sel2, cutoff=>$cutoff));
	  } elsif ($a eq "mindist") {
	    push(@results,$charmm->analyzeMinimumDistance(selection=>$sel1, xselection=>$sel2));
	  } elsif ($a eq "maxdist") {
	    push(@results,$charmm->analyzeMaximumDistance(selection=>$sel1, xselection=>$sel2));
	  } elsif ($a eq "volume") {
	    push(@results,$charmm->analyzeMolecularVolume(selection=>$sel1));
	  } elsif ($a eq "surface") {
	    push(@results,$charmm->analyzeAccessibleSurface(context=>$sel1, calculation=>$sel2));
	  } elsif ($a eq "ientropy") {
	    push(@results,$charmm->analyzeInertiaEntropy(selection=>$sel1));
          } elsif ($a eq "dist") {
            if ($#selections>=0) {
              foreach my $s ( @selections ) {
  	        push(@results,$charmm->analyzeDistance(selection1=>$s->{sel1},selection2=>$s->{sel2},
                                                       mass=>$mass, comp=>$comp));
              } 
            } else { 
  	      push(@results,$charmm->analyzeDistance(selection1=>$sel1,selection2=>$sel2,
		 				     mass=>$mass, comp=>$comp));
            }
          } elsif ($a eq "angle") {
            if ($#selections>=0) {
              foreach my $s ( @selections ) {
  	        push(@results,$charmm->analyzeAngle(selection1=>$s->{sel1},selection2=>$s->{sel2},
                                                    selection3=>$s->{sel3},mass=>$mass));
              } 
            } else { 
  	      push(@results,$charmm->analyzeAngle(selection1=>$sel1,selection2=>$sel2,
		 				  selection3=>$sel3,mass=>$mass));
            }
          } elsif ($a eq "dihedral") {
            if ($#selections>=0) {
              foreach my $s ( @selections ) {
  	        push(@results,$charmm->analyzeDihedral(selection1=>$s->{sel1},selection2=>$s->{sel2},
	  	 				       selection3=>$s->{sel3},selection4=>$s->{sel4},
                                                       mass=>$mass));
              } 
            } else { 
  	      push(@results,$charmm->analyzeDihedral(selection1=>$sel1,selection2=>$sel2,
		 				     selection3=>$sel3,selection4=>$sel4,mass=>$mass));
            }
	  } elsif ($a eq "pucker" && $#puckerlist>=0) {
            foreach my $p ( @puckerlist ) {
  	      push(@results,$charmm->analyzePucker(segid=>$p->{segid},resid=>$p->{resid}));
            }
	  }
	}
	if ($timeonly) {
	  printf "%f",$timestep;
	} else {
	  printf "%d %d %f",$iframe,$charmm->{_readframes},$timestep;
	}
	foreach my $r (@results) {
	  printf " %s",$r; 
	}
	printf "\n";
      }
    }
  }
  $charmm->closeTrajectory();
 }
}

$charmm->finish();

exit 0;

1;
 
sub _getBoxSizeFromDCD() {
  my $dcd=shift;

  my $dcdfile=&GenUtil::getInputFile($dcd);
  binmode $dcdfile;

  my $buffer;
  my $len;
  ($buffer,$len)=&GenUtil::readFortran($dcdfile);
  my ($tag,@icontrol)=unpack("A4L*",$buffer);
  
  ($buffer,$len)=&GenUtil::readFortran($dcdfile);
  ($buffer,$len)=&GenUtil::readFortran($dcdfile);

  my $fixed=$icontrol[8];
  if ($fixed) {
    ($buffer,$len)=&GenUtil::readFortran($dcdfile);
  }

  my $crystal=$icontrol[10];

  if ($crystal) {
    my ($tbuf,$tlen)=&GenUtil::readFortran($dcdfile);
    my @cdat=unpack("d*",$tbuf);
    undef $dcdfile;

    my $a=sqrt($cdat[0]*$cdat[0]+$cdat[1]*$cdat[1]+$cdat[3]*$cdat[3]);
    my $b=sqrt($cdat[1]*$cdat[1]+$cdat[2]*$cdat[2]+$cdat[4]*$cdat[4]);
    my $c=sqrt($cdat[3]*$cdat[3]+$cdat[4]*$cdat[4]+$cdat[5]*$cdat[5]);

    my $ab=$cdat[1]*($cdat[0]+$cdat[2])+$cdat[3]*$cdat[4];
    my $bc=$cdat[4]*($cdat[2]+$cdat[5])+$cdat[1]*$cdat[3];
    my $ca=$cdat[3]*($cdat[0]+$cdat[5])+$cdat[1]*$cdat[4];

    my $a1=&GenUtil::acos($bc/($b*$c))*180.0/$GenUtil::pi;
    my $a2=&GenUtil::acos($ca/($c*$a))*180.0/$GenUtil::pi;
    my $a3=&GenUtil::acos($ab/($a*$b))*180.0/$GenUtil::pi;
    
    my $shape;
    if (abs($a-$b)<0.001 && abs($a-$c)<0.001) {
      if (abs($a1-109.4712)<0.001 && abs($a2-109.4712)<0.001 && abs($a3-109.4712)<0.001) {
	$shape="octahedron";
      } else {
	$shape="cubic";
      }
    } else {
      $shape="ortho";
    }

    return ($a,$b,$c,$shape);
  } else {
    undef $dcdfile;
    return (undef,undef,undef);
  }
}

