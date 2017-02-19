#!/usr/bin/env perl

# builds homology model using various tools
#
# 2006, Michael Feig, Michigan State University

sub usage {
  printf STDERR "usage:   buildModel.pl [options] [fastaFile]\n";
  printf STDERR "options:               [-models value]\n";
  printf STDERR "                       [-maxloop value]\n";
  printf STDERR "                       [-verbose]\n";
  printf STDERR "                       [-nocenter]\n";
  printf STDERR "                       [-stats file]\n";
  printf STDERR "                       [-pdbdir name]\n";
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
use Analyze;
use Sequence;
use SICHO;

my $seqname;
my $verbose=0;
my $maxloop=12;
my $models=5; 
my $center=1;
my $statsfile;

my $pdbdir;
$pdbdir=$ENV{'PDBDIR'} if (defined $ENV{'PDBDIR'});

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-models") {
    shift @ARGV;
    $models=shift @ARGV;
  } elsif ($ARGV[0] eq "-maxloop") {
    shift @ARGV;
    $maxloop=shift @ARGV;
  } elsif ($ARGV[0] eq "-verbose") {
    shift @ARGV;
    $verbose=1;
  } elsif ($ARGV[0] eq "-nocenter") {
    shift @ARGV;
    $center=0;
  } elsif ($ARGV[0] eq "-stats") {
    shift @ARGV;
    $statsfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-pdbdir") {
    shift @ARGV;
    $pdbdir=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $seqname = shift @ARGV;
    $done=1;
  }
}

$pdbdir="." if (!defined $pdbdir);

my $tag=$$;

my $inp=&GenUtil::getInputFile($seqname);

my $pdbid;
my $pdbname;
my $chain;

<$inp>;
my $target=<$inp>;
chomp $target;

printf STDERR "%s\n",$target if ($verbose);

my $line=<$inp>;

if ($line=~/^>([\S]+\.pdb)[\_|]([A-Za-z0-9]?)/) {
  $pdbname=$1;
  $chain=$2;
} elsif ($line=~/^>([\S]+\.pdb)\s*/) {
  $pdbname=$1;
  $chain="";
} elsif ($line=~/^>([A-Za-z0-9]+)[\_|]*([A-Za-z0-9]?)/) {
  $pdbid=uc $1;
  $chain=$2;
} else {
  die "input format error";
}

my $alignment=<$inp>;
chomp $alignment;
undef $inp;

printf STDERR "%s\n",$alignment if ($verbose);

my $fname;

if (defined $pdbname) {
  $fname=$pdbname;
} elsif (-d $pdbdir) {
  my $tname=sprintf("%s/pdb%s.ent",$pdbdir,lc $pdbid);
  if (-r $tname) {
    $fname=$tname;
  } else {
    $tname=sprintf("%s/%s.pdb",$pdbdir,uc $pdbid);
    if (-r $tname) {
      $fname=$tname;
    }
  }
}

die "cannot read template PDB file" if (!defined $fname || !-r $fname);

my $mol=Molecule::new();
$mol->readPDB($fname);

if (defined $chain && $chain ne "") {
  $mol->setValidChain($chain);
}

my $tmol=$mol->clone(1);

if ($tmol->empty()) {
  $mol->readPDB($fname);
  $mol->setValidSelection("solute");
  $tmol=$mol->clone(1);
}

$tmol->setChain("");
$tmol->center() unless (!$center);
my $slen=length($alignment);

my ($pdbseq,$pdbmin,$pdbmax)=&getSeq($tmol);

my $alignseq=$alignment;
$alignseq=~s/[^A-Z]//g;

my ($amatch,$avseq)=&match($alignseq,$pdbseq);

printf STDERR "PDB match: %d\n",$amatch+$pdbmin if ($verbose);
printf STDERR "match seq: %s\n",$avseq if ($verbose);

my $valignment=$alignment;
my $iam=0;
for (my $ia=0; $ia<$slen; $ia++) {
  if (substr($valignment,$ia,1)=~/[A-Z]/) {
    substr($valignment,$ia,1)=substr($avseq,$iam++,1);
  }
}

printf STDERR "%s\n",$valignment if ($verbose);

my $last=-1;

my @frags=();

my $alen=length($valignment);

my $last=-1;
my $lastm=-1;
my $lastt=-1;
my $im=0;
my $it=0;

$valignment=~s/X/8/g;

for (my $is=0; $is<$alen; $is++) {
  my $ca=substr($valignment,$is,1);
  my $ct=substr($target,$is,1);
#  printf STDERR "is %d %s %s\n",$is,$ca,$ct;

  if ($last<0 &&
      $ca=~/[A-Z]/ && $ct=~/[A-Z]/ &&
      ($is==0 || substr($target,$is-1,1)=~/[A-Z]/)) { 
    $last=$is;
    $lastm=$im;
    $lastt=$it;
#    printf STDERR "last %d %d %d\n",$is,$im,$it;
  }
  
  if ($last>=0 && 
      ($ca!~/[A-Z]/ || $ct!~/[A-Z]/ || $is==$alen-1 || 
       substr($target,$is+1,1)!~/[A-Z]/)) {
#    printf STDERR "new\n";
    my $frec={};
    
    $frec->{afrom}=$last;
    $frec->{pfrom}=$lastm+$amatch+$pdbmin;
    $frec->{from}=$lastt;
    $frec->{tfrom}=$lastm;
#    if ($is==$alen-1) {
#      $frec->{aseq}=substr($valignment,$last,$is-$last+1);
#      $frec->{tseq}=substr($target,$last,$is-$last+1);
#      $frec->{ato}=$is;
#      $frec->{to}=$it;
#      $frec->{pto}=$im+$amatch+$pdbmin;
#      $frec->{tto}=$im;
#    } else {
      $frec->{aseq}=uc substr($valignment,$last,$is-$last);
      $frec->{tseq}=uc substr($target,$last,$is-$last);
      $frec->{ato}=$is-1;
      $frec->{to}=$it-1;
      $frec->{pto}=$im+$amatch+$pdbmin-1;
      $frec->{tto}=$im-1;
#    }

#    printf STDERR "added %d %d %d %d\n",$frec->{afrom},$frec->{pfrom},$frec->{from},$frec->{tfrom};
#    printf STDERR "added %d %d %d %d\n",$frec->{ato},$frec->{pto},$frec->{to},$frec->{tto};
    push (@frags,$frec);

    $last=-1;
  }
  $im++ if ($ca=~/[A-Z8]/);
  $it++ if ($ct=~/[A-Z]/);
}

my @have;
foreach my $f (@frags) {
#  printf STDERR "frag %d %d\n",$f->{from},$f->{to};
  if ($f->{to}-$f->{from}>0) {
   my $srec={};
   $srec->{from}=$f->{pfrom};
   $srec->{to}=$f->{pto};
   my $sarr=();
   push(@{$sarr},$srec); 

   $tmol->setValidResidues($sarr);
   my $smol=$tmol->clone(1);

   printf STDERR "%s (%d-%d) -> %s (%d-%d)\n",
    $f->{aseq},$srec->{from},$srec->{to},$f->{tseq},$f->{from}+1,$f->{to}+1 
   if ($verbose);

   my $sicho=SICHO::new(resolution=>1.0,intflag=>0);
   $sicho->genSimpleFromAllAtom($smol,ca=>1);

   my $rmol=Molecule::new();
   my $tseq=Sequence::new($f->{tseq});
   $rmol->rebuildFromSICHO($tseq,$sicho,undef,undef,1);

   $rmol->renumber($f->{from}+1);

   $f->{mol}=$rmol;
   push(@have,$f);
  }
}

die "no fragments found" if ($#have<0);

printf STDERR "******** assembling structure\n" if ($verbose);
my $nmol=Molecule::new();

my $natom=1;

$nmol->{chain}=();
$nmol->{chainlookup}={};
$nmol->{defchain}=undef;
$nmol->{ssbond}=();

my $chainrec=$nmol->_newChain("");

$target=~s/[^A-Z]//g;
for (my $i=0; $i<length($target); $i++) {
  my $resname=$Sequence::_seqlong{substr($target,$i,1)};
  my $resnum=$natom;
  my $resrec={ name => $resname, num => $resnum, start => $natom-1, chain=>"", valid=>1 };
  push (@{$chainrec->{res}},$resrec);
 
  my $asc={ atominx => $natom++, atomname => "CA", 
            resname => $resname, resnum  => $resnum, chain=> "",
            xcoor   => 0.0, ycoor => 0.0, zcoor => 0.0, valid=>1 };
  push (@{$chainrec->{atom}},$asc);

  $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
}

$nmol->_coorCache();

my @missing=();
my $firstpdb=1;
my $lastpdb=length($target);

my $last=0;
my $tlast=0;
my $lastmol;
foreach my $h (sort {$a->{from}<=>$b->{from}} @have) {
  $nmol->merge($h->{mol});
  if ($h->{from}>$last) {
    my $mrec={};
    $mrec->{from}=$last+1;
    $mrec->{to}=$h->{from};

    if ($mrec->{from}==1) {
      $firstpdb=$mrec->{to}+1 
    } else {
      push(@missing,$mrec);
    }

    my $d=0.0;
    if ($last>0) {
      my $c1=$h->{mol}->getChain();
      my $c2=$lastmol->getChain();
      my $ca1=$h->{mol}->getAtomInResidue($c1->{res}->[0],"CA",$c1);
      my $ca2=$lastmol->getAtomInResidue($c2->{res}->[$#{$c2->{res}}],"CA",$c2);
      my $dx=$ca1->{xcoor}-$ca2->{xcoor};
      my $dy=$ca1->{ycoor}-$ca2->{ycoor};
      my $dz=$ca1->{zcoor}-$ca2->{zcoor};
      $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);

      if ($tlast>$h->{tfrom}-1) {
        printf STDERR "%d:%d %f -:-\n",$last+1,$h->{from},$d if ($verbose);
      } else {
        printf STDERR "%d:%d %f %d:%d\n",$last+1,$h->{from},$d,$tlast+1,$h->{tfrom} if ($verbose);
      }
    } else {
      printf STDERR "%d:%d\n",$last+1,$h->{from} if ($verbose);
    }
  }
  $last=$h->{to}+1;
  $tlast=$h->{tto}+1;
  $lastmol=$h->{mol};
}

if ($last<length($target)) {
  printf STDERR "%d:%d\n",$last,length($target),0.0 if ($verbose);
  my $mrec={};
  $mrec->{from}=$last;
  $mrec->{to}=length($target);
#  push(@missing,$mrec);
  $lastpdb=$last-1;
}

my $targetpdb=sprintf("%s.pdb",$tag);

my @slist=();
my $srec={};
$srec->{from}=$firstpdb;
$srec->{to}=$lastpdb;
push(@slist,$srec);
$nmol->setValidResidues(\@slist);
my $sbmol=$nmol->clone(1);



printf STDERR "******** building loops\n" if ($verbose);
my @mfrags=();
foreach my $m ( @missing ) {
  if ($m->{to}-$m->{from}+3<=$maxloop) {
    printf STDERR "loop %d:%d will be built\n",$m->{from},$m->{to} if ($verbose);
    my $rec={};
    $rec->{from}=$m->{from}-1;
    $rec->{to}=$m->{to}+1;
    push(@mfrags,$rec);
    $m->{modeled}=1;
  } else {
    printf STDERR "loop %d:%d too large\n",$m->{from},$m->{to} if ($verbose);
    $m->{modeled}=0;
  } 
} 

my @notbuilt=();
foreach my $m ( @missing ) {
  if (!$m->{modeled}) {
    push(@notbuilt,$m);
  }
}

my $dir="$$.loops";
if ($#mfrags>=0) {
  system "mkdir $dir";

  $sbmol->writePDB("$dir/$$.pdb");

  open MOD,">$dir/$$.top";
  printf MOD "INCLUDE\n";
  printf MOD "SET SEQUENCE='$$'\n";
  printf MOD "SET LOOP_MODEL='$$.pdb'\n";
  printf MOD "SET LOOP_STARTING_MODEL=1\n";
  printf MOD "SET LOOP_ENDING_MODEL=$models\n";
  printf MOD "CALL ROUTINE='loop'\n";
  printf MOD "\n";
  printf MOD "SUBROUTINE ROUTINE='select_loop_atoms'\n";
  for (my $i=0; $i<=$#mfrags; $i++) {
    if ($i==0) {
      printf MOD "  PICK_ATOMS SELECTION_SEGMENT = '%d:' '%d:', SELECTION_STATUS = 'initialize'\n",
	$mfrags[$i]->{from},$mfrags[$i]->{to};
    } else {
      printf MOD "  PICK_ATOMS SELECTION_SEGMENT = '%d:' '%d:', SELECTION_STATUS = 'add'\n",
	$mfrags[$i]->{from},$mfrags[$i]->{to};
    }
  }
  printf MOD "RETURN\nEND_SUBROUTINE\n";
  close MOD;

  my $modexec;
  if ($ENV{'MODELLEREXEC'} ne "") {
     $modexec=$ENV{'MODELLEREXEC'};
  } else {
     $modexec=&GenUtil::findExecutable("mod9v7");
     $modexec=&GenUtil::findExecutable("mod8v2") if (!defined $modexec);
     $modexec=&GenUtil::findExecutable("mod8v1") if (!defined $modexec);
     $modexec=&GenUtil::findExecutable("mod9v1") if (!defined $modexec);
     $modexec=&GenUtil::findExecutable("mod9v2") if (!defined $modexec);
  }

  die "cannot find MODELLER exectuable" if (!defined $modexec);
 
  system "cd $dir; $modexec $$.top";

  my @model=();
  for (my $i=1; $i<=$models; $i++) {
    system sprintf("mv $dir/$$.BL%04d0001.pdb $dir/model.%d.pdb\n",$i,$i);
    push(@model,0);
  }
  
  my $mi=0;
  open INP,"$dir/$$.log";
  while (<INP>) {
    chomp;
    if (/Differences between/) {
      $mi++;
    } elsif (/Current energy\s+:\s+(\S+)/) {
      $model[$mi-1]=$1;
    }
  }
  close INP;

  my $minscore=999999999;
  my $mininx=-1;
  for (my $i=1; $i<=$models; $i++) {
    printf STDERR "model %d: %f\n",$i,$model[$i-1] if ($verbose);
    if ($model[$i-1]<$minscore && ($model[$i-1]>0.0001 || $model[$i-1]<-0.0001) &&
	-r "$dir/model.$i.pdb") {
      $minscore=$model[$i-1];
      $mininx=$i;
    }
  }

  my $rmol=Molecule::new();
  $rmol->readPDB(sprintf("%s/model.%d.pdb",$dir,$mininx));
  if ($#notbuilt>=0) {
    $rmol->setValidResidues(\@notbuilt,1,1);
    my $nmol=$rmol->clone(1);  
    $nmol->writePDB(\*STDOUT);
  } else {
    $rmol->writePDB(\*STDOUT);
  }
} else {
  $sbmol->writePDB(\*STDOUT);
}

my $statsout=(defined $statsfile)?&GenUtil::getOutputFile($statsfile):\*STDERR;

printf $statsout "built %d:%d\n",$firstpdb,$lastpdb;
printf $statsout "missing N-term %d:%d\n",1,$firstpdb-1 if ($firstpdb>1);
printf $statsout "missing C-term %d:%d\n",$lastpdb+1,length($target) 
       if ($lastpdb<length($target));
foreach my $m ( @missing ) {
  if (!$m->{modeled}) {
    printf $statsout "missing frag %d-%d\n",$m->{from},$m->{to};
  }
}
undef $statsout if (defined $statsfile);

&GenUtil::remove("$$.pdb") if (-r "$$.pdb");
system "rm -rf $dir" if (-d $dir);

1;

sub getSeq {
  my $m=shift;

  my $r=$m->{chain}->[0]->{res};
  my $a=$m->{chain}->[0]->{atom};

  my $nmin=$r->[0]->{num};
  my $nmax=$r->[$#{$r}]->{num};

  my $nrange=$nmax-$nmin+1;
 
  my $arr=(); 
  for (my $i=0; $i<$nrange; $i++) {
    push(@{$arr},"X");
  }

  foreach my $ir ( @{$r} ) {
    my $c=$Sequence::_seqabbrev{$ir->{name}};

    if (defined $c && $c ne "") {
      my $haveca=0;
      for (my $ia=$ir->{start}; $ia<=$ir->{end}; $ia++) {
        if ($a->[$ia]->{atomname} eq "CA") {
          $haveca=1;
        }
      }
      if (!$haveca) {
        $c="X";       
      }
    } else {
      $c="X";
    }
    $arr->[$ir->{num}-$nmin]=$c;
  }

  return (join("",@{$arr}),$nmin,$nmax);
}

sub match {
  my $a=shift;
  my $p=shift;

  my $alen=length($a);
  my $plen=length($p);

  my $maxmatch=0;
  my $maxinx=0;

  for (my $s=-$alen+1; $s<$plen+$alen-1; $s++) {
    my $matched=0;
    for (my $i=0; $i<$alen; $i++) {
      if ($i+$s>=0 && $i+$s<$plen && substr($a,$i,1) eq substr($p,$i+$s,1)) {
        $matched++;
      }
    } 
    if ($matched>$maxmatch) {
      $maxmatch=$matched;
      $maxinx=$s;
    } 
  }
 
  my $vseq=$a;
  
  for (my $i=0; $i<$alen; $i++) {
    if ($i+$maxinx<0 || $i+$maxinx>=$plen || substr($a,$i,1) ne substr($p,$i+$maxinx,1)) {
      substr($vseq,$i,1)="8";
    }
  }

  return ($maxinx,$vseq);
}
