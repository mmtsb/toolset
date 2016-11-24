#!/usr/bin/env perl

# mutates nucleic acid residues in a given PDB files
#
# 2006, Michael Feig, MSU

sub usage {
  printf STDERR "usage:   mutateNA.pl [options] [PDBfile]\n";
  printf STDERR "options: [-seq [chain]index:[complementary:]sequence[=index:sequence]]\n";
  printf STDERR "         [-rna]\n";
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

my %nalookup = ( "A" => "ADE",
                 "T" => "THY",
                 "U" => "URA",
                 "C" => "CYT",
                 "G" => "GUA" );
my %base;

my %complDNA = ( "A" => "T",
                 "T" => "A",
                 "C" => "G",
                 "G" => "C" );

my %complRNA = ( "A" => "U",
                 "U" => "A",
                 "C" => "G",
                 "G" => "C" );

foreach my $b ( qw ( N9 C5 N7 C8 H8 N1 C2 H2 N3 C4 C6 N6 H61 H62 ) ) {
   $base{"ADE:".$b}=1;
}

foreach my $b ( qw ( N1 C6 H6 C2 O2 N3 H3 C4 O4 C5 C5M H51 H52 H53 ) ) {
   $base{"THY:".$b}=1;
}

foreach my $b ( qw ( N1 C6 H6 C2 O2 N3 H3 C4 O4 C5 H5 ) ) {
   $base{"URA:".$b}=1;
}

foreach my $b ( qw ( N9 C4 N2 H21 H22 N3 C2 N1 H1 C6 O6 C5 N7 C8 H8 ) ) {
   $base{"GUA:".$b}=1;
}

foreach my $b ( qw ( N1 C6 H6 C5 H5 C2 O2 N3 C4 N4 H41 H42 ) ) {
   $base{"CYT:".$b}=1;
}

my $rna=0;

my %mutate;
my $fname="-";

my $bp=();

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-rna") {
    shift @ARGV;
    $rna=1;
  } elsif ($ARGV[0] eq "-seq") {
    shift @ARGV;
    foreach my $t ( split(/=/,shift @ARGV) ) {
      my @f=split(/:/,$t);
      if ($#f==1) {
	my $tinx=$f[0];
	my $tseq=$f[1];
	my $tchain="";
	if ($tinx=~/([A-Z])([0-9]+)/) {
	  $tchain=$1;
	  $tinx=$2;
	}
	for (my $i=0; $i<length($tseq); $i++) {
	  my $key=sprintf("%s%d",$tchain,$i+$tinx); 
	  $mutate{$key}=substr($tseq,$i,1);
	}
      } elsif ($#f==2) {
	my $tinx=$f[0];
	my $cinx=$f[1];
	my $tseq=$f[2];
	my $tchain="";
        my $cchain="";
	if ($tinx=~/([A-Z])([0-9]+)/) {
	  $tchain=$1;
	  $tinx=$2;
	}
	if ($cinx=~/([A-Z])([0-9]+)/) {
	  $cchain=$1;
	  $cinx=$2;
	}
	for (my $i=0; $i<length($tseq); $i++) {
	  my $key=sprintf("%s%d",$tchain,$i+$tinx);
	  $mutate{$key}=substr($tseq,$i,1);
	  $key=sprintf("%s%d",$cchain,$cinx-$i);
	  $mutate{$key}=($rna)?$complRNA{substr($tseq,$i,1)}:$complDNA{substr($tseq,$i,1)};
	  my $bprec={};
	  $bprec->{chain1}=$tchain;
	  $bprec->{chain2}=$cchain;
	  $bprec->{inx1}=$tinx+$i;
	  $bprec->{inx2}=$cinx-$i;
	  push(@{$bp},$bprec);
	}
      }
    }
  } else {
    $fname = shift @ARGV;
  }
}

my $mol=Molecule::new();
$mol->readPDB($fname);

my $nmol=Molecule::new();

foreach my $c ( @{$mol->activeChains()} ) {
  my $nc=undef;

  foreach my $r (@{$c->{res}}) {
    $nc=$nmol->_newChain($c->{id}) if (!defined $nc);

    my $rrec={};
    %{$rrec}=%{$r};
    $rrec->{start}=$#{$nc->{atom}}+1;

    my $fromres;
    my $tores;
    my $tkey=sprintf("%s%d",$c->{id},$rrec->{num});
    if (exists $mutate{$tkey}) {
       if ($rrec->{name}=~/ADE|GUA|THY|URA|CYT/) {
          $fromres=$rrec->{name};
          $tores=$nalookup{$mutate{$tkey}};
          $rrec->{name}=$tores;
       } else {
          printf STDERR "wrong residue %s\n",$tkey;
       }
    }
        
    push(@{$nc->{res}},$rrec);


    if (defined $fromres && defined $tores && $fromres ne $tores) {
     my $firstbase;
     my $lastbase;
     for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
       my $at=$c->{atom}->[$ia];
       if (exists $base{$at->{resname}.":".$at->{atomname}}) {
         $firstbase=$ia if (!defined $firstbase);
         $lastbase=$ia;
       }
     }
     
     for (my $ia=$r->{start}; $ia<$firstbase; $ia++) {
      my $arec={};
      %{$arec}=%{$c->{atom}->[$ia]};
      $arec->{resname}=$tores;
      push(@{$nc->{atom}},$arec);
     }

     my $newbase=&makeBase($tores,$c->{atom}->[$firstbase]);

     &fitBase($newbase,$c->{atom},$firstbase,$lastbase,$fromres,$tores); 
     
     for (my $ia=0; $ia<=$#{$newbase}; $ia++) {
      my $arec={};
      %{$arec}=%{$newbase->[$ia]};
      $arec->{resname}=$tores;
      push(@{$nc->{atom}},$arec);
     }

     for (my $ia=$lastbase+1; $ia<=$r->{end}; $ia++) {
      my $arec={};
      %{$arec}=%{$c->{atom}->[$ia]};
      $arec->{resname}=$tores;
      push(@{$nc->{atom}},$arec);
     }
    } else {
     for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
      my $arec={};
      %{$arec}=%{$c->{atom}->[$ia]};
      push(@{$nc->{atom}},$arec);
     }
    }
    $rrec->{end}=$#{$nc->{atom}};
  }
}
  
$nmol->_coorCache();
#$nmol->writePDB(\*STDOUT,translate=>"GENERIC");

foreach my $b ( @{$bp} ) {
  &fitBasePair($nmol,$b);
}

$nmol->writePDB(\*STDOUT,translate=>"GENERIC");

1;

sub fitBasePair {
  my $m=shift;
  my $b=shift;

  my $c1=$m->getChain($b->{chain1});
  my $c2=$m->getChain($b->{chain2});
  my $r1=$m->getResidueInChain($b->{inx1},$c1);
  my $r2=$m->getResidueInChain($b->{inx2},$c2);
  my $a1=$c1->{atom};
  my $a2=$c2->{atom};

  my $cmx=0.0;
  my $cmy=0.0;
  my $cmz=0.0;
  my $mn=0;

  my $cpx=0.0;
  my $cpy=0.0;
  my $cpz=0.0;
  my @parr=();

  my @marr=();

  my $newbase=&makeBase($r1->{name},$a1->[$r1->{start}]);
  for (my $ia=$r1->{start}; $ia<=$r1->{end}; $ia++) {
    my $at=$a1->[$ia];
    if (exists $base{$at->{resname}.":".$at->{atomname}}) {
      $cmx+=$at->{xcoor};
      $cmy+=$at->{ycoor};
      $cmz+=$at->{zcoor};
      push(@marr,$at);
      $mn++;

      my $found=0;
      for my $bb ( @{$newbase} ) {
	if ($bb->{atomname} eq $at->{atomname}) {
	  $cpx+=$bb->{xcoor};
	  $cpy+=$bb->{ycoor};
	  $cpz+=$bb->{zcoor};
	  push(@parr,$bb);
	  $found++;
	}
      }
      
      die "base atom $at->{atomname} not found" if ($found==0);
      die "too many matching atoms found for $at->{atomname}" if ($found>1);
    }
  }

  $newbase=&makeBase($r2->{name},$a2->[$r2->{start}]);
  for (my $ia=$r2->{start}; $ia<=$r2->{end}; $ia++) {
    my $at=$a2->[$ia];
    if (exists $base{$at->{resname}.":".$at->{atomname}}) {
      $cmx+=$at->{xcoor};
      $cmy+=$at->{ycoor};
      $cmz+=$at->{zcoor};
      push(@marr,$at);
      $mn++;

      my $found=0;
      for my $bb ( @{$newbase} ) {
	if ($bb->{atomname} eq $at->{atomname}) {
	  $cpx+=$bb->{xcoor};
	  $cpy+=$bb->{ycoor};
	  $cpz+=$bb->{zcoor};
	  push(@parr,$bb);
	  $found++;
	}
      }
      
      die "base atom $at->{atomname} not found" if ($found==0);
      die "too many matching atoms found for $at->{atomname}" if ($found>1);
    }
  }
  $cmx/=$mn;
  $cmy/=$mn;
  $cmz/=$mn;

  $cpx/=$mn;
  $cpy/=$mn;
  $cpz/=$mn;

  my $r=();
  for (my $i=1; $i<=3; $i++) {
    $r->[$i]=();
    for (my $j=1; $j<=3; $j++) {
      $r->[$i]->[$j]=0.0;
    }
  }

  for (my $i=0; $i<=$#marr; $i++) {
    my $tc=$parr[$i];
    my $tr=$marr[$i];

    my $cx=$tc->{xcoor}-$cpx;
    my $cy=$tc->{ycoor}-$cpy;
    my $cz=$tc->{zcoor}-$cpz;

    my $rx=$tr->{xcoor}-$cmx;
    my $ry=$tr->{ycoor}-$cmy;
    my $rz=$tr->{zcoor}-$cmz;

    $r->[1]->[1]+=$cx*$rx;
    $r->[2]->[1]+=$cx*$ry;
    $r->[3]->[1]+=$cx*$rz;
    
    $r->[1]->[2]+=$cy*$rx;
    $r->[2]->[2]+=$cy*$ry;
    $r->[3]->[2]+=$cy*$rz;

    $r->[1]->[3]+=$cz*$rx;
    $r->[2]->[3]+=$cz*$ry;
    $r->[3]->[3]+=$cz*$rz;
  }

  my $u=&Analyze::_frotu($r);


  for (my $in=0; $in<=$#parr; $in++) {
    my $x=$parr[$in]->{xcoor}-$cpx;
    my $y=$parr[$in]->{ycoor}-$cpy;
    my $z=$parr[$in]->{zcoor}-$cpz;

    my $tx=$u->[1]->[1]*$x+$u->[1]->[2]*$y+$u->[1]->[3]*$z;
    my $ty=$u->[2]->[1]*$x+$u->[2]->[2]*$y+$u->[2]->[3]*$z;
    my $tz=$u->[3]->[1]*$x+$u->[3]->[2]*$y+$u->[3]->[3]*$z;

    $marr[$in]->{xcoor}=$tx+$cmx;
    $marr[$in]->{ycoor}=$ty+$cmy;
    $marr[$in]->{zcoor}=$tz+$cmz;
  }
}

sub makeBase {
  my $base=shift;
  my $tat=shift;

  my $arr=();
  
  if ($base eq "CYT") {
    push(@{$arr},&makeBaseRec("N1",1.098,-21.696,2.799,$tat));
    push(@{$arr},&makeBaseRec("C6",0.914,-20.388,2.566,$tat));
    push(@{$arr},&makeBaseRec("H6",1.795,-19.753,2.372,$tat));
    push(@{$arr},&makeBaseRec("C5",-0.295,-19.896,2.399,$tat));
    push(@{$arr},&makeBaseRec("H5",-0.460,-18.824,2.167,$tat));
    push(@{$arr},&makeBaseRec("C2",0.021,-22.566,3.028,$tat));
    push(@{$arr},&makeBaseRec("O2",0.207,-23.741,3.430,$tat));
    push(@{$arr},&makeBaseRec("N3",-1.246,-22.138,2.880,$tat));
    push(@{$arr},&makeBaseRec("C4",-1.471,-20.832,2.541,$tat));
    push(@{$arr},&makeBaseRec("N4",-2.771,-20.418,2.518,$tat));
    push(@{$arr},&makeBaseRec("H41",-3.320,-21.098,2.991,$tat));
    push(@{$arr},&makeBaseRec("H42",-3.016,-19.441,2.550,$tat));
  } elsif ($base eq "GUA") {
    push(@{$arr},&makeBaseRec("N9",-6.193,-26.773,3.721,$tat));
    push(@{$arr},&makeBaseRec("C4",-5.131,-26.026,3.305,$tat));
    push(@{$arr},&makeBaseRec("N2",-1.678,-25.666,2.980,$tat));
    push(@{$arr},&makeBaseRec("H21",-1.246,-26.543,2.871,$tat));
    push(@{$arr},&makeBaseRec("H22",-1.118,-24.841,2.984,$tat));
    push(@{$arr},&makeBaseRec("N3",-3.896,-26.455,3.095,$tat));
    push(@{$arr},&makeBaseRec("C2",-2.985,-25.463,3.017,$tat));
    push(@{$arr},&makeBaseRec("N1",-3.368,-24.087,3.032,$tat));
    push(@{$arr},&makeBaseRec("H1",-2.714,-23.292,2.885,$tat));
    push(@{$arr},&makeBaseRec("C6",-4.705,-23.642,3.046,$tat));
    push(@{$arr},&makeBaseRec("O6",-4.933,-22.470,2.919,$tat));
    push(@{$arr},&makeBaseRec("C5",-5.617,-24.724,3.259,$tat));
    push(@{$arr},&makeBaseRec("N7",-6.923,-24.670,3.469,$tat));
    push(@{$arr},&makeBaseRec("C8",-7.216,-25.897,3.700,$tat));
    push(@{$arr},&makeBaseRec("H8",-8.217,-26.207,4.014,$tat));
  } elsif ($base eq "ADE") {
    push(@{$arr},&makeBaseRec("N9",-3.909,-27.336,-6.573,$tat));
    push(@{$arr},&makeBaseRec("C5",-3.059,-25.332,-7.088,$tat));
    push(@{$arr},&makeBaseRec("N7",-2.040,-26.231,-7.341,$tat));
    push(@{$arr},&makeBaseRec("C8",-2.591,-27.372,-6.973,$tat));
    push(@{$arr},&makeBaseRec("H8",-1.998,-28.305,-6.916,$tat));
    push(@{$arr},&makeBaseRec("N1",-4.277,-23.367,-7.088,$tat));
    push(@{$arr},&makeBaseRec("C2",-5.294,-24.153,-6.792,$tat));
    push(@{$arr},&makeBaseRec("H2",-6.210,-23.563,-6.652,$tat));
    push(@{$arr},&makeBaseRec("N3",-5.388,-25.479,-6.492,$tat));
    push(@{$arr},&makeBaseRec("C4",-4.216,-26.010,-6.700,$tat));
    push(@{$arr},&makeBaseRec("C6",-3.105,-23.916,-7.206,$tat));
    push(@{$arr},&makeBaseRec("N6",-2.070,-23.108,-7.544,$tat));
    push(@{$arr},&makeBaseRec("H61",-1.187,-23.562,-7.571,$tat));
    push(@{$arr},&makeBaseRec("H62",-2.074,-22.174,-7.176,$tat));
  } elsif ($base eq "THY") {
    push(@{$arr},&makeBaseRec("N1",-6.370,-19.065,-7.508,$tat));
    push(@{$arr},&makeBaseRec("C6",-5.485,-18.081,-7.322,$tat));
    push(@{$arr},&makeBaseRec("H6",-5.828,-17.081,-7.467,$tat));
    push(@{$arr},&makeBaseRec("C2",-6.037,-20.402,-7.352,$tat));
    push(@{$arr},&makeBaseRec("O2",-6.856,-21.298,-7.401,$tat));
    push(@{$arr},&makeBaseRec("N3",-4.711,-20.569,-7.117,$tat));
    push(@{$arr},&makeBaseRec("H3",-4.371,-21.482,-6.917,$tat));
    push(@{$arr},&makeBaseRec("C4",-3.755,-19.613,-6.879,$tat));
    push(@{$arr},&makeBaseRec("O4",-2.620,-19.994,-6.654,$tat));
    push(@{$arr},&makeBaseRec("C5",-4.164,-18.268,-7.049,$tat));
    push(@{$arr},&makeBaseRec("C5M",-3.152,-17.077,-6.992,$tat));
    push(@{$arr},&makeBaseRec("H51",-2.160,-17.254,-6.457,$tat));
    push(@{$arr},&makeBaseRec("H52",-3.616,-16.386,-6.273,$tat));
    push(@{$arr},&makeBaseRec("H53",-2.976,-16.525,-7.956,$tat));
  } elsif ($base eq "URA") {
    push(@{$arr},&makeBaseRec("N1",-6.416,-19.053,-7.506,$tat));
    push(@{$arr},&makeBaseRec("C6",-5.477,-18.057,-7.334,$tat));
    push(@{$arr},&makeBaseRec("H6",-5.839,-17.046,-7.485,$tat));
    push(@{$arr},&makeBaseRec("C2",-6.065,-20.390,-7.338,$tat));
    push(@{$arr},&makeBaseRec("O2",-6.877,-21.301,-7.469,$tat));
    push(@{$arr},&makeBaseRec("N3",-4.721,-20.607,-7.016,$tat));
    push(@{$arr},&makeBaseRec("H3",-4.438,-21.569,-6.902,$tat));
    push(@{$arr},&makeBaseRec("C4",-3.715,-19.639,-6.892,$tat));
    push(@{$arr},&makeBaseRec("O4",-2.544,-19.955,-6.677,$tat));
    push(@{$arr},&makeBaseRec("C5",-4.185,-18.288,-7.049,$tat));
    push(@{$arr},&makeBaseRec("H5",-3.471,-17.469,-6.969,$tat));
  } else {
    die "unknown base $base";
  }
  return $arr;
}

sub makeBaseRec {
  my $aname=shift;
  my $x=shift;
  my $y=shift;
  my $z=shift;
  my $at=shift;
  
  my $nat={};
  
  %{$nat}=%{$at};
  $nat->{atomname}=$aname;
  $nat->{xcoor}=$x;
  $nat->{ycoor}=$y;
  $nat->{zcoor}=$z;

  return $nat;
}

sub fitBase {
  my $nbase=shift;
  my $at=shift;
  my $afrom=shift;
  my $ato=shift;
  my $rfrom=shift;
  my $rto=shift;

  my @narr=();
  my @oarr=();

  if (($rfrom eq "ADE" || $rfrom eq "GUA")  && ($rto eq "GUA" || $rto eq "ADE")) {
    push(@oarr,qw ( N9 C4 C2 C6 C5 N7 ));
    push(@narr,qw ( N9 C4 C2 C6 C5 N7 ));
  } elsif (($rfrom eq "CYT" || $rfrom eq "THY" || $rfrom eq "URA") &&
	   ($rto eq "CYT" || $rto eq "THY" || $rto eq "URA")) {
    push(@oarr,qw ( N1 C6 C5 C4 N3 C2 ));
    push(@narr,qw ( N1 C6 C5 C4 N3 C2 ));
  } elsif (($rfrom eq "ADE" || $rfrom eq "GUA") &&
	   ($rto eq "CYT" || $rto eq "THY" || $rto eq "URA")) {
    push(@oarr,qw ( N9 C2 C5 C6 C8 C2 ));
    push(@narr,qw ( N1 N3 C5 C4 C6 N3 ));
  } elsif (($rfrom eq "CYT" || $rfrom eq "THY" || $rfrom eq "URA") &&
	   ($rto eq "GUA" || $rto eq "ADE")) {
    push(@narr,qw ( N9 C2 C5 C6 C8 C2 ));
    push(@oarr,qw ( N1 N3 C5 C4 C6 N3 ));
  } else {
    die "unknown base conversion $rfrom to $rto";
  }

  my %cmp;
  my %ref;

  foreach my $a ( @{$nbase} ) {
    $cmp{$a->{atomname}}=$a;
  }

  for (my $i=$afrom; $i<=$ato; $i++) {
    my $a=$at->[$i];
    $ref{$a->{atomname}}=$a;
  }

  my $cnx=0.0;
  my $cny=0.0;
  my $cnz=0.0;

  my $cox=0.0;
  my $coy=0.0;
  my $coz=0.0;

  my $nn=0;
  my $on=0;
  
  foreach my $t (@narr) {
    die "cmp atom $t does not exist" if (!exists $cmp{$t});
    $cnx+=$cmp{$t}->{xcoor};
    $cny+=$cmp{$t}->{ycoor};
    $cnz+=$cmp{$t}->{zcoor};
    $nn++;
  }

  foreach my $t (@oarr) {
    die "ref atom $t does not exist" if (!exists $ref{$t});
    $cox+=$ref{$t}->{xcoor};
    $coy+=$ref{$t}->{ycoor};
    $coz+=$ref{$t}->{zcoor};
    $on++;
  } 

  $cnx/=$nn;
  $cny/=$nn;
  $cnz/=$nn;

  $cox/=$on;
  $coy/=$on;
  $coz/=$on;

  my $r=();
  for (my $i=1; $i<=3; $i++) {
    $r->[$i]=();
    for (my $j=1; $j<=3; $j++) {
      $r->[$i]->[$j]=0.0;
    }
  }

  for (my $i=0; $i<=$#narr; $i++) {
    my $tc=$cmp{$narr[$i]};
    my $tr=$ref{$oarr[$i]};

    my $cx=$tc->{xcoor}-$cnx;
    my $cy=$tc->{ycoor}-$cny;
    my $cz=$tc->{zcoor}-$cnz;

    my $rx=$tr->{xcoor}-$cox;
    my $ry=$tr->{ycoor}-$coy;
    my $rz=$tr->{zcoor}-$coz;

    $r->[1]->[1]+=$cx*$rx;
    $r->[2]->[1]+=$cx*$ry;
    $r->[3]->[1]+=$cx*$rz;
    
    $r->[1]->[2]+=$cy*$rx;
    $r->[2]->[2]+=$cy*$ry;
    $r->[3]->[2]+=$cy*$rz;

    $r->[1]->[3]+=$cz*$rx;
    $r->[2]->[3]+=$cz*$ry;
    $r->[3]->[3]+=$cz*$rz;
  }

  my $u=&Analyze::_frotu($r);

  my $glycoatname=($rto eq "ADE" || $rto eq "GUA")?"N9":"N1";
  my $rglycoatname=($rfrom eq "ADE" || $rfrom eq "GUA")?"N9":"N1";
  my $glycoat;

  for (my $in=0; $in<=$#{$nbase}; $in++) {
    my $x=$nbase->[$in]->{xcoor}-$cnx;
    my $y=$nbase->[$in]->{ycoor}-$cny;
    my $z=$nbase->[$in]->{zcoor}-$cnz;

    my $tx=$u->[1]->[1]*$x+$u->[1]->[2]*$y+$u->[1]->[3]*$z;
    my $ty=$u->[2]->[1]*$x+$u->[2]->[2]*$y+$u->[2]->[3]*$z;
    my $tz=$u->[3]->[1]*$x+$u->[3]->[2]*$y+$u->[3]->[3]*$z;

    $nbase->[$in]->{xcoor}=$tx+$cox;
    $nbase->[$in]->{ycoor}=$ty+$coy;
    $nbase->[$in]->{zcoor}=$tz+$coz;
  }
}
