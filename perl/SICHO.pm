# SICHO model package
# read/write/generate SICHO models
#
# http://mmtsb.scripps.edu/doc/SICHO.pm.html
# 2000, Michael Feig, Brooks group, TSRI

package SICHO;

require 5.004;

use strict;

use FileHandle;
use IPC::Open2;

use GenUtil;
use Molecule;

## data: sidechain[] -> {xcoor, ycoor, zcoor}
## side chain coordinates

## data: ca[] -> {xcoor, ycoor, zcoor}
## C-alpha coordinates (optional)

## data: resolution
## lattice resolution

## data: offset -> {xcoor, ycoor, zcoor}
## lattice offset

## data: gridsize
## lattice gridsize

## data: integer
## flag selecting integer or floating point coordinates

## constructor: new([parameters])
## creates a new SICHO object, the following parameters
## may be set through hash-style key=>value pairs:
## <mark>resolution</mark>, <mark>offsetx</mark>, <mark>offsety</mark>, <mark>offsetz</mark>,
## <mark>gridsize</mark>, <mark>intflag</mark>

sub new {
  my %par=@_;
  my $self={};

  $self->{sidechain}=();
  $self->{ca}=();
  $self->{resolution}  = (defined $par{resolution}) ? $par{resolution} : 1.45;
  $self->{offset}->{xcoor} = (defined $par{offsetx})    ? $par{offsetx}    : 50;
  $self->{offset}->{ycoor} = (defined $par{offsety})    ? $par{offsety}    : 50;
  $self->{offset}->{zcoor} = (defined $par{offsetz})    ? $par{offsetz}    : 50;
  $self->{gridsize}    = (defined $par{gridsize})   ? $par{gridsize}   : 100;
  $self->{integer}     = (defined $par{intflag})    ? $par{intflag}    : 1;

  bless($self);
  return $self;
}

## method: copy(fromSICHO)
## copies the data from another SICHO object

sub copy {
  my $self=shift;
  my $from=shift;
  @{$self->{sidechain}}=@{$from->{sidechain}};
  @{$self->{ca}}=@{$from->{ca}} if (defined $from->{ca}->[0]);
  $self->{resolution}=$from->{resolution};
  %{$self->{offset}}=%{$from->{offset}};
  $self->{integer}=$from->{integer};
  $self->{gridsize}=$from->{gridsize};
}

## method: readChain(file)
## reads a SICHO chain file

sub readChain {
  my $self=shift;
  my $chainfile=&GenUtil::getInputFile(shift);

  $self->{sidechain}=();
  $self->{ca}=();

  my $num=<$chainfile>;

  return if ($num<=0);

  while (<$chainfile>) {
    s/^ +//;
    my @f=split(/ +/);
    my $rec={ xcoor => $f[0], ycoor => $f[1], zcoor => $f[2] };
    push (@{$self->{sidechain}},$rec);
  }

  undef $chainfile;
}

## method: writeChain(file)
## writes a SICHO chain file

sub writeChain {
  my $self=shift;
  my $chainfile=&GenUtil::getOutputFile(shift);

  printf $chainfile " %5d\n",$#{$self->{sidechain}}+1;
  for (my $i=0; $i<=$#{$self->{sidechain}}; $i++) {
    if ($self->{integer}) {
      my $vx=$self->{sidechain}->[$i]->{xcoor};
      my $vy=$self->{sidechain}->[$i]->{ycoor};
      my $vz=$self->{sidechain}->[$i]->{zcoor};
      $vx+=($vx>0.0)?0.5:-0.5;
      $vy+=($vy>0.0)?0.5:-0.5;
      $vz+=($vz>0.0)?0.5:-0.5;
      printf $chainfile " %5d%5d%5d\n",int($vx),int($vy),int($vz);
    } else {
      printf $chainfile " %10.5f%10.5f%10.5f\n",
      $self->{sidechain}->[$i]->{xcoor},$self->{sidechain}->[$i]->{ycoor},
      $self->{sidechain}->[$i]->{zcoor};
    }
  }

  undef $chainfile;
}

## method: convert(parameters)
## translates the lattice chain representation of the current
## structure to a new lattice resolution and/or with a new offset.
## Parameters of the new representation are given as in the
## constructor.

sub convert {
  my $self=shift;
  my %par=@_;

  my $dx = (defined $par{offsetx} ) ? $par{offsetx} : $self->{offset}->{xcoor};
  my $dy = (defined $par{offsety} ) ? $par{offsety} : $self->{offset}->{ycoor};
  my $dz = (defined $par{offsetz} ) ? $par{offsetz} : $self->{offset}->{zcoor};

  my $resolution = (defined $par{resolution} ) ? 
    $par{resolution} : $self->{resolution};

  my (@x,@y,@z);
  my (@cax,@cay,@caz);

  for (my $i=0; $i<=$#{$self->{sidechain}}; $i++) {
    ($x[$i],$y[$i],$z[$i])=$self->fromProjection($i);
    ($cax[$i],$cay[$i],$caz[$i])=$self->fromProjection($i,1) 
      if (defined $self->{ca}->[$i]);
  }

  $self->{resolution}=$resolution;
  $self->{offset}->{xcoor}=$dx;
  $self->{offset}->{ycoor}=$dy;
  $self->{offset}->{zcoor}=$dz;

  for (my $i=0; $i<=$#{$self->{sidechain}}; $i++) {
    $self->toProjection($x[$i],$y[$i],$z[$i],$i);
    $self->toProjection($cax[$i],$cay[$i],$caz[$i],$i,1) 
      if (defined $self->{ca}->[$i]);
  }
}

## method: getSimpleFromAllAtom(mol[,parameters])
## generates a simple lattice chain from an all atom model
## represented as a Molecule object
## Parameters for the new chain can be given as in the
## constructor.

sub genSimpleFromAllAtom {
  my $self=shift;
  my $mol=shift;
  my %par=@_;
  
  $self->{offset}->{xcoor}=$par{offsetx}    if (defined $par{offsetx});
  $self->{offset}->{ycoor}=$par{offsety}    if (defined $par{offsety});
  $self->{offset}->{zcoor}=$par{offsetz}    if (defined $par{offsetz});
  $self->{resolution} =$par{resolution}     if (defined $par{resolution});

  $self->{sidechain}=();
  $self->{ca}=();

  my $wantca=(defined $par{ca} && $par{ca});

  my $cmol=$mol->activeChains()->[0];

  for (my $ires=0; $ires<=$#{$cmol->{res}}; $ires++) {
    my $from=$cmol->{res}->[$ires]->{start};
    my $to=$cmol->{res}->[$ires]->{end};
    my ($x,$y,$z,$cax,$cay,$caz);
    my $nc=0;
    $x=$y=$z=0.0;
    for (my $ia=$from; $ia<=$to; $ia++) {
      my $aname=$cmol->{atom}->[$ia]->{atomname};
      if ($aname ne "N" && $aname ne "C" && $aname ne "O" && 
	  $aname!~/^[0-9]*H/) {
# && ($aname ne "CA" || !$wantca)) {
	$x+=$cmol->{atom}->[$ia]->{xcoor};
	$y+=$cmol->{atom}->[$ia]->{ycoor};
	$z+=$cmol->{atom}->[$ia]->{zcoor};
	$nc++;
      }
      if ($aname eq "CA" && $wantca) {
	$cax=$cmol->{atom}->[$ia]->{xcoor};
	$cay=$cmol->{atom}->[$ia]->{ycoor};
	$caz=$cmol->{atom}->[$ia]->{zcoor};
      }
    }

    my $rec={xcoor=>0, ycoor=>0, zcoor=>0};
    push (@{$self->{sidechain}},$rec);

    if ($nc>0) {
      $x/=$nc;
      $y/=$nc;
      $z/=$nc;
      $self->toProjection($x,$y,$z,$ires);
    } 

    if ($wantca) {
      my $carec={xcoor=>0, ycoor=>0, zcoor=>0};
      push (@{$self->{ca}},$carec);
      $self->toProjection($cax,$cay,$caz,$ires,1)
    }
  }
}

## method: getRandomMONSSTER(num[,parameters])
## generates a random MONSSTER chain
## Parameters for the new chain can be given as in the
## constructor.

sub genRandomMONSSTER {
  my $self=shift;
  my $num=shift;
  my %par=@_;
  
  $num=10 if (!defined $num);
  
  $self->{gridsize}    = $par{gridsize} if (defined $par{gridsize});
  $self->{offset}->{xcoor} = $par{offsetx}  if (defined $par{offsetx});
  $self->{offset}->{ycoor} = $par{offsety}  if (defined $par{offsety});
  $self->{offset}->{zcoor} = $par{offsetz}  if (defined $par{offsetz});

  $self->{integer}=1;
  $self->{resolution}=1.45;

  my $mchainbin=&GenUtil::findExecutable("mchain");
  die "cannot find mchain executable" 
    if (!defined $mchainbin);

  my $seed=time ^ $$ ^ unpack "%L*", `ps -ela | gzip`;
  my $option="-o $self->{offset}->{xcoor} $self->{offset}->{ycoor} $self->{offset}->{zcoor}";
  $option.=" -g $self->{gridsize} -n $num -seed $seed";

  open MCHAIN,"$mchainbin $option |";
  $self->readChain(\*MCHAIN);
  close MCHAIN;
}

## method: genMONSSTERFromAllAtom(mol[,parameters])
## generates a MONSSTER chain from an all-atom structure
## represented as a Molecule object.
## Parameters for the new chain can be given as in the
## constructor. In addition <mark>fraglist</mark> may be
## used to provide a residue list for loop/fragment modeling.
## The chain for those residues will then be generated
## randomly.

sub genMONSSTERFromAllAtom {
  my $self=shift;
  my $mol=shift;
  my %par=@_;


  $self->{gridsize} = $par{gridsize} if (defined $par{gridsize});
  $self->{offset}->{xcoor} = $par{offsetx}  if (defined $par{offsetx});
  $self->{offset}->{ycoor} = $par{offsety}  if (defined $par{offsety});
  $self->{offset}->{zcoor} = $par{offsetz}  if (defined $par{offsetz});
  
  my $fraglist;
  $fraglist = $par{fraglist} if (defined $par{fraglist});

  $self->{resolution}=1.45;
  $self->{integer}=1;

  my $option="-o $self->{offset}->{xcoor} $self->{offset}->{ycoor} $self->{offset}->{zcoor}";
  $option.=" -g $self->{gridsize}";
  $option.=" -l ".&GenUtil::fragOptionFromList($fraglist) 
    if (defined $fraglist);

  my $seed=time ^ $$ ^ unpack "%L*", `ps -ela | gzip`;
  $option.=" -seed $seed ";

  my $mchainbin=&GenUtil::findExecutable("mchain");
  die "cannot find mchain executable" 
    if (!defined $mchainbin);

  local (*READ,*WRITE);
  my $pid=open2(*READ,*WRITE,"$mchainbin $option -p -");
  $mol->writePDB(\*WRITE,ssbond=>0);  close WRITE;
  $self->readChain(\*READ); close READ;
  waitpid($pid,0);
}

## method: readMONSSTERTraj(file,timestep,cycle)
## reads a single chian configuration from a MONSSTER 
## trajectory. The configuration is identified by
## the time step and cycle indices.

sub readMONSSTERTraj {
  my $self=shift;
  my $fname=&GenUtil::getInputFile(shift);
  my $tstep=shift;
  my $cycle=shift;

  my $lastcycle=99999;
  my $cyc;
  my $t=0;

  $self->{sidechain}=();
  $self->{ca}=();
  $self->{resolution}=1.45;
  $self->{integer}=1;

 TRAJIO:
  while(<$fname>) {
    chomp;
    s/^ +//;
    my @f=split(/ +/);
    if ($#f+1==4 ) {
      last TRAJIO if ($#{$self->{sidechain}}>=0);
      $cyc=$f[2]+0;
      $t++ if ($cyc<$lastcycle);
      $lastcycle=$cyc;
    } else {
      if ($cyc==$cycle && $t==$tstep) {
	while (@f) {
	  my $rec={};
	  $rec->{xcoor}=shift @f;
	  $rec->{ycoor}=shift @f;
	  $rec->{zcoor}=shift @f;
	  push(@{$self->{sidechain}},$rec);
	}
      }
    }
  }
  
  undef $fname;
}

## method: $success = nextMONSSTERTraj(handle)
## reads the next chain configuration from a previously
## opened MONSSTER trajectory file. It returns 1 if 
## successful, 0 otherwise.

sub nextMONSSTERTraj {
  my $self=shift;
  my $fname=shift;

  $self->{sidechain}=();
  $self->{ca}=();
  $self->{resolution}=1.45;
  $self->{integer}=1;

 TRAJIO:
  while(<$fname>) {
    chomp;
    s/^ +//;
    my @f=split(/ +/);
    if ($#f==4 || $_=~/\./) {
      last TRAJIO if ($#{$self->{sidechain}}>=0);
    } else {
      while (@f) {
	my $rec={};
	$rec->{xcoor}=shift @f;
	$rec->{ycoor}=shift @f;
	$rec->{zcoor}=shift @f;
	push(@{$self->{sidechain}},$rec);
      }
    }
  }
  
  undef $fname;

  return ($#{$self->{sidechain}}>=0)?1:0;
}


## method: checkMONSSTER(show)
## checks whether the current chain is a valid
## MONSSTER chain. It returns a non-zero value
## if it is not a valid MONSSTER chain. Additional
## information is printed if the argument
## <mark>show</mark> is set to 1

sub checkMONSSTER {
  my $self=shift;
  my $show=shift;

  die "not a MONSSTER chain" 
    if (abs($self->{resolution}-1.45)>0.001);

  return -1 if ($#{$self->{sidechain}}<0);

  $show=0 if (!defined $show);

  my $violate=0;

  my $id;
  for (my $i=0; $i<$#{$self->{sidechain}}; $i++) {
    $id=&_idiff($self->{sidechain},$i,$i+1);
    if ($id>30) {
      printf "residues %3d <-> %3d are too far (%d, max: 30)\n",
      $i+1,$i+2,$id if ($show);
      $violate++;
    }
    if ($i>0) {
      $id=&_idiff($self->{sidechain},$i-1,$i+1);
      if ($id>68) {
	printf "residues %3d <-> %3d are too far (%d, max: 68)\n",
	$i,$i+2,$id if ($show);
	$violate++;
      }
    }
    if ($i>0) {
      if (&_colinear($self->{sidechain},$i-1,$i,$i+1)) {
	printf "residues %3d <-> %3d <-> %3d are colinear\n",
	$i,$i+1,$i+2 if ($show);
	$violate++;
      }
    }

    for (my $j=$i+1; $j<=$#{$self->{sidechain}}; $j++) {
      $id=&_idiff($self->{sidechain},$i,$j);
      if ($id<9) {
	printf "residues %3d <-> %3d are too close (%d, min 9)\n",
	$i+1,$j+1,$id if ($show);
	$violate++;
      }
    }
  }

  return $violate;
}

## method: fromMolecule(mol)
## generates a lattice chain from a chain representation
## as a Molecule object.

sub fromMolecule {
  my $self=shift;
  my $mol=shift;

  $self->{sidechain}=();
  $self->{ca}=();
  $self->{resolution}=0;
  $self->{integer}=0;
  $self->{offset}->{xcoor}=0;
  $self->{offset}->{ycoor}=0;
  $self->{offset}->{zcoor}=0;

  my $c=$mol->activeChains()->[0];
  my $atom=$c->{atom};
  my $res=$c->{res};

  for (my $ir=0; $ir<=$#{$res}; $ir++) {
    my $from=$res->[$ir]->{start};
    my $to=($ir<$#{$res})?$res->[$ir+1]->{start}-1:$#{$atom};
    my ($carec, $screc);
    for (my $i=$from; $i<=$to; $i++) {
      my $a=$atom->[$i];
      if ($a->{atomname} eq "CA") {
	$carec={ xcoor => $a->{xcoor}, ycoor => $a->{ycoor}, zcoor => $a->{zcoor} };
      } elsif ($a->{atomname} eq "SC") {
	$screc={ xcoor => $a->{xcoor}, ycoor => $a->{ycoor}, zcoor => $a->{zcoor} };
      }
    }
    die "need to have CA or side chain for each residue" 
      if (!defined $screc && !defined !$carec);

    $screc={ xcoor => 0.0, ycoor => 0.0, zcoor => 0.0 } 
      unless (defined $screc);
  
    push (@{$self->{sidechain}},$screc);
    push (@{$self->{ca}},$carec) if (defined $carec);
  }
  
  die "number of CA and side chain centers do not match"
    if ($#{$self->{sidechain}} != $#{$self->{ca}} && $#{$self->{ca}}>=0);
}

sub _idiff {
  my ($sidechain,$i,$j)=@_;
  my ($dx,$dy,$dz);
  
  $dx=$sidechain->[$i]->{xcoor}-$sidechain->[$j]->{xcoor};
  $dy=$sidechain->[$i]->{ycoor}-$sidechain->[$j]->{ycoor};
  $dz=$sidechain->[$i]->{zcoor}-$sidechain->[$j]->{zcoor};

  return ($dx*$dx+$dy*$dy+$dz*$dz);
}
  
sub _colinear {
  my ($sidechain,$i,$j,$k)=@_;
  my ($ix,$iy,$iz,$jx,$jy,$jz,$kx,$ky,$kz);

  $ix=$sidechain->[$i]->{xcoor}-$sidechain->[$j]->{xcoor};
  $iy=$sidechain->[$i]->{ycoor}-$sidechain->[$j]->{ycoor};
  $iz=$sidechain->[$i]->{zcoor}-$sidechain->[$j]->{zcoor};

  $jx=$sidechain->[$j]->{xcoor}-$sidechain->[$k]->{xcoor};
  $jy=$sidechain->[$j]->{ycoor}-$sidechain->[$k]->{ycoor};
  $jz=$sidechain->[$j]->{zcoor}-$sidechain->[$k]->{zcoor};

  $kx=$iy*$jz-$iz*$jy;
  $ky=$iz*$jx-$ix*$jz;
  $kz=$ix*$jy-$iy*$jx;
  
  return ($kx*$kx+$ky*$ky+$kz*$kz)==0;
}

sub toProjection {
  my $self=shift;
  my $x=shift;
  my $y=shift;
  my $z=shift;
  my $inx=shift;
  my $ca=shift;

  my $s=(defined $ca && $ca)?$self->{ca}->[$inx]:$self->{sidechain}->[$inx];

  if ($self->{resolution}>0) {
    $s->{xcoor}=$x/$self->{resolution};
    $s->{ycoor}=$y/$self->{resolution};
    $s->{zcoor}=$z/$self->{resolution};
  } else {
    $s->{xcoor}=$x;
    $s->{ycoor}=$y;
    $s->{zcoor}=$z;
  }
  
  $s->{xcoor}+=$self->{offset}->{xcoor};
  $s->{ycoor}+=$self->{offset}->{ycoor};
  $s->{zcoor}+=$self->{offset}->{zcoor};
}

sub fromProjection {
  my $self=shift;
  my $inx=shift;
  my $ca=shift;

  my $s=(defined $ca && $ca)?$self->{ca}->[$inx]:$self->{sidechain}->[$inx];

  my $newx=$s->{xcoor}-$self->{offset}->{xcoor};
  my $newy=$s->{ycoor}-$self->{offset}->{ycoor};
  my $newz=$s->{zcoor}-$self->{offset}->{zcoor};

  if ($self->{resolution}>0) {
    $newx*=$self->{resolution};
    $newy*=$self->{resolution};
    $newz*=$self->{resolution};
  } 

  return ($newx,$newy,$newz);
}


1;
