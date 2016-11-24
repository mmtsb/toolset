# Sequence package
# read/write/convert sequence information
#
# http://mmtsb.scripps.edu/doc/Sequence.pm.html
# 2000, Michael Feig, Brooks group, TSRI

package Sequence;

require 5.004;

use strict;

use FileHandle;
use IPC::Open2;

use GenUtil;
use Molecule;

## data: sequence[] -> { residue secondary index valid }
## sequence data

## data: resinx[]
## lookup table for residue index

use vars qw ( %_sectrans %_seqabbrev %_seqlong %_secnum %_secletter );

BEGIN {
  %_sectrans  = ( ALPHA   => 'H', A => 'H',
  	          BETA    => 'E', B => 'E', SHEET => 'E', EXTENDED => 'E',
	          COIL    => 'C', 
                  OTHER   => 'U', UNKNOWN => 'U' );

  %_secnum    = ( H => 2, E => 4, C => 1, U => 1 ); 

  %_secletter = ( 1 => 'U', 2 => 'H', 4 => 'E' );

  %_seqabbrev = ('GLY', 'G', 'PRO', 'P', 'ALA', 'A', 'VAL', 'V',
	         'LEU', 'L', 'ILE', 'I', 'MET', 'M', 'CYS', 'C',
	         'PHE', 'F', 'TYR', 'Y', 'TRP', 'W', 'HIS', 'H',
	         'HSD', 'H', 'HSE', 'H', 'LYS', 'K', 'ARG', 'R',
	         'GLN', 'Q', 'ASN', 'N', 'GLU', 'E', 'ASP', 'D',
	         'SER', 'S', 'THR', 'T', 'UNK', 'X', 'HSP', 'H');

  %_seqlong   = ('G', 'GLY', 'P', 'PRO', 'A', 'ALA', 'V', 'VAL', 
	         'L', 'LEU', 'I', 'ILE', 'M', 'MET', 'C', 'CYS',
	         'F', 'PHE', 'Y', 'TYR', 'W', 'TRP', 'H', 'HIS',
	         'K', 'LYS', 'R', 'ARG', 'Q', 'GLN', 'N', 'ASN',
	         'E', 'GLU', 'D', 'ASP', 'S', 'SER', 'T', 'THR',
                 'X', 'UNK' );
}

## constructor: new([seqstring|molecule][,fraglist])
## creates a new Sequence object. An abbreviated
## sequence string or a Molecule object can be
## given as arguments to extract the sequence.
## A fragment list can be given as a second argument
## to complete fragments missing in a molecule structure
## from abbreviated sequence strings.

sub new {
  my $arg=shift;
  my $slist=shift;
  my $fill=shift;

  my $self={};

  if (!defined $arg) {
    $self->{sequence}=();
  } elsif (!ref($arg)) {
    $arg=~s/[ \n\t]+//;
    for (my $i=0; $i<length($arg); $i++) {
      my $seqrec={};
      my $char=uc substr($arg,$i,1);
      die "unknown residue name abbreviation $char found" 
	if (!exists $_seqlong{$char});
      $seqrec->{index}=$i+1;
      $seqrec->{residue}=$_seqlong{$char};
      $seqrec->{secondary}='U';
      $seqrec->{valid}=1;
      push(@{$self->{sequence}},$seqrec);
    }
  } else {
    my $c=$arg->activeChains()->[0];

    my $ts={};
    foreach my $r ( @{$c->{res}} ) {
      $ts->{$r->{num}}=$r->{name};
    }

    foreach my $s ( @{$slist} ) {
      for (my $i=0; $i<length($s->{seq}); $i++) {
	my $char=uc substr($s->{seq},$i,1);
	die "unknown residue name abbreviation $char found" 
	  if (!exists $_seqlong{$char});
	my $index=$i+$s->{inx};
	$ts->{$index}=$_seqlong{$char};
      }
    }

    my $lastinx;
    foreach my $t ( sort {$a<=>$b} keys %{$ts} ) {
      my $seqrec={};
      $seqrec->{index}=$t;
      $seqrec->{residue}=$ts->{$t};
      $seqrec->{secondary}='U';
      $seqrec->{valid}=1;
     
      if (defined $lastinx && $seqrec->{index}!=$lastinx+1) {
        if (defined $fill && $fill) {
          for (my $ifill=$lastinx+1; $ifill<$seqrec->{index}; $ifill++) {
            my $tsrec={};
            $tsrec->{index}=$ifill;
            $tsrec->{residue}="UNK";
            $tsrec->{secondary}="U";
            $tsrec->{valid}=1;
            push(@{$self->{sequence}},$tsrec);
          }
        } else { 
          printf STDOUT "warning: non-continuous sequence at $lastinx\n"
        }
      }

      push(@{$self->{sequence}},$seqrec);
      $lastinx=$seqrec->{index};
    }
  }

  bless($self);
  return $self;
}

## method: readMONSSTER(file)
## read a sequence from a MONSSTER sequence file

sub readMONSSTER {
  my $self=shift;
  my $seqfile=&GenUtil::getInputFile(shift);

  $self->{sequence}=();

  while(<$seqfile>) {
    s/^[ \t]+//;
    my @f=split(/ +/);
    my $seqrec={};

    die "unknown residue name $f[1]" 
      if (!exists $_seqabbrev{$f[1]});
    die "unknown secondary index $f[2]"
      if (!exists $_secletter{$f[2]});

    $seqrec->{index}=$f[0];
    $seqrec->{residue}=$f[1];
    $seqrec->{secondary}=$_secletter{$f[2]};
    $seqrec->{valid}=1;
    
    push(@{$self->{sequence}},$seqrec);
  }

  undef $seqfile;
}

## method: writeMONSSTER(file)
## write the current sequence to a MONSSTER sequence file

sub writeMONSSTER {
  my $self=shift;
  my $seqfile=&GenUtil::getOutputFile(shift);
  my $usevalid=shift;

  for (my $i=0; $i<=$#{$self->{sequence}}; $i++) {
    if ($self->{sequence}->[$i]->{valid} || !defined $usevalid || !$usevalid) {
      my $resname=uc $self->{sequence}->[$i]->{residue};
      $resname=~s/HSD|HSE|HSP/HIS/;
      my $sec=$_secnum{$self->{sequence}->[$i]->{secondary}};
      $sec=1 if ($sec==0);
      printf $seqfile "%5d   %3s%5d%5d\n",
      $self->{sequence}->[$i]->{index}, $resname, $sec, 1;
    }
  }

  undef $seqfile;
}

## method: $string = abbrevSeq()
## generates an abbreviated sequence string

sub abbrevSeq {
  my $self=shift;
  my $usevalid=shift;
  my $tstr="";

  for (my $i=0; $i<=$#{$self->{sequence}}; $i++) {
    if ($self->{sequence}->[$i]->{valid} || !defined $usevalid || !$usevalid) {
      my $resname=uc $self->{sequence}->[$i]->{residue};
      die "unknown residue name $resname" 
        if (!exists $_seqabbrev{$resname});
      $tstr.=$_seqabbrev{$resname};
    }
  }
  return $tstr;
}

## method: $string = abbrevSec()
## generates an abbreviated string for the secondary
## structure information

sub abbrevSec {
  my $self=shift;
  my $usevalid=shift;
  my $tstr="";

  for (my $i=0; $i<=$#{$self->{sequence}}; $i++) {
    if ($self->{sequence}->[$i]->{valid} || !defined $usevalid || !$usevalid) {
    $tstr.=$self->{sequence}->[$i]->{secondary};
    }
  }
  return $tstr;
}

## method: modifyResidueName(newreslist) 
## sets residue names to new names in the
## argument

sub modifyResidueName {
  my $self=shift;
  my $newreslist=shift;

  return if (!defined $newreslist);
  
  foreach my $n ( split(/,/,$newreslist) )  {
    my ($num,$newname)=split(/:/,$n);
    if (defined $num && defined $newname) {
      my $r=$self->getResidue($num);
      if (defined $r->{residue}) {
	my $name=(length($newname)==1)?$_seqlong{uc $newname}:uc $newname;
	$r->{residue}=$name;
      }
    }
  }
}

## method: setSecondary(index,type)
## sets the secondary structure information
## for a specific residue

sub setSecondary {
  my $self=shift;
  my $num=shift;
  my $type=uc shift;

  my $sec = (exists $_secnum{$type}) ? $type : $_sectrans{$type};

  die "invalid secondary structure specification" 
    if (!defined $sec || $sec eq "");
  
  die "index out of range" 
    if (!defined $self->{sequence}->[$num]);

  $self->{sequence}->[$num]->{secondary}=$sec;
}

## method: secFromOneFile(file)
## extracts secondary structure information from
## a file containing single letter abbreviations

sub secFromOneFile {
  my $self=shift;
  my $fname=shift;

  my $secstr="";
  my $inp=&GenUtil::getInputFile($fname);
  while(<$inp>) { chomp $_; $secstr.=$_; }
  $secstr=~s/ +//g;
  close $inp;

  for (my $i=0; $i<=$#{$self->{sequence}}; $i++) {  
    my $t="U";
    if (substr($secstr,$i,1) eq "E") {
      $t="E";
    } elsif(substr($secstr,$i,1) eq "H") {
      $t="H";
    }
    $self->{sequence}->[$i]->{secondary}=$t;
  }
}

## method: secFromPredict(filelist)
## extracts secondary structure information from
## a list of file names containing output from
## common prediction servers

sub secFromPredict {
  my $self=shift;
  my @fileArr=@_;
  
  my @allsec=();
  foreach my $filename (@fileArr) {
    my $sec=&_readPredictOutput($filename,$self->{sequence});
    push (@allsec,$sec);
  }
  
  for (my $i=0; $i<=$#{$self->{sequence}}; $i++) {
    my $cntE=0;
    my $cntH=0;
    my $cntU=0;

    for (my $j=0; $j<=$#allsec; $j++) {
      my $c=substr($allsec[$j],$i,1);
      if ($c eq "E") {
	$cntE++;
      } elsif ($c eq "H") {
	$cntH++;
      } else {
	$cntU++;
      }
    }

    if ($cntE>$cntH && $cntE>=$cntU) {
      $self->{sequence}->[$i]->{secondary}="E";
    } elsif ($cntH>$cntE && $cntH>=$cntU) {
      $self->{sequence}->[$i]->{secondary}="H";
    } else {
      $self->{sequence}->[$i]->{secondary}="U";
    }
  }
}

## method: secFromDSSP(molecule)
## sets secondary structure information from the
## output of the external DSSP program for the
## protein structure given as Molecule object argument 

sub secFromDSSP {
  my $self=shift;
  my $mol=shift;

  &GenUtil::log("Sequence::secFromDSSP");

#  my $dsspbin=&GenUtil::findExecutable("dssp");
  my $dsspbin=&GenUtil::findExecutable("dsspcmbi");
  die "Cannot find dssp executable"
    if (!defined $dsspbin);
  
  local (*READ,*WRITE);
  my $pid=open2(*READ,*WRITE,"$dsspbin -- 2>/dev/null");

  $mol->writePDB(\*WRITE,translate=>"GENERIC",ssbond=>0);  close WRITE;

  my $readdata=0;
  while(<READ>) {
    if ($readdata) {
      my $num=substr($_,5,5);
      $num=~s/ +//g;
      if ($num ne "") {
      my $ab=uc substr($_,13,1);
      my $sec=uc substr($_,16,1);

      my $eab=$_seqabbrev{$self->getResidue($num)->{residue}};
      if ( $eab ne $ab && $ab ne "X") {
	printf STDERR "Residues do not match (dssp: $ab, expected: $eab)\n";
	close READ;
	waitpid($pid,0);
	return 0;
      }
	
	$self->getResidue($num)->{secondary}=($sec eq "H")?"H":(($sec eq "E")?"E":"U");
      }
    } elsif (/^ +\#  RESIDUE AA.*/) {
      $readdata=1;
    }
  }
  close READ;
  waitpid($pid,0);
  return 1;
}

## method: number = firstResNum()
## number of first residue

sub firstResNum {
  my $self=shift;
  return $self->{sequence}->[0]->{index};
}

## method: number = lastResNum()
## number of last residue

sub lastResNum {
  my $self=shift;
  return $self->{sequence}->[$#{$self->{sequence}}]->{index};
}

## method: setValidResidues(fraglist)
## sets the <mark>valid</mark> flag to 1
## for residues in the fragment list and to
## 0 for residues outside the list

sub setValidResidues {
  my $self=shift;
  my $fraglist=shift;
  my $exclmode=shift;

  $exclmode=0 if (!defined $exclmode);

  $self->resetValidResidues($exclmode?1:0);

  foreach my $l ( @{$fraglist} ) {
    if (!defined $l->{from}) {
      $l->{from}=$self->firstResNum();
      $l->{to}=$self->lastResNum();
    } elsif (!defined $l->{to}) {
      $l->{to}=$l->{from};
    }

    $l->{from}=$self->firstResNum()
      if ($l->{from}<$self->firstResNum());
    $l->{to}=$self->lastResNum()
      if ($l->{to}>$self->lastResNum());

    for (my $i=$l->{from}; $i<=$l->{to}; $i++) {
      my $r=$self->getResidue($i);
      $r->{valid}=($exclmode?0:1) if (defined $r);
    }    
  }
}

## method: resetValidResidues([value])
## sets the <mark>valid</mark> flag of all
## residues to the given value (default: 1)

sub resetValidResidues {
  my $self=shift;
  my $value=shift;

  $value=1 if (!defined $value);
  
  foreach my $r ( @{$self->{sequence}} ) {
    $r->{valid}=$value;
  }
}

## method: $list = listFromValid(force)
## returns a list of residues from the residues 
## previously set with <mark>setValidResidues</mark>.

sub listFromValid {
  my $self=shift;

  my $retlist=();
  my @arr;
  foreach my $r ( @{$self->{sequence}} ) {
    push(@arr,$r->{index}) if ($r->{valid});
  }
  if ($#arr>=0) {
    foreach my $tlist ( @{&GenUtil::fragListFromArray(\@arr)} ) {
      push(@{$retlist},$tlist);
    }
  } else {
    my $rec={};
    $rec->{from}=$self->firstResNum();
    $rec->{to}=$self->lastResNum();
    push(@{$retlist},$rec);
  }

  return $retlist;
}

## method: $index = getResidue(resnum)
## returns the residue index for a residue number. 
## A reference to a chain structure may be given as 
## additional argument for multi-domain structures

sub getResidue {
  my $self=shift;
  my $inx=shift;

  if (!defined $self->{resinx}) {
    $self->{resinx}={};
    for (my $ir=0; $ir<=$#{$self->{sequence}}; $ir++) {
      my $key=$self->{sequence}->[$ir]->{index};
      $self->{resinx}->{$key}=$ir;
    }
  }

  return $self->{sequence}->[$self->{resinx}->{$inx}];
}  


sub _readPredictOutput {
  my $filename=shift;
  my $seq=shift;
  my $secstr="";
  
  die "cannot read prediction output filename $filename" 
    if (!-r $filename);

  my $current=0;
  my $sspro=0;
  my $phd=0;
  my $phdprof=0;

  my $readsa=0;
  my @sasave;

  my $pirmode=0;
  my $pirinx=0;

  my $lastline;

  open PREDINP,"$filename";
 PREDIO:
  while (<PREDINP>) {
    chomp;

    my $sa="";
    my $sec="";

    if (/SSpro/) {
      $sspro=1;
    } elsif ($sspro==1 && /^Prediction:/) {
      $sspro=2;
    } elsif ($sspro==2 && /^[A-Z]+/) {                         # sspro
      ($sa=$_)=~s/ +//;
      $sec=<PREDINP>;
      chomp $sec;
      $sec=~s/ +//;
    } elsif (/PROF results .normal./) {
      $phdprof=1;
    } elsif ($phdprof==1 && /^[ \t]+AA[ \t]+\|/) {
      ($sa=$_)=~s/^[ \t]+AA[ \t]+\|//;
      $sa=~s/\|.*$//;
      $sec=<PREDINP>;
      $sec=<PREDINP>
	if ($sec=~/OBS/);
      chomp $sec;
      $sec=~s/^[ \t]+PROF_[ \ta-z]+\|//;
      $sec=~s/\|.*$//;
    } elsif (/PHD output/) {
      $phd=1;
    } elsif ($phd==1 && /^[ \t]+protein:[ \t]+predict[ \t]+length/) {
      $phd=2;
    } elsif ($phd==2 && /^[ \t]+AA[ \t]+\|/) {                 # phd
      ($sa=$_)=~s/^[ \t]+AA[ \t]+\|//;
      $sa=~s/\|.*$//;
      $sec=<PREDINP>;
      chomp $sec;
      $sec=~s/^[ \t]+PHD[ \ta-z]+\|//;
      $sec=~s/\|.*$//;
    } elsif (/^Pred: [CEH]+/) {                                # psipred
      ($sec=$_)=~s/^Pred: //;
      $sec=~s/ +//;
      $sa=<PREDINP>;
      chomp $sa;
      $sa=~s/^  AA: //;
      $sa=~s/ +//;
    } elsif (/^[ \t]*[A-Z]+[ \t]+[CEHU][ \t]+[01]\.[0-9]/) {   # jpred2, pssp, prof2, pred2ary
      s/^[ \t]+//;
      my @f=split(/[ \t]+/);
      
      $sa=(exists $_seqabbrev{$f[0]}) ? $_seqabbrev{$f[0]} : $f[0];
      $sec=$f[1];
    } elsif (/^ *[0-9]+ [A-Z] . [HE-] [0-9]/) {                # pred2ary
      s/^[ \t]+//;
      my @f=split(/[ \t]+/);
      
      $sa=(exists $_seqabbrev{$f[1]}) ? $_seqabbrev{$f[1]} : $f[1];
      $sec=$f[3];
    } elsif (/^>P.; .+/) {
      <PREDINP>;
      $pirmode++;
      $pirinx=0;
    } elsif (/[A-Z]/ && $pirmode>0) {
      chomp;
      s/ +//g;
      s/[ \*]+//g;
      if ($pirmode==1) {
         $sasave[$pirinx++]=$_;
      } else {
         $sa=$sasave[$pirinx++];
         $sec=$_;
      }
    } else {
      $lastline=$_;
    }

    for (my $i=0; $i<length($sa); $i++) {
      $readsa=1;
      my $sai=substr($sa,$i,1);
      my $seci=substr($sec,$i,1);
      if ($sai eq $_seqabbrev{$seq->[$current]->{residue}}) {
	$secstr.=$seci;
	last PREDIO if (++$current>$#{$seq}); 
      }
    }
  }

  if (!$readsa) {
    $lastline=~s/^ +//g;
    $secstr=substr($lastline,0,$#{$seq}+1);
    $current=length($secstr);
  } 

  die "did not read secondary prediction for complete sequence in $filename"
    if ($current<=$#{$seq});

  close PREDINP;

  return $secstr;
}


## method: readMoleculeChain(file)
## read a sequence from a molecule chain 

sub readMoleculeChain {
  my $self=shift;
  my $c=shift;

  $self->{sequence}=();

  my $lastinx;
  foreach my $r ( @{$c->{res}} ) {
    foreach  (my $i=$lastinx+1; $i<$r->{num}; $i++) {
      push(@{$self->{sequence}},_newRec($i,"UNK",'U',1));
    }
    push(@{$self->{sequence}},_newRec($r->{num},$r->{name},'U',1));
    $lastinx=$r->{num};
  }
}

sub _newRec {
  my $inx=shift;
  my $res=shift;
  my $sec=shift;
  my $val=shift;

  my $seqrec={};
  $seqrec->{index}=$inx;
  $seqrec->{residue}=(exists $_seqabbrev{$res})?$res:"UNK";
  $seqrec->{secondary}=$sec;
  $seqrec->{valid}=$val;
  return $seqrec;
}

1;
