# Molecule package
# read/write/convert structure info
#
# http://mmtsb.scripps.edu/doc/Molecule.pm.html
# 2000, Michael Feig, Brooks group, TSRI
#

package Molecule;

require 5.004;

use strict;

use FileHandle;
use IPC::Open2;
use IPC::Open3;

use GenUtil;
use Sequence;
use SICHO;
use Analyze;

## data: chain[] -> { id atom[] res[] resinx[]
## data:              xcoor[] ycoor[] zcoor[] }
## data:   atom[] -> { atominx atomname resname resnum
## data:               chain xcoor ycoor zcoor hyd seg aux1 aux2 hetero }
## data:   res[]  -> { name num start end seg valid chain hetero }
## molecule information with substructures containing atom and residue
## information as well as a residue lookup table and coordinate cache
## arrays.

## data: chainlookup -> { chainid ... }
## lookup hash table for multiple chains

## data: defchain 
## currently selected default chain (chain ID)

## data: selchain
## currently selected chain (object)

## data: segmentlist[] -> { name first last }
## list of segment IDs

## data: ssbond[] -> { chain1 resnum1 chain2 resnum2 }
## disulfide bonds

## data: havesegments

## constructor: new([PDBfile])
## creates a new Molecule object and reads a PDB structures
## if a file name is given

sub new {
  my $farg = shift;

  my $self = {};
  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;
  $self->{segmentlist}=undef;
  $self->{ssbond}=();
  $self->{havesegments}=0;
  $self->{have}={};

  bless $self;

  $self->readPDB($farg) 
    if (defined $farg);

  return $self;
}

## constructor: newInherit([PDBfile])
## creates a new Molecule object and reads a PDB structures
## if a file name is given
## (this form can act as a base class for derived classes)

sub newInherit {
  my $proto=shift;
  my $class=ref($proto) || $proto;
  my $farg = shift;

  my $self = {};
  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;
  $self->{segmentlist}=undef;
  $self->{ssbond}=();

  bless($self,$class);

  $self->readPDB($farg) 
    if (defined $farg);

  return $self;
}

## method: readPDB(file[,translate[,ignoreseg]])
## reads a protein structures from a PDB file.
## <mark>translate</mark> may be set to <mark>CHARMM19</mark>
## for proper recognition of histidine residues.
## If <mark>ignoreseg</mark> is set segment IDs from the
## PDB are not read.

sub readPDB {
  my $self=shift;
  my $fname=&GenUtil::getInputFile(shift);
  my %par=@_;

  my $translate=(defined $par{translate})?$par{translate}:"";
  my $ignoreseg=(defined $par{ignoreseg} && $par{ignoreseg})?1:0;
  my $chainfromseg=(defined $par{chainfromseg} && $par{chainfromseg})?1:0;
  my $activemodel;

  my $firstmodelonly=(defined $par{firstmodel} && $par{firstmodel})?1:0;
  my $splitseg=(defined $par{splitseg} && $par{splitseg})?1:0;
  my $readalternate=(defined $par{alternate} && $par{alternate})?1:0;
  my $readmodels=0;

  my $lastchain=".";
  my $chainrec;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;
  $self->{ssbond}=();

  my $lastseg="XXX";
  my $lastnum=-999;
  my $lastinum=-999;
  my $ignore;

  my $addrnum=0;

  my $newchain=0;
 
  my $lastainx=0;

  my $hetero;
 
  my %sshave; 
 READPDB:
  while(<$fname>) {
    if (/^CRYST1/) {
      chomp;
      my @f=split(/\s+/);
      my $trec={};
      $trec->{a}=$f[1];
      $trec->{b}=$f[2];
      $trec->{c}=$f[3];
      $trec->{aa}=$f[4];
      $trec->{ab}=$f[5];
      $trec->{ac}=$f[6];
      $self->{cryst}=$trec;
    } elsif (/^SSBOND/) {
      my $trec={};
      ($trec->{chain1}=substr($_,15,1))=~s/ //g;
      ($trec->{resnum1}=substr($_,16,5))=~s/[A-Z\s]+//g;
      ($trec->{chain2}=substr($_,29,1))=~s/ //g;
      ($trec->{resnum2}=substr($_,30,5))=~s/[A-Z\s]+//g;

      my $key1=sprintf("%s:%d",$trec->{chain1},$trec->{resnum1});
      my $key2=sprintf("%s:%d",$trec->{chain2},$trec->{resnum2});

      push(@{$self->{ssbond}},$trec) 
        unless ($sshave{$key1} || $sshave{$key2} ||
                ($trec->{chain1} eq $trec->{chain2}  
              && $trec->{resnum1} == $trec->{resnum2}));
      $sshave{$key1}=1;
      $sshave{$key2}=1;
    } elsif (/^MODEL +([0-9]+)/) {
      if ((!defined $par{model} && (!$firstmodelonly || $readmodels==0)) || $1 == $par{model}) {
	$activemodel=1;
	$readmodels++;
      } else {
	$activemodel=0;
      }
    }
    last READPDB if (/^(END|End|\#End)/ && (!defined $activemodel || $activemodel));

    if(/HETAT/) {
      $hetero=1;
    } else { 
      $hetero=0;
    }

    #s/HETATM/ATOM  / if (/MSE/ || /CSD/ || /ABA/ || /CGU/ || /CME/ || /MLY/ || /PCA/ || /PTR/ || /SEP/ || /TPO/);
    #s/HETAT/ATOM / if (/MSE/ || /CSD/ || /ABA/ || /CGU/ || /CME/ || /MLY/ || /PCA/ || /PTR/ || /SEP/ || /TPO/);
    #s/HETATM/ATOM  / if (substr($_,21,1) eq $lastchain && $lastchain=~/[A-Z0-9]/ );

    if (/MSE/ || /CSD/ || /ABA/ || /CGU/ || /CME/ || /MLY/ || /PCA/ || /PTR/ || /SEP/ || /TPO/) {
      s/HETATM/ATOM  /;
      s/HETAT/ATOM /;
      $hetero=0;
    }
    if (substr($_,21,1) eq $lastchain && $lastchain=~/[A-Za-z0-9\-\_\=]/ ) {
      s/HETATM/ATOM  /;
    }
    
    s/^MSE$/MET/;
    s/^CSD$/ALA/;
    s/^ABA$/ALA/;
    s/^CGU$/GLU/;
    s/^CME$/CYS/;
    s/^MLY$/LYS/;
    s/^PCA$/ALA/;
    s/^PTR$/TYR/;
    s/^SEP$/SER/;
    s/^TPO$/THR/;
    if (/^ATOM/ && (!defined $activemodel || $activemodel)) {
      my ($atomname, $resname, $resnum, $iresnum, $alt, $chain, $seg);
      
      if (($alt=substr($_,16,1))=~/[ A0-9]/) {
	($atomname=substr($_,12,4))=~s/ //g;
	$atomname.=$alt if ($alt=~/[0-9]/);
	($resname=substr($_,17,4))=~s/ //g;
	$resnum=substr($_,22,6);
        $resnum=~s/\s//g;
	($iresnum=$resnum)=~s/[A-Z]+//g;
	$iresnum+=0;
	$chain=substr($_,21,1);
	($seg=substr($_,72,4))=~s/[ \n]//g;

	if ($chainfromseg && !$ignoreseg) {
	  if ($seg=~/HETA/) {
	    $chain="+";
	  } elsif ($seg=~/...([A-Za-z0-9\+\-\=\_])/) {
   	    $chain=$1;
          }
	}

	if ($resnum ne $lastnum) {
	  $ignore=($iresnum == $lastinum && !$readalternate);
	}

	if (!$ignore) {
# J
          $resname=~s/^WAT$/TIP3/;
          $atomname="OH2" if ($resname eq "TIP3" && $atomname eq "O");
          $atomname="OH2" if ($resname eq "TIP4" && $atomname eq "O");

          $resname=~s/Na\+/SOD/;
          $atomname=~s/Na\+/SOD/;

          $resname="GUA"
            if ($resname eq "G" || $resname =~/^[DR]G[35]*$/);
          $resname="ADE"
            if ($resname eq "A" || $resname =~/^[DR]A[35]*$/);
          $resname="CYT"
            if ($resname eq "C" || $resname =~/^[DR]C[35]*$/);
          $resname="THY"
            if ($resname eq "T" || $resname =~/^[DR]T[35]*$/);
          $resname="URA"
            if ($resname eq "U" || $resname =~/^[DR]U[35]*$/);
          $atomname="O5'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "O5*");
          $atomname="C5'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "C5*");
          $atomname="C4'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "C4*");
          $atomname="O4'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "O4*");
          $atomname="C3'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "C3*");
          $atomname="O3'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "O3*");
          $atomname="C2'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "C2*");
          $atomname="O2'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "O2*");
          $atomname="C1'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "C1*");
          $atomname="H5'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "1H5*");
          $atomname="H5''"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "2H5*");
          $atomname="H4'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "H4*");
          $atomname="H3'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "H3*");
          $atomname="H1'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "H1*");
          $atomname="H41"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "1H4*");
          $atomname="H42"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "2H4*");
          $atomname="H2''"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "1H2*");
          $atomname="H2'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "2H2*");
          $atomname="H2'"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "2HO*");
          $atomname="H21"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "1H2*");
          $atomname="H22"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "2H2*");
          $atomname="H51"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "1H5M");
          $atomname="H52"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "2H5M");
          $atomname="H53"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "3H5M");
          $atomname="H61"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "1H6");
          $atomname="H62"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "2H6");
          $atomname="O1P"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "OP1");
          $atomname="O2P"
            if ($resname =~/^CYT$|^GUA$|^ADE$|^URA$|^THY$|^GTP$|^ATP$|^UTP$|^CTP$/ && $atomname eq "OP2");

#PNA names
          $atomname="N"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "N1*");
          $atomname="H1'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && ($atomname eq "H1*2" || $atomname eq "2H1*"));
          $atomname="N2'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "N4*");
          $atomname="C5'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "C3*");
          $atomname="H5'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && ($atomname eq "H3*1" || $atomname eq "1H3*"));
          $atomname="H5''"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && ($atomname eq "H3*2" || $atomname eq "2H3*"));
          $atomname="C6'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "C2*");
          $atomname="H6'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && ($atomname eq "H2*1" || $atomname eq "1H2*"));
          $atomname="H6''"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && ($atomname eq "H2*2" || $atomname eq "2H2*"));
          $atomname="C3'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "C7*");
          $atomname="O3'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "O7*");
          $atomname="C4'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "C8*");
          $atomname="H4'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && ($atomname eq "H8*1" || $atomname eq "1H8*"));
          $atomname="H4''"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && ($atomname eq "H8*2" || $atomname eq "2H8*"));
          $atomname="C2'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "C5*");
          $atomname="H2'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && ($atomname eq "H5*1" || $atomname eq "1H5*"));
          $atomname="H2''"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && ($atomname eq "H5*2" || $atomname eq "2H5*"));
          $atomname="C"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "C*");
          $atomname="O1'"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "O1*");
          $atomname="H51"
            if ($resname =~/^TPN$/ && $atomname eq "1H7");
          $atomname="H52"
            if ($resname =~/^TPN$/ && $atomname eq "2H7");
          $atomname="H53"
            if ($resname =~/^TPN$/ && $atomname eq "3H7");
          $atomname="C5M"
            if ($resname =~/^TPN$/ && $atomname eq "C7");
          $atomname="H21"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "1H2");
          $atomname="H22"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "2H2");
          $atomname="H61"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "1H6");
          $atomname="H62"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "2H6");
          $atomname="H41"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "1H4");
          $atomname="H42"
            if ($resname =~/^TPN$|^GPN$|^CPN$|^APN$|^UPN$/ && $atomname eq "2H4");

          $atomname="C5M"
            if ($resname =~/THY/ && $atomname eq "C7");
#J end
	  $atomname="CD1" 
	    if ($resname eq "ILE" && $atomname eq "CD");
	  $atomname="SD"
	    if ($resname eq "MET" && $atomname eq "SE");
	  $atomname="O"
	    if (($atomname eq "OT1" || $atomname eq "O1" || $atomname eq "OCT1")
		&& $resname =~/^ALA$|^ARG$|^ASN$|^ASP$|^CYS$|^GLN$|^GLU$|^GLY$|^HSD$|^HSE$|^HSP$|^HIS$|^HID$|^HIE$|^HIP$|^HSP$|^ILE$|^LEU$|^LYS$|^MET$|^PHE$|^PRO$|^SER$|^THR$|^TRP$|^TYR$|^VAL$|^CYX$/);
	  $atomname="OXT"
	    if (($atomname eq "OT2" || $atomname eq "O2" || $atomname eq "OCT2")
		&& $resname =~/^ALA$|^ARG$|^ASN$|^ASP$|^CYS$|^GLN$|^GLU$|^GLY$|^HSD$|^HSE$|^HSP$|^HIS$|^HID$|^HIE$|^HIP$|^HSP$|^ILE$|^LEU$|^LYS$|^MET$|^PHE$|^PRO$|^SER$|^THR$|^TRP$|^TYR$|^VAL$|^CYX$/);
	  
	  if ($translate eq "CHARMM19") {
	    $resname=~s/^HSD$/HSE/;
	    $resname=~s/^HIS$/HSD/;
          } elsif ($translate eq "AMBER") {
	    $resname=~s/^CYX$/CYS/;
	    $resname=~s/^HID$/HSD/;
	    $resname=~s/^HIE$/HSE/;
	    $resname=~s/^HIP$/HSP/;
            $resname=~s/^ASH$/ASP/;
            $resname=~s/^GLH$/GLU/;
            $resname=~s/^CTG$/GLY/;
        
            if ($atomname =~ /([0-9])(H.+)/) {
               $atomname=$2.$1;
            }
	    $atomname="HT1" if ($atomname eq "H1");
	    $atomname="HT2" if ($atomname eq "H2");
	    $atomname="HT3" if ($atomname eq "H3");
            $atomname="HN1" if ($atomname eq "HT2" && $resname=~/^PRO$/);
            $atomname="HN2" if ($atomname eq "HT3" && $resname=~/^PRO$/);
	    $atomname="HN" if ($atomname eq "H");
	    $atomname="HB1" if ($atomname eq "HB2" && $resname=~/^MET$|^ASP$|^ASN$|^GLU$|^GLN$|^TRP$|^PHE$|^TYR$|^LYS$|^ARG$|^LEU$|^PRO$|^SER$|^CYS$|^HSD$|^HSE$|^HSP$/);
	    $atomname="HB2" if ($atomname eq "HB3" && $resname=~/^MET$|^ASP$|^ASN$|^GLU$|^GLN$|^TRP$|^PHE$|^TYR$|^LYS$|^ARG$|^LEU$|^PRO$|^SER$|^CYS$|^HSD$|^HSE$|^HSP$/);
	    $atomname="HG1" if ($atomname eq "HG2" && $resname=~/^MET$|^GLU$|^GLN$|^LYS$|^ARG$|^PRO$/);
	    $atomname="HG2" if ($atomname eq "HG3" && $resname=~/^MET$|^GLU$|^GLN$|^LYS$|^ARG$|^PRO$/);
	    $atomname="HG1" if ($atomname eq "HSG" && $resname=~/^CYS$/);
	    $atomname="HG1" if ($atomname eq "HG" && $resname=~/^CYS$/);

	    $atomname="HD1" if ($atomname eq "HD2" && $resname=~/^LYS$|^ARG$|^PRO$/);
	    $atomname="HD2" if ($atomname eq "HD3" && $resname=~/^LYS$|^ARG$|^PRO$/);
	    $atomname="HE1" if ($atomname eq "HE2" && $resname=~/^LYS$/);
	    $atomname="HE2" if ($atomname eq "HE3" && $resname=~/^LYS$/);
	    $atomname="HA1" if ($atomname eq "HA2" && $resname=~/^GLY$/);
	    $atomname="HA2" if ($atomname eq "HA3" && $resname=~/^GLY$/);
	    $atomname="HG1"  if ($atomname eq "HG" && $resname=~/^SER$/);
	    $atomname="HD1" if ($atomname eq "HD11" && $resname=~/^ILE$/);
	    $atomname="HD2" if ($atomname eq "HD12" && $resname=~/^ILE$/);
	    $atomname="HD3" if ($atomname eq "HD13" && $resname=~/^ILE$/);
	    $atomname="HG11" if ($atomname eq "HG12" && $resname=~/^ILE$/);
	    $atomname="HG12" if ($atomname eq "HG13" && $resname=~/^ILE$/);
	    $atomname="HN1" if ($atomname eq "H2" && $resname=~/^PRO$/);
	    $atomname="HN2" if ($atomname eq "H3" && $resname=~/^PRO$/);
            
          } elsif ($translate eq "IMPACT") {
	    $resname="HSD" if ($resname eq "HID");
	    $resname="HSE" if ($resname eq "HIE");
	    $atomname="HT1" if ($atomname eq "H1");
	    $atomname="HT2" if ($atomname eq "H2");
	    $atomname="HT3" if ($atomname eq "H3");
	    $atomname="HN" if ($atomname eq "H");
            $atomname="HN1" if ($atomname eq "2H" && $resname eq "PRO");
            $atomname="HN2" if ($atomname eq "3H" && $resname eq "PRO");
	    $atomname="HA1" if ($atomname eq "1HA");
	    $atomname="HA2" if ($atomname eq "2HA");
	    $atomname="HB1" if ($atomname eq "1HB");
	    $atomname="HB2" if ($atomname eq "2HB");
	    $atomname="HB3" if ($atomname eq "3HB");
	    $atomname="HG1" if ($atomname eq "1HG");
	    $atomname="HG2" if ($atomname eq "2HG");
	    $atomname="HE1" if ($atomname eq "1HE");
	    $atomname="HE2" if ($atomname eq "2HE");
	    $atomname="HE3" if ($atomname eq "2HE");

	    $atomname="HD1" if ($atomname eq "1HD1" && $resname eq "ILE");
	    $atomname="HD2" if ($atomname eq "2HD1" && $resname eq "ILE");
	    $atomname="HD3" if ($atomname eq "3HD1" && $resname eq "ILE");

	    $atomname="HD11" if ($atomname eq "1HD1");
	    $atomname="HD12" if ($atomname eq "2HD1");
	    $atomname="HD13" if ($atomname eq "3HD1");
	    $atomname="HD21" if ($atomname eq "1HD2");
	    $atomname="HD22" if ($atomname eq "2HD2");
	    $atomname="HD23" if ($atomname eq "3HD2");
	    $atomname="HD1" if ($atomname eq "1HD");
	    $atomname="HD2" if ($atomname eq "2HD");
	    $atomname="HD3" if ($atomname eq "2HD");
	    $atomname="HZ1" if ($atomname eq "1HZ");
	    $atomname="HZ2" if ($atomname eq "2HZ");
	    $atomname="HZ3" if ($atomname eq "2HZ");
	    $atomname="HG11" if ($atomname eq "1HG1");
	    $atomname="HG12" if ($atomname eq "2HG1");
	    $atomname="HG13" if ($atomname eq "3HG1");
	    $atomname="HG21" if ($atomname eq "1HG2");
	    $atomname="HG22" if ($atomname eq "2HG2");
	    $atomname="HG23" if ($atomname eq "3HG2");
	    $atomname="HH11" if ($atomname eq "1HH1");
	    $atomname="HH12" if ($atomname eq "2HH1");
	    $atomname="HH21" if ($atomname eq "1HH2");
	    $atomname="HH22" if ($atomname eq "2HH2");
	    $atomname="HG1"  if ($atomname eq "HG" && $resname=~/SER|CYS/);
	  } else { 
#	    $resname=~s/CYX/CYS/;
#	    $resname=~s/HID/HSD/;
#	    $resname=~s/HIE/HSE/;
#	    $resname=~s/HIP/HSP/;
	  }

	  if ($chain ne $lastchain || ($seg ne $lastseg && $splitseg)) {
            my $crec=$self->{chainlookup}->{$chain};
	    $chainrec=(defined $crec)?$crec:$self->_newChain($chain);
	    $newchain=1;
            $addrnum=0;
	  } else {
	    $newchain=0;
	    $addrnum+=10000 if (($lastinum-$addrnum)==9999 && $iresnum<9999);
	  }

          $iresnum+=$addrnum;
	  
	  my $pdbrec={};
	 
          my $ainx=substr($_,4,7);
          if ($ainx!~/^[0-9\s]+$/) {
            $pdbrec->{atominx}=$lastainx+1;
          } else { 
  	    $pdbrec->{atominx}=$ainx+0;
          }
          $lastainx=$pdbrec->{atominx};
          
	  $pdbrec->{atomname}=$atomname;
	  $pdbrec->{resname}=$resname;
	  $pdbrec->{resnum}=$iresnum;
          $pdbrec->{aresnum}=$resnum if ($readalternate);

	  $pdbrec->{chain}=$chainrec->{id};
	  $pdbrec->{xcoor}=substr($_,30,8)+0.0;
	  $pdbrec->{ycoor}=substr($_,38,8)+0.0;
	  $pdbrec->{zcoor}=substr($_,46,8)+0.0;
	  $pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
	  $pdbrec->{aux1}=substr($_,54,6)+0.0;
	  $pdbrec->{aux2}=substr($_,60,6)+0.0;

	  $pdbrec->{seg}=$seg
	    unless (defined $ignoreseg && $ignoreseg);
          $pdbrec->{valid}=1;
          $pdbrec->{hetero}=$hetero;
	  
	  push (@{$chainrec->{atom}}, $pdbrec);

#	  printf "%s %s %d %s %s %s %s\n",$chainrec->{id},$resnum,$iresnum,$lastnum,$resname,$seg,$atomname;
	  
	  if ($resnum ne $lastnum || $newchain || ($seg ne $lastseg && (!defined $ignoreseg || !$ignoreseg))) {
	    my $resrec={};
	    $resrec->{name}=$resname;
	    $resrec->{num}=$iresnum;
            $resrec->{anum}=$resnum if ($readalternate);
	    $resrec->{chain}=$chainrec->{id};
	    $resrec->{start}=$#{$chainrec->{atom}};
	    $resrec->{end}=$resrec->{start};
	    $resrec->{valid}=1;
            $resrec->{hetero}=$hetero;
	    $resrec->{seg}=$seg
	      unless (defined $ignoreseg && $ignoreseg);
	    push(@{$chainrec->{res}},$resrec);
	  } else {
	    $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
	  }

	}
	$lastseg=$seg;
	$lastnum=$resnum;
	$lastinum=$iresnum;
	$lastchain=$chain;
      }
    } elsif (/^HETAT/ && (!defined $activemodel || $activemodel)) {
      my ($atomname, $resname, $resnum, $iresnum, $alt, $chain, $seg);
      
      ($atomname=substr($_,12,4))=~s/ //g;
      $atomname.=$alt if ($alt=~/[0-9]/);
      ($resname=substr($_,17,4))=~s/ //g;
      $resnum=substr($_,22,6);
      ($iresnum=$resnum)=~s/[A-Z]+//g;
      $iresnum+=0;
      ($seg=substr($_,72,4))=~s/[ \n]//g;

      $chain="+";

      my $crec=$self->{chainlookup}->{$chain};
      if (defined $crec) {
        $chainrec=$crec;
        $newchain=0;
      } else {
        $chainrec=$self->_newChain($chain);
	$newchain=1;
      }

      my $pdbrec={};

      my $ainx;
      if (/^HETATM/) {	
        $ainx=substr($_,6,5);
      } else {
        $ainx=substr($_,5,6);
      }
      if ($ainx!~/[0-9]/) {
        $pdbrec->{atominx}=$lastainx+1;
      } else { 
        $pdbrec->{atominx}=$ainx+0;
      }
      $lastainx=$pdbrec->{atominx};

      $pdbrec->{atomname}=$atomname;
      $pdbrec->{resname}=$resname;
      $pdbrec->{resnum}=$iresnum;

      $pdbrec->{chain}=$chainrec->{id};
      $pdbrec->{xcoor}=substr($_,30,8)+0.0;
      $pdbrec->{ycoor}=substr($_,38,8)+0.0;
      $pdbrec->{zcoor}=substr($_,46,8)+0.0;
      $pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
      $pdbrec->{aux1}=substr($_,55,6)+0.0;
      $pdbrec->{aux2}=substr($_,61,6)+0.0;
      $pdbrec->{seg}=$seg
	unless (defined $ignoreseg && $ignoreseg);
      $pdbrec->{valid}=1;
      $pdbrec->{hetero}=$hetero;
      
      push (@{$chainrec->{atom}}, $pdbrec);

#      printf "HA %s %d %d %s %s %s\n",$chainrec->{id},$resnum,$lastnum,$resname,$seg,$atomname;
      
      if ($resnum ne $lastnum || $newchain) {
	my $resrec={};
	$resrec->{name}=$resname;
	$resrec->{num}=$iresnum;
	$resrec->{chain}=$chainrec->{id};
	$resrec->{start}=$#{$chainrec->{atom}};
	$resrec->{end}=$resrec->{start};
	$resrec->{valid}=1;
        $resrec->{hetero}=$hetero;
	$resrec->{seg}=$seg
	  unless (defined $ignoreseg && $ignoreseg);
	push(@{$chainrec->{res}},$resrec);
      } else {
	$chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
      }
      $lastnum=$resnum;
      $lastinum=$iresnum;
    }
  }

  my $nchain=();
  my $keephet;
  foreach my $c ( @{$self->{chain}} ) {
    if ($c->{id} ne "+") {
      push(@{$nchain},$c);
    } else {
      $keephet=$c;
    }
  }

  if (defined $keephet) {
    push(@{$nchain},$keephet);
  }
  
  $self->{chain}=$nchain;

  $self->{segmentlist}=undef;

  $self->_coorCache();

  undef $fname;
}

## method: readCIF(file)
## reads a protein structures from an mmCIF file.

sub readCIF {
  my $self=shift;
  my $fname=&GenUtil::getInputFile(shift);

  my $datatag;
  my %data;
  my $looptag=undef;
  my $loop=0;
  my @keys=();
  while (<$fname>) {
    chomp;
    s/^\s+//;
    s/\s+$//;
    if (/^data_(.+)$/) {
      $datatag=$1;
    } elsif (/^#/) {
      $looptag=undef;
      $loop=0;
    } elsif (/loop_/) {
      $loop=1;
      @keys=();
    } elsif (!$loop && /^_(\S+)\.(\S+)\s+(\S+)/) {
      if (!exists $data{$1}) {
        $data{$1}=();
      }
      my $trec={};
      $trec->{$2}=$3;
      push(@{$data{$1}},$trec);
    } elsif (!$loop && /^_(\S+)\s+(\S+)/) {
      if (!exists $data{$1}) {
        $data{$1}=();
      }
      my $trec={};
      $trec->{value}=$2;
      push(@{$data{$1}},$trec);
    } elsif ($loop && /^_(\S+)\.(\S+)/) {
      if (!exists $data{$1}) {
        $looptag=$1;
        $data{$looptag}=();
      }
      push(@keys,$2);
    } elsif ($loop) {
      my @f=split(/\s+/);
      my $trec={};
      for (my $if=0; $if<=$#f; $if++) {
        $trec->{$keys[$if]}=$f[$if];
      }
      push(@{$data{$looptag}},$trec);
    }
  }

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;

  my $lastchain=".";
  my $lastseg=".";
  my $chainrec;

  my $lastnum=-999;
  my $ignore;

  my $newchain=0;

  if (!exists $data{"atom_site"}) {
    printf STDERR "cannot find atom information in mmCIF file\n";
    exit 1;
  }

  my $extchain=0;
  foreach my $t ( @{$data{"atom_site"}} ) {
    my $chain=$t->{auth_asym_id};
    if (length($chain)>1) {
      $extchain=1;
    }
  }

  foreach my $t ( @{$data{"atom_site"}} ) {
    my $atomname=$t->{label_atom_id};
    my $resnum=$t->{auth_seq_id};
    my $resname=$t->{auth_comp_id};
    my $chain=($extchain)?" ":$t->{auth_asym_id};
    my $seg=$t->{auth_asym_id};
    my $ainx=$t->{id};
    my $xcoor=$t->{Cartn_x};
    my $ycoor=$t->{Cartn_y};
    my $zcoor=$t->{Cartn_z};
    my $occupancy=$t->{occupancy};
    my $bvalue=$t->{B_iso_or_equiv};
    my $iresnum=$resnum+0;

    $atomname=~s/\"//g;
    $resname=~s/\"//g;
	  
    $atomname="CD1" 
      if ($resname eq "ILE" && $atomname eq "CD");
    $atomname="O"
      if (($atomname eq "OT1" || $atomname eq "O1" || $atomname eq "OCT1")
          && $resname =~/^ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX$/);
    $atomname="OXT"
      if (($atomname eq "OT2" || $atomname eq "O2" || $atomname eq "OCT2")
          && $resname =~/^ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX$/);
	
    if ($t->{group_PDB} eq "HETATM") {
      $chain="+";
    }
	
    if ($chain ne $lastchain || $seg ne $lastseg) { 
      my $crec=$self->{chainlookup}->{$chain};
      $chainrec=(defined $crec)?$crec:$self->_newChain($chain);
      $newchain=1;
    } else {
      $newchain=0;
    }
	
    my $pdbrec={};
	
    $pdbrec->{atominx}=$ainx; 
    $pdbrec->{atomname}=$atomname;
    $pdbrec->{resname}=$resname;
    $pdbrec->{resnum}=$resnum+0;
	
    $pdbrec->{chain}=$chainrec->{id};
    $pdbrec->{xcoor}=$xcoor;
    $pdbrec->{ycoor}=$ycoor;
    $pdbrec->{zcoor}=$zcoor;
    $pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
    $pdbrec->{seg}=$seg;

    $pdbrec->{aux1}=$occupancy;
    $pdbrec->{aux2}=$bvalue;

    $pdbrec->{valid}=1;
    	  
    push (@{$chainrec->{atom}}, $pdbrec);

#	printf "%s %s %d %s %s %s %s\n",$chainrec->{id},$resnum,$iresnum,$lastnum,$resname,$seg,$atomname;
	
    if ($iresnum != $lastnum || $newchain) {
      my $resrec={};
      $resrec->{name}=$resname;
      $resrec->{num}=$iresnum;
      $resrec->{chain}=$chainrec->{id};
      $resrec->{start}=$#{$chainrec->{atom}};
      $resrec->{end}=$resrec->{start};
      $resrec->{valid}=1;
      $resrec->{seg}=$seg;
      push(@{$chainrec->{res}},$resrec);
    } else {
      $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
    }
    $lastnum=$iresnum;
    $lastchain=$chain;
    $lastseg=$seg;
  }

  $self->{segmentlist}=undef;

  $self->_coorCache();

  undef $fname;
}


## method: readCRD(file)
## reads a protein structures from a CHARMM CRD file.

sub readCRD {
  my $self=shift;
  my $fname=&GenUtil::getInputFile(shift);

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;

  my $lastchain=".";
  my $chainrec;

  my $lastnum=-999;
  my $ignore;

  my $newchain=0;
  
  my $first=1;

  my $extendedformat=0;

  while(<$fname>) {
    if (!/^\*/) {
      if ($first) {
        if (/EXT/) {
          $extendedformat=1;
        } 
	$first=0;
      } else {
	my ($atomname, $resname, $resnum, $iresnum, $chain,$seg,$ainx,$xcoor,$ycoor,$zcoor);
     
        if ($extendedformat) { 
   	  ($atomname=substr($_,32,4))=~s/ //g;
	  ($resname=substr($_,22,4))=~s/ //g;
	  $resnum=substr($_,15,5);
	  ($seg=substr($_,102,4))=~s/[ \n]//g;
          $ainx=substr($_,0,10)+0;
	  $xcoor=substr($_,40,20)+0.0;
	  $ycoor=substr($_,60,20)+0.0;
	  $zcoor=substr($_,80,20)+0.0;
        } else {
   	  ($atomname=substr($_,16,4))=~s/ //g;
	  ($resname=substr($_,11,4))=~s/ //g;
	  $resnum=substr($_,5,5);
	  ($seg=substr($_,51,4))=~s/[ \n]//g;
          $ainx=substr($_,0,5)+0;
	  $xcoor=substr($_,20,10)+0.0;
	  $ycoor=substr($_,30,10)+0.0;
	  $zcoor=substr($_,40,10)+0.0;
        }
        
	$iresnum=$resnum+0;

	$chain=($resname eq "TIP3" || $resname eq "TIP4" || $resname eq "HOH" || $resname eq "SPC")?"W":" ";
        $chain=" ";
	  
	$atomname="CD1" 
	  if ($resname eq "ILE" && $atomname eq "CD");
	$atomname="O"
	  if (($atomname eq "OT1" || $atomname eq "O1" || $atomname eq "OCT1")
	      && $resname =~/^ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX$/);
	$atomname="OXT"
	  if (($atomname eq "OT2" || $atomname eq "O2" || $atomname eq "OCT2")
	      && $resname =~/^ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX$/);
	
	$resname=~s/^HSD$/HSE/;
	$resname=~s/^HIS$/HSD/;
	
	if ($chain ne $lastchain) {
	  my $crec=$self->{chainlookup}->{$chain};
	  $chainrec=(defined $crec)?$crec:$self->_newChain($chain);
	  $newchain=1;
	} else {
	  $newchain=0;
	}
	
	my $pdbrec={};
	
	$pdbrec->{atominx}=$ainx; 
	$pdbrec->{atomname}=$atomname;
	$pdbrec->{resname}=$resname;
	$pdbrec->{resnum}=$resnum+0;
	
	$pdbrec->{chain}=$chainrec->{id};
	$pdbrec->{xcoor}=$xcoor;
	$pdbrec->{ycoor}=$ycoor;
	$pdbrec->{zcoor}=$zcoor;
	$pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
	$pdbrec->{seg}=$seg;
        $pdbrec->{valid}=1;
	  
	push (@{$chainrec->{atom}}, $pdbrec);

#	printf "%s %s %d %s %s %s %s\n",$chainrec->{id},$resnum,$iresnum,$lastnum,$resname,$seg,$atomname;
	
	if ($iresnum != $lastnum || $newchain) {
	  my $resrec={};
	  $resrec->{name}=$resname;
	  $resrec->{num}=$iresnum;
	  $resrec->{chain}=$chainrec->{id};
	  $resrec->{start}=$#{$chainrec->{atom}};
	  $resrec->{end}=$resrec->{start};
	  $resrec->{valid}=1;
	  $resrec->{seg}=$seg;
	  push(@{$chainrec->{res}},$resrec);
	} else {
	  $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
	}
	$lastnum=$iresnum;
	$lastchain=$chain;
      }
    }
  }

  $self->{segmentlist}=undef;

  $self->_coorCache();

  undef $fname;
}


## method: readMol2(file)
## reads a molecular structure from a MOL2 file

sub readMol2 {
  my $self=shift;
  my $fname=&GenUtil::getInputFile(shift);

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;

  my $lastchain=".";
  my $chainrec;

  my $lastnum=-999;
  my $ignore;

  my $newchain=0;

  my $readatoms=0;
  my ($atomname, $resname, $resnum, $iresnum, $chain,$seg, $atominx, $xcoor,$ycoor,$zcoor);
  while(<$fname>) {
    if (/\@\<TRIPOS\>ATOM/) {
      $readatoms=1;
    } elsif (/\@\<TRIPOS\>/) {
      $readatoms=0;
    } elsif ($readatoms && /^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s[A-Z]\S*\s+(\d+)\s+(\S+)/) {
      $atominx=$1;
      $seg="";
      $resnum=$6;
      $resname=$7;
      $atomname=$2;
      $iresnum=$resnum+0;
      $xcoor=$3;
      $ycoor=$4;
      $zcoor=$5;

      $chain=($resname eq "TIP3" || $resname eq "TIP4" || $resname eq "HOH" || $resname eq "SPC")?"+":" ";
      
      if ($chain ne $lastchain) {
	my $crec=$self->{chainlookup}->{$chain};
	$chainrec=(defined $crec)?$crec:$self->_newChain($chain);
	$newchain=1;
      } else {
	$newchain=0;
      }
      
      my $pdbrec={};
      
      $pdbrec->{atominx}=$atominx;
      $pdbrec->{atomname}=$atomname;
      $pdbrec->{resname}=$resname;
      $pdbrec->{resnum}=$resnum+0;
      
      $pdbrec->{chain}=$chainrec->{id};
      $pdbrec->{xcoor}=$xcoor;
      $pdbrec->{ycoor}=$ycoor;
      $pdbrec->{zcoor}=$zcoor;
      $pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
      $pdbrec->{seg}=$seg;
      $pdbrec->{valid}=1;
      
      push (@{$chainrec->{atom}}, $pdbrec);

#	printf STDERR "%s %s %d %s %s %s %s\n",$chainrec->{id},$resnum,$iresnum,$lastnum,$resname,$seg,$atomname;
	
      if ($iresnum != $lastnum || $newchain) {
	my $resrec={};
	$resrec->{name}=$resname;
	$resrec->{num}=$iresnum;
	$resrec->{chain}=$chainrec->{id};
	$resrec->{start}=$#{$chainrec->{atom}};
	$resrec->{end}=$resrec->{start};
	$resrec->{valid}=1;
	$resrec->{seg}=$seg;
	push(@{$chainrec->{res}},$resrec);
      } else {
	$chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
      }
      $lastnum=$iresnum;
      $lastchain=$chain;
    }
  }

  $self->{segmentlist}=undef;

  $self->_coorCache();

  undef $fname;
} 

## method: readPSF(file)
## reads a molecular structure from a CHARMM PSF file.

sub readPSF {
  my $self=shift;
  my $fname=&GenUtil::getInputFile(shift);

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;

  my $lastchain=".";
  my $chainrec;

  my $lastnum=-999;
  my $ignore;

  my $newchain=0;

  my ($atomname, $atomtype,$resname, $resnum, $iresnum, $chain,$seg, $atominx,$charge,$mass);
  while(<$fname>) {
    if (/^\s*([\d]+) ([\S]...) +([\d]+) +([\S]+) +([\S]+) +([\S]+) +([\S]+) +([\d\.Ee\+\-]+) .+$/) {
        $atominx=$1;
        $seg=$2;
	$resnum=$3;
        $resname=$4;
	$atomname=$5;
	$atomtype=$6;
	$charge=$7;
	$mass=$8;
        $iresnum=$resnum+0;

        if ($seg=~/...([A-Za-z0-9\+\-\=\_])/) {
          $chain=$1;
	} elsif ($resname eq "TIP3" || $resname eq "TIP4" || $resname eq "HOH" || $resname eq "SPC") {
          $chain="+";
        } else {
          $chain=" ";
        }
	#$chain=($resname eq "TIP3" || $resname eq "HOH" || $resname eq "SPC")?"+":" ";
	  
	if ($chain ne $lastchain) {
	  my $crec=$self->{chainlookup}->{$chain};
	  $chainrec=(defined $crec)?$crec:$self->_newChain($chain);
	  $newchain=1;
	} else {
	  $newchain=0;
	}
	
	my $pdbrec={};
	
	$pdbrec->{atominx}=$atominx;
	$pdbrec->{atomname}=$atomname;
	$pdbrec->{atomtype}=$atomtype;
	$pdbrec->{resname}=$resname;
	$pdbrec->{resnum}=$resnum+0;
	
	$pdbrec->{chain}=$chainrec->{id};
	$pdbrec->{xcoor}=0.0;
	$pdbrec->{ycoor}=0.0;
	$pdbrec->{zcoor}=0.0;
	$pdbrec->{hyd}=($atomname=~/^[0-9]*H.*/)?1:0;
	$pdbrec->{seg}=$seg;
        $pdbrec->{valid}=1;
	$pdbrec->{mass}=$mass;
	$pdbrec->{charge}=$charge;
	  
	push (@{$chainrec->{atom}}, $pdbrec);

#	printf STDERR "%s %s %d %s %s %s %s\n",$chainrec->{id},$resnum,$iresnum,$lastnum,$resname,$seg,$atomname;
	
	if ($iresnum != $lastnum || $newchain) {
	  my $resrec={};
	  $resrec->{name}=$resname;
	  $resrec->{num}=$iresnum;
	  $resrec->{chain}=$chainrec->{id};
	  $resrec->{start}=$#{$chainrec->{atom}};
	  $resrec->{end}=$resrec->{start};
	  $resrec->{valid}=1;
	  $resrec->{seg}=$seg;
	  push(@{$chainrec->{res}},$resrec);
	} else {
	  $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
	}
	$lastnum=$iresnum;
	$lastchain=$chain;
    }
  }

  $self->{segmentlist}=undef;

  $self->_coorCache();

  undef $fname;
}

## method: writePDB(file[,translate=>format, longaux2=>1, ssbond=>0, cleanaux=>0])
## writes out the current structure in PDB format
## the output format may be specified through <mark>translate</mark>.
## Possible formats are <mark>CHARMM19</mark>, <mark>CHARMM22</mark>,
## <mark>AMBER</mark>, and <mark>GENERIC</mark>. A segment ID 
## as required by CHARMM may be given as the third argument

sub writePDB {
  my $self=shift;
  my $fname=&GenUtil::getOutputFile(shift);
  my %wpar=@_;

  my $translate=$wpar{translate};
  my $longaux2=$wpar{longaux2};
  my $ssbond=$wpar{ssbond};
  my $cleanaux=$wpar{cleanaux};
  my $dohetero=$wpar{dohetero};
  my $genresno=$wpar{genresno};
  my $writeend=$wpar{writeend};

  $dohetero=1 if (!defined $dohetero);
  $translate="" if (!defined $translate);
  $genresno=0 if (!defined $genresno);
  $writeend=1 if (!defined $writeend);

  die "empty molecule, nothing to write"
    if (!defined $self->{chain} || $#{$self->{chain}}<0);

  if (!defined $ssbond || $ssbond) {
    my %havess;
    my $sinx=0;
    foreach my $s ( @{$self->{ssbond}}) {
      my $c1=$self->getChain($s->{chain1});
      my $c2=$self->getChain($s->{chain2});
      if (defined $c1 && defined $c2) {
	my $r1=$self->getResidueInChain($s->{resnum1},$c1);
	my $r2=$self->getResidueInChain($s->{resnum2},$c2);

	if (defined $r1 && defined $r2 && 
	    $r1->{name} eq "CYS" && $r2->{name} eq "CYS") {
	  printf $fname "SSBOND%4d CYS %1s%5d    CYS %1s%5d\n",
	    ++$sinx,$s->{chain1},$s->{resnum1},$s->{chain2},$s->{resnum2}
	      unless ($havess{"$s->{chain1}:$s->{resnum1}:$s->{chain2}:$s->{resnum2}"});

	  $havess{"$s->{chain1}:$s->{resnum1}:$s->{chain2}:$s->{resnum2}"}=1;
	  $havess{"$s->{chain2}:$s->{resnum2}:$s->{chain1}:$s->{resnum1}"}=1;
	}
      }
    }
  }

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{atom}}>=0) {
      my $cterm=0;
      my $lastres=$c->{res}->[$#{$c->{res}}];
      for (my $i=$lastres->{start}; $i<=$lastres->{end}; $i++) {
	if ($c->{atom}->[$i]->{atomname} eq "OXT") {
	  $cterm=1;
	}
      }

      my $prevres;
      my $HisType="";
      foreach my $a ( @{$c->{atom}} ) {
	my $ta;
	%{$ta}=%{$a};

	if ($translate=~/CHA/) {
	  $ta->{atomname}="CD" 
	    if ($ta->{resname} eq "ILE" && $ta->{atomname} eq "CD1");
	  $ta->{atomname}="OT1"
	    if ($ta->{atomname} eq "O" && $ta->{resnum} eq $lastres->{num} && $cterm);
	  $ta->{atomname}="OT2"
	    if ($ta->{atomname} eq "OXT" && $ta->{resnum} eq $lastres->{num});
#	  $ta->{resname}=~s/HOH/TIP3/;

          if ($translate eq "CHARMM19") {
            $ta->{atomname}="H" if ($ta->{atomname} eq "HN");
            $ta->{atomname}="CLM" if ($ta->{atomname} eq "CAY");
            $ta->{atomname}="CAM" if ($ta->{atomname} eq "CY");
            $ta->{atomname}="OAM" if ($ta->{atomname} eq "OY"); 
            $ta->{atomname}="NCB" if ($ta->{atomname} eq "NT");
            $ta->{atomname}="HCB" if ($ta->{atomname} eq "HNT");
            $ta->{atomname}="CAC" if ($ta->{atomname} eq "CAT");
          }
	} 

	if (($ta->{resname}=~/^HSD$|^HSE$|^HIS$|^HSP$/) && ($translate =~ /CHARMM/)) {
	    if ((! defined $prevres) || ($ta->{resnum} != $prevres)) {
		$HisType=$self->getHisType($c,$ta->{resnum});
		if ($translate =~ /CHARMM19/) {
		    if ($HisType eq "HSD") {
			$HisType="HIS";
		    } elsif ($HisType eq "HSE") {
			$HisType="HSD";
		    } elsif ($HisType eq "HSP") {
			$HisType="HSC";
		    }
		} elsif (($HisType eq "HIS") && ($translate =~ /CHARMM22/)) {
		    $HisType="HSD";
		}
	    }
	    $ta->{resname}=$HisType;

	} elsif ($translate =~ /AMBER/) {
	  $ta->{resname}=~s/^HSD$/HID/;
	  $ta->{resname}=~s/^HSE$/HIE/;
	  $ta->{resname}=~s/^HSP$/HIP/;
#	  $ta->{atomname}="CD1" 
#	    if ($ta->{resname} eq "ILE" && $ta->{atomname} eq "CD");

	  if ($translate !~ /CHAMBER/) {
	    $ta->{atomname}="H1" if ($ta->{atomname} eq "HT1");
	    $ta->{atomname}="H2" if ($ta->{atomname} eq "HT2");
	    $ta->{atomname}="H3" if ($ta->{atomname} eq "HT3");
	    $ta->{atomname}="H" if ($ta->{atomname} eq "HN");
	    $ta->{atomname}="HB3" if ($ta->{atomname} eq "HB2" && $ta->{resname}=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HID|HIE/);
	    $ta->{atomname}="HB2" if ($ta->{atomname} eq "HB1" && $ta->{resname}=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HID|HIE/);
	    $ta->{atomname}="HG3" if ($ta->{atomname} eq "HG2" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $ta->{atomname}="HG2" if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $ta->{atomname}=(($ta->{resnum} eq $c->{atom}->[0]->{resnum}) || ($ta->{resnum} eq $lastres->{num} && $cterm))
	      ?"HSG":"HG" if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/CYS/);

	    $ta->{atomname}="HD3" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/LYS|ARG|PRO/);
	    $ta->{atomname}="HD2" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/LYS|ARG|PRO/);
	    $ta->{atomname}="HE3" if ($ta->{atomname} eq "HE2" && $ta->{resname}=~/LYS/);
	    $ta->{atomname}="HE2" if ($ta->{atomname} eq "HE1" && $ta->{resname}=~/LYS/);
	    $ta->{atomname}="HA3" if ($ta->{atomname} eq "HA2" && $ta->{resname}=~/GLY/);
	    $ta->{atomname}="HA2" if ($ta->{atomname} eq "HA1" && $ta->{resname}=~/GLY/);
	    $ta->{atomname}="HG"  if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/SER/);
	    $ta->{atomname}="HD11" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="HD12" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="HD13" if ($ta->{atomname} eq "HD3" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="HG13" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="HG12" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/ILE/);
	    $ta->{atomname}="H2" if ($ta->{atomname} eq "HN1" && $ta->{resname}=~/PRO/);
	    $ta->{atomname}="H3" if ($ta->{atomname} eq "HN2" && $ta->{resname}=~/PRO/);
	  } 
	} elsif ($translate =~ /IMPACT/) {
	  $ta->{resname}="HID" if ($ta->{resname} =~/HSD/);
	  $ta->{resname}="HIE" if ($ta->{resname} =~/HSE/);
	  $ta->{atomname}="1H" if ($ta->{atomname} eq "HT1");
	  $ta->{atomname}="2H" if ($ta->{atomname} eq "HT2");
	  $ta->{atomname}="3H" if ($ta->{atomname} eq "HT3");
	  $ta->{atomname}="H" if ($ta->{atomname} eq "HN");
	  $ta->{atomname}="2H" if ($ta->{atomname} eq "HN1");
	  $ta->{atomname}="3H" if ($ta->{atomname} eq "HN2");
	  $ta->{atomname}="3HB" if ($ta->{atomname} eq "HB3" && $ta->{resname}=~/ALA/);
	  $ta->{atomname}="2HB" if ($ta->{atomname} eq "HB2" && $ta->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $ta->{atomname}="1HB" if ($ta->{atomname} eq "HB1" && $ta->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $ta->{atomname}="2HG" if ($ta->{atomname} eq "HG2" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $ta->{atomname}="1HG" if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $ta->{atomname}="2HA" if ($ta->{atomname} eq "HA2" && $ta->{resname}=~/GLY/);
	  $ta->{atomname}="1HA" if ($ta->{atomname} eq "HA1" && $ta->{resname}=~/GLY/);
	  $ta->{atomname}="3HE" if ($ta->{atomname} eq "HE3" && $ta->{resname}=~/MET/);
	  $ta->{atomname}="2HE" if ($ta->{atomname} eq "HE2" && $ta->{resname}=~/MET|LYS/);
	  $ta->{atomname}="1HE" if ($ta->{atomname} eq "HE1" && $ta->{resname}=~/MET|LYS/);
	  $ta->{atomname}="2HD" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ARG|PRO/);
	  $ta->{atomname}="1HD" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ARG|PRO/);
	  $ta->{atomname}="1HH1" if ($ta->{atomname} eq "HH11" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="2HH1" if ($ta->{atomname} eq "HH12" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="1HH2" if ($ta->{atomname} eq "HH21" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="2HH2" if ($ta->{atomname} eq "HH22" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="1HE2" if ($ta->{atomname} eq "HE21" && $ta->{resname}=~/GLN/);
	  $ta->{atomname}="2HE2" if ($ta->{atomname} eq "HE22" && $ta->{resname}=~/GLN/);
	  $ta->{atomname}="1HD1" if ($ta->{atomname} eq "HD11" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="2HD1" if ($ta->{atomname} eq "HD12" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="3HD1" if ($ta->{atomname} eq "HD13" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="1HD2" if ($ta->{atomname} eq "HD21" && $ta->{resname}=~/LEU|ASN/);
	  $ta->{atomname}="2HD2" if ($ta->{atomname} eq "HD22" && $ta->{resname}=~/LEU|ASN/);
	  $ta->{atomname}="3HD2" if ($ta->{atomname} eq "HD23" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="1HD" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/LYS|ARG|PRO/);
	  $ta->{atomname}="2HD" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/LYS|ARG|PRO/);
	  $ta->{atomname}="3HZ" if ($ta->{atomname} eq "HZ3" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="2HZ" if ($ta->{atomname} eq "HZ2" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="1HZ" if ($ta->{atomname} eq "HZ1" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="1HG1" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="2HG1" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="3HG1" if ($ta->{atomname} eq "HG13" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="1HG2" if ($ta->{atomname} eq "HG21" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="2HG2" if ($ta->{atomname} eq "HG22" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="3HG2" if ($ta->{atomname} eq "HG23" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="HG"  if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/SER|CYS/);
	  $ta->{atomname}="1HD1" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HD1" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="3HD1" if ($ta->{atomname} eq "HD3" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="1HG1" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HG1" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="1HG2" if ($ta->{atomname} eq "HG21" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HG2" if ($ta->{atomname} eq "HG22" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="3HG2" if ($ta->{atomname} eq "HG23" && $ta->{resname}=~/ILE/);
	} elsif ($translate =~ /ICM/) {
	  $ta->{resname}="HIS" if ($ta->{resname} =~/HSD/);
	  $ta->{resname}="HIS" if ($ta->{resname} =~/HID/);
	  $ta->{resname}="HIE" if ($ta->{resname} =~/HSE/);
	  $ta->{atomname}="1H" if ($ta->{atomname} eq "HT1");
	  $ta->{atomname}="2H" if ($ta->{atomname} eq "HT2");
	  $ta->{atomname}="3H" if ($ta->{atomname} eq "HT3");
	  $ta->{atomname}="HN" if ($ta->{atomname} eq "HN");
	  $ta->{atomname}="2H" if ($ta->{atomname} eq "HN1");
	  $ta->{atomname}="3H" if ($ta->{atomname} eq "HN2");
	  $ta->{atomname}="3HB" if ($ta->{atomname} eq "HB3" && $ta->{resname}=~/ALA/);
	  $ta->{atomname}="2HB" if ($ta->{atomname} eq "HB2" && $ta->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $ta->{atomname}="1HB" if ($ta->{atomname} eq "HB1" && $ta->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $ta->{atomname}="2HG" if ($ta->{atomname} eq "HG2" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $ta->{atomname}="1HG" if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $ta->{atomname}="2HA" if ($ta->{atomname} eq "HA2" && $ta->{resname}=~/GLY/);
	  $ta->{atomname}="1HA" if ($ta->{atomname} eq "HA1" && $ta->{resname}=~/GLY/);
	  $ta->{atomname}="3HE" if ($ta->{atomname} eq "HE3" && $ta->{resname}=~/MET/);
	  $ta->{atomname}="2HE" if ($ta->{atomname} eq "HE2" && $ta->{resname}=~/MET|LYS/);
	  $ta->{atomname}="1HE" if ($ta->{atomname} eq "HE1" && $ta->{resname}=~/MET|LYS/);
	  $ta->{atomname}="2HD" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ARG|PRO/);
	  $ta->{atomname}="1HD" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ARG|PRO/);
	  $ta->{atomname}="1HH1" if ($ta->{atomname} eq "HH11" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="2HH1" if ($ta->{atomname} eq "HH12" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="1HH2" if ($ta->{atomname} eq "HH21" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="2HH2" if ($ta->{atomname} eq "HH22" && $ta->{resname}=~/ARG/);
	  $ta->{atomname}="1HE2" if ($ta->{atomname} eq "HE21" && $ta->{resname}=~/GLN/);
	  $ta->{atomname}="2HE2" if ($ta->{atomname} eq "HE22" && $ta->{resname}=~/GLN/);
	  $ta->{atomname}="1HD1" if ($ta->{atomname} eq "HD11" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="2HD1" if ($ta->{atomname} eq "HD12" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="3HD1" if ($ta->{atomname} eq "HD13" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="1HD2" if ($ta->{atomname} eq "HD21" && $ta->{resname}=~/LEU|ASN/);
	  $ta->{atomname}="2HD2" if ($ta->{atomname} eq "HD22" && $ta->{resname}=~/LEU|ASN/);
	  $ta->{atomname}="3HD2" if ($ta->{atomname} eq "HD23" && $ta->{resname}=~/LEU/);
	  $ta->{atomname}="1HD" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/LYS|ARG|PRO/);
	  $ta->{atomname}="2HD" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/LYS|ARG|PRO/);
	  $ta->{atomname}="3HZ" if ($ta->{atomname} eq "HZ3" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="2HZ" if ($ta->{atomname} eq "HZ2" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="1HZ" if ($ta->{atomname} eq "HZ1" && $ta->{resname}=~/LYS/);
	  $ta->{atomname}="1HG1" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="2HG1" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="3HG1" if ($ta->{atomname} eq "HG13" && $ta->{resname}=~/VAL/);
	  $ta->{atomname}="1HG2" if ($ta->{atomname} eq "HG21" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="2HG2" if ($ta->{atomname} eq "HG22" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="3HG2" if ($ta->{atomname} eq "HG23" && $ta->{resname}=~/VAL|THR/);
	  $ta->{atomname}="HG"  if ($ta->{atomname} eq "HG1" && $ta->{resname}=~/SER|CYS/);
	  $ta->{atomname}="1HD1" if ($ta->{atomname} eq "HD1" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HD1" if ($ta->{atomname} eq "HD2" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="3HD1" if ($ta->{atomname} eq "HD3" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="1HG1" if ($ta->{atomname} eq "HG11" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HG1" if ($ta->{atomname} eq "HG12" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="1HG2" if ($ta->{atomname} eq "HG21" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="2HG2" if ($ta->{atomname} eq "HG22" && $ta->{resname}=~/ILE/);
	  $ta->{atomname}="3HG2" if ($ta->{atomname} eq "HG23" && $ta->{resname}=~/ILE/);
	} elsif ($translate =~ /GENERIC/) {
	  $ta->{resname}=~s/^HSD$/HIS/;
	  $ta->{resname}=~s/^HSE$/HIS/;
	  $ta->{resname}=~s/^HSP$/HIS/;	  
	  $ta->{atomname}="CD1" 
	    if ($ta->{resname} eq "ILE" && $ta->{atomname} eq "CD");
	}
	
	if ($cleanaux) {
	  $ta->{aux1}=1.0;
	  $ta->{aux2}=0.0;
	}

	printf $fname "%s\n",&_pdbLine($ta,($translate=~/CHA/ && !$genresno),$longaux2)
	  if (($translate !~ /NOH/ || $ta->{atomname}!~/^[0-9]*H.*/) && !(!$dohetero && $ta->{hetero}));
	$prevres=$ta->{resnum};
      }
      printf $fname "TER\n";
    }
  }
  printf $fname "END\n" if ($writeend);
  
  undef $fname;
}

## method: writeCRD(file)
## writes out the current structure in PDB format

sub writeCRD {
  my $self=shift;
  my $fname=&GenUtil::getOutputFile(shift);
  my %wpar=@_;

  my $translate=$wpar{translate};
  my $extendedformat=$wpar{extend};

  $translate="" if (!defined $translate);

  die "empty molecule, nothing to write"
    if (!defined $self->{chain} || $#{$self->{chain}}<0);

  printf $fname "* CHARMM CRD file\n";
  printf $fname "* generated by convpdb.pl (MMTSB Tool Set)\n";
  printf $fname "*\n";

  my $natoms=0;
  foreach my $c ( @{$self->activeChains()} ) {
    $natoms+=$#{$c->{atom}}+1;
  }

  if ($extendedformat) {
    printf $fname "%10d  EXT\n",$natoms;
  } else {
    printf $fname "%d\n",$natoms;
  }
  
  my $ainx=0;
  my $rinx=0;
  my $prevseg="";
  my $prevres;
  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{atom}}>=0) {
      my %cterm;
      foreach my $r ( @{$c->{res}} ) {
	for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	  if ($c->{atom}->[$i]->{atomname} eq "OXT") {
	    $cterm{$r->{num}}=1;
	  }
	}
      }

      my $HisType="";
      foreach my $a ( @{$c->{atom}} ) {
	my $ta;
	%{$ta}=%{$a};

	if ($translate=~/CHA/) {
	  $ta->{atomname}="CD" 
	    if ($ta->{resname} eq "ILE" && $ta->{atomname} eq "CD1");
	  $ta->{atomname}="OT1"
	    if ($ta->{atomname} eq "O" && $cterm{$ta->{resnum}});
	  $ta->{atomname}="OT2"
	    if ($ta->{atomname} eq "OXT" && $cterm{$ta->{resnum}});
#	  $ta->{resname}=~s/HOH/TIP3/;

          if ($translate eq "CHARMM19") {
            $ta->{atomname}="H" if ($ta->{atomname} eq "HN");
            $ta->{atomname}="CLM" if ($ta->{atomname} eq "CAY");
            $ta->{atomname}="CAM" if ($ta->{atomname} eq "CY");
            $ta->{atomname}="OAM" if ($ta->{atomname} eq "OY"); 
            $ta->{atomname}="NCB" if ($ta->{atomname} eq "NT");
            $ta->{atomname}="HCB" if ($ta->{atomname} eq "HNT");
            $ta->{atomname}="CAC" if ($ta->{atomname} eq "CAT");
          }
	} 

	if (($ta->{resname}=~/^HSD$|^HSE$|^HIS$|^HSP$/) && ($translate =~ /CHARMM/)) {
	    if ((! defined $prevres) || ($ta->{resnum} != $prevres)) {
		$HisType=$self->getHisType($c,$ta->{resnum});
		if ($translate =~ /CHARMM19/) {
		    if ($HisType eq "HSD") {
			$HisType="HIS";
		    } elsif ($HisType eq "HSE") {
			$HisType="HSD";
		    } elsif ($HisType eq "HSP") {
			$HisType="HSC";
		    }
		} elsif (($HisType eq "HIS") && ($translate =~ /CHARMM22/)) {
		    $HisType="HSD";
		}
	    }
	    $ta->{resname}=$HisType;
	} elsif ($translate =~ /GENERIC/) {
	  $ta->{resname}=~s/HSD/HIS/;
	  $ta->{resname}=~s/HSE/HIS/;
	  $ta->{resname}=~s/HSP/HIS/;	  
	  $ta->{atomname}="CD1" 
	    if ($ta->{resname} eq "ILE" && $ta->{atomname} eq "CD");
	}

        $rinx++ if ($prevres!=$ta->{resnum} || $prevseg ne $ta->{seg});
        $ainx++;

        if ($extendedformat) {
  	  printf $fname "%10d%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-10d  %16.10f\n",
	    $ainx,$rinx,$ta->{resname},$ta->{atomname},$ta->{xcoor},$ta->{ycoor},$ta->{zcoor},$ta->{seg},$ta->{resnum},$ta->{aux1};
        } else {
  	  printf $fname "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %4s %-4d %9.5f\n",
	    $ainx,$rinx,$ta->{resname},$ta->{atomname},$ta->{xcoor},$ta->{ycoor},$ta->{zcoor},$ta->{seg},$ta->{resnum},$ta->{aux1};
        }

	$prevres=$ta->{resnum};
        $prevseg=$ta->{seg};
      }
    }
  }
  
  undef $fname;
}

## method: readAmber(partopfile,coorfile)
## reads protein structure from Amber topology
## and coordinate files

sub readAmberPre6 {
  my $self=shift;
  my $partop=&GenUtil::getInputFile(shift);
  my $coorfile=&GenUtil::getInputFile(shift);

  my $mode="";

  my $natom;
  my $ntype;
  my $nres;

  my @aname;
  my @lbres;
  my @ipres;

  while (<$partop>) {
    chomp;
    if (/\%FLAG (.+)/) {
      ($mode=$1)=~s/ +//g;
    } elsif ($mode eq "POINTERS") {
      my @num=&GenUtil::readFormat($partop,8,30);
      $natom=$num[0];
      $ntype=$num[1];
      $nres=$num[11];
      $mode="";
    } elsif ($mode eq "ATOM_NAME") {
      @aname=&GenUtil::readFormat($partop,4,$natom);
      $mode="";
    } elsif ($mode eq "RESIDUE_LABEL") {
      @lbres=&GenUtil::readFormat($partop,4,$nres);
      $mode="";
    } elsif ($mode eq "RESIDUE_POINTER") {
      @ipres=&GenUtil::readFormat($partop,8,$nres);
      $mode="";
    }
  }
  undef $partop;

  <$coorfile>;
  my $coornatom=<$coorfile>;
  chomp $coornatom;
  
  die "number of atoms does not match: expected: $natom, found: $coornatom"
    if ($coornatom != $natom);

  my @x=&GenUtil::readFormat($coorfile,12,$coornatom*3);
  undef $coorfile;

  $ipres[$nres]=$natom+1;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;
  
  my $chainrec=$self->_newChain("");

  my $ainx=1;
  for (my $ir=0; $ir<$nres; $ir++) {
    my $resrec={};
    ($resrec->{name}=$lbres[$ir])=~s/ +//g;
    $resrec->{num}=$ir+1;
    $resrec->{chain}="";
    $resrec->{start}=$ipres[$ir]-1;    
    $resrec->{end}=$ipres[$ir+1]-2;
    $resrec->{valid}=1;

    push (@{$chainrec->{res}},$resrec);
    
    for (my $ia=$ipres[$ir]; $ia<$ipres[$ir+1]; $ia++) {
      my $atomrec={};
      $atomrec->{atominx}=$ainx++;
      ($atomrec->{atomname}=$aname[$ia-1])=~s/ +//g;
      $atomrec->{resname}=$resrec->{name};
      $atomrec->{resnum}=$ir+1;
      $atomrec->{chain}="";
      $atomrec->{xcoor}=$x[($ia-1)*3];
      $atomrec->{ycoor}=$x[($ia-1)*3+1];
      $atomrec->{zcoor}=$x[($ia-1)*3+2];
      push (@{$chainrec->{atom}},$atomrec);
    }
  }

  $self->_coorCache();
}

## method: readAmber(partopfile,coorfile)
## reads protein structure from Amber topology
## and coordinate files
## Amber6

sub readAmber6 {
  my $self=shift;
  my $partop=&GenUtil::getInputFile(shift);
  my $coorfile=&GenUtil::getInputFile(shift);

  my $mode="";

  my $natom;
  my $ntype;
  my $nres;

  <$partop>;
  my @num=&GenUtil::readFormat($partop,6,30);
  $natom=$num[0];
  $ntype=$num[1];
  $nres=$num[11];
  my @aname=&GenUtil::readFormat($partop,4,$natom);
  my @junk=&GenUtil::readFormat($partop,16,$natom);
  @junk=&GenUtil::readFormat($partop,16,$natom);
  @junk=&GenUtil::readFormat($partop,6,$natom);
  @junk=&GenUtil::readFormat($partop,6,$natom);
  @junk=&GenUtil::readFormat($partop,6,$ntype*$ntype);
  my @lbres=&GenUtil::readFormat($partop,4,$nres);
  my @ipres=&GenUtil::readFormat($partop,6,$nres);
  undef $partop;

  <$coorfile>;
  my $coornatom=<$coorfile>;
  chomp $coornatom;
  
  die "number of atoms does not match: expected: $natom, found: $coornatom"
    if ($coornatom != $natom);

  my @x=&GenUtil::readFormat($coorfile,12,$coornatom*3);
  undef $coorfile;

  $ipres[$nres]=$natom+1;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;
  
  my $chainrec=$self->_newChain("");

  my $ainx=1;
  for (my $ir=0; $ir<$nres; $ir++) {
    my $resrec={};
    ($resrec->{name}=$lbres[$ir])=~s/ +//g;
    $resrec->{num}=$ir+1;
    $resrec->{chain}="";
    $resrec->{start}=$ipres[$ir]-1;    
    $resrec->{end}=$ipres[$ir+1]-2;
    $resrec->{valid}=1;

    push (@{$chainrec->{res}},$resrec);
    
    for (my $ia=$ipres[$ir]; $ia<$ipres[$ir+1]; $ia++) {
      my $atomrec={};
      $atomrec->{atominx}=$ainx++;
      ($atomrec->{atomname}=$aname[$ia-1])=~s/ +//g;
      $atomrec->{resname}=$resrec->{name};
      $atomrec->{resnum}=$ir+1;
      $atomrec->{chain}="";
      $atomrec->{xcoor}=$x[($ia-1)*3];
      $atomrec->{ycoor}=$x[($ia-1)*3+1];
      $atomrec->{zcoor}=$x[($ia-1)*3+2];
      push (@{$chainrec->{atom}},$atomrec);
    }
  }

  $self->_coorCache();
}

## method: readAmber(partopfile,coorfile)
## reads protein structure from Amber topology
## and coordinate files
## Amber7

sub readAmber {
  my $self=shift;
  my $partop=&GenUtil::getInputFile(shift);
  my $coorfile=&GenUtil::getInputFile(shift);

  my $mode="";
  my $items=0;
  my $len=0;

  my @num=();
  my @aname=();
  my @lbres=();
  my @ipres=();

  my $format=0;

  while (<$partop>) {
    chomp;

    if (/\%FLAG ([A-Za-z0-9_]+)/) {
      $mode=$1;
    } elsif ( /\%FORMAT\(([0-9]+)[a-zA-Z]([0-9]+).*\)/) { 
      $items=$1;
      $len=$2;
     } elsif ($mode eq "POINTERS") {
      for (my $i=0; $i+$len<=length($_); $i+=$len) {
	push(@num,substr($_,$i,$len)+0);
      }      
    } elsif ($mode eq "ATOM_NAME") {
      for (my $i=0; $i+$len<=length($_); $i+=$len) {
	push(@aname,substr($_,$i,$len));
      }      
    } elsif ($mode eq "RESIDUE_LABEL") {
      for (my $i=0; $i+$len<=length($_); $i+=$len) {
	push(@lbres,substr($_,$i,$len));
      }      
    } elsif ($mode eq "RESIDUE_POINTER") {
      for (my $i=0; $i+$len<=length($_); $i+=$len) {
	push(@ipres,substr($_,$i,$len)+0);
      }      
    }
  }

  my $natom=$num[0];
  my $ntype=$num[1];
  my $nres=$num[11];

  undef $partop;

  <$coorfile>;
  my $coornatom=<$coorfile>;
  chomp $coornatom;
  
  die "number of atoms does not match: expected: $natom, found: $coornatom"
    if ($coornatom != $natom);

  my @x=&GenUtil::readFormat($coorfile,12,$coornatom*3);
  undef $coorfile;

  $ipres[$nres]=$natom+1;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{selchain}=undef;
  
  my $chainrec=$self->_newChain("");

  my $ainx=1;
  for (my $ir=0; $ir<$nres; $ir++) {
    my $resrec={};
    ($resrec->{name}=$lbres[$ir])=~s/ +//g;
    $resrec->{num}=$ir+1;
    $resrec->{chain}="";
    $resrec->{start}=$ipres[$ir]-1;    
    $resrec->{end}=$ipres[$ir+1]-2;
    $resrec->{valid}=1;

    push (@{$chainrec->{res}},$resrec);
    
    for (my $ia=$ipres[$ir]; $ia<$ipres[$ir+1]; $ia++) {
      my $atomrec={};
      $atomrec->{atominx}=$ainx++;
      ($atomrec->{atomname}=$aname[$ia-1])=~s/ +//g;
      $atomrec->{resname}=$resrec->{name};
      $atomrec->{resnum}=$ir+1;
      $atomrec->{chain}="";
      $atomrec->{xcoor}=$x[($ia-1)*3];
      $atomrec->{ycoor}=$x[($ia-1)*3+1];
      $atomrec->{zcoor}=$x[($ia-1)*3+2];
      push (@{$chainrec->{atom}},$atomrec);
    }
  }

  $self->_coorCache();
}



## method: writeAmber(file)
## writes Amber coordinate file

sub writeAmber {
  my $self=shift;
  my $file=&GenUtil::getOutputFile(shift);
  
  my $c=$self->activeChains()->[0];

  print $file "\n";
  printf $file "%5d\n",$#{$c->{atom}}+1;
  my $i;
  my $a=$c->{atom};
  for ($i=0; $i<=$#{$a}; $i++) {
    printf $file "%12.7f%12.7f%12.7f",
    $a->[$i]->{xcoor}, $a->[$i]->{ycoor},$a->[$i]->{zcoor};
    printf $file "\n" if (($i%2)==1); 
  }
  printf $file "\n" if (($i%2)==1);
}

## method: slist = generateSegNames(segname)
## generates segment names for the current structure from
## the template given as the argument. It returns a list
## of generated segment names.

sub generateSegNames {
  my $self=shift;
  my $segname=shift;
  my $lastseg;

  $segname="PRO0" if (!defined $segname);

  die "empty molecule"
    if ($#{$self->activeChains()}<0);

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{res}}>=0 && (!defined $c->{res}->[0]->{seg} || $c->{res}->[0]->{seg} eq "")) {
      my @part;
      my $ipart=0;
      my $lastnum=-99999;
      for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
	my $r=$c->{res}->[$ir];
	my $num=$r->{num};
        $ipart++ if ($num!=$lastnum+1); 
# && ($num!=$c->{res}->[0]->{num} || $lastnum!=-99999));
	$part[$ir]=$ipart;
	$lastnum=$num;
      }

      my $nwatpart=0;
      for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
	my $r=$c->{res}->[$ir];
	if ($r->{name} eq "TIP3" || $r->{name} eq "HOH" || $r->{name} eq "SPC" || $r->{name} eq "TIP4" ) {
          $nwatpart=int($r->{num}/10000);
          $r->{altnum}=$r->{num}%10000;
          if ($r->{chain}=~/[A-Za-z0-9\-\_\=]/) {
  	    $r->{seg}=sprintf("W%02d%s",$nwatpart,$r->{chain});
          } else {
  	    $r->{seg}=sprintf("WT%02d",$nwatpart);
          }
	} elsif ($r->{name} =~ /GUA|ADE|URA|THY|CYT|PSU|5MU|G7M|AET|DHU|OMC|H2U|S4U|QUO|4SU|2MA|MAD|RT/) {
	  $r->{seg}=sprintf("N%02d%1s",$part[$ir],$c->{id});
        } elsif ($r->{chain} eq "+") {
	  $r->{seg}="HETA";
	} else {
	  if ($#{$self->{chain}}>0 || $self->{chain}->[0]->{id}=~/[A-Z]/) {
	    if ($ipart>1) {
	      $r->{seg}=sprintf("%1s%02d%1s",
		 substr($segname,0,1),$part[$ir],($c->{id} eq "" || $c->{id} eq " ")?"0":$c->{id});
	    } else {
	      $r->{seg}=sprintf("%3s%1s",
		 substr($segname,0,3),($c->{id} eq "" || $c->{id} eq " ")?"0":$c->{id});
	    }
	  } else {
	    if ($ipart>1) {
	      $r->{seg}=sprintf("%1s%03d",
                 substr($segname,0,1),$part[$ir]);
	    } else {
	      $r->{seg}=$segname;
	    }
	  }
	}

	for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	  $c->{atom}->[$i]->{seg}=$r->{seg};
	}
      }
    }
  }

  $self->{segmentlist}=undef;
}

## method: slist = generateSplitSegNames()
## generates segment names, splitting where chain
## breaks are found. Names go AA00-ZZ00. Returns a list
## of generated segment names.

sub generateSplitSegNames {
  my $self=shift;

  die "empty molecule" if ($#{$self->activeChains()}<0);

  my $segnum1=65;
  my $segnum2=65;
  my $segname=chr($segnum1).chr($segnum2)."00";
  my $Sqdist_cut=2.0*2.0;

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{res}}>=0) {
	my $r=$c->{res}->[0];
	my $prevC;
	if ($r->{name} eq "TIP3" || $r->{name} eq "TIP4" || $r->{name} eq "HOH" || $r->{name} eq "WAT") {
	  $r->{seg}="WATR";
	} elsif ($r->{name} =~ /GUA|ADE|URA|THY|CYT/) {
	  $r->{seg}="NA00";
        } else {
	  $r->{seg}=$segname;
	  for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	      my $a=$c->{atom}->[$i];
	      $a->{seg}=$r->{seg};
	      $prevC=$a if ($a->{atomname} eq "C");
	  }
        }

        for (my $ir=1; $ir<=$#{$c->{res}}; $ir++) {
	   my $r=$c->{res}->[$ir];
	   if ($r->{name} eq "TIP3" || $r->{name} eq "TIP4" || $r->{name} eq "HOH" || $r->{name} eq "WAT") {
	       $r->{seg}="WATR";
	   } elsif ($r->{name} =~ /GUA|ADE|URA|THY|CYT/) {
	       $r->{seg}="NA00";
	   } else {
  	       # If the peptide bond length to the previous residue is
               # defined and is reasonable, keep the segid, otherwise
               # increment
	       my $currN;
	       my $currC;
	       for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
		   my $a=$c->{atom}->[$i];
		   $currN=$a if ($a->{atomname} eq "N");
		   $currC=$a if ($a->{atomname} eq "C");
	       }

	       if ((defined $prevC) && (defined $currN)) {
		   my $x1=$prevC->{xcoor};
		   my $x2=$currN->{xcoor};
		   my $y1=$prevC->{ycoor};
		   my $y2=$currN->{ycoor};
		   my $z1=$prevC->{zcoor};
		   my $z2=$currN->{zcoor};
 		   my $SqDist=($x1-$x2)*($x1-$x2);
		   $SqDist+=($y1-$y2)*($y1-$y2);
		   $SqDist+=($z1-$z2)*($z1-$z2);
		   if ($SqDist > $Sqdist_cut) {
		       $segnum2++;
		       if ($segnum2 == 91) {
			   $segnum1++;
			   $segnum2=65;
		       }
		       $segname=chr($segnum1).chr($segnum2)."00";
		   }
	       }
	       $prevC=$currC;
	       $r->{seg}=$segname;
	       for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
		   my $a=$c->{atom}->[$i];
		   $a->{seg}=$r->{seg};
		   $prevC=$a if ($a->{atomname} eq "C");
	       }

	   }
	}
    }

    # Increment the segid for the new chain
    $segnum2++;
    if ($segnum2 == 91) {
	$segnum1++;
	$segnum2=65;
    }
    $segname=chr($segnum1).chr($segnum2)."00";

  }

  $self->{segmentlist}=undef;
}

## method: slist = getSegNames()
## generates and returns segment name list for the current structure 

sub getSegNames {
  my $self=shift;

  my $lastseg;

  if (!defined $self->{segmentlist} || $#{$self->{segmentlist}}<0) {
    $self->{segmentlist}=();

    my $rec;

    foreach my $c ( @{$self->activeChains()} ) {
      if ($#{$c->{res}}>=0) {
	$lastseg=undef;
        for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
          my $r=$c->{res}->[$ir];
#	foreach my $r ( @{$c->{res}} ) { 
          my $rseg=(defined $r->{seg})?$r->{seg}:"";
	  if (!defined $lastseg || $rseg ne $lastseg) {
	    $rec={};
	    $rec->{name}=$rseg;
	    $rec->{first}=(!defined $lastseg);
	    $rec->{last}=0;
	    $rec->{chain}=$c->{id};
	    $rec->{chainrec}=$c;
	    $rec->{from}=$r->{num};
	    $rec->{frominx}=$ir;
	    push(@{$self->{segmentlist}},$rec);
	    $lastseg=$rseg;
	  } 
	  $rec->{to}=$r->{num};
	  $rec->{toinx}=$ir;
	}
	$rec->{last}=1;
      }
    }
  } 

  return $self->{segmentlist};
}


## method: copySegNames(mol)
## sets the segment names from another Molecule object

sub copySegNames {
  my $self=shift;
  my $mol=shift;

  my %segs;
  foreach my $c ( @{$mol->{chain}} ) {
    foreach my $r ( @{$c->{res}} ) {
      my $key1=sprintf("%s:%s:%d",$c->{id},$r->{name},$r->{num});
      my $key2=sprintf(":%s:%d",$r->{name},$r->{num});
      $segs{$key1}=() if (!exists $segs{key1});
      push(@{$segs{$key1}},$r);
      $segs{$key2}=() if (!exists $segs{key2});
      push(@{$segs{$key2}},$r);
    }
  }

  foreach my $c ( @{$self->{chain}} ) {
    foreach my $r ( @{$c->{res}} ) {
      my $key1=sprintf("%s:%s:%d",$c->{id},$r->{name},$r->{num});
      my $key2=sprintf(":%s:%d",$r->{name},$r->{num});
      
      my $mr;
      if (exists $segs{$key1}) {
	if ($#{$segs{$key1}}>=0) {
	  $mr=shift @{$segs{$key1}};
	} 
      } elsif (exists $segs{$key2}) {
	if ($#{$segs{$key2}}>=0) {
	  $mr=shift @{$segs{$key2}};
	} 
      }
      
      if (defined $mr) {
	$r->{seg}=$mr->{seg} if (!defined $r->{seg});
	for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	  $c->{atom}->[$i]->{seg}=$r->{seg} 
	}
      }
    }
  }

  $self->{segmentlist}=undef;
}

## method: getChain([chainid])
## returns chain for given ID or first chain if
## no argument is given

sub getChain {
  my $self=shift;
  my $chainid=shift;

  return $self->{chain}->[0] 
    if (!defined $chainid || $chainid eq "");

  return $self->{chainlookup}->{$chainid};
}

## method: activeChains([chainid])
## returns a list of active chains

sub activeChains {
  my $self=shift;
  my $chainid=shift;

  my $list=();

  if (defined $self->{selchain}) {
    push(@{$list},$self->{selchain});
    return $list;
  } elsif (defined $chainid) {
    push(@{$list},$self->getChain($chainid));
    return $list;
  } elsif (defined $self->{defchain}) {
    push(@{$list},$self->getChain($self->{defchain}));
    return $list;
  } else {
    return $self->{chain};
  }
}

## method: selectChain(chainid)
## set a chain to be used as default

sub selectChain {
  my $self=shift;
  my $chainid=shift;
  
  $self->{defchain}=$chainid;
}

## method: numberAtoms()
## numbers all atoms continously across chain/segment boundaries

sub numberAtoms {
  my $self=shift;

  my $n=0;
  foreach my $c ( @{$self->activeChains()} ) {
    foreach my $a ( @{$c->{atom}} ) {
      $a->{tnum}=$n++;
    }
  }
}



## method: removeChain()

sub removeChain {
  my $self=shift;
  my $cid=shift;

  for (my $i=0; $i<=$#{$self->{chain}}; $i++) {
      if ($self->{chain}->[$i]->{id} eq $cid) {
	  splice(@{$self->{chain}},$i,1);
      }
  }
  undef $self->{chainlookup}->{$cid};

  return;
}

## method: removeHetero()
## if there's a heteroatom chain, remove it

sub removeHetero {
  my $self=shift;

  $self->removeChain("+");

  return;
}

## method: countInterfaceContacts([distCut][,chainA][,chainB][,CACBonly])
## Count the number of atomic contacts across an interface
## Operates on all pairs of chains if no chain IDs are specified
## Operates on all partners for given chain if one chain ID is specified
## Operates only on one pair if two chain IDs are specified
## Returns an array of hashes with {AAcontacts}, {cidA}, {cidB}

sub countInterfaceContacts {
  my $self=shift;
  my $distCut=shift;
  my $cidA=shift;
  my $cidB=shift;
  my $CACBonly=shift;

  $distCut=5 if (! defined $distCut);
  my $SqDistCut=$distCut*$distCut;

  $CACBonly=0 if (! defined $CACBonly);

  my @ac=@{$self->activeChains()};
  my $AAlist=();

  for (my $i=0; $i<=$#ac; $i++) {
    my $c1=$ac[$i];
    for (my $j=$i+1; $j<=$#ac; $j++) {
      my $c2=$ac[$j];
      if (($#{$c1->{res}}>=0) && ($#{$c2->{res}}>=0) &&
           ((! defined $cidA) || ($cidA eq $c1->{id}) ||
           ($cidA eq $c2->{id})) && ((! defined $cidB) ||
           ($cidB eq $c1->{id}) || ($cidB eq $c2->{id}))) {

	  my $AAcontacts=0;
	  foreach my $a1 ( @{$c1->{atom}} ) {
	    if (($a1->{aux1} >= 0.0) && (! $a1->{hyd}) && ((! $CACBonly) || 
                ($a1->{atomname} eq "CA") || ($a1->{atomname} eq "CB"))) {
	      my $x1=$a1->{xcoor};
	      my $y1=$a1->{ycoor};
	      my $z1=$a1->{zcoor};
	      foreach my $a2 ( @{$c2->{atom}} ) {
		if (($a2->{aux1} >= 0.0) && (! $a2->{hyd}) && ((! $CACBonly) ||
                    ($a2->{atomname} eq "CA") || ($a2->{atomname} eq "CB"))) {

		  my $x2=$a2->{xcoor};
		  my $y2=$a2->{ycoor};
		  my $z2=$a2->{zcoor};
		  my $SqDist=($x1-$x2)*($x1-$x2);
		  $SqDist+=($y1-$y2)*($y1-$y2);
		  $SqDist+=($z1-$z2)*($z1-$z2);
		  $AAcontacts++ if ($SqDist < $SqDistCut);

	        }
	      }
	    }
	  }
	  my $h={};
	  $h->{AAcontacts}=$AAcontacts;
	  $h->{cidA}=$c1->{id};
	  $h->{cidB}=$c2->{id};
	  push(@{$AAlist},$h);

      }
    }
  }

  return $AAlist;
}

## method: setValidSequence(sequencestring)
## sets the <mark>valid</mark> flag to 1
## for residues that match the given sequence string
## this method recognizes a previous residue subselection

sub setValidSequence {
  my $self=shift;
  my $seqstr=shift;

  my $foundany=0;

  foreach my $c ( @{$self->activeChains()} ) {
    my $ir;
    for ($ir=0; $ir<=$#{$c->{res}}+1-length($seqstr); $ir++) {
      my $tir=0;
      while ($c->{res}->[$tir+$ir]->{valid} && $tir<length($seqstr) && 
	     &_cmpResName($c->{res}->[$tir+$ir]->{name},$Sequence::_seqlong{substr($seqstr,$tir,1)})) {
	$tir++;
      }
      if ($tir<length($seqstr)) {
	$c->{res}->[$ir]->{valid}=0;
      } else {
	$foundany=1;
	$ir+=length($seqstr)-1;
      }
    }
    for ( ;$ir<=$#{$c->{res}}; $ir++) {
      $c->{res}->[$ir]->{valid}=0;
    }
  } 
  return $foundany;
}

## method: setValidChain(chain[,exclude])
## sets the <mark>valid</mark> flag to 1 for
## all residues in the given chain 

sub setValidChain {
  my $self=shift;
  my $chain=shift;
  my $exclmode=shift;
  my $noreset=shift;

  return if (!defined $chain);

  $self->resetValidResidues($exclmode?1:0) unless (defined $noreset && $noreset);

  $exclmode=0 if (!defined $exclmode);

  my $c=$self->getChain($chain);
  if (defined $c) {
    foreach my $r ( @{$c->{res}} ) {
      $r->{valid}=($exclmode)?0:1;
    }    
  }
}

## method: setValidResidues(fraglist)
## sets the <mark>valid</mark> flag to 1
## for residues in the fragment list and to
## 0 for residues outside the list

sub setValidResidues {
  my $self=shift;
  my $fraglist=shift;
  my $exclmode=shift;
  my $noreset=shift;

  return if (!defined $fraglist);

  $exclmode=0 if (!defined $exclmode);

  $self->resetValidResidues($exclmode?1:0) unless (defined $noreset && $noreset);

  foreach my $l ( @{$fraglist} ) {
    my $c=$self->getChain($l->{chain});
    if (defined $c) {
      if (!defined $l->{from}) {
	$l->{from}=$self->firstResNum($c);
	$l->{to}=$self->lastResNum($c);
      } elsif (!defined $l->{to}) {
	$l->{to}=$l->{from};
      }
      for (my $i=$l->{from}; $i<=$l->{to}; $i++) {
	my $r=$self->getResidueInChain($i,$c);
	$r->{valid}=($exclmode?0:1) if (defined $r);
      }    
    }
  }
}


## method: setValidSegment(segname[,exclmode])
## sets the <mark>valid</mark> flag to 1
## for residues with the given segment name 

sub setValidSegment {
  my $self=shift;
  my $segname=shift;
  my $exclmode=shift;

  return if (!defined $segname);
  $exclmode=0 if (!defined $exclmode);

  $self->resetValidResidues($exclmode?1:0);

  foreach my $c ( @{$self->activeChains()} ) {
    foreach my $r ( @{$c->{res}} ) {
      $r->{valid}=($exclmode?0:1) if ($r->{seg} eq $segname);
    }
  }
}


## method: markClashes() 
## finds and marks atom-atom clashes

sub markClashes {
  my $self=shift;
  my $clashes=shift;
  my $thresh=shift;
   
  $thresh=0.1 if (!defined $thresh);

  my @ax=();
  my @ay=();
  my @az=();
  my @alist=();

  my $validval=(defined $clashes && $clashes)?0:1;
  my $notvalidval=(defined $clashes && $clashes)?1:0;
  foreach my $c ( @{$self->{chain}} ) {
    my $atom=$c->{atom};
    foreach my $r ( @{$c->{res}} ) {
      $r->{valid}=$validval;
      for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	$atom->[$ia]->{valid}=$validval;
	push(@alist,$atom->[$ia]);
	push(@ax,$atom->[$ia]->{xcoor});
	push(@ay,$atom->[$ia]->{ycoor});
	push(@az,$atom->[$ia]->{zcoor});
      }
    }
  }

  my $nn=0;
  for (my $i=0; $i<=$#alist; $i++) {
    if ($alist[$i]->{valid}==$validval) {
      my $axi=$ax[$i];
      my $ayi=$ay[$i];
      my $azi=$az[$i];
      for (my $j=$i+1; $j<=$#alist; $j++) {
	my $dx=$axi-$ax[$j];
	my $dy=$ayi-$ay[$j];
	my $dz=$azi-$az[$j];
	my $d2=$dx*$dx+$dy*$dy+$dz*$dz;
	if ($d2<$thresh) {
	  $alist[$i]->{valid}=$notvalidval;
	  $alist[$j]->{valid}=$notvalidval;
	  $nn++;
	}
      }
    }
  }
}


## method: resetValidResidues([value])
## sets the <mark>valid</mark> flag of all
## residues to the given value (default: 1)

sub resetValidResidues {
  my $self=shift;
  my $value=shift;
  my $atoms=shift;

  $atoms=0 if (!defined $atoms);
  $value=1 if (!defined $value);
  
  foreach my $c ( @{$self->activeChains()} ) {
    foreach my $r ( @{$c->{res}} ) {
      $r->{valid}=$value;
    }
    if ($atoms) {
      foreach my $a ( @{$c->{atom}} ) {
	$a->{valid}=$value;
      }
    }
  }
}

## method: $list = listFromValid(force)
## returns a list of residues from the residues 
## previously set with <mark>setValidResidues</mark>.

sub listFromValid {
  my $self=shift;

  my $retlist=();
  foreach my $c ( @{$self->activeChains()} ) {
    my @arr;
    foreach my $r ( @{$c->{res}} ) {
      push(@arr,$r->{num}) if ($r->{valid});
    }
    if ($#arr>=0) {
      foreach my $tlist ( @{&GenUtil::fragListFromArray(\@arr,$c->{id})} ) {
	push(@{$retlist},$tlist);
      }
    } else {
      my $rec={};
      $rec->{from}=$self->firstResNum($c);
      $rec->{to}=$self->lastResNum($c);
      $rec->{chain}=$c->{id};
      push(@{$retlist},$rec);
    }
  }
  return $retlist;
}

## method: zapCoordinates([chain])
## sets all coordinates in the current structure
## to zero 

sub zapCoordinates {
  my $self=shift;
  my $chainid=shift;
  my $start=shift;
  my $stop=shift;

  foreach my $c ( @{$self->activeChains($chainid)} ) {
    foreach my $a ( @{$c->{atom}} ) {
      if (((! defined $start) || ($a->{resnum} >= $start)) &&
	  ((! defined $stop) || ($a->{resnum} <= $stop))) {
	$a->{xcoor}=0.0;
	$a->{ycoor}=0.0;
	$a->{zcoor}=0.0;
      }
    }
    $self->_coorCache($c->{id});
  }
}

## method: fillCoorFromPDB(file)
## updates only coordinates from a PDB file

sub fillCoorFromPDB {
  my $self=shift;
  my $pdb=Molecule::new(shift);

  foreach my $c ( @{$pdb->{chain}} ) {
    my $lc=$self->getChain($c->{id});

    if (defined $lc) {
      my %have;
      foreach my $a ( @{$lc->{atom}} ) {
	my $key=sprintf("%s%d:%s",
			$a->{resname},$a->{resnum},$a->{atomname});
	$have{$key}=$a;
      }

      foreach my $a ( @{$c->{atom}} ) {
	my $key=sprintf("%s%d:%s",
			$a->{resname},$a->{resnum},$a->{atomname});
	if (exists $have{$key}) {
	  $have{$key}->{xcoor}=$a->{xcoor};
	  $have{$key}->{ycoor}=$a->{ycoor};
	  $have{$key}->{zcoor}=$a->{zcoor};
	}
      }
      $self->_coorCache($c->{id});
    }
  }

  undef $pdb;
}

## method: merge(refmol)
## merges the current structure with the structure
## in <mark>refmol</mark> on a per residue basis

sub merge {
  my $self=shift;
  my $ref=shift;

  my $have={};
  my $cid=();

  foreach my $c ( @{$ref->activeChains()} ) {
    my $sc=$self->getChain($c->{id});
    if (!defined $sc) {
      my $crec=$self->_newChain($c->{id});
      my $ainx=1;
      foreach my $r ( @{$c->{res}} ) {
	if ($r->{valid}) {
	  my $rrec={};
	  %{$rrec}=%{$r};
	  push(@{$crec->{res}},$rrec);
	  $rrec->{start}=$#{$crec->{atom}}+1;
	  for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	    my $arec={};
	    %{$arec}=%{$c->{atom}->[$ia]};
	    $arec->{atominx}=$ainx++;
	    push(@{$crec->{atom}},$arec);
	  }
	  $rrec->{end}=$#{$crec->{atom}};
	}
      }
    } else {
      my %have;
      foreach my $r ( @{$sc->{res}} ) {
	my $trec={res=>$r, atom=>$sc->{atom}};
	$have{$r->{num}}=$trec;
      }
      foreach my $r ( @{$c->{res}} ) {
	if ($r->{valid}) {
	  my $trec={res=>$r, atom=>$c->{atom}};
	  $have{$r->{num}}=$trec;
	}
      }
      
      my $newatom=();
      my $newres=();
      my $ainx=1;
      foreach my $k ( sort { $a<=>$b } keys %have ) {
	my $r=$have{$k}->{res};
	my $a=$have{$k}->{atom};
	my $rrec={};
	%{$rrec}=%{$r};
	push(@{$newres},$rrec);
	$rrec->{start}=$#{$newatom}+1;
	for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	  my $arec={};
	  %{$arec}=%{$a->[$ia]};
	  $arec->{atominx}=$ainx++;
	  push(@{$newatom},$arec);
	}
	$rrec->{end}=$#{$newatom};
      }
      $sc->{atom}=$newatom;
      $sc->{res}=$newres;
    }
  }

  $self->{segmentlist}=undef;
  $self->_coorCache();
}

## method: $mol = clone([valid])
## generates a copy of the current molecule
## and returns it. If the valid flag is set, only 
## valid residues are copied

sub clone {
  my $self=shift;
  my $valid=shift;

  my $n={};
  $n->{chain}=();
  $n->{chainlookup}={};
  $n->{defchain}=undef;
  $n->{segmentlist}=undef;
  $n->{ssbond}=();

  bless $n;

  foreach my $s ( @{$self->{ssbond}} ) {
    my $srec={};
    %{$srec}=%{$s};
    push (@{$n->{ssbond}},$srec);
  }

  foreach my $c ( @{$self->activeChains()} ) {
    my $nc=undef;

    foreach my $r (@{$c->{res}}) {
      if ($r->{valid} || !defined $valid || !$valid) {
	$nc=$n->_newChain($c->{id}) if (!defined $nc);

	my $rrec={};
	%{$rrec}=%{$r};
	$rrec->{start}=$#{$nc->{atom}}+1;
	for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
         if ($c->{atom}->[$ia]->{valid} || !defined $valid || !$valid) {
	  my $arec={};
	  %{$arec}=%{$c->{atom}->[$ia]};
	  push(@{$nc->{atom}},$arec);
         }
	}
        if ($#{$nc->{atom}}>=$rrec->{start}) {
  	   $rrec->{end}=$#{$nc->{atom}};
	   push(@{$nc->{res}},$rrec);
        }
      }
    }
  }

  $n->_coorCache();

  return $n;
}

sub copycoor {
  my $self=shift;
  my $from=shift;
  my $to=shift;

  foreach my $c ( @{$self->activeChains()} ) {
    foreach my $a ( @{$c->{atom}}) {
      $a->{"x".$to}=$a->{"x".$from};
      $a->{"y".$to}=$a->{"y".$from};
      $a->{"z".$to}=$a->{"z".$from};
    }
  }
}

## method: translate(mode)
## translates the current structure according to the mode argument.
## As with <mark>writePDB</mark> the following modes are recognized:
## <mark>CHARMM19</mark>, <mark>CHARMM22</mark>,
## <mark>AMBER</mark>, and <mark>GENERIC</mark>.

sub translate {
  my $self=shift;
  my $translate=shift;

  $translate="" if (!defined $translate);

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{atom}}>=0) {
      my $cterm=0;
      my $lastres=$c->{res}->[$#{$c->{res}}];
      for (my $i=$lastres->{start}; $i<=$lastres->{end}; $i++) {
	if ($c->{atom}->[$i]->{atomname} eq "OXT") {
	  $cterm=1;
	}
      }

      foreach my $r ( @{$c->{res}} ) {
	if ($translate=~/CHARMM/) {
#	  $r->{name}=~s/HOH/TIP3/;
	  if ($r->{name}=~/^HSD$|^HSE$|^HIS$|^HSP$/) {
	    my $HisType=$self->getHisType($c,$r->{num});
	    if ($translate =~ /CHARMM19/) {
	      if ($HisType eq "HSD") {
		$HisType="HIS";
	      } elsif ($HisType eq "HSE") {
		$HisType="HSD";
	      } elsif ($HisType eq "HSP") {
		$HisType="HSC";
	      }
	    } elsif (($HisType eq "HIS") && ($translate =~ /CHARMM22/)) {
	      $HisType="HSD";
	    }
            for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	      $c->{atom}->[$i]->{resname}=$HisType;
	    }
	  }

        } elsif ($translate =~ /AMBER/) {
	  $r->{name}=~s/^HSD$/HID/;
	  $r->{name}=~s/^HSE$/HIE/;
	} elsif ($translate =~ /GENERIC/) {
	  $r->{name}=~s/^HSD$/HIS/;
	  $r->{name}=~s/^HSE$/HIS/;
	}
      }

      foreach my $a ( @{$c->{atom}} ) {

	if ($translate=~/CHARMM/) {
	  $a->{atomname}="CD" 
	    if ($a->{resname} eq "ILE" && $a->{atomname} eq "CD1");
	  $a->{atomname}="OT1"
	    if ($a->{atomname} eq "O" && $a->{resnum} eq $lastres->{num} && $cterm &&
		$a->{resname} =~/ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HIS|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX/);
	  $a->{atomname}="OT2"
	    if ($a->{atomname} eq "OXT" && $a->{resnum} eq $lastres->{num});
#	  $a->{resname}=~s/HOH/TIP3/;
	} 

	if ($translate =~ /AMBER/) {
	  $a->{resname}=~s/^HSD$/HID/;
	  $a->{resname}=~s/^HSE$/HIE/;
	  $a->{atomname}="CD1" 
	    if ($a->{resname} eq "ILE" && $a->{atomname} eq "CD");


	  if ($translate !~/CHAMBER/) {
	    $a->{atomname}="H1" if ($a->{atomname} eq "HT1");
	    $a->{atomname}="H2" if ($a->{atomname} eq "HT2");
	    $a->{atomname}="H3" if ($a->{atomname} eq "HT3");
	    $a->{atomname}="H" if ($a->{atomname} eq "HN");
	    $a->{atomname}="HB3" if ($a->{atomname} eq "HB2" && $a->{resname}=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HID|HIE/);
	    $a->{atomname}="HB2" if ($a->{atomname} eq "HB1" && $a->{resname}=~/MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HID|HIE/);
	    $a->{atomname}="HG3" if ($a->{atomname} eq "HG2" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $a->{atomname}="HG2" if ($a->{atomname} eq "HG1" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	    $a->{atomname}=(($a->{resnum} eq $c->{atom}->[0]->{resnum}) || ($a->{resnum} eq $lastres->{num} && $cterm))
	      ?"HSG":"HG" if ($a->{atomname} eq "HG1" && $a->{resname}=~/CYS/);
	    $a->{atomname}="HD3" if ($a->{atomname} eq "HD2" && $a->{resname}=~/LYS|ARG|PRO/);
	    $a->{atomname}="HD2" if ($a->{atomname} eq "HD1" && $a->{resname}=~/LYS|ARG|PRO/);
	    $a->{atomname}="HE3" if ($a->{atomname} eq "HE2" && $a->{resname}=~/LYS/);
	    $a->{atomname}="HE2" if ($a->{atomname} eq "HE1" && $a->{resname}=~/LYS/);
	    $a->{atomname}="HA3" if ($a->{atomname} eq "HA2" && $a->{resname}=~/GLY/);
	    $a->{atomname}="HA2" if ($a->{atomname} eq "HA1" && $a->{resname}=~/GLY/);
	    $a->{atomname}="HG"  if ($a->{atomname} eq "HG1" && $a->{resname}=~/SER/);
	    $a->{atomname}="HD11" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ILE/);
	    $a->{atomname}="HD12" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ILE/);
	    $a->{atomname}="HD13" if ($a->{atomname} eq "HD3" && $a->{resname}=~/ILE/);
	    $a->{atomname}="HG13" if ($a->{atomname} eq "HG12" && $a->{resname}=~/ILE/);
	    $a->{atomname}="HG12" if ($a->{atomname} eq "HG11" && $a->{resname}=~/ILE/);
	    $a->{atomname}="H2" if ($a->{atomname} eq "HN1" && $a->{resname}=~/PRO/);
	    $a->{atomname}="H3" if ($a->{atomname} eq "HN2" && $a->{resname}=~/PRO/);
	  }
	} elsif ($translate =~ /IMPACT/) {
	  $a->{resname}="HID" if ($a->{resname} =~/HSD/);
	  $a->{resname}="HIE" if ($a->{resname} =~/HSE/);
	  $a->{atomname}="1H" if ($a->{atomname} eq "HT1");
	  $a->{atomname}="2H" if ($a->{atomname} eq "HT2");
	  $a->{atomname}="3H" if ($a->{atomname} eq "HT3");
	  $a->{atomname}="H" if ($a->{atomname} eq "HN");
	  $a->{atomname}="2H" if ($a->{atomname} eq "HN1");
	  $a->{atomname}="3H" if ($a->{atomname} eq "HN2");
	  $a->{atomname}="3HB" if ($a->{atomname} eq "HB3" && $a->{resname}=~/ALA/);
	  $a->{atomname}="2HB" if ($a->{atomname} eq "HB2" && $a->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $a->{atomname}="1HB" if ($a->{atomname} eq "HB1" && $a->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $a->{atomname}="2HG" if ($a->{atomname} eq "HG2" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $a->{atomname}="1HG" if ($a->{atomname} eq "HG1" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $a->{atomname}="2HA" if ($a->{atomname} eq "HA2" && $a->{resname}=~/GLY/);
	  $a->{atomname}="1HA" if ($a->{atomname} eq "HA1" && $a->{resname}=~/GLY/);
	  $a->{atomname}="3HE" if ($a->{atomname} eq "HE3" && $a->{resname}=~/MET/);
	  $a->{atomname}="2HE" if ($a->{atomname} eq "HE2" && $a->{resname}=~/MET|LYS/);
	  $a->{atomname}="1HE" if ($a->{atomname} eq "HE1" && $a->{resname}=~/MET|LYS/);
	  $a->{atomname}="2HD" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ARG|PRO/);
	  $a->{atomname}="1HD" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ARG|PRO/);
	  $a->{atomname}="1HH1" if ($a->{atomname} eq "HH11" && $a->{resname}=~/ARG/);
	  $a->{atomname}="2HH1" if ($a->{atomname} eq "HH12" && $a->{resname}=~/ARG/);
	  $a->{atomname}="1HH2" if ($a->{atomname} eq "HH21" && $a->{resname}=~/ARG/);
	  $a->{atomname}="2HH2" if ($a->{atomname} eq "HH22" && $a->{resname}=~/ARG/);
	  $a->{atomname}="1HE2" if ($a->{atomname} eq "HE21" && $a->{resname}=~/GLN/);
	  $a->{atomname}="2HE2" if ($a->{atomname} eq "HE22" && $a->{resname}=~/GLN/);
	  $a->{atomname}="1HD1" if ($a->{atomname} eq "HD11" && $a->{resname}=~/LEU/);
	  $a->{atomname}="2HD1" if ($a->{atomname} eq "HD12" && $a->{resname}=~/LEU/);
	  $a->{atomname}="3HD1" if ($a->{atomname} eq "HD13" && $a->{resname}=~/LEU/);
	  $a->{atomname}="1HD2" if ($a->{atomname} eq "HD21" && $a->{resname}=~/LEU|ASN/);
	  $a->{atomname}="2HD2" if ($a->{atomname} eq "HD22" && $a->{resname}=~/LEU|ASN/);
	  $a->{atomname}="3HD2" if ($a->{atomname} eq "HD23" && $a->{resname}=~/LEU/);
	  $a->{atomname}="1HD" if ($a->{atomname} eq "HD1" && $a->{resname}=~/LYS|ARG|PRO/);
	  $a->{atomname}="2HD" if ($a->{atomname} eq "HD2" && $a->{resname}=~/LYS|ARG|PRO/);
	  $a->{atomname}="3HZ" if ($a->{atomname} eq "HZ3" && $a->{resname}=~/LYS/);
	  $a->{atomname}="2HZ" if ($a->{atomname} eq "HZ2" && $a->{resname}=~/LYS/);
	  $a->{atomname}="1HZ" if ($a->{atomname} eq "HZ1" && $a->{resname}=~/LYS/);
	  $a->{atomname}="1HG1" if ($a->{atomname} eq "HG11" && $a->{resname}=~/VAL/);
	  $a->{atomname}="2HG1" if ($a->{atomname} eq "HG12" && $a->{resname}=~/VAL/);
	  $a->{atomname}="3HG1" if ($a->{atomname} eq "HG13" && $a->{resname}=~/VAL/);
	  $a->{atomname}="1HG2" if ($a->{atomname} eq "HG21" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="2HG2" if ($a->{atomname} eq "HG22" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="3HG2" if ($a->{atomname} eq "HG23" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="HG"  if ($a->{atomname} eq "HG1" && $a->{resname}=~/SER|CYS/);
	  $a->{atomname}="1HD1" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HD1" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ILE/);
	  $a->{atomname}="3HD1" if ($a->{atomname} eq "HD3" && $a->{resname}=~/ILE/);
	  $a->{atomname}="1HG1" if ($a->{atomname} eq "HG11" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HG1" if ($a->{atomname} eq "HG12" && $a->{resname}=~/ILE/);
	  $a->{atomname}="1HG2" if ($a->{atomname} eq "HG21" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HG2" if ($a->{atomname} eq "HG22" && $a->{resname}=~/ILE/);
	  $a->{atomname}="3HG2" if ($a->{atomname} eq "HG23" && $a->{resname}=~/ILE/);
	} elsif ($translate =~ /ICM/) {
	  $a->{resname}="HIS" if ($a->{resname} =~/HSD/);
	  $a->{resname}="HIS" if ($a->{resname} =~/HID/);
	  $a->{resname}="HIE" if ($a->{resname} =~/HSE/);
	  $a->{atomname}="1H" if ($a->{atomname} eq "HT1");
	  $a->{atomname}="2H" if ($a->{atomname} eq "HT2");
	  $a->{atomname}="3H" if ($a->{atomname} eq "HT3");
	  $a->{atomname}="HN" if ($a->{atomname} eq "HN");
	  $a->{atomname}="2H" if ($a->{atomname} eq "HN1");
	  $a->{atomname}="3H" if ($a->{atomname} eq "HN2");
	  $a->{atomname}="3HB" if ($a->{atomname} eq "HB3" && $a->{resname}=~/ALA/);
	  $a->{atomname}="2HB" if ($a->{atomname} eq "HB2" && $a->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $a->{atomname}="1HB" if ($a->{atomname} eq "HB1" && $a->{resname}=~/ALA|MET|ASP|ASN|GLU|GLN|TRP|PHE|TYR|LYS|ARG|LEU|PRO|SER|CYS|HIS|HID|HIE/);
	  $a->{atomname}="2HG" if ($a->{atomname} eq "HG2" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $a->{atomname}="1HG" if ($a->{atomname} eq "HG1" && $a->{resname}=~/MET|GLU|GLN|LYS|ARG|PRO/);
	  $a->{atomname}="2HA" if ($a->{atomname} eq "HA2" && $a->{resname}=~/GLY/);
	  $a->{atomname}="1HA" if ($a->{atomname} eq "HA1" && $a->{resname}=~/GLY/);
	  $a->{atomname}="3HE" if ($a->{atomname} eq "HE3" && $a->{resname}=~/MET/);
	  $a->{atomname}="2HE" if ($a->{atomname} eq "HE2" && $a->{resname}=~/MET|LYS/);
	  $a->{atomname}="1HE" if ($a->{atomname} eq "HE1" && $a->{resname}=~/MET|LYS/);
	  $a->{atomname}="2HD" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ARG|PRO/);
	  $a->{atomname}="1HD" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ARG|PRO/);
	  $a->{atomname}="1HH1" if ($a->{atomname} eq "HH11" && $a->{resname}=~/ARG/);
	  $a->{atomname}="2HH1" if ($a->{atomname} eq "HH12" && $a->{resname}=~/ARG/);
	  $a->{atomname}="1HH2" if ($a->{atomname} eq "HH21" && $a->{resname}=~/ARG/);
	  $a->{atomname}="2HH2" if ($a->{atomname} eq "HH22" && $a->{resname}=~/ARG/);
	  $a->{atomname}="1HE2" if ($a->{atomname} eq "HE21" && $a->{resname}=~/GLN/);
	  $a->{atomname}="2HE2" if ($a->{atomname} eq "HE22" && $a->{resname}=~/GLN/);
	  $a->{atomname}="1HD1" if ($a->{atomname} eq "HD11" && $a->{resname}=~/LEU/);
	  $a->{atomname}="2HD1" if ($a->{atomname} eq "HD12" && $a->{resname}=~/LEU/);
	  $a->{atomname}="3HD1" if ($a->{atomname} eq "HD13" && $a->{resname}=~/LEU/);
	  $a->{atomname}="1HD2" if ($a->{atomname} eq "HD21" && $a->{resname}=~/LEU|ASN/);
	  $a->{atomname}="2HD2" if ($a->{atomname} eq "HD22" && $a->{resname}=~/LEU|ASN/);
	  $a->{atomname}="3HD2" if ($a->{atomname} eq "HD23" && $a->{resname}=~/LEU/);
	  $a->{atomname}="1HD" if ($a->{atomname} eq "HD1" && $a->{resname}=~/LYS|ARG|PRO/);
	  $a->{atomname}="2HD" if ($a->{atomname} eq "HD2" && $a->{resname}=~/LYS|ARG|PRO/);
	  $a->{atomname}="3HZ" if ($a->{atomname} eq "HZ3" && $a->{resname}=~/LYS/);
	  $a->{atomname}="2HZ" if ($a->{atomname} eq "HZ2" && $a->{resname}=~/LYS/);
	  $a->{atomname}="1HZ" if ($a->{atomname} eq "HZ1" && $a->{resname}=~/LYS/);
	  $a->{atomname}="1HG1" if ($a->{atomname} eq "HG11" && $a->{resname}=~/VAL/);
	  $a->{atomname}="2HG1" if ($a->{atomname} eq "HG12" && $a->{resname}=~/VAL/);
	  $a->{atomname}="3HG1" if ($a->{atomname} eq "HG13" && $a->{resname}=~/VAL/);
	  $a->{atomname}="1HG2" if ($a->{atomname} eq "HG21" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="2HG2" if ($a->{atomname} eq "HG22" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="3HG2" if ($a->{atomname} eq "HG23" && $a->{resname}=~/VAL|THR/);
	  $a->{atomname}="HG"  if ($a->{atomname} eq "HG1" && $a->{resname}=~/SER|CYS/);
	  $a->{atomname}="1HD1" if ($a->{atomname} eq "HD1" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HD1" if ($a->{atomname} eq "HD2" && $a->{resname}=~/ILE/);
	  $a->{atomname}="3HD1" if ($a->{atomname} eq "HD3" && $a->{resname}=~/ILE/);
	  $a->{atomname}="1HG1" if ($a->{atomname} eq "HG11" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HG1" if ($a->{atomname} eq "HG12" && $a->{resname}=~/ILE/);
	  $a->{atomname}="1HG2" if ($a->{atomname} eq "HG21" && $a->{resname}=~/ILE/);
	  $a->{atomname}="2HG2" if ($a->{atomname} eq "HG22" && $a->{resname}=~/ILE/);
	  $a->{atomname}="3HG2" if ($a->{atomname} eq "HG23" && $a->{resname}=~/ILE/);
	} elsif ($translate =~ /GENERIC/) {
	  $a->{resname}=~s/^HSD$/HIS/;
	  $a->{resname}=~s/^HSE$/HIS/;
	  $a->{resname}=~s/^HSP$/HIS/;
	  $a->{atomname}="CD1" 
	    if ($a->{resname} eq "ILE" && $a->{atomname} eq "CD");
	  $a->{atomname}="O"
	    if ($a->{atomname} eq "OT1");
	  $a->{atomname}="OXT"
	    if ($a->{atomname} eq "OT2");
	}
      }
    }
  }
}

## method: $map = renumber([map|startindex][,chain])
## renumbers the residues either according
## to a map given as an argument or continuously
## beginning from a start index. The method returns
## the translation map which may be applied
## to other structures or used to revert the
## residue numbering to its original form using
## <mark>numberReset</mark>

sub renumber {
  my $self=shift;

  my $start=1;
  my $map;
  my $t=shift;
  if (defined $t) {
    if (ref $t) {
      $map=$t;
    } else {
      $start=$t;
    }
  }

  my $c=$self->activeChains(shift)->[0];

  my $atom=$c->{atom};
  my $res=$c->{res};

  my $nres=$start-1;
  my $lastnum=-999;
  my $lastseg="XXX";
  my $lastares="-999";

  my $havemap=(defined $map);

  for (my $i=0; $i<=$#{$atom}; $i++) {  
    $atom->[$i]->{atominx}=$i+1;

    if ($atom->[$i]->{resnum} != $lastnum || 
        (defined $atom->[$i]->{aresnum} && $atom->[$i]->{aresnum} ne $lastares) ||
        $atom->[$i]->{seg} ne $lastseg) {
      $lastnum=$atom->[$i]->{resnum};
      $lastares=$atom->[$i]->{aresnum} if (defined $atom->[$i]->{aresnum});
      $lastseg=$atom->[$i]->{seg};
      $nres++;
    }

    my $newnum=$havemap ? $map->{$atom->[$i]->{resnum}} : $nres;
    $atom->[$i]->{resnum}=$newnum;
    $atom->[$i]->{aresnum}=undef if (defined $atom->[$i]->{aresnum});
  }

  my @nn;
  for (my $i=0; $i<=$#{$res}; $i++) {
    my $newnum;
    if ($havemap) {
      $newnum=$map->{$res->[$i]->{num}};
    } else {
      $map->{$res->[$i]->{num}}=$i+$start;
      $newnum=$i+$start;
    }
    $nn[$i]=$newnum;
  }

  foreach my $s ( @{$self->{ssbond}} ) {
  ADD1:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum1}==$res->[$i]->{num} &&
	  $s->{chain1}==$res->[$i]->{chain}) {
	$s->{resnum1}=$nn[$i];
	last ADD1;
      }
    }
  ADD2:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum2}==$res->[$i]->{num} &&
	  $s->{chain2}==$res->[$i]->{chain}) {
	$s->{resnum2}=$nn[$i];
	last ADD2;
      }
    }
  }

  for (my $i=0; $i<=$#{$res}; $i++) {
    $res->[$i]->{num}=$nn[$i];
    $res->[$i]->{anum}=undef if (defined $res->[$i]->{anum});
  }

  $c->{resinx}=undef;

  return $map;
}

## method: renumberWaterSegments()
## renumber water segments to range from 1-10000 and avoid overflow
## water segments are identifed by WT* labels

sub renumberWaterSegments {
  my $self=shift;

  my $renum=0;

  my $slist=$self->getSegNames();
  foreach my $s ( @{$slist} ) {
    if ($s->{name}=~/^WT/) {
      my $c=$s->{chainrec};
      my $r=$c->{res};
      my $a=$c->{atom};
      my $num=1;
      for (my $ir=$s->{frominx}; $ir<=$s->{toinx}; $ir++) {
        my $tr=$r->[$ir];
        $tr->{num}=$num;
        for (my $ia=$tr->{start}; $ia<=$tr->{end}; $ia++) {
          $a->[$ia]->{resnum}=$num;
        } 
        $num++;
      }
      $renum=1;
    } 
  }

  if (!$renum) {
    my $num=9999;
    my $nseg=0;
    my $segname;
    foreach my $c (@{$self->{chain}}) {
      foreach my $a ( @{$c->{atom}} ) {
        if ($a->{resname} eq "TIP3" || $a->{resname} eq "HOH") {
          if ($a->{atomname} eq "OH2") {
            $num++;
          }  
          if ($num>9999) {
            $num=1;
            $nseg++;
            $segname=sprintf("WT%02d",$nseg);
          }
          $a->{resnum}=$num;
          $a->{seg}=$segname;
        }       
      }
    }  
  }
}

## method: numberReset(map[,chain])
## applies the reverse translation from the one
## given in the map argument

sub numberReset {
  my $self=shift;
  my $map=shift;

  my $c=$self->activeChains(shift)->[0];

  my $reverseMap={};
  
  foreach my $m ( keys %{$map} ) {
    $reverseMap->{$map->{$m}}=$m;
  }

  my $atom=$c->{atom};
  my $res=$c->{res};

  for (my $i=0; $i<=$#{$atom}; $i++) {  
    my $newnum=$reverseMap->{$atom->[$i]->{resnum}};
    $atom->[$i]->{resnum}=$newnum;
  }

  foreach my $s ( @{$self->{ssbond}} ) {
  ADD1:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum1}==$res->[$i]->{num} &&
	  $s->{chain1}==$res->[$i]->{chain}) {
	$s->{resnum1}=$reverseMap->{$res->[$i]->{num}};
	last ADD1;
      }
    }
  ADD2:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum2}==$res->[$i]->{num} &&
	  $s->{chain2}==$res->[$i]->{chain}) {
	$s->{resnum2}=$reverseMap->{$res->[$i]->{num}};
	last ADD2;
      }
    }
  }

  for (my $i=0; $i<=$#{$res}; $i++) {
    my $newnum=$reverseMap->{$res->[$i]->{num}};
    $res->[$i]->{num}=$newnum;
  }

  $c->{resinx}=undef;
}

## method: shiftResNumber(add[,chain])
## shifts all residue numbers by a constant. A chain ID
## may be given to shift residues in that chain

sub shiftResNumber {
  my $self=shift;
  my $add=shift;
  my $cid=shift;

  my $c=$self->activeChains($cid)->[0];

  my $atom=$c->{atom};
  my $res=$c->{res};

  for (my $i=0; $i<=$#{$atom}; $i++) {  
    $atom->[$i]->{resnum}+=$add
  }

  foreach my $s ( @{$self->{ssbond}} ) {
  ADD1:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum1}==$res->[$i]->{num} &&
	  $s->{chain1}==$res->[$i]->{chain}) {
	$s->{resnum1}+=$add;
	last ADD1;
      }
    }
  ADD2:
    for (my $i=0; $i<=$#{$res}; $i++) {
      if ($s->{resnum2}==$res->[$i]->{num} &&
	  $s->{chain2}==$res->[$i]->{chain}) {
	$s->{resnum2}+=$add;
	last ADD2;
      }
    }
  }

  for (my $i=0; $i<=$#{$res}; $i++) {
    $res->[$i]->{num}+=$add
  }
  
  $c->{resinx}=undef;
}

## method: ($cx,$cy,$cz) = centerOfMass()
## calculates the center of mass for the current structure

sub centerOfMass {
  my $self=shift;
  my $chain=shift;
  
  my $n;
  my ($cx,$cy,$cz)=(0.0,0.0,0.0);

  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    $n+=$#{$atom}+1;
    
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $cx+=$atom->[$i]->{xcoor};
      $cy+=$atom->[$i]->{ycoor};
      $cz+=$atom->[$i]->{zcoor};
    }
  }

  $cx/=$n;
  $cy/=$n;
  $cz/=$n;

  return ($cx,$cy,$cz);
}  

## method: wrap(by,boxx,boxy,boxz)
## wraps the current structure with respect to the origin
## by: atom, 
##     chain (based on center of mass for each system), 
##     system (entire system)
##     reimage (needs scx,scy,scz)

sub wrap {
  my $self=shift;
  my $by=shift;
  my $boxx=shift;
  my $boxy=shift;
  my $boxz=shift;
  my $scx=shift;
  my $scy=shift;
  my $scz=shift;

  $boxx=10 if (!defined $boxx || $boxx<=0.001);
  $boxy=10 if (!defined $boxy || $boxy<=0.001);
  $boxz=10 if (!defined $boxz || $boxz<=0.001);
  
  if ($by eq "system") {
    my ($cx,$cy,$cz)=$self->centerOfMass();
    my $dx=$cx;
    my $dy=$cy;
    my $dz=$cz;
    $dx-=$boxx*&GenUtil::nint($dx/$boxx); 
    $dy-=$boxy*&GenUtil::nint($dy/$boxy); 
    $dz-=$boxz*&GenUtil::nint($dz/$boxz); 
    $self->move($dx-$cx,$dy-$cy,$dz-$cz); 
  } elsif ($by eq "chain") {
    foreach my $c ( @{$self->{chain}} ) {
      my ($cx,$cy,$cz)=(0.0,0.0,0.0);
      my $atom=$c->{atom};
      my $n+=$#{$atom}+1;
    
      for (my $i=0; $i<=$#{$atom}; $i++) {
        $cx+=$atom->[$i]->{xcoor};
        $cy+=$atom->[$i]->{ycoor};
        $cz+=$atom->[$i]->{zcoor};
      }
      $cx/=$n;
      $cy/=$n;
      $cz/=$n;

      my $dx=$cx;
      my $dy=$cy;
      my $dz=$cz;
      $dx-=$boxx*&GenUtil::nint($dx/$boxx); 
      $dy-=$boxy*&GenUtil::nint($dy/$boxy); 
      $dz-=$boxz*&GenUtil::nint($dz/$boxz); 
      for (my $i=0; $i<=$#{$atom}; $i++) {
        $atom->[$i]->{xcoor}+=$dx-$cx;
        $atom->[$i]->{ycoor}+=$dy-$cy;
        $atom->[$i]->{zcoor}+=$dz-$cz;
      }
    } 
  } elsif ($by eq "atom") {
    foreach my $c ( @{$self->{chain}} ) {
      my $atom=$c->{atom};
      for (my $i=0; $i<=$#{$atom}; $i++) {
        my $cx=$atom->[$i]->{xcoor};
        my $cy=$atom->[$i]->{ycoor};
        my $cz=$atom->[$i]->{zcoor};
        my $dx=$cx;
        my $dy=$cy;
        my $dz=$cz;
        $dx-=$boxx*&GenUtil::nint($dx/$boxx); 
        $dy-=$boxy*&GenUtil::nint($dy/$boxy); 
        $dz-=$boxz*&GenUtil::nint($dz/$boxz); 
        $atom->[$i]->{xcoor}+=$dx-$cx;
        $atom->[$i]->{ycoor}+=$dy-$cy;
        $atom->[$i]->{zcoor}+=$dz-$cz;
      }
    }
  } elsif ($by eq "reimage") {
    foreach my $c ( @{$self->{chain}} ) {
      my $atom=$c->{atom};
      for (my $i=0; $i<=$#{$atom}; $i++) {
        my $cx=$atom->[$i]->{xcoor};
        my $cy=$atom->[$i]->{ycoor};
        my $cz=$atom->[$i]->{zcoor};
        my $dx=($scx-$cx);
        my $dy=($scy-$cy);
        my $dz=($scz-$cz);
        $dx=$boxx*&GenUtil::nint($dx/$boxx); 
        $dy=$boxy*&GenUtil::nint($dy/$boxy); 
        $dz=$boxz*&GenUtil::nint($dz/$boxz); 
        $atom->[$i]->{xcoor}+=$dx;
        $atom->[$i]->{ycoor}+=$dy;
        $atom->[$i]->{zcoor}+=$dz;
      }
    }
  }   
}

## method: center()
## centers the current structure by subtracting the center 
## of mass from all coordinates

sub center {
  my $self=shift;
  my $chain=shift;
  
  my ($cx,$cy,$cz)=$self->centerOfMass($chain);

  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $atom->[$i]->{xcoor}-=$cx;
      $atom->[$i]->{ycoor}-=$cy;
      $atom->[$i]->{zcoor}-=$cz;
    }
  }
}

## method: orient()
## aligns principal axes to x-y-z (largest moment of interia)

sub orient {
  my $self=shift;

  my ($cx,$cy,$cz)=$self->centerOfMass();
 
  my ($xx,$xy,$xz,$yy,$yz,$zz)=(0,0,0,0,0,0); 
  foreach my $c ( @{$self->{chain}} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $atom->[$i]->{xcoor}-=$cx;
      $atom->[$i]->{ycoor}-=$cy;
      $atom->[$i]->{zcoor}-=$cz;
  
      my $tx=$atom->[$i]->{xcoor};
      my $ty=$atom->[$i]->{ycoor};
      my $tz=$atom->[$i]->{zcoor};

      $xx+=$tx*$tx;
      $xy+=$tx*$ty;
      $xz+=$tx*$tz;
      $yy+=$ty*$ty;
      $yz+=$ty*$tz;
      $zz+=$tz*$tz;
    }
  }

  my $w=();
  $w->[1]=$zz;
  $w->[2]=$yz;
  $w->[3]=$xz;
  $w->[4]=$yy;
  $w->[5]=$xy;
  $w->[6]=$xx;

  my ($u,$ev)=&Analyze::_diagq(3,3,$w);

  # generate transpose matrix
  my $det;
  for (my $i=1; $i<=3; $i++) {
    $det=$u->[$i];
    $u->[$i]=$u->[$i+6];
    $u->[$i+6]=$det;
  }
  for (my $i=1; $i<=3; $i++) {
    my $ipt=($i-1)*3;
    $det=$u->[$ipt+1];
    $u->[$ipt+1]=$u->[$ipt+3];
    $u->[$ipt+3]=$det;
    if ($u->[$ipt+$i]<0.0) {
      for (my $j=1; $j<=3; $j++) {
        $ipt++;
        $u->[$ipt]=-$u->[$ipt];
      }
    }
  }

  $det=$u->[1]*($u->[5]*$u->[9]-$u->[6]*$u->[8])+$u->[2]*($u->[6]*$u->[7]-$u->[4]*$u->[9])+$u->[3]*($u->[4]*$u->[8]-$u->[5]*$u->[7]);
  if ($det<0.0) {
    $u->[7]=-$u->[7];
    $u->[8]=-$u->[8];
    $u->[9]=-$u->[9];
    $det=-$det;
  }
  if (abs($det-1.0) > 0.0001) {
    printf STDERR "Molecule::orient: rotation matrix not unitary\n";
  }

#  for (my $i=1; $i<=9; $i++) {
#    printf STDERR "%lf ",$u->[$i];
#  }
#  printf STDERR "\n";

  my ($r,$phi)=&Analyze::_fndrot($u);

#  for (my $i=1; $i<=9; $i++) {
#    printf STDERR "%lf ",$r->[$i];
#  }
#  printf STDERR "\n";
#  printf STDERR "%lf\n",$phi;
  
  foreach my $c ( @{$self->{chain}} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      my $tx=$atom->[$i]->{xcoor};
      my $ty=$atom->[$i]->{ycoor};
      my $tz=$atom->[$i]->{zcoor};

      $atom->[$i]->{xcoor}=$u->[1]*$tx+$u->[2]*$ty+$u->[3]*$tz;
      $atom->[$i]->{ycoor}=$u->[4]*$tx+$u->[5]*$ty+$u->[6]*$tz;
      $atom->[$i]->{zcoor}=$u->[7]*$tx+$u->[8]*$ty+$u->[9]*$tz;
    }
  }

  foreach my $c ( @{$self->{chain}} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $atom->[$i]->{xcoor}+=$cx;
      $atom->[$i]->{ycoor}+=$cy;
      $atom->[$i]->{zcoor}+=$cz;
    }
  }
}

## method: setaux1()
## set value of aux1 (occupancy column)
 
sub setaux1 {
  my $self=shift;
  my $val=shift;
  my $valid=shift;

  foreach my $c ( @{$self->{chain}} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $atom->[$i]->{aux1}=$val if (!defined $valid || !$valid || $atom->[$i]->{valid});
    }
  }
}

## method: setaux2()
## set value of aux2 (occupancy column)
 
sub setaux2 {
  my $self=shift;
  my $val=shift;
  my $valid=shift;

  foreach my $c ( @{$self->{chain}} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $atom->[$i]->{aux2}=$val if (!defined $valid || !$valid || $atom->[$i]->{valid});
    }
  }
}

## method: move(dx,dy,dz,chain)
## shifts the current structure by the given distances
## in x/y/z direction

sub move {
  my $self=shift;
  my $dx=shift;
  my $dy=shift;
  my $dz=shift;
  my $chain=shift;
  
  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $atom->[$i]->{xcoor}+=$dx;
      $atom->[$i]->{ycoor}+=$dy;
      $atom->[$i]->{zcoor}+=$dz;
    }
  }
}

## method: add (mol)
## adds $cmp (a vector in the form of a pdb object) to mol
## (similar to but simpler than "displace")

sub add {
  my $self=shift;
  my $cmp=shift;
  
  for (my $ic=0; $ic<=$#{$self->{chain}}; $ic++) {
    my $atom=$self->{chain}->[$ic]->{atom};
    my $catom=$cmp->{chain}->[$ic]->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
        die "molecule mismatch" unless ($atom->[$i]->{atomname} eq $catom->[$i]->{atomname});
	$atom->[$i]->{xcoor}+=$catom->[$i]->{xcoor};
	$atom->[$i]->{ycoor}+=$catom->[$i]->{ycoor};
	$atom->[$i]->{zcoor}+=$catom->[$i]->{zcoor};
    }
  }
}

## method: subtract(mol)
## calculates the difference vector with respect to
## given molecule

sub subtract {
  my $self=shift;
  my $cmp=shift;
  
  for (my $ic=0; $ic<=$#{$self->{chain}}; $ic++) {
    my $atom=$self->{chain}->[$ic]->{atom};
    my $catom=$cmp->{chain}->[$ic]->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
        die "molecule mismatch" unless ($atom->[$i]->{atomname} eq $catom->[$i]->{atomname});
	$atom->[$i]->{xcoor}-=$catom->[$i]->{xcoor};
	$atom->[$i]->{ycoor}-=$catom->[$i]->{ycoor};
	$atom->[$i]->{zcoor}-=$catom->[$i]->{zcoor};
    }
  }
}

## method: scale(factor,[chain])
## scales all coordinates by given factor

sub scale {
  my $self=shift;
  my $factor=shift;
  my $chain=shift;
  
  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      $atom->[$i]->{xcoor}*=$factor;
      $atom->[$i]->{ycoor}*=$factor;
      $atom->[$i]->{zcoor}*=$factor;
    }
  }
}

## method: displace(vectorarr,amplitude,massweighting)
## displaces atoms along given vectors with amplitude
## mass weighting is done if flag is set

sub displace {
  my $self=shift;
  my $vectors=shift;
  my $amplitude=shift;
  my $weight=shift;
   
  my $n=0; 
  foreach my $c ( @{$self->activeChains()} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      my $mweight=1.0;
    
      if (defined $weight && $weight) {
       my $mass=1.0;
       my $aname=$atom->[$i]->{atomname};
       $mass=1.008 if ($aname=~/^[0-9]*H/);
       $mass=12.011 if ($aname=~/^[0-9]*C/);
       $mass=15.999 if ($aname=~/^[0-9]*O/);
       $mass=14.007 if ($aname=~/^[0-9]*N/);
       $mass=32.06 if ($aname=~/^[0-9]*S/);
       $mass=30.974 if ($aname=~/^[0-9]*P/);
       $mass=40.08 if ($aname=~/^[0-9]*CAL/);
       $mass=65.37 if ($aname=~/^[0-9]*ZN/);
       $mass=55.847 if ($aname=~/^[0-9]*FE/);
       $mass=24.305 if ($aname=~/^[0-9]*MG/);
       $mass=22.98977 if ($aname=~/^[0-9]*SOD/);
       $mass=22.98977 if ($aname=~/^[0-9]*NA/);
       $mass=39.102 if ($aname=~/^[0-9]*POT/);
       $mass=35.45 if ($aname=~/^[0-9]*CL/);

       $mweight=1.0/sqrt($mass);
      }

      $atom->[$i]->{xcoor}+=$vectors->[$n]->{x}*$amplitude*$mweight;
      $atom->[$i]->{ycoor}+=$vectors->[$n]->{y}*$amplitude*$mweight;
      $atom->[$i]->{zcoor}+=$vectors->[$n]->{z}*$amplitude*$mweight;
      $n++;
    }
  }
  
}

## method: rotate(m11,m12,m13,m21,m22,m23,m31,m32,m33,chain)
## rotates the current structure according to
## the given rotation matrix

sub rotate {
  my $self=shift;
  my $m11=shift;
  my $m12=shift;
  my $m13=shift;
  my $m21=shift;
  my $m22=shift;
  my $m23=shift;
  my $m31=shift;
  my $m32=shift;
  my $m33=shift;
  my $chain=shift;
  
  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      my $tx=$atom->[$i]->{xcoor};
      my $ty=$atom->[$i]->{ycoor};
      my $tz=$atom->[$i]->{zcoor};
      $atom->[$i]->{xcoor}=$m11*$tx+$m12*$ty+$m13*$tz;
      $atom->[$i]->{ycoor}=$m21*$tx+$m22*$ty+$m23*$tz;
      $atom->[$i]->{zcoor}=$m31*$tx+$m32*$ty+$m33*$tz;
    }
  }
}

sub rotatex {
  my $self=shift;
  my $phi=shift;
  $phi=$phi/180.0*$GenUtil::pi;
  $self->rotate(1,0,0,0,cos($phi),sin($phi),0,-sin($phi),cos($phi));
}

sub rotatey {
  my $self=shift;
  my $phi=shift;
  $phi=$phi/180.0*$GenUtil::pi;
  $self->rotate(cos($phi),0,-sin($phi),0,1,0,sin($phi),0,cos($phi));
}
 
sub rotatez {
  my $self=shift;
  my $phi=shift;
  $phi=$phi/180.0*$GenUtil::pi;
  $self->rotate(cos($phi),sin($phi),0,-sin($phi),cos($phi),0,0,0,1);
}

## method: matrixOperation(A,B,C,D,E,F,G,H,I [,chain])
## apply the matrix operator:
##  A  B  C
##  D  E  F
##  G  H  I

sub matrixOperation {
  my $self=shift;
  my $A=shift;
  my $B=shift;
  my $C=shift;
  my $D=shift;
  my $E=shift;
  my $F=shift;
  my $G=shift;
  my $H=shift;
  my $I=shift;
  my $chain=shift;

  foreach my $c ( @{$self->activeChains($chain)} ) {
    my $atom=$c->{atom};
    for (my $i=0; $i<=$#{$atom}; $i++) {
      my $x=$atom->[$i]->{xcoor};
      my $y=$atom->[$i]->{ycoor};
      my $z=$atom->[$i]->{zcoor};

      $atom->[$i]->{xcoor}=($A*$x)+($B*$y)+($C*$z);
      $atom->[$i]->{ycoor}=($D*$x)+($E*$y)+($F*$z);
      $atom->[$i]->{zcoor}=($G*$x)+($H*$y)+($I*$z);
    }
  }
}

## method: resetResidueName(resname,resnum[,chain])
## sets the residue name for a given residue number

sub resetResidueName {
  my $self=shift;
  my $resname=shift;
  my $resnum=shift;
  my $chain=shift;

  my $c=$self->activeChains($chain)->[0];
  my $r=$self->getResidueInChain($resnum,$c);
  if (defined $r) {
    for (my $ai=$r->{start}; $ai<=$r->{end}; $ai++) {
      $c->{atom}->[$ai]->{resname}=$resname;
    }
    $r->{name}=$resname;
    $r->{lockname}=1;
  }
} 

## method: fromSequence(inx,atomname,sequence)
## generates a molecule structure from a list of amino acids 

sub fromSequence {
  my $self=shift;
  my $inx=shift;
  my $aname=shift;
  my $sequence=shift;

  my $natom;
  my $chainid;

  my $natom=1;
  my $atinx;
  if ($inx=~/([A-Z])([0-9]+)/) {
    $chainid=$1;
    $atinx=$2;
  } else {
    $chainid="";
    $atinx=$inx;
  }

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{ssbond}=();

  my $chainrec=$self->_newChain($chainid);

  for (my $i=0; $i<length($sequence); $i++) {
    my $char=uc substr($sequence,$i,1);
    die "unknown residue name" if (!exists $Sequence::_seqlong{$char});
    my $resname=$Sequence::_seqlong{$char};
    my $resnum=$atinx++;

    my $resrec={ name => $resname, num => $resnum, start => $natom-1, 
                 chain=>$chainid, valid=>1 };
    push (@{$chainrec->{res}},$resrec);
   
    my $asc={ atominx => $natom++, atomname => $aname, 
	      resname => $resname, resnum  => $resnum, chain=> $chainid,
              xcoor   => 0.0, ycoor => 0.0, zcoor => 0.0 , valid=>1};
    push (@{$chainrec->{atom}},$asc);

    $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
  }

  $self->_coorCache();
}



## method: rebuildFromSICHO(sequence,chain[,fraglistoption,refPDB])
## rebuilds an all-atom structure from a SICHO chain and
## sequence file. For loop modeling loop/fragment residues
## and a protein template in PDB form may be given as
## method: fromSICHO(sequence,chain)
## generates a molecule structure from a SICHO chain file

sub fromSICHO {
  my $self=shift;
  my $seq=shift;
  my $sicho=shift;

  my $nsicho=$#{$sicho->{sidechain}};
  my $nseq=$#{$seq->{sequence}};

  my $coff=0;
  $coff=1 if ($nsicho == $nseq+2);

  my $natom=1;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{ssbond}=();

  my $chainrec=$self->_newChain("");

  for (my $i=0; $i<=$nseq; $i++) {
    my $resname=$seq->{sequence}->[$i]->{residue};
    my $resnum=$seq->{sequence}->[$i]->{index};

    my $resrec={ name => $resname, num => $resnum, start => $natom-1, chain=>"", valid=>1 };
    push (@{$chainrec->{res}},$resrec);
   
    my ($x,$y,$z)=$sicho->fromProjection($i+$coff);
    my $asc={ atominx => $natom++, atomname => "SC", 
	      resname => $resname, resnum  => $resnum, chain=> "",
              xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
    push (@{$chainrec->{atom}},$asc);

    if (defined $sicho->{ca}->[$i+$coff]) {
      my ($cax,$cay,$caz)=$sicho->fromProjection($i+$coff,1);
      my $aca={ atominx => $natom++, atomname => "CA", 
		resname => $resname, resnum  => $resnum, chain=> "",
		xcoor => $cax, ycoor => $cay, zcoor => $caz};
      push (@{$chainrec->{atom}},$aca);
    }

    $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
  }

  $self->_coorCache();
}

## method: fromXYZ(xyz,atomname)
## generates a molecule structure from a list of XYZ coordinates

sub fromXYZ {
  my $self=shift;
  my $xyz=shift;
  my $aname=shift;

  my $natom=1;

  $self->{chain}=();
  $self->{chainlookup}={};
  $self->{defchain}=undef;
  $self->{ssbond}=();

  my $chainrec=$self->_newChain("");

  for (my $i=0; $i<=$#{$xyz}; $i++) {
    my $resname=$aname;
    my $resnum=$natom;

    my $resrec={ name => $resname, num => $resnum, start => $natom-1, chain=>"", valid=>1 };
    push (@{$chainrec->{res}},$resrec);
   
    my $asc={ atominx => $natom++, atomname => $aname, 
	      resname => $resname, resnum  => $resnum, chain=> "",
              xcoor   => $xyz->[$i]->{x}, ycoor => $xyz->[$i]->{y}, zcoor => $xyz->[$i]->{z}, valid=>1 };
    push (@{$chainrec->{atom}},$asc);

    $chainrec->{res}->[$#{$chainrec->{res}}]->{end}=$#{$chainrec->{atom}};
  }

  $self->_coorCache();
}



## method: rebuildFromSICHO(sequence,chain[,fraglistoption,refPDB])
## rebuilds an all-atom structure from a SICHO chain and
## sequence file. For loop modeling loop/fragment residues
## and a protein template in PDB form may be given as
## extra arguments

sub rebuildFromSICHO {
  my $self=shift;
  my $seq=shift;
  my $sicho=shift;
  my $fraglistopt=shift;
  my $refpdb=shift;
  my $fixca=shift;

  my $nsicho=$#{$sicho->{sidechain}};
  my $nseq=$#{$seq->{sequence}};

  my $option="";

  die "empty SICHO chain" if ($nsicho<0);
  die "missing sequence" if ($nseq<0);

  if (defined $fraglistopt) {  
    die "Need reference PDB file" if (!defined $refpdb);
    die "Chain ID in fragment list is not supported for all-atom rebuilding"
      if ($fraglistopt =~ /[A-Za-z]/);
    $option.=" -l ".$fraglistopt;
    $option.=" ".$refpdb;
  }

  if (defined $fixca && $fixca) {
    $option.=" -fixca";
  }

  my $rebuildbin=&GenUtil::findExecutable("rebuild");
  die "cannot find rebuild executable"
    if (!defined $rebuildbin);

  local (*READ,*WRITE);
  my $pid=open2(*READ,*WRITE,"$rebuildbin $option");

  die "cannot open2" if (!defined $pid || $pid<0);

  my $coff=0;
  $coff=1 if ($nsicho == $nseq+2);

  for (my $i=0; $i<=$nseq; $i++) {
    my ($x,$y,$z)=$sicho->fromProjection($i+$coff);
    my $str=sprintf("%d %s %f %f %f",
		    $seq->{sequence}->[$i]->{index},$seq->{sequence}->[$i]->{residue},
		    $x,$y,$z);
    if (defined $sicho->{ca}->[$i+$coff]) {
      my ($cax,$cay,$caz)=$sicho->fromProjection($i+$coff,1);
      $str.=sprintf(" %f %f %f",$cax,$cay,$caz);
    }
    print WRITE $str,"\n";
  }
  close WRITE;

  $self->readPDB(\*READ); close READ;
  waitpid($pid,0);
}

sub listSegments {
  my $self=shift;

  my $segs=$self->getSegNames();
  foreach my $s ( @{$segs}) {
    printf "%s\n",$s->{name};
  }
}  

sub info {
  my $self=shift;
 
  printf "%d chains total\n",$#{$self->activeChains()}+1;
  foreach my $c ( @{$self->activeChains()} ) {
    printf "== chain %1s\n",$c->{id};
    if ($#{$c->{res}}>0) {
      my $lastnum=-99999;
      my $first=$c->{res}->[0];
      for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
	my $r=$c->{res}->[$ir];
	my $num=$r->{num};
        if ($num!=$lastnum+1 && $lastnum>=0) {
          if ($lastnum==$first->{num}) {
           printf "%d\n",$first->{num};
          } else {
           printf "%d-%d ",$first->{num},$lastnum;
          }
          $first=$r; 
        }
	$lastnum=$num;
      }
      if ($lastnum==$first->{num}) {
        printf "%d\n",$first->{num};
      } else {
        printf "%d-%d\n",$first->{num},$lastnum;
      }

      my $minx=999999999;
      my $miny=999999999;
      my $minz=999999999;
 
      my $maxx=-999999999;
      my $maxy=-999999999;
      my $maxz=-999999999;

      foreach my $a ( @{$c->{atom}}) {
        $minx=$a->{xcoor} if ($a->{xcoor}<$minx);
        $miny=$a->{ycoor} if ($a->{ycoor}<$miny);
        $minz=$a->{zcoor} if ($a->{zcoor}<$minz);

        $maxx=$a->{xcoor} if ($a->{xcoor}>$maxx);
        $maxy=$a->{ycoor} if ($a->{ycoor}>$maxy);
        $maxz=$a->{zcoor} if ($a->{zcoor}>$maxz);
      }
      printf "xrange %lf %lf : %lf\n",$minx,$maxx,$maxx-$minx;
      printf "yrange %lf %lf : %lf\n",$miny,$maxy,$maxy-$miny;
      printf "zrange %lf %lf : %lf\n",$minz,$maxz,$maxz-$minz;
    } 
  }
}

## method: completeWater()
## adds water hydrogens

sub completeWater {
  my $self=shift;

  srand();

  my $n={};
  $n->{chain}=();
  $n->{chainlookup}={};
  $n->{defchain}=undef;
  $n->{segmentlist}=undef;
  $n->{ssbond}=();

  bless $n;

  foreach my $s ( @{$self->{ssbond}} ) {
    my $srec={};
    %{$srec}=%{$s};
    push (@{$n->{ssbond}},$srec);
  }

  foreach my $c ( @{$self->activeChains()} ) {
    my $nc=undef;
    
    foreach my $r (@{$c->{res}}) {
      $nc=$n->_newChain($c->{id}) if (!defined $nc);
      
      my $rrec={};
      %{$rrec}=%{$r};
      $rrec->{start}=$#{$nc->{atom}}+1;
      push(@{$nc->{res}},$rrec);

      for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	my $arec={};
	%{$arec}=%{$c->{atom}->[$ia]};
	$arec->{atomname}="OH2" 
	  if ($arec->{atomname} eq "O" && ($arec->{resname} eq "TIP3" || $arec->{resname} eq "TIP4" || $arec->{resname} eq "HOH"));
	push(@{$nc->{atom}},$arec);
      }
      if (($r->{name} eq "TIP3" || $r->{name} eq "TIP4" || $r->{name} eq "HOH") && $r->{start} == $r->{end}) {
	my $ox=$c->{atom}->[$r->{start}]->{xcoor};
	my $oy=$c->{atom}->[$r->{start}]->{ycoor};
	my $oz=$c->{atom}->[$r->{start}]->{zcoor};
	
	my $arec={};
	%{$arec}=%{$c->{atom}->[$r->{start}]};
	$arec->{atomname}="H1";
	
	my $phi=rand($GenUtil::pi);
	my $theta=rand(2.0*$GenUtil::pi);
	$arec->{xcoor}=$ox+0.9572*sin($phi)*cos($theta);
	$arec->{ycoor}=$oy+0.9572*sin($phi)*sin($theta);
	$arec->{zcoor}=$oz+0.9572*cos($phi);
	push(@{$nc->{atom}},$arec);
	
	$arec={};
	%{$arec}=%{$c->{atom}->[$r->{start}]};
	$arec->{atomname}="H2";
	$arec->{xcoor}=$ox+0.9572*sin($phi)*cos($theta+1.8239);
	$arec->{ycoor}=$oy+0.9572*sin($phi)*sin($theta+1.8239);
	$arec->{zcoor}=$oz+0.9572*cos($phi);
	push(@{$nc->{atom}},$arec);
      }
      $rrec->{end}=$#{$nc->{atom}};
    }
  }

  $n->_coorCache();

  return $n;
}


## method: completeResidue()
## completes missing peptide backbone
## and side chain atoms 

sub completeResidue {
  my $self=shift;
  my $sicho=shift;
  my $fixca=shift;

  $fixca=1 if (!defined $fixca);

  my $defchain=$self->{defchain};

  $self->translate("GENERIC");

  foreach my $c ( @{$self->activeChains()} ) {
    my @sicholist;
    my @backlist;
    my @scwrllist;
    my $sichonogly=0;
    my $scwrlnogly=0;
    my %resseg;
#    if ($c->{id} ne "+") {
    foreach my $r ( @{$c->{res}} ) {
      $resseg{$r->{num}}=$r->{seg};
      if ($r->{valid} && 
	  $r->{name}=~/^(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HSD|HSE|HSP|HIS|HID|HIE|HIP|HSP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|CYX)$/) {
	my %have;
	for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
	  $have{$c->{atom}->[$i]->{atomname}}=1;
	}

	my $rn=$r->{name};
	my $side=(&_cside($rn,\%have,qw( ALA CB )) ||
	    &_cside($rn,\%have,qw( CYS CB SG )) ||
	    &_cside($rn,\%have,qw( LEU CB CG CD1 CD2)) ||
	    &_cside($rn,\%have,qw( ILE CB CG2 CD1 )) ||
	    &_cside($rn,\%have,qw( ARG CB CG CD NE CZ NH1 NH2)) ||
	    &_cside($rn,\%have,qw( MET CB CG SD CE )) ||
	    &_cside($rn,\%have,qw( LYS CB CG CD CE NZ )) ||
	    &_cside($rn,\%have,qw( GLY CA )) ||
	    &_cside($rn,\%have,qw( SER CB OG )) ||
	    &_cside($rn,\%have,qw( TRP CB CG CD2 CE2 CE3 CD1 NE1 CZ2 CZ3 CH2 )) ||
	    &_cside($rn,\%have,qw( PRO CA CD CB CG)) ||
	    &_cside($rn,\%have,qw( HIS CB CG ND1 CD2 NE2 CE1)) ||
	    &_cside($rn,\%have,qw( HSE CB CG ND1 CD2 NE2 CE1)) ||
	    &_cside($rn,\%have,qw( HSD CB CG ND1 CD2 NE2 CE1)) ||
	    &_cside($rn,\%have,qw( HSP CB CG ND1 CD2 NE2 CE1)) ||
	    &_cside($rn,\%have,qw( GLU CB CG CD OE1 OE2)) || 
	    &_cside($rn,\%have,qw( GLN CB CG CD OE1 NE2)) || 
	    &_cside($rn,\%have,qw( ASP CB CG OD1 OD2)) || 
	    &_cside($rn,\%have,qw( ASN CB CG OD1 ND2)) || 
	    &_cside($rn,\%have,qw( TYR CB CG CD1 CE1 CD2 CE2 CZ OH)) || 
	    &_cside($rn,\%have,qw( PHE CB CG CD1 CE1 CD2 CE2 CZ)) || 
	    &_cside($rn,\%have,qw( THR CB OG1 CG2)) ||
	    &_cside($rn,\%have,qw( VAL CB CG1 CG2)) );
	my $back=(defined $have{CA} && defined $have{C} &&
		  defined $have{N} && 
		  (defined $have{O} || defined $have{OT1}));

	if (!$back) {
	  push(@backlist,$r->{num});
	}
	if (!$side && $rn ne "GLY") { #|| $rn eq "GLY") {
	  if (defined $have{CB}) { #|| $rn eq "GLY") {
	    push(@sicholist,$r->{num});
	    $sichonogly=1 if ($rn ne "GLY");
	  } 
	  if (!defined $have{CB}) {
	    push(@scwrllist,$r->{num});
	    $scwrlnogly=1 if ($rn ne "GLY");
	  } 
	}
      }
    }

#    printf STDERR "backlist %s\n",join(" ",@backlist);
#    printf STDERR "sicholist %s\n",join(" ",@sicholist);
#    printf STDERR "scwrllist %s\n",join(" ",@scwrllist);
   
#    printf STDERR "sichonogly: %d\n",$sichonogly; 
#    printf STDERR "scwrlnogly: %d\n",$scwrlnogly; 

    if ($#backlist>=0 && join("::",@backlist) ne join("::",@sicholist)) {
      my $option="-backonly ";  
      $option.="-fixca " if ($fixca);
      my $refpdb;
      my $fraglist;
      if ($#backlist+1 < $#{$c->{res}}) {
	$self->selectChain($c->{id});
	$refpdb="ref-".int(rand(1000000)).".pdb";
	$self->writePDB($refpdb,ssbond=>0);
	$fraglist=&GenUtil::fragListFromArray(\@backlist);
	$option.="-l ".&GenUtil::fragOptionFromList($fraglist)." ".$refpdb;
      }

      my $rebuildbin=&GenUtil::findExecutable("rebuild");
      die "cannot find rebuild executable"
	if (!defined $rebuildbin);

#      printf STDERR "running 1 $rebuildbin $option\n";

      local (*READ,*WRITE);
      my $pid=open2(*READ,*WRITE,"$rebuildbin $option");

      die "cannot open2" if (!defined $pid || $pid<0);

      foreach my $a ( @{$c->{atom}} ) {
	if ($a->{atomname} eq "CA") {
	  printf WRITE "%d %s %f %f %f %f %f %f\n",
	  $a->{resnum},$a->{resname},$a->{xcoor},$a->{ycoor},$a->{zcoor},
	  $a->{xcoor},$a->{ycoor},$a->{zcoor};
	}
      }

      close WRITE;

      my $tmol=Molecule::new();
      $tmol->readPDB(\*READ); close READ;
      $tmol->setValidResidues($fraglist);
      $tmol->setChain($c->{id});

      foreach my $tr ( @{$tmol->{chain}->[0]->{res}} ) {
        $tr->{seg}=$resseg{$tr->{num}};
        for (my $ia=$tr->{start}; $ia<=$tr->{end}; $ia++) {
          $tmol->{chain}->[0]->{atom}->[$ia]->{seg}=$resseg{$tr->{num}};
        }
      }

      $self->merge($tmol);
      &GenUtil::remove($refpdb);
      waitpid($pid,0);
    }
    if ($#sicholist>=0 && $sichonogly && (!defined $sicho || $sicho)) {
      my $option="";
      $option.="-fixca " if ($fixca);
      my $refpdb;
      my $fraglist;
      if ($#sicholist+1 < $#{$c->{res}}) {
	$self->selectChain($c->{id});
	$refpdb="ref-".int(rand(1000000)).".pdb";
	$self->writePDB($refpdb,ssbond=>0);
	$fraglist=&GenUtil::fragListFromArray(\@sicholist);
	$option.="-l ".&GenUtil::fragOptionFromList($fraglist)." ".$refpdb;
      }

      my $rebuildbin=&GenUtil::findExecutable("rebuild");
      die "cannot find rebuild executable"
	if (!defined $rebuildbin);

#     printf STDERR "running 2 $rebuildbin $option\n";

      local (*READ,*WRITE);
      my $pid=open2(*READ,*WRITE,"$rebuildbin $option");

      die "cannot open2" if (!defined $pid || $pid<0);

      foreach my $r ( @{$c->{res}} ) {
	my $ca=undef;
	my $cb=undef;
	for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	  if ($c->{atom}->[$ia]->{atomname} eq "CA") {
	    $ca=$c->{atom}->[$ia];
	  } elsif ($c->{atom}->[$ia]->{atomname} eq "CB") {
	    $cb=$c->{atom}->[$ia];
	  }
	}

	if (defined $cb) {
	  printf WRITE "%d %s %f %f %f %f %f %f\n",
	  $ca->{resnum},$ca->{resname},
	  $ca->{xcoor}+$cb->{xcoor},$ca->{ycoor}+$cb->{ycoor},
	  $ca->{zcoor}+$cb->{zcoor},
	  $ca->{xcoor},$ca->{ycoor},$ca->{zcoor};
	} else {
	  printf WRITE "%d %s %f %f %f %f %f %f\n",
	  $ca->{resnum},$ca->{resname},$ca->{xcoor},$ca->{ycoor},$ca->{zcoor},
	  $ca->{xcoor},$ca->{ycoor},$ca->{zcoor};
	}
      }
      close WRITE;

      my $tmol=Molecule::new();
      $tmol->readPDB(\*READ); close READ;
      $tmol->setValidResidues($fraglist);
      $tmol->setChain($c->{id});

      foreach my $tr ( @{$tmol->{chain}->[0]->{res}} ) {
        $tr->{seg}=$resseg{$tr->{num}};
        for (my $ia=$tr->{start}; $ia<=$tr->{end}; $ia++) {
          $tmol->{chain}->[0]->{atom}->[$ia]->{seg}=$resseg{$tr->{num}};
        }
      }

      $self->merge($tmol);
      &GenUtil::remove($refpdb);
      waitpid($pid,0);
    } 

    if ($#scwrllist>=0 && $scwrlnogly) {
      my $refseq="seq-".int(rand(1000000));
      my $option="-s $refseq";
      my $fraglist;
# Needed for scwrl3
      my $initpdb;
      my $outpdb;
      my $tmpdat;
      my $pid;

      my %haves;
      foreach my $s ( @scwrllist ) {
	$haves{$s}=1;
      }
      my $seqstr="";
      foreach my $r (@{$c->{res}}) {
	if (defined $haves{$r->{num}}) {
	  $seqstr.=$Sequence::_seqabbrev{$r->{name}};
	} else {
	  $seqstr.="x";
	}
      }
      
      open OUT,">$refseq";
      print OUT $seqstr;
      close OUT;
### Added Mayako Michino procedure to use scwrl 3.x or swrl 2.x      
      my $scwrlbin=&GenUtil::findExecutable("scwrl");
      if(defined $scwrlbin) {

#	  printf STDERR "running $scwrlbin $option\n";

	  local (*READ,*WRITE);
	  my $pid=open2(*READ,*WRITE,"$scwrlbin $option");

	  die "cannot open2" if (!defined $pid || $pid<0);

	  $self->selectChain($c->{id});
	  $self->writePDB(\*WRITE,ssbond=>0);
	  close WRITE;
      } else {
	  $scwrlbin=&GenUtil::findExecutable("scwrl3");
	  die "cannot find scwrl executable"
	      if (!defined $scwrlbin);
	  $initpdb = "init.pdb";
	  $outpdb = "out.pdb";
	  $tmpdat = "tmp.dat";

	  $self->selectChain($c->{id});
	  $self->writePDB($initpdb,ssbond=>0);
	  $option .= " -i $initpdb -o $outpdb > $tmpdat";

#	  printf STDERR "running $scwrlbin $option\n";
	  my $pid = system("$scwrlbin $option");

	  open READ,"<$outpdb"; 
      }

      my $tmol=Molecule::new();
      $tmol->readPDB(\*READ); close READ;
      $tmol->setValidResidues(&GenUtil::fragListFromArray(\@scwrllist));
      $tmol->setChain($c->{id});

      foreach my $tr ( @{$tmol->{chain}->[0]->{res}} ) {
        $tr->{seg}=$resseg{$tr->{num}};
        for (my $ia=$tr->{start}; $ia<=$tr->{end}; $ia++) {
          $tmol->{chain}->[0]->{atom}->[$ia]->{seg}=$resseg{$tr->{num}};
        }
      }

      $self->merge($tmol);
      &GenUtil::remove($refseq);
      if(defined $initpdb){
      &GenUtil::remove($initpdb);
      &GenUtil::remove($outpdb);
      &GenUtil::remove($tmpdat);
      }
      waitpid($pid,0);
    }
    #$self->writePDB("lat.$c->{id}.pdb");
  }
  $self->{defchain}=$defchain;
  return $self;
}  

sub _cside {
  my $name=shift;
  my $have=shift;
  my $cname=shift;
  my @atoms=@_;
  return 0 if ($name ne $cname);
  foreach my $a (@atoms) {
    return 0 
      if (!defined $have->{$a});
  }
  return 1;
}

## method: fixCOO()
## corrects the position of OT2/OXT if necessary

sub fixCOO {
  my $self=shift;
  my $nocterm=shift;
 
  foreach my $ch ( @{$self->activeChains()} ) {
    my $lastres=$ch->{res}->[$#{$ch->{res}}];

    my $ca;
    my $c;
    my $ot1;
    my $ot2;
    
    for (my $i=$lastres->{start}; $i<=$lastres->{end}; $i++) {
      my $aname=$ch->{atom}->[$i]->{atomname};

      $ch->{atom}->[$i]->{atomname}="O"
	if ($nocterm && $aname eq "OT1");
      $ch->{atom}->[$i]->{atomname}="OXX"
	if ($nocterm && $aname eq "OXT");

      if ($aname eq "CA") {
	$ca=$ch->{atom}->[$i];
      } elsif ($aname eq "C") {
	$c=$ch->{atom}->[$i];
      } elsif ($aname eq "O" || $aname eq "OT1") {
	$ot1=$ch->{atom}->[$i];
      } elsif ($aname eq "OXT" || $aname eq "OT2") {
	$ot2=$ch->{atom}->[$i];
      }
    }

    return if ($nocterm || !defined $ca || !defined $c || !defined $ot1);

    my $dot=0.0;
    if (defined $ot2) {
      my $dotx=$ot1->{xcoor}-$ot2->{xcoor};
      my $doty=$ot1->{ycoor}-$ot2->{ycoor};
      my $dotz=$ot1->{zcoor}-$ot2->{zcoor};
      $dot=sqrt($dotx*$dotx+$doty*$doty+$dotz*$dotz);
    } else {
      my $arec={};
      %{$arec}=%{$ch->{atom}->[$#{$ch->{atom}}]};
      $arec->{atominx}++;
      $arec->{atomname}="OXT";
      push(@{$ch->{atom}},$arec);
      $ot2=$arec;
    }

    if ($dot<1.5) {
      my $tcx=$c->{xcoor}-$ca->{xcoor};
      my $tcy=$c->{ycoor}-$ca->{ycoor};
      my $tcz=$c->{zcoor}-$ca->{zcoor};

      my $ltc=sqrt($tcx*$tcx+$tcy*$tcy+$tcz*$tcz);
      
      $tcx/=$ltc;
      $tcy/=$ltc;
      $tcz/=$ltc;
      
      my $t1x=$ot1->{xcoor}-$c->{xcoor};
      my $t1y=$ot1->{ycoor}-$c->{ycoor};
      my $t1z=$ot1->{zcoor}-$c->{zcoor};

      my $c1=($tcx*$t1x+$tcy*$t1y+$tcz*$t1z);
      
      $ot2->{xcoor}=$ot1->{xcoor}+2.0*($c1*$tcx-$t1x);
      $ot2->{ycoor}=$ot1->{ycoor}+2.0*($c1*$tcy-$t1y);
      $ot2->{zcoor}=$ot1->{zcoor}+2.0*($c1*$tcz-$t1z);
    }
  }
  $self->_coorCache();
}

## method: solvate([cutoff,[shape]])
## solvates a PDB from a pre-equilibrated water box
## returns error output from solvate program
sub solvate {
  my $self=shift;
  my $cutoff=shift;
  my $shape=shift;
  my $fraglist=shift;
  my $solvcut=shift;
  my $tip4p=shift;
  my $center=shift;
  my $splitseg=shift;

  #print STDERR "cutoff= $cutoff shape= $shape tip4p= $tip4p center= $center\n";
  $cutoff=9.0 if (!defined $cutoff);
  my $option=" -cutoff $cutoff ";
  $option.=" -box $ENV{MMTSBDIR}/data/water_tip4p.pdb " if ($tip4p==1);
  $option.=" -box $ENV{MMTSBDIR}/data/water.pdb " if ($tip4p!=1);
  $option.="-$shape " if (defined $shape);
  #$option.="-tip4p " if (defined $tip4p && $tip4p);
  $option.="-solvcut $solvcut " if (defined $solvcut);
  $option.="-tip4p " if (defined $tip4p && $tip4p);
  $option.="-center " if (defined $center && $center);
  $option.="-nocenter " if (defined $center && !$center);
  my $solvatebin=&GenUtil::findExecutable("solvate");
  die "cannot find solvate executable"
    if (!defined $solvatebin);

#  printf STDERR "running $solvatebin $option -\n";

  local (*READ,*WRITE,*ERR);
  my $pid=open3(*WRITE,*READ,*ERR,"$solvatebin $option -");
  
  die "cannot open3" if (!defined $pid || $pid <0);

  my $ss=$self->getSSBonds();
  $self->writePDB(\*WRITE,ssbond=>0);
  close WRITE;

  $self->readPDB(\*READ,splitseg=>$splitseg); close READ;
  
  my $err="";
  while (<ERR>) {
    $err.=$_;
  }
  close ERR;

  waitpid($pid,0);

  if ($shape eq "octahedron") {
    $self->rotate(2.0/3.0,-1.0/3.0,-2.0/3.0,
                 -1.0/3.0,2.0/3.0,-2.0/3.0, 
                  2.0/3.0,2.0/3.0,1.0/3.0); 
  } 
  $self->setSSBonds($ss);

  foreach my $c ( @{$self->{chain}} ) {
    foreach my $a ( @{$c->{atom}} ) {
      $a->{resname}=~s/HOH/TIP3/;
    }
    foreach my $r ( @{$c->{res}} ) {
      $r->{name}=~s/HOH/TIP3/;
    }
  }

  return $err;
}

## method: replaceIons(list)
## replaces water molecules randomly with ions

sub replaceIons {
  my $self=shift;
  my $ions=shift;

  return if (!defined $ions || $#{$ions}<0);

  my @wat=();
  foreach my $c ( @{$self->{chain}} ) {
    foreach my $r ( @{$c->{res}} ) {
       if ($r->{name} eq "TIP3" || $r->{name} eq "TIP4" || $r->{name} eq "HOH") {
         $r->{atom}=$c->{atom}->[$r->{start}];
         push(@wat,$r);
       }
    }
  }

  my $ainx=1; 
  my $rinx=1;
  my $chain="+";
  my $chainrec=$self->{chainlookup}->{$chain};
  if (defined $chainrec) {
    $ainx=$chainrec->{atom}->[$#{$chainrec->{atom}}]->{atominx}+1;
    $rinx=$chainrec->{res}->[$#{$chainrec->{res}}]->{num}+1;
  } else {
    $chainrec=$self->_newChain($chain);
  }

  foreach my $i ( @{$ions} ) {
    for (my $ii=0; $ii<$i->{num}; $ii++) {
      my $iw;
      my $rw;

      do {
        $iw=int(rand($#wat+1));
        $rw=$wat[$iw];
      } while ($rw->{valid}==0);
  
      $rw->{valid}=0;

      my $pdbrec={};
	
      $pdbrec->{atominx}=$ainx++;
      $pdbrec->{atomname}=$i->{name};
      $pdbrec->{resname}=$i->{name};
      $pdbrec->{resnum}=$rinx;
      $pdbrec->{chain}=$chainrec->{id};
      $pdbrec->{xcoor}=$rw->{atom}->{xcoor};
      $pdbrec->{ycoor}=$rw->{atom}->{ycoor};
      $pdbrec->{zcoor}=$rw->{atom}->{zcoor};
      $pdbrec->{hyd}=0;
      $pdbrec->{aux1}=0.0;
      $pdbrec->{aux2}=0.0;
      $pdbrec->{seg}="";
      $pdbrec->{valid}=1;
      
      push (@{$chainrec->{atom}}, $pdbrec);

      my $resrec={};

      $resrec->{name}=$i->{name};
      $resrec->{num}=$rinx++;
      $resrec->{chain}=$chainrec->{id};
      $resrec->{start}=$#{$chainrec->{atom}};
      $resrec->{end}=$resrec->{start};
      $resrec->{valid}=1;
      $resrec->{seg}="";
 
      push(@{$chainrec->{res}},$resrec);
    } 
  } 
}

## method: replaceWaterWithMolecule(mol, num)
## replaces water molecules randomly with small molecules

sub replaceWaterWithMolecule {
  my $self=shift;
  my $rmol=shift;
  my $cutoff=shift;
  my $rnum=shift;

  my @wat=();
  my @other=();
  foreach my $c ( @{$self->{chain}} ) {
    foreach my $r ( @{$c->{res}} ) {
       if ($r->{name} eq "TIP3" || $r->{name} eq "TIP4" || $r->{name} eq "HOH") {
         $r->{atom}=$c->{atom}->[$r->{start}];
         push(@wat,$r);
       } 
    }
    foreach my $a ( @{$c->{atom}} ) {
       if ($a->{resname} ne "TIP3" && $a->{resname} ne "TIP4" && $a->{resname} ne "HOH") {
         push(@other,$a);
       }
    }
  }

  my $ainx=1; 
  my $rinx=1;
  my $chain="X";
  my $chainrec=$self->{chainlookup}->{$chain};
  my $newchain=1;
  if (defined $chainrec) {
    $chain="Y";
    $chainrec=$self->{chainlookup}->{$chain};
    if (defined $chainrec) {
      $chain="Z";
      $chainrec=$self->{chainlookup}->{$chain};
      if (defined $chainrec) {
        $chain="+"; 
        $chainrec=$self->{chainlookup}->{$chain};
        if (defined $chainrec) {
          $ainx=$chainrec->{atom}->[$#{$chainrec->{atom}}]->{atominx}+1;
          $rinx=$chainrec->{res}->[$#{$chainrec->{res}}]->{num}+1;
          $newchain=0;
        } 
      }
    }
  }
  if ($newchain) {
    $chainrec=$self->_newChain($chain);
  }

  for (my $ii=0; $ii<$rnum; $ii++) {
    my $iw;
    my $rw;
    my $wx;
    my $wy;
    my $wz;

    my $clash=0;
    do {
     do {
      $iw=int(rand($#wat+1));
      $rw=$wat[$iw];
     } while ($rw->{valid}==0);

     $wx=$rw->{atom}->{xcoor};
     $wy=$rw->{atom}->{ycoor};
     $wz=$rw->{atom}->{zcoor};

     $clash=0; 
     foreach my $a ( @other ) {
       my $dx=$wx-$a->{xcoor};
       my $dy=$wy-$a->{ycoor};
       my $dz=$wz-$a->{zcoor};
       my $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);
       $clash=1 if ($d<$cutoff);
     } 
    } while ($clash);
     
    $rw->{valid}=0;

    $rmol->center();
    $rmol->rotatex(rand(360));
    $rmol->rotatez(rand(360));
    $rmol->rotatex(rand(360));

    foreach my $tw ( @wat ) {
      foreach my $ra ( @{$rmol->{chain}->[0]->{atom}} ) {
        my $dx=$wx+$ra->{xcoor}-$tw->{atom}->{xcoor};
        my $dy=$wy+$ra->{ycoor}-$tw->{atom}->{ycoor};
        my $dz=$wz+$ra->{zcoor}-$tw->{atom}->{zcoor};
        my $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);
        if ($d<1.5) {
          $tw->{valid}=0;
        }
      }
    } 

    my $start=$#{$chainrec->{atom}}+1;

    foreach my $ra ( @{$rmol->{chain}->[0]->{atom}} ) {
      my $pdbrec={};
      $pdbrec->{atominx}=$ainx++;
      $pdbrec->{atomname}=$ra->{atomname};
      $pdbrec->{resname}=$ra->{resname};
      $pdbrec->{resnum}=$rinx;
      $pdbrec->{chain}=$chainrec->{id};
      $pdbrec->{xcoor}=$wx+$ra->{xcoor};
      $pdbrec->{ycoor}=$wy+$ra->{ycoor};
      $pdbrec->{zcoor}=$wz+$ra->{zcoor};
      $pdbrec->{hyd}=0;
      $pdbrec->{aux1}=0.0;
      $pdbrec->{aux2}=0.0;
      $pdbrec->{seg}="";
      $pdbrec->{valid}=1;
      push (@{$chainrec->{atom}}, $pdbrec);
    }

    my $resrec={};

    $resrec->{name}=$rmol->{chain}->[0]->{res}->[0]->{name};
    $resrec->{num}=$rinx++;
    $resrec->{chain}=$chainrec->{id};
    $resrec->{start}=$start;
    $resrec->{end}=$#{$chainrec->{atom}};
    $resrec->{valid}=1;
    $resrec->{seg}="";
 
    push(@{$chainrec->{res}},$resrec);
  } 
}

## method: changeResName(list)
## changes the residue names according to the list

sub changeResName {
  my $self=shift;
  my $list=shift;

  if (defined $list && $list ne "") {
    foreach my $n ( split(/_/,$list)) {
      my ($newname,$rnum)=split(/:/,$n);
      my ($chainid,$num)=($rnum=~/([A-Za-z]*)([0-9]+)/);
      $self->resetResidueName($newname,$num,$chainid);
    }
  }
}

## method: fixHistidine(hsdlist,hselist,hsplist)
## changes the residue names to HSD and HSE respectively
## according to the lists of residues given as arguments

sub fixHistidine {
  my $self=shift;
  my $hsd=shift;
  my $hse=shift;
  my $hsp=shift;

  if (defined $hsd && $hsd ne "") {
    foreach my $n ( split(/:/,$hsd)) {
      my ($chainid,$num)=($n=~/([A-Za-z]*)([0-9]+)/);
      $self->resetResidueName("HSD",$num,$chainid);
    }
  }

  if (defined $hse && $hse ne "") {
    foreach my $n ( split(/:/,$hse)) {
      my ($chainid,$num)=($n=~/([A-Za-z]*)([0-9]+)/);
      $self->resetResidueName("HSE",$num,$chainid);
    }
  }

  if (defined $hsp && $hsp ne "") {
    foreach my $n ( split(/:/,$hsp)) {
      my ($chainid,$num)=($n=~/([A-Za-z]*)([0-9]+)/);
      $self->resetResidueName("HSP",$num,$chainid);
    }
  }
}


## method: setChain(chainid)
## sets the chain ID for the current structure

sub setChain {
  my $self=shift;
  my $chainid=shift;
  my $setall=shift;

  $setall=0 if (!defined $setall);

  my $achains=$self->activeChains();
  my $maxinx=($setall)?$#{$achains}:0;

  for (my $ic=0; $ic<=$maxinx; $ic++) {
    my $c=$achains->[$ic];

    $self->{chainlookup}=undef;
    $self->{defchain}=undef;
    $self->{segmentlist}=undef;
  
    $c->{id}=$chainid;
  
    foreach my $a ( @{$c->{atom}} ) {
      $a->{chain}=$chainid;
    }

    foreach my $r ( @{$c->{res}} ) {
     foreach my $s ( @{$self->{ssbond}} ) {
      if ($s->{resnum1}==$r->{num} &&
	  $s->{chain1}==$r->{chain}) {
	$s->{chain1}=$chainid;
      }
      if ($s->{resnum2}==$r->{num} &&
	  $s->{chain2}==$r->{chain}) {
	$s->{chain2}=$chainid;
      }
     }
    }

    foreach my $r ( @{$c->{res}} ) {
      $r->{chain}=$chainid;
    }

    $self->{chainlookup}->{$chainid}=$c;
  }
}

sub setSegment {
  my $self=shift;
  my $segid=shift;
  my $setall=shift;

  $setall=0 if (!defined $setall);

  my $achains=$self->activeChains();
  my $maxinx=($setall)?$#{$achains}:0;

  for (my $ic=0; $ic<=$maxinx; $ic++) {
    my $c=$achains->[$ic];

    foreach my $a ( @{$c->{atom}} ) {
      $a->{seg}=$segid;
    }

    foreach my $r ( @{$c->{res}} ) {
      $r->{seg}=$segid;
    }
  }
}

## 
## method: match(refMolecule)
## matches residues between two structures
## with different residue numbering and/or
## missing residues. This method changes
## the residue numbering of the current structure
## to match the reference structure given
## as the argument

sub match {
  my $self=shift;
  my $ref=shift;

  my $sc=$self->activeChains(shift)->[0];

  die "cannot find chain"
    if (!defined $sc);

  my $rc=$ref->activeChains($sc->{id})->[0];
  $rc=$ref->activeChains()->[0] if (!defined $rc);
  
  die "cannot match chain"
    if (!defined $rc);

  my $nself=$#{$sc->{res}}+1;
  my $nref=$#{$rc->{res}}+1; 

  my $maxshift=0;
  my $maxfirst=0;
  my $maxlast=0;
  my $maxmatch=0;
  
  for (my $shift=-$nref; $shift<=$nself+$nref; $shift++) {
    my $nmatch=0;
    my $first=undef;
    my $last=undef;

    for (my $j=0; $j<$nself; $j++) {
      my $sr=$sc->{res}->[$j];
      my $rr=($j>=$shift)?$rc->{res}->[$j-$shift]:undef;
      if (defined $rr && &_cmpResName($sr->{name},$rr->{name})) {
	$first=$j if (!defined $first);
	$last=$j;
	$nmatch++;
      } 
    } 
    if ($nmatch>$maxmatch) {
      $maxmatch=$nmatch;
      $maxshift=$shift;
      $maxfirst=$first;
      $maxlast=$last;
    }
  }    

  my $delta=$rc->{res}->[$maxfirst-$maxshift]->{num}-$sc->{res}->[$maxfirst]->{num};
  $self->shiftResNumber($delta,$sc->{id});

  return ($sc->{res}->[$maxfirst]->{num},$sc->{res}->[$maxlast]->{num});
}

## method: $index = getResidue(resnum[,chain])
## returns the residue index for a residue number. 
## A chain ID may be given as additional argument for
## multi-domain structures

sub getResidue {
  my $self=shift;
  my $inx=shift;
  my $chain=shift;

  my $c=$self->activeChains($chain)->[0];
  return $self->getResidueInChain($inx,$c);
}  

## method: $index = getResidueInChain(resnum[,chain])
## returns the residue index for a residue number. 
## A reference to a chain structure may be given as 
## additional argument for multi-domain structures

sub getResidueInChain {
  my $self=shift;
  my $inx=shift;
  my $c=shift;
  my $seg=shift;

  return undef if (!defined $c);

  if (!defined $c->{resinx}) {
    $c->{resinx}={};
    for (my $ir=0; $ir<=$#{$c->{res}}; $ir++) {
      my $key=$c->{res}->[$ir]->{num};
      my $segkey=$c->{res}->[$ir]->{num}."s".$c->{res}->[$ir]->{seg};
      $c->{resinx}->{$key}=$ir;
      $c->{resinx}->{$segkey}=$ir;
    }
  }

  my $ikey=$inx;
  $ikey.="s".$seg if (defined $seg);
  return undef if (!defined $c->{resinx}->{$ikey});

  return $c->{res}->[$c->{resinx}->{$ikey}];
}  

## method: $atom = getAtomInResidue(residue,name,chain)
## returns the atom index for an atom name

sub getAtomInResidue {
  my $self=shift;
  my $res=shift;
  my $aname=shift;
  my $chain=shift;

  $chain=$self->activeChains($chain)->[0] if (!defined $chain);
  my $a=$chain->{atom};

  my $found=-1;
  for (my $i=$res->{start}; $i<=$res->{end} && $found<0; $i++) {
    if ($a->[$i]->{atomname} eq $aname) {
      $found=$i;
    }
  }
  if ($found>=0) {
    return $a->[$found];
  } else {
    return undef;
  }
}


## method: $distance = minDistance(ir,jr)
## returns the minimum distance between heavy atoms of residue
## structures ir and jr.

sub minDistance {
  my $self=shift;
  my $ir=shift;
  my $jr=shift;

  my $ic=$self->getChain($ir->{chain});
  my $jc=$self->getChain($jr->{chain});

  my $ia=$ic->{atom};
  my $ja=$jc->{atom};

  my $ix=$ic->{xcoor};
  my $iy=$ic->{ycoor};
  my $iz=$ic->{zcoor};

  my $jx=$jc->{xcoor};
  my $jy=$jc->{ycoor};
  my $jz=$jc->{zcoor};

  my $mind=1E99;
  for (my $i=$ir->{start}; $i<=$ir->{end}; $i++) {
    if (!$ia->[$i]->{hyd}) {
      for (my $j=$jr->{start}; $j<=$jr->{end}; $j++) {
	if (!$ja->[$j]->{hyd}) {
	  my $dx=$ix->[$i]-$jx->[$j];
	  my $dy=$iy->[$i]-$jy->[$j];
	  my $dz=$iz->[$i]-$jz->[$j];
	  my $d=($dx*$dx+$dy*$dy+$dz*$dz);
	  $mind=$d if ($d<$mind);
	}
      }
    }
  }
  return sqrt($mind);
}

## method: number = firstResNum([chain])
## returns the number of the first residue. If
## a chain reference is given as argument it returns the
## number of the first residue of that chain

sub firstResNum {
  my $self=shift;
  my $chain=shift;

  $chain=$self->activeChains()->[0]
    if (!defined $chain);

  return undef if (!defined $chain);

  return $chain->{res}->[0]->{num};
}

## method: number = lastResNum([chain])
## returns the number of the last residue. If
## a chain reference is given as argument it returns the
## number of the last residue of that chain.

sub lastResNum {
  my $self=shift;
  my $chain=shift;

  $chain=$self->activeChains()->[0]
    if (!defined $chain);

  return undef if (!defined $chain);

  return $chain->{res}->[$#{$chain->{res}}]->{num};
}

## method: ret = empty()
## returns 1 or 0 whether the current molecule
## is empty

sub empty {
  my $self=shift;
  
  if (!defined $self->{chain} || $#{$self->{chain}}<0) {
    return 1;
  }

  my $nat=0;
  foreach my $c ( @{$self->{chain}} ) {
    $nat+=$#{$c->{atom}}+1;
  }

  if ($nat==0) {
    return 1;
  }

  return 0;
}

sub _coorCache {
  my $self=shift;
  my $chainid=shift;

  my $chains=$self->activeChains($chainid);
  if (defined $chains) {
    foreach my $c ( @{$self->activeChains($chainid)} ) {
      my $a=$c->{atom};
      for (my $i=0; $i<=$#{$a}; $i++) {
	$c->{xcoor}->[$i]=$a->[$i]->{xcoor};
	$c->{ycoor}->[$i]=$a->[$i]->{ycoor};
	$c->{zcoor}->[$i]=$a->[$i]->{zcoor};
      }
    }
  }

  $self->_generateLookupTable();
}

sub _pdbLine {
  my $pdbrec=shift;
  my $chmode=shift;
  my $longaux2=shift;
  $chmode=0 if (!defined $chmode);

  my $chainid=$pdbrec->{chain};
  $chainid=" " if (!defined $chainid || $chainid eq "" || $chainid eq "+");

  my $resnumstr;
  if ($pdbrec->{resnum}>999 && $chmode) {
    if (defined $pdbrec->{aresnum} && $pdbrec->{aresnum}=~/[A-Z]/) {
      $resnumstr=sprintf(" %4s ",$pdbrec->{aresnum});
    } elsif ($pdbrec->{resnum}>9999) {
      $resnumstr=sprintf("%5d",$pdbrec->{resnum});
    } else {
      $resnumstr=sprintf(" %4d ",$pdbrec->{resnum});
    }
  } else {
    if (defined $pdbrec->{aresnum} && $pdbrec->{aresnum}=~/[A-Z]/) {
      $resnumstr=sprintf(" %4s ",$pdbrec->{aresnum});
    } elsif ($pdbrec->{resnum}>9999) {
      $resnumstr=sprintf("%5d",$pdbrec->{resnum});
    } else {
      $resnumstr=sprintf("%4d  ",$pdbrec->{resnum});
    }
  }

  my $aux1=1.0;
  my $aux2=0.0;

  $aux1=$pdbrec->{aux1} if (defined $pdbrec->{aux1});
  $aux2=$pdbrec->{aux2} if (defined $pdbrec->{aux2});

  my $atomnum;
  if ($pdbrec->{chain} eq "+") {
    if ($pdbrec->{atominx}<100000) {
      $atomnum=sprintf("HETATM%5d",$pdbrec->{atominx});
    } else {
      $atomnum=sprintf("HETAT%6d",$pdbrec->{atominx});
    }
  } else {
    $atomnum=sprintf("ATOM%7d",$pdbrec->{atominx});
  }
  if (length($pdbrec->{atomname})>3 || ($pdbrec->{atomname}=~/[0-9]H.*/)) {
    if (defined $longaux2 && $longaux2) {
      return sprintf "%11s %-4s %-4s%1s%-6s  %8.3f%8.3f%8.3f %5.2f %6.3f     %-4s",
	$atomnum,$pdbrec->{atomname}, 
	  $pdbrec->{resname},$chainid,$resnumstr, 
	    $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $aux1, $aux2, $pdbrec->{seg};
    } else {
      return sprintf "%11s %-4s %-4s%1s%-6s  %8.3f%8.3f%8.3f %5.2f%6.2f      %-4s",
	$atomnum,$pdbrec->{atomname}, 
	  $pdbrec->{resname},$chainid,$resnumstr, 
	    $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $aux1, $aux2, $pdbrec->{seg};
    }
  } else {
    if (defined $longaux2 && $longaux2) {
      return sprintf "%11s  %-3s %-4s%1s%-6s  %8.3f%8.3f%8.3f %5.2f %6.3f     %-4s",
	$atomnum,$pdbrec->{atomname}, 
	  $pdbrec->{resname},$chainid,$resnumstr, 
	    $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $aux1, $aux2, $pdbrec->{seg};
    } else {
      return sprintf "%11s  %-3s %-4s%1s%-6s  %8.3f%8.3f%8.3f %5.2f%6.2f      %-4s",
	$atomnum,$pdbrec->{atomname}, 
	  $pdbrec->{resname},$chainid,$resnumstr, 
	    $pdbrec->{xcoor}, $pdbrec->{ycoor}, $pdbrec->{zcoor}, $aux1, $aux2, $pdbrec->{seg};
    }      
   } 
}

sub _newChain {
  my $self=shift;
  my $cid=shift;

#  printf "new Chain >%s<\n",$cid;

  my $chainrec={};
  $chainrec->{id}=($cid eq " " || $cid eq "")?"":$cid;
  $chainrec->{atom}=();
  $chainrec->{res}=();
  $chainrec->{xcoor}=();
  $chainrec->{ycoor}=();
  $chainrec->{zcoor}=();
  push(@{$self->{chain}},$chainrec);

  $self->{chainlookup}->{$chainrec->{id}}=$chainrec
    if ($chainrec->{id} ne "");

  return $chainrec;
}

sub _cmpResName {
  my $name1=shift;
  my $name2=shift;
  return ($name1 eq $name2 || 
	  ($name1 =~ /HIS|HSD|HSE|HSP/ && $name2 =~ /HIS|HSD|HSE|HSP/));
}

# Return the name a His should have by looking for proton(s)
sub getHisType {
  my $self=shift;
  my $c=shift;
  my $resnum=shift;

  my $type="HIS";

  my $r=$self->getResidueInChain($resnum,$c);
   if (exists $r->{lockname} || $r->{name} eq "HSE" || $r->{name} eq "HSD" || $r->{name} eq "HSD2" || $r->{name} eq "HSE2") {
       return $r->{name};
  }

  for (my $i=$r->{start}; $i<=$r->{end}; $i++) {
      my $aname=$c->{atom}->[$i]->{atomname};
      my $atype=$c->{atom}->[$i]->{atomtype};

      if ($aname eq "HD1") {
	  if ($type eq "HSE") {
	      $type="HSP";
	  } else {
	      $type="HSD";
	  }
      } elsif ($aname eq "HE2") {
	  if ($type eq "HSD") {
	      $type="HSP";
	  } else {
	      $type="HSE";
	  }
      }

  }

  return $type;
}


sub getSSBonds {
  my $self=shift;
  
  my $ns=();
  if ($#{$self->{ssbond}}>=0) {
    foreach my $s ( @{$self->{ssbond}} ) {
      my $nrec={};
      %{$nrec}=%{$s};
      push(@{$ns},$nrec);
    }
  }
  return $ns;
}

sub setSSBonds {
  my $self=shift;
  my $ss=shift;
  
#  $self->{ssbond}=();
  if (defined $ss && $#{$ss}>=0) {
    foreach my $s ( @{$ss} ) { 
      my $nrec={};
      %{$nrec}=%{$s};
      push(@{$self->{ssbond}},$nrec);
    }
  }
}  

sub findSSBonds {
  my $self=shift;

  my $distCut=3.0;
  my $sqDistCut=$distCut*$distCut;
  $self->{ssbond}=();

  # First collect all the Cys S atoms
  my $Slist=();

  foreach my $c ( @{$self->activeChains()} ) {
    if ($#{$c->{atom}}>=0) {
      my $aa=$c->{atom};
      foreach my $a ( @{$aa} ) {
	  push(@{$Slist},$a) if ($a->{atomname} eq "SG");
      }
    }
  }

  # Make nearby S atoms SS bonded
  for (my $i=0; $i<$#{$Slist}; $i++) {
    my $ia=$Slist->[$i];
    my $ix=$ia->{xcoor};
    my $iy=$ia->{ycoor};
    my $iz=$ia->{zcoor};

    for (my $j=$i+1; $j<=$#{$Slist}; $j++) {
       my $ja=$Slist->[$j];
       my $jx=$ja->{xcoor};
       my $jy=$ja->{ycoor};
       my $jz=$ja->{zcoor};
       my $dx=$ix-$jx;
       my $dy=$iy-$jy;
       my $dz=$iz-$jz;
       my $d=($dx*$dx+$dy*$dy+$dz*$dz);

       if ($d<$sqDistCut) {
	 my $trec={};
	 $trec->{chain1}=$ia->{chain};
	 $trec->{resnum1}=$ia->{resnum};
	 $trec->{chain2}=$ja->{chain};
	 $trec->{resnum2}=$ja->{resnum};
	 push(@{$self->{ssbond}},$trec);
       }
    }
  }


  # Ensure none have multiple partners
  # (required because of the naive method used for pairing)
  my $usedRes={};
  foreach my $SS ( @{$self->{ssbond}} ) {
      my $Sa=$SS->{chain1}.$SS->{resnum1};
      if (! defined $usedRes->{$Sa}) {
	  $usedRes->{$Sa}=1;
      } else {
	  die "Unable to perform SS pairing";
      }
      $Sa=$SS->{chain2}.$SS->{resnum2};
      if (! defined $usedRes->{$Sa}) {
	  $usedRes->{$Sa}=1;
      } else {
	  die "Unable to perform SS pairing";
      }
  }

  # Return the number of SS bonds found
  return ($#{$self->{ssbond}}+1);
}

## function: $mol = genHTRNA()
## generate Hori & Takada CG model for RNA
## Hori, N. & Takada, S. (2012), JCTC 8, 3384-3394
##

sub genHTRNA {
  my $self=shift;

  my $n={};
  $n->{chain}=();
  $n->{chainlookup}={};
  $n->{defchain}=undef;
  $n->{segmentlist}=undef;
  $n->{ssbond}=();

  bless $n;

  foreach my $s ( @{$self->{ssbond}} ) {
    my $srec={};
    %{$srec}=%{$s};
    push (@{$n->{ssbond}},$srec);
  }

  foreach my $c ( @{$self->activeChains()} ) {
    my $nc=undef;

    my $natom=1;
    foreach my $r (@{$c->{res}}) {
      $nc=$n->_newChain($c->{id}) if (!defined $nc);

      my $rrec={};
      %{$rrec}=%{$r};
      $rrec->{start}=$#{$nc->{atom}}+1;
      push(@{$nc->{res}},$rrec);
      my %alook;
      for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	my $a=$c->{atom}->[$ia];
	$alook{$a->{atomname}}=$a;
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{'P'}) {
	my $x=$alook{'P'}->{xcoor};
	my $y=$alook{'P'}->{ycoor};
	my $z=$alook{'P'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "P", seg => $r->{seg},
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C1'"} && exists $alook{"C2'"} && exists $alook{"C3'"} && exists $alook{"C4'"} && exists $alook{"O4'"}) {
	my $x=($alook{"C1'"}->{xcoor}+$alook{"C2'"}->{xcoor}+$alook{"C3'"}->{xcoor}+$alook{"C4'"}->{xcoor}+$alook{"O4'"}->{xcoor})/5.0;
	my $y=($alook{"C1'"}->{ycoor}+$alook{"C2'"}->{ycoor}+$alook{"C3'"}->{ycoor}+$alook{"C4'"}->{ycoor}+$alook{"O4'"}->{ycoor})/5.0;
	my $z=($alook{"C1'"}->{zcoor}+$alook{"C2'"}->{zcoor}+$alook{"C3'"}->{zcoor}+$alook{"C4'"}->{zcoor}+$alook{"O4'"}->{zcoor})/5.0;
	my $arec={ atominx => $natom++, atomname => "S", seg => $r->{seg},
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && exists $alook{"N3"}) {
	my $x=$alook{"N3"}->{xcoor};
	my $y=$alook{"N3"}->{ycoor};
	my $z=$alook{"N3"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "B", seg => $r->{seg},
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" ) && exists $alook{"N1"}) {
	my $x=$alook{"N1"}->{xcoor};
	my $y=$alook{"N1"}->{ycoor};
	my $z=$alook{"N1"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "B", seg => $r->{seg},
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      $rrec->{end}=$#{$nc->{atom}};
    }
  }

  $n->_coorCache();

  return $n;
}


## function: $mol = genPRIMO()
## generate reduced PRIMO (protein intermediate model) representation
## from all-atom structure (handles both proteins and nucleic acids)
## 

sub genPRIMO {
  my $self=shift;

  my $n={};
  $n->{chain}=();
  $n->{chainlookup}={};
  $n->{defchain}=undef;
  $n->{segmentlist}=undef;
  $n->{ssbond}=();

  bless $n;

  foreach my $s ( @{$self->{ssbond}} ) {
    my $srec={};
    %{$srec}=%{$s};
    push (@{$n->{ssbond}},$srec);
  }

  foreach my $c ( @{$self->activeChains()} ) {
    my $nc=undef;

    my $natom=1;
    foreach my $r (@{$c->{res}}) {
      $nc=$n->_newChain($c->{id}) if (!defined $nc);

      my $rrec={};
      %{$rrec}=%{$r};
      $rrec->{start}=$#{$nc->{atom}}+1;
      push(@{$nc->{res}},$rrec);
      my %alook;
      for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	my $a=$c->{atom}->[$ia];
	$alook{$a->{atomname}}=$a;
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{"N"}) {
	my $x=$alook{"N"}->{xcoor};
	my $y=$alook{"N"}->{ycoor};
	my $z=$alook{"N"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "N", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      } elsif ($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") {   
        # do nothing if the residues belong to DNA or RNA
      } else {
	printf STDERR "cannot find atom N for residue %s:%s:%d\n",$r->{chain},$r->{name},$r->{num};
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{'CA'}) {
	my $x=$alook{'CA'}->{xcoor};
	my $y=$alook{'CA'}->{ycoor};
	my $z=$alook{'CA'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "CA", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      } elsif ($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") {
        # do nothing if the residues belong to DNA or RNA
      } else {
	printf STDERR "cannot find atom CA for residue %s:%s:%d\n",$r->{chain},$r->{name},$r->{resnum};
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{'C'} && exists $alook{'O'}) {
	my $x=($alook{'C'}->{xcoor}+$alook{'O'}->{xcoor})/2.0;
	my $y=($alook{'C'}->{ycoor}+$alook{'O'}->{ycoor})/2.0;
	my $z=($alook{'C'}->{zcoor}+$alook{'O'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "CO", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      } elsif (exists $alook{'C'} && exists $alook{'OT1'}) {
	my $x=($alook{'C'}->{xcoor}+$alook{'OT1'}->{xcoor})/2.0;
	my $y=($alook{'C'}->{ycoor}+$alook{'OT1'}->{ycoor})/2.0;
	my $z=($alook{'C'}->{zcoor}+$alook{'OT1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "CO", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      } elsif ($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") {
        # do nothing if the residues belong to DNA or RNA
      } else {
	printf STDERR "cannot find atom CO for residue %s:%s:%d\n",$r->{chain},$r->{name},$r->{resnum};
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP") && 
	  exists $alook{'CB'}) {
	my $x=$alook{'CB'}->{xcoor};
	my $y=$alook{'CB'}->{ycoor};
	my $z=$alook{'CB'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

       if (($r->{name} eq "VAL")  &&
         exists $alook{'CB'} && exists $alook{'CG1'}) {
       my $x=($alook{'CB'}->{xcoor}+$alook{'CG1'}->{xcoor})/2.0;
       my $y=($alook{'CB'}->{ycoor}+$alook{'CG1'}->{ycoor})/2.0;
       my $z=($alook{'CB'}->{zcoor}+$alook{'CG1'}->{zcoor})/2.0;
       my $arec={ atominx => $natom++, atomname => "SC1", 
                  resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                  xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
       push(@{$nc->{atom}},$arec);
     }

       if (($r->{name} eq "VAL")  &&
         exists $alook{'CB'} && exists $alook{'CG2'}) {
       my $x=($alook{'CB'}->{xcoor}+$alook{'CG2'}->{xcoor})/2.0;
       my $y=($alook{'CB'}->{ycoor}+$alook{'CG2'}->{ycoor})/2.0;
       my $z=($alook{'CB'}->{zcoor}+$alook{'CG2'}->{zcoor})/2.0;
       my $arec={ atominx => $natom++, atomname => "SC2", 
                  resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                  xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
       push(@{$nc->{atom}},$arec);
     }


      if ($r->{name} eq "LEU" && exists $alook{'CG'} && exists $alook{'CD1'} && exists $alook{'CD2'}) {
	my $x=($alook{'CG'}->{xcoor}+$alook{'CD1'}->{xcoor}+$alook{'CD2'}->{xcoor})/3.0;
	my $y=($alook{'CG'}->{ycoor}+$alook{'CD1'}->{ycoor}+$alook{'CD2'}->{ycoor})/3.0;
	my $z=($alook{'CG'}->{zcoor}+$alook{'CD1'}->{zcoor}+$alook{'CD2'}->{zcoor})/3.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ILE") &&
          exists $alook{'CG2'}) {
        my $x=$alook{'CG2'}->{xcoor};
        my $y=$alook{'CG2'}->{ycoor};
        my $z=$alook{'CG2'}->{zcoor};
        my $arec={ atominx => $natom++, atomname => "SC1",
                   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
        push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ILE" && exists $alook{'CD'} && exists $alook{'CG1'}) {
        my $x=($alook{'CD'}->{xcoor}+$alook{'CG1'}->{xcoor})/2.0;
        my $y=($alook{'CD'}->{ycoor}+$alook{'CG1'}->{ycoor})/2.0;
        my $z=($alook{'CD'}->{zcoor}+$alook{'CG1'}->{zcoor})/2.0;
        my $arec={ atominx => $natom++, atomname => "SC2",
                   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
        push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ILE" && exists $alook{'CD1'} && exists $alook{'CG1'}) {
        my $x=($alook{'CD1'}->{xcoor}+$alook{'CG1'}->{xcoor})/2.0;
        my $y=($alook{'CD1'}->{ycoor}+$alook{'CG1'}->{ycoor})/2.0;
        my $z=($alook{'CD1'}->{zcoor}+$alook{'CG1'}->{zcoor})/2.0;
        my $arec={ atominx => $natom++, atomname => "SC2",
                   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
        push(@{$nc->{atom}},$arec);
      }


      if ($r->{name} eq "PRO" && exists $alook{'CB'} && exists $alook{'CG'} && exists $alook{'CD'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'CG'}->{xcoor}+$alook{'CD'}->{xcoor})/3.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'CG'}->{ycoor}+$alook{'CD'}->{ycoor})/3.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'CG'}->{zcoor}+$alook{'CD'}->{zcoor})/3.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "CYS" && exists $alook{'CB'} && exists $alook{'SG'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'SG'}->{xcoor})/2.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'SG'}->{ycoor})/2.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'SG'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || 
	   $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" ||
           $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" ||
           $r->{name} eq "LEU" )  && 
	  exists $alook{'CB'} && exists $alook{'CG'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'CG'}->{xcoor})/2.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'CG'}->{ycoor})/2.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'CG'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "MET" && exists $alook{'SD'}) {
	my $x=$alook{'SD'}->{xcoor};
	my $y=$alook{'SD'}->{ycoor};
	my $z=$alook{'SD'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "MET" && exists $alook{'CE'}) {
	my $x=$alook{'CE'}->{xcoor};
	my $y=$alook{'CE'}->{ycoor};
	my $z=$alook{'CE'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }


      if (($r->{name} eq "PHE" || $r->{name} eq "TYR") && exists $alook{'CE1'}) {
	my $x=$alook{'CE1'}->{xcoor};
	my $y=$alook{'CE1'}->{ycoor};
	my $z=$alook{'CE1'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "PHE" || $r->{name} eq "TYR") && exists $alook{'CE2'}) {
	my $x=$alook{'CE2'}->{xcoor};
	my $y=$alook{'CE2'}->{ycoor};
	my $z=$alook{'CE2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "TRP" && exists $alook{'NE1'}) {
	my $x=$alook{'NE1'}->{xcoor};
	my $y=$alook{'NE1'}->{ycoor};
	my $z=$alook{'NE1'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "TRP" && exists $alook{'CE3'}) {
	my $x=$alook{'CE3'}->{xcoor};
	my $y=$alook{'CE3'}->{ycoor};
	my $z=$alook{'CE3'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "TRP" && exists $alook{'CZ2'}) {
	my $x=$alook{'CZ2'}->{xcoor};
	my $y=$alook{'CZ2'}->{ycoor};
	my $z=$alook{'CZ2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }
      
      if ($r->{name} eq "SER" && exists $alook{'CB'} && exists $alook{'OG'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'OG'}->{xcoor})/2.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'OG'}->{ycoor})/2.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'OG'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "THR" && exists $alook{'CB'} && exists $alook{'OG1'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'OG1'}->{xcoor})/2.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'OG1'}->{ycoor})/2.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'OG1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }
      
      if ($r->{name} eq "THR" && exists $alook{'CG2'}) {
	my $x=$alook{'CG2'}->{xcoor};
	my $y=$alook{'CG2'}->{ycoor};
	my $z=$alook{'CG2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "TYR" && exists $alook{'OH'}) {
	my $x=$alook{'OH'}->{xcoor};
	my $y=$alook{'OH'}->{ycoor};
	my $z=$alook{'OH'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ASN" && exists $alook{'CG'} && exists $alook{'OD1'}) {
	my $x=($alook{'CG'}->{xcoor}+$alook{'OD1'}->{xcoor})/2.0;
	my $y=($alook{'CG'}->{ycoor}+$alook{'OD1'}->{ycoor})/2.0;
	my $z=($alook{'CG'}->{zcoor}+$alook{'OD1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ASN" && exists $alook{'ND2'}) {
	my $x=$alook{'ND2'}->{xcoor};
	my $y=$alook{'ND2'}->{ycoor};
	my $z=$alook{'ND2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "GLN" && exists $alook{'CD'} && exists $alook{'OE1'}) {
	my $x=($alook{'CD'}->{xcoor}+$alook{'OE1'}->{xcoor})/2.0;
	my $y=($alook{'CD'}->{ycoor}+$alook{'OE1'}->{ycoor})/2.0;
	my $z=($alook{'CD'}->{zcoor}+$alook{'OE1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "GLN" && exists $alook{'NE2'}) {
	my $x=$alook{'NE2'}->{xcoor};
	my $y=$alook{'NE2'}->{ycoor};
	my $z=$alook{'NE2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }
      
      if ($r->{name} eq "LYS" && exists $alook{'CD'}) {
	my $x=$alook{'CD'}->{xcoor};
	my $y=$alook{'CD'}->{ycoor};
	my $z=$alook{'CD'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "LYS" && exists $alook{'CE'}) {
	my $x=$alook{'CE'}->{xcoor};
	my $y=$alook{'CE'}->{ycoor};
	my $z=$alook{'CE'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }


      if ($r->{name} eq "LYS" && exists $alook{'NZ'}) {
	my $x=$alook{'NZ'}->{xcoor};
	my $y=$alook{'NZ'}->{ycoor};
	my $z=$alook{'NZ'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }


      if ($r->{name} eq "ARG" && exists $alook{'CD'}) {
	my $x=$alook{'CD'}->{xcoor};
	my $y=$alook{'CD'}->{ycoor};
	my $z=$alook{'CD'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ARG" && exists $alook{'NE'} && exists $alook{'CZ'}) {
	my $x=$alook{'NE'}->{xcoor};
	my $y=$alook{'NE'}->{ycoor};
	my $z=$alook{'NE'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ARG" && exists $alook{'NH1'}) {
	my $x=$alook{'NH1'}->{xcoor};
	my $y=$alook{'NH1'}->{ycoor};
	my $z=$alook{'NH1'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ARG" && exists $alook{'NH2'}) {
	my $x=$alook{'NH2'}->{xcoor};
	my $y=$alook{'NH2'}->{ycoor};
	my $z=$alook{'NH2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "HSD" || $r->{name} eq "HSE" || 
	   $r->{name} eq "HIS" || $r->{name} eq "HSP") && exists $alook{'ND1'}) {
	my $x=$alook{'ND1'}->{xcoor};
	my $y=$alook{'ND1'}->{ycoor};
	my $z=$alook{'ND1'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "HSD" || $r->{name} eq "HSE" || 
	   $r->{name} eq "HIS" || $r->{name} eq "HSP") && exists $alook{'NE2'}) {
	my $x=$alook{'NE2'}->{xcoor};
	my $y=$alook{'NE2'}->{ycoor};
	my $z=$alook{'NE2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ASP" && exists $alook{'CG'} && exists $alook{'OD1'}) {
	my $x=($alook{'CG'}->{xcoor}+$alook{'OD1'}->{xcoor})/2.0;
	my $y=($alook{'CG'}->{ycoor}+$alook{'OD1'}->{ycoor})/2.0;
	my $z=($alook{'CG'}->{zcoor}+$alook{'OD1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ASP" && exists $alook{'CG'} && exists $alook{'OD2'}) {
	my $x=($alook{'CG'}->{xcoor}+$alook{'OD2'}->{xcoor})/2.0;
	my $y=($alook{'CG'}->{ycoor}+$alook{'OD2'}->{ycoor})/2.0;
	my $z=($alook{'CG'}->{zcoor}+$alook{'OD2'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "GLU" && exists $alook{'CD'} && exists $alook{'OE1'}) {
	my $x=($alook{'CD'}->{xcoor}+$alook{'OE1'}->{xcoor})/2.0;
	my $y=($alook{'CD'}->{ycoor}+$alook{'OE1'}->{ycoor})/2.0;
	my $z=($alook{'CD'}->{zcoor}+$alook{'OE1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "GLU" && exists $alook{'CD'} && exists $alook{'OE2'}) {
	my $x=($alook{'CD'}->{xcoor}+$alook{'OE2'}->{xcoor})/2.0;
	my $y=($alook{'CD'}->{ycoor}+$alook{'OE2'}->{ycoor})/2.0;
	my $z=($alook{'CD'}->{zcoor}+$alook{'OE2'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

## section for nucleic acids

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{'O1P'}) {
	my $x=$alook{'O1P'}->{xcoor};
	my $y=$alook{'O1P'}->{ycoor};
	my $z=$alook{'O1P'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BB8", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{'O2P'}) {
	my $x=$alook{'O2P'}->{xcoor};
	my $y=$alook{'O2P'}->{ycoor};
	my $z=$alook{'O2P'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BB7", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"O5'"} && exists $alook{"C5'"}) {
	my $x=($alook{"O5'"}->{xcoor}+$alook{"C5'"}->{xcoor})/2.0;
	my $y=($alook{"O5'"}->{ycoor}+$alook{"C5'"}->{ycoor})/2.0;
	my $z=($alook{"O5'"}->{zcoor}+$alook{"C5'"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BB6", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }


## for 5' terminal  

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{'O5T'}) {
	my $x=$alook{'O5T'}->{xcoor};
	my $y=$alook{'O5T'}->{ycoor};
	my $z=$alook{'O5T'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "B5T", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C5'"} && exists $alook{"O5'"} && exists $alook{"O1P"} && exists $alook{"O2P"}) {
# do nothing
      }
      elsif (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C5'"} && exists $alook{"O5'"}) {
	my $x=$alook{"O5'"}->{xcoor};
	my $y=$alook{"O5'"}->{ycoor};
	my $z=$alook{"O5'"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BT5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"O4'"}) {
	my $x=$alook{"O4'"}->{xcoor};
	my $y=$alook{"O4'"}->{ycoor};
	my $z=$alook{"O4'"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BB5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C4'"} && exists $alook{"C3'"}) {
	my $x=($alook{"C4'"}->{xcoor}+$alook{"C3'"}->{xcoor})/2.0;
	my $y=($alook{"C4'"}->{ycoor}+$alook{"C3'"}->{ycoor})/2.0;
	my $z=($alook{"C4'"}->{zcoor}+$alook{"C3'"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BB4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

## for RNA O2'-C2' else for DNA C2'

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C2'"} && exists $alook{"O2'"}) {
	my $x=($alook{"C2'"}->{xcoor}+$alook{"O2'"}->{xcoor})/2.0;
	my $y=($alook{"C2'"}->{ycoor}+$alook{"O2'"}->{ycoor})/2.0;
	my $z=($alook{"C2'"}->{zcoor}+$alook{"O2'"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BR2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }
      elsif (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C2'"}) {
	my $x=$alook{"C2'"}->{xcoor};
	my $y=$alook{"C2'"}->{ycoor};
	my $z=$alook{"C2'"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BD2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA") && 
	  exists $alook{"C1'"} && exists $alook{"N9"}) {
	my $x=($alook{"C1'"}->{xcoor}+$alook{"N9"}->{xcoor})/2.0;
	my $y=($alook{"C1'"}->{ycoor}+$alook{"N9"}->{ycoor})/2.0;
	my $z=($alook{"C1'"}->{zcoor}+$alook{"N9"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BB1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") &&           exists $alook{"C1'"} && exists $alook{"N1"}) {
	my $x=($alook{"C1'"}->{xcoor}+$alook{"N1"}->{xcoor})/2.0;
	my $y=($alook{"C1'"}->{ycoor}+$alook{"N1"}->{ycoor})/2.0;
	my $z=($alook{"C1'"}->{zcoor}+$alook{"N1"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BB1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") &&           exists $alook{"C2"} && exists $alook{"O2"}) {
	my $x=($alook{"C2"}->{xcoor}+$alook{"O2"}->{xcoor})/2.0;
	my $y=($alook{"C2"}->{ycoor}+$alook{"O2"}->{ycoor})/2.0;
	my $z=($alook{"C2"}->{zcoor}+$alook{"O2"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BS1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") &&           exists $alook{"N3"}) {
	my $x=$alook{"N3"}->{xcoor};
	my $y=$alook{"N3"}->{ycoor};
	my $z=$alook{"N3"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "THY" || $r->{name} eq "URA") && exists $alook{"C4"} && exists $alook{"O4"}) {
	my $x=($alook{"C4"}->{xcoor}+$alook{"O4"}->{xcoor})/2.0;
	my $y=($alook{"C4"}->{ycoor}+$alook{"O4"}->{ycoor})/2.0;
	my $z=($alook{"C4"}->{zcoor}+$alook{"O4"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BS3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT") && exists $alook{"C4"} && exists $alook{"N4"}) {
	my $x=$alook{"N4"}->{xcoor};
	my $y=$alook{"N4"}->{ycoor};
	my $z=$alook{"N4"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "THY") && exists $alook{"C5M"}) {
	my $x=$alook{"C5M"}->{xcoor};
	my $y=$alook{"C5M"}->{ycoor};
	my $z=$alook{"C5M"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT" || $r->{name} eq "URA") && exists $alook{"C5"} && exists $alook{"C6"}) {
	my $x=($alook{"C5"}->{xcoor}+$alook{"C6"}->{xcoor})/2.0;
	my $y=($alook{"C5"}->{ycoor}+$alook{"C6"}->{ycoor})/2.0;
	my $z=($alook{"C5"}->{zcoor}+$alook{"C6"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BS4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "THY") && exists $alook{"C6"}) {
	my $x=$alook{"C6"}->{xcoor};
	my $y=$alook{"C6"}->{ycoor};
	my $z=$alook{"C6"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA") && exists $alook{"N3"}) {
	my $x=$alook{"N3"}->{xcoor};
	my $y=$alook{"N3"}->{ycoor};
	my $z=$alook{"N3"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE") && exists $alook{"N1"}) {
	my $x=$alook{"N1"}->{xcoor};
	my $y=$alook{"N1"}->{ycoor};
	my $z=$alook{"N1"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "GUA") && exists $alook{"N2"}) {
	my $x=$alook{"N2"}->{xcoor};
	my $y=$alook{"N2"}->{ycoor};
	my $z=$alook{"N2"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE") && exists $alook{"N6"}) {
	my $x=$alook{"N6"}->{xcoor};
	my $y=$alook{"N6"}->{ycoor};
	my $z=$alook{"N6"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "GUA") && exists $alook{"N1"}) {
	my $x=$alook{"N1"}->{xcoor};
	my $y=$alook{"N1"}->{ycoor};
	my $z=$alook{"N1"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE") && exists $alook{"N7"}) {
	my $x=$alook{"N7"}->{xcoor};
	my $y=$alook{"N7"}->{ycoor};
	my $z=$alook{"N7"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "GUA") && exists $alook{"C6"} && exists $alook{"O6"}) {
	my $x=($alook{"C6"}->{xcoor}+$alook{"O6"}->{xcoor})/2.0;
	my $y=($alook{"C6"}->{ycoor}+$alook{"O6"}->{ycoor})/2.0;
	my $z=($alook{"C6"}->{zcoor}+$alook{"O6"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BS4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "GUA") && exists $alook{"N7"}) {
	my $x=$alook{"N7"}->{xcoor};
	my $y=$alook{"N7"}->{ycoor};
	my $z=$alook{"N7"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"O3'"}) {
	my $x=$alook{"O3'"}->{xcoor};
	my $y=$alook{"O3'"}->{ycoor};
	my $z=$alook{"O3'"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BB3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      $rrec->{end}=$#{$nc->{atom}};
    }
  }

  $n->_coorCache();

  return $n;
}

## function: $mol = genPRIMO2()
## generate reduced PRIMO (protein intermediate model) representation
## from all-atom structure (handles both proteins and nucleic acids)
## 

sub genPRIMO2 {
  my $self=shift;

  my $n={};
  $n->{chain}=();
  $n->{chainlookup}={};
  $n->{defchain}=undef;
  $n->{segmentlist}=undef;
  $n->{ssbond}=();

  bless $n;

  foreach my $s ( @{$self->{ssbond}} ) {
    my $srec={};
    %{$srec}=%{$s};
    push (@{$n->{ssbond}},$srec);
  }


  foreach my $c ( @{$self->activeChains()} ) {
    my $nc=undef;

    my $natom=1;
    foreach my $r (@{$c->{res}}) {
      $nc=$n->_newChain($c->{id}) if (!defined $nc);

      my $rrec={};
      %{$rrec}=%{$r};
      $rrec->{start}=$#{$nc->{atom}}+1;
      push(@{$nc->{res}},$rrec);
      my %alook;
      for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
	my $a=$c->{atom}->[$ia];
	$alook{$a->{atomname}}=$a;
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{'CY'} && exists $alook{'OY'}) {
	my $x=($alook{'CY'}->{xcoor}+$alook{'OY'}->{xcoor})/2.0;
	my $y=($alook{'CY'}->{ycoor}+$alook{'OY'}->{ycoor})/2.0;
	my $z=($alook{'CY'}->{zcoor}+$alook{'OY'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "COY", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{'CAY'}) {
	my $x=$alook{'CAY'}->{xcoor};
	my $y=$alook{'CAY'}->{ycoor};
	my $z=$alook{'CAY'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SNT", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{"N"}) {
	my $x=$alook{"N"}->{xcoor};
	my $y=$alook{"N"}->{ycoor};
	my $z=$alook{"N"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "N1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      } elsif ($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") {   
        # do nothing if the residues belong to DNA or RNA
      } else {
	printf STDERR "cannot find atom N for residue %s:%s:%d\n",$r->{chain},$r->{name},$r->{num};
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{'CA'}) {
	my $x=$alook{'CA'}->{xcoor};
	my $y=$alook{'CA'}->{ycoor};
	my $z=$alook{'CA'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "CA1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      } elsif ($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") {
        # do nothing if the residues belong to DNA or RNA
      } else {
	printf STDERR "cannot find atom CA for residue %s:%s:%d\n",$r->{chain},$r->{name},$r->{resnum};
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{'C'} && exists $alook{'O'}) {
	my $x=($alook{'C'}->{xcoor}+$alook{'O'}->{xcoor})/2.0;
	my $y=($alook{'C'}->{ycoor}+$alook{'O'}->{ycoor})/2.0;
	my $z=($alook{'C'}->{zcoor}+$alook{'O'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "CO", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      } elsif (exists $alook{'C'} && exists $alook{'OT1'}) {
	my $x=($alook{'C'}->{xcoor}+$alook{'OT1'}->{xcoor})/2.0;
	my $y=($alook{'C'}->{ycoor}+$alook{'OT1'}->{ycoor})/2.0;
	my $z=($alook{'C'}->{zcoor}+$alook{'OT1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "CO", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      } elsif ($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") {
        # do nothing if the residues belong to DNA or RNA
      } else {
	printf STDERR "cannot find atom CO for residue %s:%s:%d\n",$r->{chain},$r->{name},$r->{resnum};
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY")) {
	if (exists $alook{'OT2'} || exists $alook{'OXT'} ) {
           my $x=0.0;
           my $y=0.0;
           my $z=0.0;

           if (exists $alook{'OT2'}) {
  	     $x=$alook{'OT2'}->{xcoor};
	     $y=$alook{'OT2'}->{ycoor};
	     $z=$alook{'OT2'}->{zcoor};
           } elsif (exists $alook{'OXT'}) {
  	     $x=$alook{'OXT'}->{xcoor};
	     $y=$alook{'OXT'}->{ycoor};
	     $z=$alook{'OXT'}->{zcoor};
           } 
	   my $arec={ atominx => $natom++, atomname => "OX", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	   push(@{$nc->{atom}},$arec);
        }
      }


      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "LEU"  ) && 
	  exists $alook{'CB'}) {
	my $x=$alook{'CB'}->{xcoor};
	my $y=$alook{'CB'}->{ycoor};
	my $z=$alook{'CB'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

       if (($r->{name} eq "VAL")  &&
         exists $alook{'CB'} && exists $alook{'CG1'}) {
       my $x=($alook{'CB'}->{xcoor}+$alook{'CG1'}->{xcoor})/2.0;
       my $y=($alook{'CB'}->{ycoor}+$alook{'CG1'}->{ycoor})/2.0;
       my $z=($alook{'CB'}->{zcoor}+$alook{'CG1'}->{zcoor})/2.0;
       my $arec={ atominx => $natom++, atomname => "SC1", 
                  resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                  xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
       push(@{$nc->{atom}},$arec);
     }

       if (($r->{name} eq "VAL")  &&
         exists $alook{'CB'} && exists $alook{'CG2'}) {
       my $x=($alook{'CB'}->{xcoor}+$alook{'CG2'}->{xcoor})/2.0;
       my $y=($alook{'CB'}->{ycoor}+$alook{'CG2'}->{ycoor})/2.0;
       my $z=($alook{'CB'}->{zcoor}+$alook{'CG2'}->{zcoor})/2.0;
       my $arec={ atominx => $natom++, atomname => "SC2", 
                  resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                  xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
       push(@{$nc->{atom}},$arec);
     }
      
      if ($r->{name} eq "LEU" && exists $alook{'CG'} && exists $alook{'CD1'}) {
	my $x=($alook{'CG'}->{xcoor}+$alook{'CD1'}->{xcoor})/2.0;
	my $y=($alook{'CG'}->{ycoor}+$alook{'CD1'}->{ycoor})/2.0;
	my $z=($alook{'CG'}->{zcoor}+$alook{'CD1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "LEU" && exists $alook{'CG'} && exists $alook{'CD2'}) {
        my $x=($alook{'CG'}->{xcoor}+$alook{'CD2'}->{xcoor})/2.0;
        my $y=($alook{'CG'}->{ycoor}+$alook{'CD2'}->{ycoor})/2.0;
        my $z=($alook{'CG'}->{zcoor}+$alook{'CD2'}->{zcoor})/2.0;
        my $arec={ atominx => $natom++, atomname => "SC3",
                   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
        push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ILE") &&
          exists $alook{'CG2'}) {
        my $x=$alook{'CG2'}->{xcoor};
        my $y=$alook{'CG2'}->{ycoor};
        my $z=$alook{'CG2'}->{zcoor};
        my $arec={ atominx => $natom++, atomname => "SC1",
                   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
        push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ILE" && exists $alook{'CD'} && exists $alook{'CG1'}) {
        my $x=($alook{'CD'}->{xcoor}+$alook{'CG1'}->{xcoor})/2.0;
        my $y=($alook{'CD'}->{ycoor}+$alook{'CG1'}->{ycoor})/2.0;
        my $z=($alook{'CD'}->{zcoor}+$alook{'CG1'}->{zcoor})/2.0;
        my $arec={ atominx => $natom++, atomname => "SC2",
                   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
        push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ILE" && exists $alook{'CD1'} && exists $alook{'CG1'}) {
        my $x=($alook{'CD1'}->{xcoor}+$alook{'CG1'}->{xcoor})/2.0;
        my $y=($alook{'CD1'}->{ycoor}+$alook{'CG1'}->{ycoor})/2.0;
        my $z=($alook{'CD1'}->{zcoor}+$alook{'CG1'}->{zcoor})/2.0;
        my $arec={ atominx => $natom++, atomname => "SC2",
                   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
                   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
        push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "PRO" && exists $alook{'CB'} && exists $alook{'CG'} && exists $alook{'CD'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'CG'}->{xcoor}+$alook{'CD'}->{xcoor})/3.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'CG'}->{ycoor}+$alook{'CD'}->{ycoor})/3.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'CG'}->{zcoor}+$alook{'CD'}->{zcoor})/3.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "CYS" && exists $alook{'CB'} && exists $alook{'SG'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'SG'}->{xcoor})/2.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'SG'}->{ycoor})/2.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'SG'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || 
	   $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" ||
           $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" ) &&
	  exists $alook{'CB'} && exists $alook{'CG'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'CG'}->{xcoor})/2.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'CG'}->{ycoor})/2.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'CG'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "MET" && exists $alook{'SD'}) {
	my $x=$alook{'SD'}->{xcoor};
	my $y=$alook{'SD'}->{ycoor};
	my $z=$alook{'SD'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "MET" && exists $alook{'CE'}) {
	my $x=$alook{'CE'}->{xcoor};
	my $y=$alook{'CE'}->{ycoor};
	my $z=$alook{'CE'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }


      if (($r->{name} eq "PHE" || $r->{name} eq "TYR") && exists $alook{'CE1'}) {
	my $x=$alook{'CE1'}->{xcoor};
	my $y=$alook{'CE1'}->{ycoor};
	my $z=$alook{'CE1'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "PHE" || $r->{name} eq "TYR") && exists $alook{'CE2'}) {
	my $x=$alook{'CE2'}->{xcoor};
	my $y=$alook{'CE2'}->{ycoor};
	my $z=$alook{'CE2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "TRP" && exists $alook{'NE1'}) {
	my $x=$alook{'NE1'}->{xcoor};
	my $y=$alook{'NE1'}->{ycoor};
	my $z=$alook{'NE1'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "TRP" && exists $alook{'CE3'}) {
	my $x=$alook{'CE3'}->{xcoor};
	my $y=$alook{'CE3'}->{ycoor};
	my $z=$alook{'CE3'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "TRP" && exists $alook{'CZ2'}) {
	my $x=$alook{'CZ2'}->{xcoor};
	my $y=$alook{'CZ2'}->{ycoor};
	my $z=$alook{'CZ2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }
      
      if ($r->{name} eq "SER" && exists $alook{'CB'} && exists $alook{'OG'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'OG'}->{xcoor})/2.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'OG'}->{ycoor})/2.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'OG'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "THR" && exists $alook{'CB'} && exists $alook{'OG1'}) {
	my $x=($alook{'CB'}->{xcoor}+$alook{'OG1'}->{xcoor})/2.0;
	my $y=($alook{'CB'}->{ycoor}+$alook{'OG1'}->{ycoor})/2.0;
	my $z=($alook{'CB'}->{zcoor}+$alook{'OG1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }
      
      if ($r->{name} eq "THR" && exists $alook{'CG2'}) {
	my $x=$alook{'CG2'}->{xcoor};
	my $y=$alook{'CG2'}->{ycoor};
	my $z=$alook{'CG2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "TYR" && exists $alook{'OH'}) {
	my $x=$alook{'OH'}->{xcoor};
	my $y=$alook{'OH'}->{ycoor};
	my $z=$alook{'OH'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ASN" && exists $alook{'CG'} && exists $alook{'OD1'}) {
	my $x=($alook{'CG'}->{xcoor}+$alook{'OD1'}->{xcoor})/2.0;
	my $y=($alook{'CG'}->{ycoor}+$alook{'OD1'}->{ycoor})/2.0;
	my $z=($alook{'CG'}->{zcoor}+$alook{'OD1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ASN" && exists $alook{'ND2'}) {
	my $x=$alook{'ND2'}->{xcoor};
	my $y=$alook{'ND2'}->{ycoor};
	my $z=$alook{'ND2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "GLN" && exists $alook{'CD'} && exists $alook{'OE1'}) {
	my $x=($alook{'CD'}->{xcoor}+$alook{'OE1'}->{xcoor})/2.0;
	my $y=($alook{'CD'}->{ycoor}+$alook{'OE1'}->{ycoor})/2.0;
	my $z=($alook{'CD'}->{zcoor}+$alook{'OE1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "GLN" && exists $alook{'NE2'}) {
	my $x=$alook{'NE2'}->{xcoor};
	my $y=$alook{'NE2'}->{ycoor};
	my $z=$alook{'NE2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }
      
      if ($r->{name} eq "LYS" && exists $alook{'CD'}) {
	my $x=$alook{'CD'}->{xcoor};
	my $y=$alook{'CD'}->{ycoor};
	my $z=$alook{'CD'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "LYS" && exists $alook{'CE'}) {
	my $x=$alook{'CE'}->{xcoor};
	my $y=$alook{'CE'}->{ycoor};
	my $z=$alook{'CE'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }


      if ($r->{name} eq "LYS" && exists $alook{'NZ'}) {
	my $x=$alook{'NZ'}->{xcoor};
	my $y=$alook{'NZ'}->{ycoor};
	my $z=$alook{'NZ'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }


      if ($r->{name} eq "ARG" && exists $alook{'CD'}) {
	my $x=$alook{'CD'}->{xcoor};
	my $y=$alook{'CD'}->{ycoor};
	my $z=$alook{'CD'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ARG" && exists $alook{'NE'} && exists $alook{'CZ'}) {
	my $x=$alook{'NE'}->{xcoor};
	my $y=$alook{'NE'}->{ycoor};
	my $z=$alook{'NE'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ARG" && exists $alook{'NH1'}) {
	my $x=$alook{'NH1'}->{xcoor};
	my $y=$alook{'NH1'}->{ycoor};
	my $z=$alook{'NH1'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ARG" && exists $alook{'NH2'}) {
	my $x=$alook{'NH2'}->{xcoor};
	my $y=$alook{'NH2'}->{ycoor};
	my $z=$alook{'NH2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "HSD" || $r->{name} eq "HSE" || 
	   $r->{name} eq "HIS" || $r->{name} eq "HSP") && exists $alook{'ND1'}) {
	my $x=$alook{'ND1'}->{xcoor};
	my $y=$alook{'ND1'}->{ycoor};
	my $z=$alook{'ND1'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "HSD" || $r->{name} eq "HSE" || 
	   $r->{name} eq "HIS" || $r->{name} eq "HSP") && exists $alook{'NE2'}) {
	my $x=$alook{'NE2'}->{xcoor};
	my $y=$alook{'NE2'}->{ycoor};
	my $z=$alook{'NE2'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ASP" && exists $alook{'CG'} && exists $alook{'OD1'}) {
	my $x=($alook{'CG'}->{xcoor}+$alook{'OD1'}->{xcoor})/2.0;
	my $y=($alook{'CG'}->{ycoor}+$alook{'OD1'}->{ycoor})/2.0;
	my $z=($alook{'CG'}->{zcoor}+$alook{'OD1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "ASP" && exists $alook{'CG'} && exists $alook{'OD2'}) {
	my $x=($alook{'CG'}->{xcoor}+$alook{'OD2'}->{xcoor})/2.0;
	my $y=($alook{'CG'}->{ycoor}+$alook{'OD2'}->{ycoor})/2.0;
	my $z=($alook{'CG'}->{zcoor}+$alook{'OD2'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "GLU" && exists $alook{'CD'} && exists $alook{'OE1'}) {
	my $x=($alook{'CD'}->{xcoor}+$alook{'OE1'}->{xcoor})/2.0;
	my $y=($alook{'CD'}->{ycoor}+$alook{'OE1'}->{ycoor})/2.0;
	my $z=($alook{'CD'}->{zcoor}+$alook{'OE1'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if ($r->{name} eq "GLU" && exists $alook{'CD'} && exists $alook{'OE2'}) {
	my $x=($alook{'CD'}->{xcoor}+$alook{'OE2'}->{xcoor})/2.0;
	my $y=($alook{'CD'}->{ycoor}+$alook{'OE2'}->{ycoor})/2.0;
	my $z=($alook{'CD'}->{zcoor}+$alook{'OE2'}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "SC3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{"NT"}) {
	my $x=$alook{"NT"}->{xcoor};
	my $y=$alook{"NT"}->{ycoor};
	my $z=$alook{"NT"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "NTT", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ALA" || $r->{name} eq "ASN" || $r->{name} eq "ASP" || $r->{name} eq "VAL" || $r->{name} eq "ILE" || $r->{name} eq "LEU" || $r->{name} eq "PRO" || $r->{name} eq "CYS" || $r->{name} eq "PHE" || $r->{name} eq "TYR" || $r->{name} eq "TRP" || $r->{name} eq "GLN" || $r->{name} eq "GLU" || $r->{name} eq "LYS" || $r->{name} eq "ARG" || $r->{name} eq "HSD" || $r->{name} eq "HSP" || $r->{name} eq "HIS" || $r->{name} eq "HSE" || $r->{name} eq "MET" || $r->{name} eq "THR" || $r->{name} eq "SER" || $r->{name} eq "GLY") && exists $alook{"CAT"}) {
	my $x=$alook{"CAT"}->{xcoor};
	my $y=$alook{"CAT"}->{ycoor};
	my $z=$alook{"CAT"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "SCT", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

## section for nucleic acids

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{'O1P'}) {
	my $x=$alook{'O1P'}->{xcoor};
	my $y=$alook{'O1P'}->{ycoor};
	my $z=$alook{'O1P'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BB8", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{'O2P'}) {
	my $x=$alook{'O2P'}->{xcoor};
	my $y=$alook{'O2P'}->{ycoor};
	my $z=$alook{'O2P'}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BB7", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"O5'"} && exists $alook{"C5'"}) {
	my $x=($alook{"O5'"}->{xcoor}+$alook{"C5'"}->{xcoor})/2.0;
	my $y=($alook{"O5'"}->{ycoor}+$alook{"C5'"}->{ycoor})/2.0;
	my $z=($alook{"O5'"}->{zcoor}+$alook{"C5'"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BB6", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }


## for 5' terminal  

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C5'"} && exists $alook{"O5'"} && exists $alook{"O1P"} && exists $alook{"O2P"}) {
# do nothing
      }
      elsif (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C5'"} && exists $alook{"O5'"}) {
	my $x=$alook{"O5'"}->{xcoor};
	my $y=$alook{"O5'"}->{ycoor};
	my $z=$alook{"O5'"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BT5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"O4'"}) {
	my $x=$alook{"O4'"}->{xcoor};
	my $y=$alook{"O4'"}->{ycoor};
	my $z=$alook{"O4'"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BB5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C4'"} && exists $alook{"C3'"}) {
	my $x=($alook{"C4'"}->{xcoor}+$alook{"C3'"}->{xcoor})/2.0;
	my $y=($alook{"C4'"}->{ycoor}+$alook{"C3'"}->{ycoor})/2.0;
	my $z=($alook{"C4'"}->{zcoor}+$alook{"C3'"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BB4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

## for RNA O2'-C2' else for DNA C2'

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C2'"} && exists $alook{"O2'"}) {
	my $x=($alook{"C2'"}->{xcoor}+$alook{"O2'"}->{xcoor})/2.0;
	my $y=($alook{"C2'"}->{ycoor}+$alook{"O2'"}->{ycoor})/2.0;
	my $z=($alook{"C2'"}->{zcoor}+$alook{"O2'"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BR2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }
      elsif (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"C2'"}) {
	my $x=$alook{"C2'"}->{xcoor};
	my $y=$alook{"C2'"}->{ycoor};
	my $z=$alook{"C2'"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BD2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA") && 
	  exists $alook{"C1'"} && exists $alook{"N9"}) {
	my $x=($alook{"C1'"}->{xcoor}+$alook{"N9"}->{xcoor})/2.0;
	my $y=($alook{"C1'"}->{ycoor}+$alook{"N9"}->{ycoor})/2.0;
	my $z=($alook{"C1'"}->{zcoor}+$alook{"N9"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BB1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") &&           exists $alook{"C1'"} && exists $alook{"N1"}) {
	my $x=($alook{"C1'"}->{xcoor}+$alook{"N1"}->{xcoor})/2.0;
	my $y=($alook{"C1'"}->{ycoor}+$alook{"N1"}->{ycoor})/2.0;
	my $z=($alook{"C1'"}->{zcoor}+$alook{"N1"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BB1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") &&           exists $alook{"C2"} && exists $alook{"O2"}) {
	my $x=($alook{"C2"}->{xcoor}+$alook{"O2"}->{xcoor})/2.0;
	my $y=($alook{"C2"}->{ycoor}+$alook{"O2"}->{ycoor})/2.0;
	my $z=($alook{"C2"}->{zcoor}+$alook{"O2"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BS1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") &&           exists $alook{"N3"}) {
	my $x=$alook{"N3"}->{xcoor};
	my $y=$alook{"N3"}->{ycoor};
	my $z=$alook{"N3"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "THY" || $r->{name} eq "URA") && exists $alook{"C4"} && exists $alook{"O4"}) {
	my $x=($alook{"C4"}->{xcoor}+$alook{"O4"}->{xcoor})/2.0;
	my $y=($alook{"C4"}->{ycoor}+$alook{"O4"}->{ycoor})/2.0;
	my $z=($alook{"C4"}->{zcoor}+$alook{"O4"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BS3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT") && exists $alook{"C4"} && exists $alook{"N4"}) {
	my $x=$alook{"N4"}->{xcoor};
	my $y=$alook{"N4"}->{ycoor};
	my $z=$alook{"N4"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "THY") && exists $alook{"C5M"}) {
	my $x=$alook{"C5M"}->{xcoor};
	my $y=$alook{"C5M"}->{ycoor};
	my $z=$alook{"C5M"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "CYT" || $r->{name} eq "URA") && exists $alook{"C5"} && exists $alook{"C6"}) {
	my $x=($alook{"C5"}->{xcoor}+$alook{"C6"}->{xcoor})/2.0;
	my $y=($alook{"C5"}->{ycoor}+$alook{"C6"}->{ycoor})/2.0;
	my $z=($alook{"C5"}->{zcoor}+$alook{"C6"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BS4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "THY") && exists $alook{"C6"}) {
	my $x=$alook{"C6"}->{xcoor};
	my $y=$alook{"C6"}->{ycoor};
	my $z=$alook{"C6"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA") && exists $alook{"N3"}) {
	my $x=$alook{"N3"}->{xcoor};
	my $y=$alook{"N3"}->{ycoor};
	my $z=$alook{"N3"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS1", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE") && exists $alook{"N1"}) {
	my $x=$alook{"N1"}->{xcoor};
	my $y=$alook{"N1"}->{ycoor};
	my $z=$alook{"N1"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "GUA") && exists $alook{"N2"}) {
	my $x=$alook{"N2"}->{xcoor};
	my $y=$alook{"N2"}->{ycoor};
	my $z=$alook{"N2"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS2", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE") && exists $alook{"N6"}) {
	my $x=$alook{"N6"}->{xcoor};
	my $y=$alook{"N6"}->{ycoor};
	my $z=$alook{"N6"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "GUA") && exists $alook{"N1"}) {
	my $x=$alook{"N1"}->{xcoor};
	my $y=$alook{"N1"}->{ycoor};
	my $z=$alook{"N1"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE") && exists $alook{"N7"}) {
	my $x=$alook{"N7"}->{xcoor};
	my $y=$alook{"N7"}->{ycoor};
	my $z=$alook{"N7"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "GUA") && exists $alook{"C6"} && exists $alook{"O6"}) {
	my $x=($alook{"C6"}->{xcoor}+$alook{"O6"}->{xcoor})/2.0;
	my $y=($alook{"C6"}->{ycoor}+$alook{"O6"}->{ycoor})/2.0;
	my $z=($alook{"C6"}->{zcoor}+$alook{"O6"}->{zcoor})/2.0;
	my $arec={ atominx => $natom++, atomname => "BS4", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "GUA") && exists $alook{"N7"}) {
	my $x=$alook{"N7"}->{xcoor};
	my $y=$alook{"N7"}->{ycoor};
	my $z=$alook{"N7"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BS5", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      if (($r->{name} eq "ADE" || $r->{name} eq "GUA" || $r->{name} eq "CYT" || $r->{name} eq "THY" || $r->{name} eq "URA") && 
	  exists $alook{"O3'"}) {
	my $x=$alook{"O3'"}->{xcoor};
	my $y=$alook{"O3'"}->{ycoor};
	my $z=$alook{"O3'"}->{zcoor};
	my $arec={ atominx => $natom++, atomname => "BB3", 
		   resname => $r->{name}, resnum  => $r->{num}, chain=> $r->{chain},
		   xcoor   => $x, ycoor => $y, zcoor => $z, valid=>1 };
	push(@{$nc->{atom}},$arec);
      }

      $rrec->{end}=$#{$nc->{atom}};
    }
  }

  foreach my $c ( @{$n->{chain}} ) {
    foreach my $r ( @{$c->{res}} ) {
      $r->{name}.="2" if ($r->{name}=~/^[A-Z][A-Z][A-Z]$/);
    }
    foreach my $a ( @{$c->{atom}} ) {
      $a->{resname}.="2" if ($a->{resname}=~/^[A-Z][A-Z][A-Z]$/);
    }
  }

  $n->_coorCache();

  return $n;
}


## function: $sel = parseSelection($selection) 
## parses selection data structure and returns 
## results in hash data structure with the following
## elements: chain, segment, residue, atom, mode, resrange

#  selection := {^}{chain|segment}{:}{residues}{.}{atoms}[,...]
#
#  chain(s)   := [A-Z]{+[A-Z]}
#  segment(s) := name{+name}
#
#  residue(s) := [0-9]{-[0-9]}
#                name
#                peptide|protein|basic|hydrophobic|acidic|charged
#                nucleic|purine|pyrimidine
#                water|hetero|metal
#                @sequence
#
#  atom(s)    := [0-9]+{-[0-9]}=...
#                name
#                backbone|sidechain
#                sugar|base
#                heavy
#                hydrogen|oxygen|carbon|sulfur|phosphorus
#                all
#                [A-Z]*

sub parseSelection {
  my $self=shift;
  my $arg=shift;

  return undef if (!defined $arg);

  my $sel=();

  foreach my $s ( split( /_/,$arg ) ) {
    my $trec={};
    $trec->{resrange}=();
    $trec->{search}=();
    $trec->{mode}="include";
    
    if (substr($s,0,1) eq "!") {
      $trec->{mode}="exclude";
      $s=substr($s,1);
    }
    
    $trec->{chain}=();
    $trec->{segment}=();

    my $rest;
    if ($s=~/:/) {
      my @fc=split(/:/,$s);
      foreach my $cs ( split(/\+/,$fc[0]) ) {
        if ($cs=~/^seg=(.+)$/) {
          push(@{$trec->{segment}},$1);
        } elsif ($cs=~/^chain=(.+)$/) {
          push(@{$trec->{chain}},$1);
        } elsif ($cs=~/^[A-Z0-9a-z\_\-\=]$/) {
	  push(@{$trec->{chain}},$cs);
	} else {
	  push(@{$trec->{segment}},$cs);
	}
      }
      $rest=$fc[1];
    } else {
      $rest=$s;
    }


    my @reslist;
    my @atomlist;
    
    my @fr=split(/\./,$rest);
    if ($#fr>0) {
      @reslist=split(/\+/,$fr[0]);
      @atomlist=split(/\+/,$fr[1]);
    } else {
      @atomlist=split(/\+/,$rest);
      @reslist=();
    }

    $trec->{atom}=();
    foreach my $atom ( @atomlist ) {
      my $ttrec=();
      my @fatom=split(/\//,$atom);
      for (my $itat=0; $itat<=$#fatom; $itat++) {
	my $tat=$fatom[$itat];
        my $negate=0;
        if ($tat=~/\^(.*)/) {
	  $tat=$1;
          $negate=1;
        }
	my @tatom=();
	if ($tat eq "backbone") {
	  # peptide backbone atoms
	  push(@tatom,qw(C N O CA HA HA1 HA2 HN1 HN OT1 OT2 OXT HT1 HT2 HT3));            
	  # nucleic acid backbone atoms
	  push(@tatom,qw(P O1P O2P O5' O5\* O3' O3\* C3' C3\* H3' H3\* C4' C4\* H4' H4\* C5' C5\* H5\* H5' H5'' H3T H5T));
	  push(@reslist,qw(protein nucleic)) if ($#reslist<0);
	} elsif ($tat eq "sidechain") {
	  # protein side chain atoms
	  push(@tatom,qw(CB CD CD1 CD2 CE CE1 CE2 CE3 CG CG1 CG2 CH2 CZ CZ2 CZ3 HB HB1 HB2 HB3 HD1 HD11 HD12 HD13 HD2 HD21 HD22 HD23 HD3 HE HE1 HE2 HE21 HE22 HE3 HG HG1 HG11 HG12 HG13 HG2 HG21 HG22 HG23 HH HH11 HH12 HH2 HH21 HH22 HZ HZ1 HZ2 HZ3 ND1 ND2 NE NE1 NE2 NH1 NH2 NZ OD1 OD2 OE1 OE2 OG OG1 OH SD SG));
	  push(@reslist,"protein") if ($#reslist<0);
	} elsif ($tat eq "sugar") {
	  # nucleic acid sugar atoms
	  push(@tatom,qw(C1' C1\* O4' O4\* H1' H1\* C2' C2\* H2' H2'' H2\* C3' C3\* H3' H3\* C4' C4\* H4' H4\* O3\* O3'));
	  push(@reslist,"nucleic") if ($#reslist<0);
	} elsif ($tat eq "base") {
	  # nucleic acid base atoms
	  push(@tatom,qw(C2 C4 C5 C5M C6 C8 H1 H2 H21 H22 H3 H41 H42 H5 H51 H52 H53 H6 H61 H62 H8 N1 N2 N3 N4 N6 N7 N9 O2 O4 O6));
	  push(@reslist,"nucleic") if ($#reslist<0);
	} elsif ($tat eq "heavy") {
	  push(@tatom,grep(!/^[H0-9].*/,keys %{$self->{have}->{atom}}));
	} elsif ($tat eq "hydrogen") {
	  push(@tatom,grep(/^[0-9]*H.*/,keys %{$self->{have}->{atom}}));
	} elsif ($tat eq "oxygen") {
	  push(@tatom,grep(/^[0-9]*O.*/,keys %{$self->{have}->{atom}}));
	} elsif ($tat eq "nitrogen") {
	  push(@tatom,grep(/^[0-9]*N.*/,keys %{$self->{have}->{atom}}));
	} elsif ($tat eq "carbon") {
	  push(@tatom,grep(/^[0-9]*C.*/,keys %{$self->{have}->{atom}}));
	} elsif ($tat eq "sulfur") {
	  push(@tatom,grep(/^[0-9]*S.*/,keys %{$self->{have}->{atom}}));
	} elsif ($tat eq "phosphorus") {
	  push(@tatom,grep(/^[0-9]*P.*/,keys %{$self->{have}->{atom}}));
	} elsif (exists $self->{have}->{atom}->{$tat}) {
	  push(@tatom,$tat);
        } elsif ($tat eq "all") {
	  push(@tatom,keys %{$self->{have}->{atom}});
	} else {
	  push(@reslist,$atom);
	  $itat=$#fatom+1;
	}
	if ($#tatom>=0) {
	  my $tttrec={};
	  $tttrec->{list}=();
          @{$tttrec->{list}}=@tatom;
	  $tttrec->{pattern}=&_makePattern($tttrec->{list});
	  $tttrec->{negate}=$negate;
	  push(@{$ttrec},$tttrec);
	}
      }
      push(@{$trec->{atom}},$ttrec) if ($#{$ttrec}>=0);
    }

    $trec->{residue}=();
    foreach my $ares ( @reslist ) {
      my $ttrec=();
      my @fres=split(/\//,$ares);
      for (my $itres=0; $itres<=$#fres; $itres++) {
	my $res=$fres[$itres];
	my $negate=0;
        if ($res=~/\^(.*)/) {
	  $res=$1;
          $negate=1;
        }
	my @tres=();
        if ($res eq "solute") {
          $res="solvent";
          $negate=!$negate;
        }
	if ($res eq "peptide" || $res eq "protein") {
	  push(@tres,qw(ALA CYS VAL LEU ILE ASP GLU GLY GLN ASN HSD HSE HSP HIE HID HIS PRO TRP MET SER THR PHE TYR LYS ARG));
	} elsif ($res eq "basic") {
	  push(@tres,qw(ARG LYS));
	} elsif ($res eq "acidic") {
	  push(@tres,qw(ASP GLU));
	} elsif ($res eq "charged") {
	  push(@tres,qw(ASP GLU ARG LYS));
	} elsif ($res eq "hydrophobic") {
	  push(@tres,qw(ALA VAL LEU ILE PHE GLY PRO CYS MET TRP));
	} elsif ($res eq "polar") {
	  push(@tres,qw(SER THR TYR ASN GLN));
	} elsif ($res eq "nucleic") {
	  push(@tres,qw(ADE THY CYT GUA URA A T C G U));
	} elsif ($res eq "purine") {
	  push(@tres,qw(ADE GUA A G));
	} elsif ($res eq "pyrimidine") {
	  push(@tres,qw(THY CYT URA T C U));
	} elsif ($res eq "water") {
	  push(@tres,qw(TIP3 TIP HOH SPC SPCE TIP4 TIP5));
	} elsif ($res eq "metal") {
	  push(@tres,qw(ZN FE NI MN CU CO CA BE));
	} elsif ($res eq "ion") {
	  push(@tres,qw(MG NA CL K SOD CLA CLM NAP POT));
        } elsif ($res eq "solvent") {
	  push(@tres,qw(MG NA CL K SOD CLA CLM NAP POT));
	  push(@tres,qw(TIP3 TIP HOH SPC SPCE TIP4 TIP5));
	} elsif ($res =~/^[\-0-9=]+$/) {
          my $tt={};
          $tt->{negate}=$negate;
          $tt->{range}=();
          foreach my $r ( split(/=/,$res) ) {
            my $ttrange={};
            if ($r =~/^([\-0-9]+)-([\-0-9]+)$/) {
  	      $ttrange->{from}=$1;
	      $ttrange->{to}=$2;
	    } elsif ($r =~/^([\-0-9]+)$/) {
  	      $ttrange->{from}=$1;
	      $ttrange->{to}=$1;
            } 
            if (exists $ttrange->{from}) {
              push(@{$tt->{range}},$ttrange);
            }
          }
          push(@{$ttrec},$tt) if ($#{$tt->{range}}>=0);
        } elsif ($res =~ /^@([A-Z]+)/ ) {
          my $tt={};
          $tt->{search}=$1;
          $tt->{negate}=$negate;
          push(@{$ttrec},$tt);
	} elsif (exists $self->{have}->{residue}->{$res}) {
	  push(@tres,$res);
	} else {
	  printf STDERR "parseSelection: cannot find residue or atom %s in current structure\n",
	    $res;
	}   
	if ($#tres>=0) {
	  my $tt={};
	  $tt->{list}=();
	  @{$tt->{list}}=@tres;
	  $tt->{pattern}=&_makePattern($tt->{list});
	  $tt->{negate}=$negate;
	  push(@{$ttrec},$tt);
	}
      }
      push(@{$trec->{residue}},$ttrec) if ($#{$ttrec}>=0);
    }
    push(@{$sel},$trec);
  }

  if (0) {
    for (my $is=0; $is<=$#{$sel}; $is++) {
      my $s=$sel->[$is];
      printf STDERR "selection #%d\n",$is;
      printf STDERR "chain: >%s<\n",&_makePattern($s->{chain});
      printf STDERR "segment: >%s<\n",&_makePattern($s->{segment});

      printf STDERR "residue: %d OR elements\n",$#{$s->{residue}}+1;
      foreach my $a ( @{$s->{residue}} ) {
	printf STDERR "    %d AND elements\n",$#{$a}+1;
	foreach my $as ( @{$a} ) {
	  printf STDERR "   pattern >%s< neg: %d\n",&_makePattern($as->{list}),$as->{negate} 
           if (exists $as->{list} && $#{$as->{list}}>=0);
	  printf STDERR "   search  >%s< neg: %d\n",$as->{search},$as->{negate} 
           if (exists $as->{search});
          if (exists $as->{range}) {
            foreach my $tas ( @{$as->{range}} ) {
	      printf STDERR "   range   >%d-%d< neg: %d\n",
               $tas->{from},$tas->{to},$as->{negate};
            }
          }
	}
      }

      printf STDERR "atom: %d OR elements\n",$#{$s->{atom}}+1;
      foreach my $a ( @{$s->{atom}} ) {
	printf STDERR "    %d AND elements\n",$#{$a}+1;
	foreach my $as ( @{$a} ) {
	  printf STDERR "    >%s< neg: %d\n",&_makePattern($as->{list}),$as->{negate};
	}
      }
    }
  }

  return $sel;
}

## method: setValidSelection(selection)
## sets the <mark>valid</mark> flag to 1
## for residues/atoms according to the given
## selection

sub setValidSelection {
  my $self=shift;
  my $selection=shift;

  my $sel=$self->parseSelection($selection);
 
  my $setval=($sel->[0]->{mode} eq "include")?0:1;
  $self->resetValidResidues($setval,1); 

  foreach my $s ( @{$sel}) {
    my $cpat=&_makePattern($s->{chain});
    my $spat=&_makePattern($s->{segment});
    my $rpat=$s->{residue};
    my $apat=$s->{atom};
    
    my $val=($s->{mode} eq "include")?1:0;
   
    my $found=0;

    foreach my $c ( @{$self->{chain}} ) {
      my $cid=$c->{id};
      if (!defined $cpat || $cid=~/$cpat/) {
        &_searchSequences($s->{residue},$c);
        if (defined $cpat || defined $spat || 
            (defined $rpat && $#{$rpat}>=0) || 
            (defined $apat && $#{$apat}>=0)) { 
	  foreach my $r ( @{$c->{res}} ) {
	    if (!defined $spat || $r->{seg}=~/$spat/) {
	      if (!defined $rpat || $#{$rpat}<0 || &_checkSelection($r->{name},$r->{num},$rpat))  {
		$r->{valid}=$val;
		for (my $ia=$r->{start}; $ia<=$r->{end}; $ia++) {
		  my $a=$c->{atom}->[$ia];
		  if (!defined $apat || $#{$apat}<0 || &_checkSelection($a->{atomname},0,$apat)) {
		    $a->{valid}=$val;
		    $found=1;
		  }
	        }
	      }
	    }
	  }
        }
      }
    }
    printf STDERR "no matching atoms found\n" if (!$found);
  }
}

sub setAtomPropFromList {
  my $self=shift;
  my $tag=shift;
  my $data=shift;

  my $i=0;
  
  foreach my $c ( @{$self->{chain}} ) {
    foreach my $a ( @{$c->{atom}} ) {
      $a->{$tag}=$data->[$i++];
    } 
  }  
}

sub _makePattern {
  my $list=shift;
  
  return undef if (!defined $list || $#{$list}<0);

  my @retarr=();
  foreach my $l ( @{$list} ) {
    push(@retarr,"^$l\$") if (defined $l);
  }
  return join("|",@retarr);
}

sub _checkSelection { 
  my $name=shift;
  my $num=shift;
  my $sel=shift;

  return 1 if (!defined $sel || $#{$sel}<0);
  foreach my $s ( @{$sel} ) {
    return 1 if (!defined $s || $#{$s}<0);
    my $found=1;
    foreach my $ss ( @{$s} ) {
      if (exists $ss->{pattern}) {
        my $ppp=$ss->{pattern};
        if ($ppp ne "") {
          if ($ss->{negate}) {
            $found=0 if ($name=~/$ppp/);
          } else {
            $found=0 if ($name!~/$ppp/);
          }
        }
      } elsif (exists $ss->{range}) {
        my $r=$ss->{range};
        if ($#{$r}>=0) {
          my $rfound=0;
          foreach my $rr ( @{$r} ) {
            $rfound=1 if ($num>=$rr->{from} && $num<=$rr->{to});
          } 
          if ($ss->{negate}) { 
            $found=0 if ($rfound);  
          } else {
            $found=0 if (!$rfound);
          }
        }  
      }
    }
    return 1 if ($found);
  }
 
  return 0;
}

sub _searchSequences {
  my $s=shift;
  my $chain=shift;

  return if ($#{$s}<0);

  my $seq=Sequence::new();
  $seq->readMoleculeChain($chain);
  my $seqstr=$seq->abbrevSeq(1);

  foreach my $ts ( @{$s} ) {
    foreach my $tts ( @{$ts} ) {
      if (exists $tts->{search}) {
        $tts->{range}=();
        my $inx;
        my $last=0;
        do {
          $inx=index($seqstr,$tts->{search},$last);
          if ($inx>=0) {
            my $trec={};
            $trec->{from}=$chain->{res}->[$inx]->{num};
            $trec->{to}=$trec->{from}+length($tts->{search})-1;
            push(@{$tts->{range}},$trec);
          }
          $last=$inx+1;
        } while ($inx>=0 && $last<length($seqstr));  
      } 
    }
  }
}

sub _generateLookupTable {
  my $self=shift;

  $self->{have}->{chain}={};
  $self->{have}->{atom}={};
  $self->{have}->{residue}={};
  $self->{have}->{segment}={};

  foreach my $c ( @{$self->{chain}} ) {
    $self->{have}->{chain}->{$c->{id}}=1;

    foreach my $r ( @{$c->{res}} ) {
      $self->{have}->{residue}->{$r->{name}}=1;
      $self->{have}->{segment}->{$r->{seg}}=1 if ($r->{seg} ne "");
    }
    
    foreach my $a ( @{$c->{atom}} ) {
      $self->{have}->{atom}->{$a->{atomname}}=1;
    }
  }
}

1;
