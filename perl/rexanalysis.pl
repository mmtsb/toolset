#!/usr/bin/env perl

# replica exchange analysis
# 2002, 2004, 2006 Michael Feig, TSRI, MSU
#

sub usage {
  printf STDERR "usage:   rexanalysis.pl [options]\n";
  printf STDERR "options: [-dir workdir]\n";
  printf STDERR "         [-inx from:to] [-step value]\n";
  printf STDERR "         [-byclient clientid]\n";
  printf STDERR "         [-bytemp temp]\n";
  printf STDERR "         [-bycond condindex]\n";
  printf STDERR "         [-order client|temp]\n";
  printf STDERR "         [-apply cmd]\n";
  printf STDERR "         [-function file]\n";
  printf STDERR "         [-nanvalue value]\n";
  printf STDERR "         [-wham prop:fname:from:to:nbins[=...]]\n";
  printf STDERR "         [-whamtemp temp[:temp2...]\n";
  printf STDERR "         [-whamener file]\n";
  printf STDERR "         [-out dcdfile]\n";
  printf STDERR "         [-atoms from:to]  [-nsel Selection]\n";
  printf STDERR "         [-verbose]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use Math::Trig;
use IPC::Open2;

use GenUtil;
use Molecule;
use Ensemble;
use Sequence;
use ReXServer;
use Analyze;

my $dir=".";
my $condinx;
my $temp;
my $clientid;
my $order="temp";

my $from=1;
my $to=999999999;
my $step=1;

my $apply;
my $ffile;

my @wham=();
my @whamtemp=();
my $whamener;

my $nanvalue=99;

my $dcdfile;

my $atomsfrom;
my $atomsto;
my $nsel;
my $verbose=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-wham") {
    shift @ARGV;
    foreach my $w ( split(/=/,shift @ARGV) ) {
      my @fw=split(/:/,$w);
      my $trec={};
      $trec->{name}=$fw[0];
      $trec->{fname}=$fw[1];
      $trec->{min}=$fw[2];
      $trec->{max}=$fw[3];
      $trec->{nbins}=$fw[4];
      push(@wham,$trec);
    }
  } elsif ($ARGV[0] eq "-whamtemp") {
    shift @ARGV;
    push(@whamtemp,split(/:/,shift @ARGV));
  } elsif ($ARGV[0] eq "-order") {
    shift @ARGV;
    $order=shift @ARGV;
  } elsif ($ARGV[0] eq "-whamener") {
    shift @ARGV;
    $whamener=shift @ARGV;
  } elsif ($ARGV[0] eq "-byclient") {
    shift @ARGV;
    $clientid=shift @ARGV;
  } elsif ($ARGV[0] eq "-bytemp") {
    shift @ARGV;
    $temp=shift @ARGV;
  } elsif ($ARGV[0] eq "-bycond") {
    shift @ARGV;
    $condinx=shift @ARGV;
  } elsif ($ARGV[0] eq "-step") {
    shift @ARGV;
    $step=shift @ARGV;
  } elsif ($ARGV[0] eq "-apply") {
    shift @ARGV;
    $apply=shift @ARGV;
  } elsif ($ARGV[0] eq "-nanvalue") {
    shift @ARGV;
    $nanvalue=shift @ARGV;
  } elsif ($ARGV[0] eq "-function") {
    shift @ARGV;
    $ffile=shift @ARGV;
  } elsif ($ARGV[0] eq "-inx") {
    shift @ARGV;
    my @f=split(/:/,shift @ARGV);
    $from=$f[0];
    $to=($#f>0)?$f[1]:$f[0];
  } elsif ($ARGV[0] eq "-out") {
    shift @ARGV;
    $dcdfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-nsel") {
    shift @ARGV;
    $nsel=shift @ARGV;
  } elsif ($ARGV[0] eq "-atoms") {
    shift @ARGV;
    ($atomsfrom,$atomsto)=split(/:/,shift @ARGV);
  } elsif ($ARGV[0] eq "-verbose") {
    shift @ARGV;
    $verbose=1;
  } elsif ($ARGV[0]=~/^-/) {
    printf STDERR "unknown option %s\n",shift @ARGV;
    &usage();
  }
}

require "$ffile" if (defined $ffile);

&start() if (defined $ffile && defined &start);

$dir="." if (!defined $dir);

my $rex=ReXServer->new(0,$dir);
my $condfile=$dir."/rexserver.cond";
$rex->setup($condfile);
$rex->readData();

if ($#wham>=0) {
  printf STDERR "Preparing input files\n";

  my @alltemps=();
  for (my $i=0; $i<=$#{$rex->{cond}}; $i++) {
    push(@alltemps,$rex->{cond}->[$i]->{temp});
  }
  my @sortedtemps=sort { $a<=>$b } @alltemps;

  my @clients=sort @{$rex->{clientid}};

  if (!defined $whamener || !-r $whamener) {
    $whamener="$$.ener";
    open OUT,">$whamener";
    for (my $ir=$from; $ir<=$rex->{trun} && $ir<=$to; $ir+=$step) {    
      my @earr=();
      for (my $i=0; $i<=$#{$rex->{cond}}; $i++) {
        my $cdat=$rex->getClientData($clients[$i],$ir);      
        push(@earr,$cdat->{ener});
      }
      printf OUT "%s\n",join(" ",@earr);
    }
    close OUT;
  }

  my $condinxfile="$$.cond";
  open OUT,">$condinxfile";
  for (my $ir=$from; $ir<=$rex->{trun} && $ir<=$to; $ir+=$step) {    
    my @cinxarr=();
    for (my $i=0; $i<=$#{$rex->{cond}}; $i++) {
      my $cdat=$rex->getClientData($clients[$i],$ir);      
      push(@cinxarr,$cdat->{cond}->{inx}+1);
    }
    printf OUT "%s\n",join(" ",@cinxarr);
  }
  close OUT;

  my @wnames=();

  my $totbins=1;
  foreach my $w ( @wham ) {
    die "cannot read $w->{fname}" if (!-r $w->{fname});

    push(@wnames,$w->{name});

    $w->{data}=();
    open INP,"$w->{fname}";
    while (<INP>) {
      chomp;
      s/^\s+//;
      my $f=();
      @{$f}=split(/\s+/);
      push(@{$w->{data}},$f);
      foreach my $df ( @{$f} ) {
        $w->{min}=$df if ($df<$w->{min});
        $w->{max}=$df if ($df>$w->{max});
      } 
    } 
    close INP;
    $w->{delta}=($w->{max}-$w->{min})/($w->{nbins}-1);
    $totbins*=$w->{nbins};

    printf STDERR "%s: [%f:%f] %d bins\n",$w->{name},$w->{min},$w->{max},$w->{nbins};
  }

  open OUT,">$$.state";
  for (my $iset=0; $iset<=$#{$wham[0]->{data}}; $iset++) {
    my @iarr=();
    for (my $irep=0; $irep<=$#{$wham[0]->{data}->[0]}; $irep++) {
      my $inx=0;
      my $mult=1;
      for (my $iw=0; $iw<=$#wham; $iw++) {
        my $val=$wham[$iw]->{data}->[$iset]->[$irep];
        my $linx=int(($val-$wham[$iw]->{min})/$wham[$iw]->{delta}+0.5);
        $inx+=$mult*$linx;
        $mult*=$wham[$iw]->{nbins};
      }
      push(@iarr,$inx+1);
    }
    printf OUT "%s\n",join(" ",@iarr);
  }
  close OUT;

  printf STDERR "Running ptWHAM\n";

  my $rexwhamexec=&GenUtil::findExecutable("ptwham");
    
  local (*READ, *WRITE);

  push(@whamtemp,$sortedtemps[0]) if ($#whamtemp<0);
    
  foreach my $wt ( @whamtemp ) {
    my $pid=open2(*READ,*WRITE,"$rexwhamexec");

#    <READ>;
    printf WRITE "%d %d %d\n",$#{$rex->{cond}}+1,$#{$rex->{cond}}+1,$#{$wham[0]->{data}}+1;
#    <READ>;
    printf WRITE "%f\n",$wt;
#    <READ>;
    foreach my $t ( @sortedtemps ) {
      printf WRITE "%f\n",$t;
    }
#    <READ>;
    printf WRITE "%s\n","$$.ener";
#    <READ>;
    printf WRITE "%s\n","$$.cond";
#    <READ>;
    printf WRITE "%d\n",$totbins;
#    <READ>;
    printf WRITE "%s\n","$$.state";
    close WRITE;

    my $buffer="";
    while (<READ>) {
      $buffer.=$_;
    } 
    close READ;
    waitpid($pid,0);

    my $outfname=sprintf("%s.%d.pmf",join("_",@wnames),int($wt+0.5));
    open OUT,">$outfname";
    foreach my $l ( split(/\n/,$buffer) ) {
      if ($l=~/PMF/) {
        my $tinx=substr($l,5,7)+0;
        my $tpmf=substr($l,12,8);
        $tpmf=$nanvalue if ($tpmf eq "Infinity" || $tpmf eq "NaN");
        $tpmf+=0.0;
          
        my @tarr=();
        $tinx--;
        for (my $iw=0; $iw<=$#wham; $iw++) {
          my $inx=($tinx%$wham[$iw]->{nbins});
          my $val=$wham[$iw]->{min}+$inx*$wham[$iw]->{delta};
          push(@tarr,sprintf("%1.8f",$val));
          $tinx-=$inx;
          $tinx/=$wham[$iw]->{nbins};
        } 
        push(@tarr,sprintf("%1.5f",$tpmf));
        printf OUT "%s\n",join(" ",@tarr);
      }
    }
    close OUT;
    printf STDERR "results written to %s\n",$outfname;
  }  
  &GenUtil::remove("$$.cond");
  &GenUtil::remove("$$.ener") if (-r "$$.ener");
  &GenUtil::remove("$$.state");
} else {
  die "need to specify command to apply" if (!defined $apply && !defined $ffile && !defined $dcdfile);
  
  if (defined $temp) {
    my $diff=1E99;
    for (my $i=0; $i <= $#{$rex->{cond}}; $i++) {
      if ((! defined $diff) || (abs($rex->{cond}->[$i]->{temp} - $temp) < $diff)) {
	$condinx=$i;
	$diff=abs($rex->{cond}->[$i]->{temp} - $temp);
      }
    }
  }
  
  my $mol=Molecule::new();
  my ($len,$count,$ft,$aux,$NATOM,$natom,$NFILE,$atoms,$dcdheader);
  if ($rex->{par}->{archive}) {
    $mol->readPDB(sprintf("%s/%s/final.pdb",$dir,$rex->{clientid}->[0]));
    
    if (defined $dcdfile){
	my $archive=sprintf("%s/%s/prod.coor.archive",$dir,$rex->{clientid}->[0]);
	($len,$count,$ft,$aux)=&GenUtil::readArchiveHeader($archive);
	$NATOM=$len/24;
	if ((!defined $atomsfrom || !defined $atomsto) && !defined $nsel){
	    $atomsfrom=1;
	    $atomsto=$NATOM;
	}
	elsif ((defined $atomsfrom || defined $atomsto) && defined $nsel){
	    die "Options \"-atoms\" and \"-nsel\" are mutually exlusive!\n";
	}
	elsif(defined $nsel){
	    $atomsfrom=-1;
	    $atomsto=-1;
	}
	$natom=$atomsto-$atomsfrom+1;
	$from=1 if ($from<0 || $from > $rex->{trun}); #Check from value
	$to=$rex->{trun} if ($to>$rex->{trun} || $to<0); #Check to value
	$NFILE=$to-$from+1;

	if (defined $nsel){
	    ($atoms,$natom)=&getSelection(sprintf("%s/%s/final.pdb",$dir,$rex->{clientid}->[0]),$nsel,$NATOM,$atoms,$verbose);
	    print "$natom atoms found\n" if ($verbose);
	}

	$dcdheader=&genDCDheader($NFILE,$natom);
	
	open (DCD, ">$dcdfile");
	binmode DCD;
	print DCD "$dcdheader";
    }
  }

  for (my $ir=$from; $ir<=$rex->{trun} && $ir<=$to; $ir+=$step) {
    my $c=$rex->getClientData(undef,$ir);

    my $cdat;
    my $scid;

    my $output=0;
    for (my $icond=0; $icond<=$#{$rex->{cond}}; $icond++) {
      if ((defined $clientid && $icond==0) || (!defined $condinx && !defined $clientid) || 
	  $icond==$condinx ) {
	
	if (defined $clientid) {
	  $scid=$clientid;
	  $cdat=$c->{$clientid};
	  die "cannot find data for client $clientid" if (!defined $cdat);
	} else { 
	  if ($order eq "temp") {
	    foreach my $cid (@{$rex->{clientid}}) {
	      if ($c->{$cid}->{cond}->{inx} == $icond) {
		$scid=$cid;
		$cdat=$c->{$cid};
	      }
	    }
	  } elsif ($order eq "client") {
	    my @sc=sort @{$rex->{clientid}};
	    $scid=$sc[$icond];
	    $cdat=$c->{$scid};
	  } else {
            die "unknown sort mode $order";
	  } 
	  die "cannot find data for condition $icond" if (!defined $cdat);
	}
    
	my $ener=$cdat->{ener};
	my $temp=$cdat->{cond}->{temp};
    
	my $dataok=1; 
	if ($rex->{par}->{archive}) {
	  my $arfile=sprintf("%s/%s/prod.coor.archive",$dir,$scid);
	  my $data=&GenUtil::readArchiveFile($arfile,$ir);
	  if (defined $data) {
	    my $start=0;
	    foreach my $c ( @{$mol->{chain}} ) {
	      foreach my $a ( @{$c->{atom}} ) {
		$a->{xcoor}=substr($data,$start,8)+0.0;
		$a->{ycoor}=substr($data,$start+8,8)+0.0;
		$a->{zcoor}=substr($data,$start+16,8)+0.0;
		$start+=24;
	      }
	    }
            $mol->_coorCache();
	  } else {
	    $dataok=0;
	  }
	} else {
	  my $pdbname=sprintf("%s/%s/prod/%s/final.pdb",$dir,$scid,&GenUtil::dataDir($ir));
	  if (-r $pdbname || -r "$pdbname.gz") {
	    $mol->readPDB($pdbname);
	  } else {
	    $dataok=0;
	  }
	}
    
	if ($dataok) {
	  if (defined $ffile) {
	    my @res=&analyze($mol);
            if (defined $res[0]) {
	      printf "%d %d %f %s ",$ir,$cdat->{cond}->{inx},$cdat->{cond}->{temp},$scid 
	       unless (!defined $clientid && !defined $condinx);
	      printf "%s ",join(" ",@res);
	      printf "\n" if (defined $clientid || defined $condinx);
              $output=1;
            }
	} elsif (defined $dcdfile) {
	    for (my $xyz=1; $xyz<=3; $xyz++){
	      my $coor.=pack("i",$natom*4);
	      foreach my $c ( @{$mol->{chain}} ) {
	        foreach my $a ( @{$c->{atom}} ) {
		  if ($a->{atominx} >= $atomsfrom && $a->{atominx} <= $atomsto
		    || (defined $atoms->{$a->{atominx}} && $atoms->{$a->{atominx}})){
		    if ($xyz == 1){
		      $coor.=pack("f",$a->{xcoor});  
		    }
		    elsif ($xyz == 2){
		      $coor.=pack("f",$a->{ycoor});
                    }
                    else{
                      $coor.=pack("f",$a->{zcoor}); 
                    }
	          }
	        }
	      }
	      $coor.=pack("i",$natom*4);
	      print DCD "$coor";
	    }
	    print "Processing frame $ir of $to...\n" if ($verbose);
	  } else {
	    die "Option \"-nsel\" not available for analysis\n" if (defined $nsel);
	    printf "%d %d %f %s ",$ir,$cdat->{cond}->{inx},$cdat->{cond}->{temp},$scid 
	      unless (!defined $clientid && !defined $condinx);
	    open OUT,"|$apply";
	    $mol->writePDB(\*OUT);
	    close OUT;
            $output=1;
	  }
	} else {
	  printf "data not found\n";
	}
      }
    }
    printf "\n" if (!defined $clientid && !defined $condinx && $output);
  }
}

close (DCD) if (defined $dcdfile);

&end() if (defined $ffile && defined &end);

exit 0;



####################################################
#                                                  #
# Process PDB file                                 #
#                                                  #
####################################################
sub getSelection {
    my $pdbfile=shift;
    my $sel=shift;
    my $NATOM=shift;
    my $atoms=shift;
    my $verbose=shift;
    my $count=0;
    my $mol=Molecule::new();
    if (defined $sel && !defined $pdbfile){
	die "Please specify a pdb file for selection\n";
    }
    elsif (defined $pdbfile && !defined $sel){
	die "Please specify a selection for the pdb file\n";
    }
    elsif (defined $sel && defined $pdbfile){
	#Both are defined
	print "Processing selection in PDB file...\n" if ($verbose);
	$mol->readPDB($pdbfile);
	$count=0;
	my $c;
	my $a;
	for $c (@{$mol->activeChains()}){
	    foreach $a (@{$c->{atom}}){ 
		$count++;
	    }
	}
	die "The number of atoms in archive file and PDB file do not match\n"
	    if ($NATOM != $count);
	$mol->setValidSelection($sel);
	$mol=$mol->clone(1);
#	my $atoms;
	$count=0;
	for $c (@{$mol->activeChains()}){
	    foreach $a (@{$c->{atom}}){ 
		$atoms->{$a->{atominx}}=1;
		$count++;
	    }
	}
	if (!$count){
	    exit 1;
	}
    }
    $mol=undef;
    return ($atoms,$count);
}

####################################################
#                                                  #
# Generate DCD header information                  #
#                                                  #
####################################################
sub genDCDheader {
    
    my $NFILE=shift;
    my $NATOM=shift;

    my $archead; #Header for archive file
    my $HDR="CORD"; #DCD Descriptor
    my $NPRIV=10; #Step number for the first frame
    my $NSAVC=10; #The step frequency
    my $NSTEP=$NFILE; #Number of steps = Number of structures!!
    my @ICNTRL;
    for (my $c=5; $c<=20; $c++){
	if ($c == 10){
	    $ICNTRL[$c]=2.045472601; #DELTA Timestep
	    #$ICNTRL[$c]=unpack("L",pack("f",$ICNTRL[$c]));
	}
	else{
	    $ICNTRL[$c]=0;
	}
    }

#   Pack DCD header
    my $dcdheadparam="iA4i9fi11";
    my $dcdhead=pack("$dcdheadparam",84,$HDR,$NFILE,$NPRIV,$NSAVC,$NSTEP,$ICNTRL[5],$ICNTRL[6],$ICNTRL[7],$ICNTRL[8],$ICNTRL[9],$ICNTRL[10],$ICNTRL[11],$ICNTRL[12],$ICNTRL[13],$ICNTRL[14],$ICNTRL[15],$ICNTRL[16],$ICNTRL[17],$ICNTRL[18],$ICNTRL[19],$ICNTRL[20],84);
    
#   Pack Titles and Remarks
    my $user=`whoami`;
    chomp $user;
    (my $sec,my $min,my $hour,my $mday,my $mon,my $year,my $wday,my $yday,my $isdst)=localtime(time);
    $mon=$mon+1;
    $year=substr($year+1900,2);
    my $remark1="REMARKS FILENAME=$dcdfile CREATED BY archive2dcd.pl";
    my $remark2="REMARKS DATE: $mon/$mday/$year CREATED BY USER: $user";
    my $date="*  DATE:     $mon/$mday/$year     $hour:$min:$sec   CREATED BY USER: $user";
    my $titleparam="i2A80A80A80i4";
    $dcdhead.=pack("$titleparam", 244,3,$remark1,$remark2,$date,244,4,$NATOM,4);
    return ($dcdhead);
}

####################################################
#                                                  #
# Trim whitespace                                  #
#                                                  #
####################################################
sub trim($){
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}
