# General utilities for internal use
# 
# http://mmtsb.scripps.edu/doc/GenUtil.pm.html
# 2000, Michael Feig, Brooks group, TSRI

package GenUtil;

require 5.004;

use strict;
use Fcntl;
use IO::File;
use IO::Handle;
use IO::Socket;
use IPC::Open2;
use Net::hostent;
use Sys::Hostname;

#use Time::HiRes 'gettimeofday';

use vars qw ( $pi $logfile $scriptloc $fortranheadersize $version );

BEGIN {
  $version="1.5 Feig group (2018)";
  $pi=3.141592654;
  $fortranheadersize=-1;
  $fortranheadersize=$ENV{FORTRANHEADERSIZE} 
    if (exists $ENV{FORTRANHEADERSIZE});

  ($scriptloc=$0)=~s/[^\/]+$//;
  
  my ($vinx)=grep( $ARGV[$_]=~/^-+version$/ || $ARGV[$_]=~/^-v$/ , 0..$#ARGV);
  if (defined $vinx) {
    printf STDERR "MMTSB ToolSet Version %s\n",$version;
    printf STDERR "M. Feig, J. Karanicolas, C. L. Brooks: J.Mol.Graph.Model (2004) 22, 377-395\n";
    exit (1);
    #splice(@ARGV,$vinx,1);
    #push(@ARGV,"-help");
  }
}

#sub htime {
#  my ($s,$usec)=gettimeofday();
#  return $s+$usec/1000000;
#}

## function: $fullpath = findExecutable(name)
## attempts to locate the full path of a executable
## with the given name

sub findExecutable {
  my $name=shift;
  
  return "$ENV{MMTSBDIR}/bin/$name" 
    if (defined $ENV{MMTSBDIR} && 
	-d "$ENV{MMTSBDIR}/bin" && 
	-x "$ENV{MMTSBDIR}/bin/$name");

  my $scriptbin=$scriptloc;
  $scriptbin=~s/perl\/$/bin/;
  return "$scriptbin/$name" 
    if (-d $scriptbin && 
	-x "$scriptbin/$name");

  foreach my $p ( split(/:/,$ENV{PATH}) ) {
    return "$p/$name"
      if (-d $p && -x "$p/$name");
  }

  return undef;
}

## function: $dir = findDataDir(name)
## attempts to locate the MMTSB data directory

sub findDataDir {
  return "$ENV{MMTSBDIR}/data" 
    if (defined $ENV{MMTSBDIR} && 
	-d "$ENV{MMTSBDIR}/data");

  my $scriptdata=$scriptloc;
  $scriptdata=~s/perl\/$/data/;
  return "$scriptdata" 
    if (-d $scriptdata);

  return ".";
}

## function: setLogFile(name)
## sets a log file for debug messages

sub setLogFile {
  $logfile=getAppendFile(shift);
  $logfile->autoflush(1);
}

## function: log(tag,text[,stderr])
## writes a debug message to a log file
## if one has been set previously with <mark>setLogFile</mark>
## The message consists of a tag and the actual text.

sub log {
  my $tag=shift;
  my $s=shift;
  my $stderr=shift;

  my $time_now=localtime;
  my $message;
  if (defined $s) {
    $message=sprintf("## %s ## %s ##  %s",$time_now,$tag,$s);
  } else {
    $message=sprintf("## %s ##  %s",$time_now,$tag);
  }
  $message.="\n" if ($s !~ /\n$/);

  printf $logfile $message
    if (defined $logfile);

  printf STDERR $message
    if (defined $stderr && $stderr);
}

## function: makeDir(directory)
## creates a directory unless it exists
## already

sub makeDir {
  my $dir=shift;
  my @f=split(/\/+/,$dir);

  my $t="";
  while (@f) { 
    $t.=(shift @f)."/";
    system "mkdir $t 2>/dev/null" if (!-d $t);
  }
}

## function: $handle = getInputFile(file)
## returns an input file handle for a given file in 
## a rather flexible manner.
## If a file name is given as argument, the file 
## is opened. Without any argument or if "-" is
## given, the handle for standard input is returned.
## If the file is compressed, <mark>gunzip</mark>
## is called automatically through a pipe to
## provide the uncompressed data.<BR>
## A previously opened file handle can also be 
## given as argument which will then be simply
## returned without further action.

sub getInputFile {
  my $fname=shift;

  $fname="-" if (!defined $fname || $fname eq "");

  if (ref $fname) {
    return $fname;
  } else {
    die "cannot open file $fname" 
      if (!-r $fname && !-r "$fname.gz" && $fname ne "-");

    my $newfile=new IO::File;

    if ((-r "$fname.gz" && !-r "$fname") || $fname=~/\.gz$/) {
      $newfile->open("gunzip -c $fname |");
    } else {
      $newfile->open("$fname");
    }

    return $newfile;
  }
}

## function: $handle = getOutputFile(file)
## returns an output file handle for a given file
## handle in a similar manner as <mark>getInputFile</mark>
## for file input.<BR>
## If a file name is given a file handle is opened
## for reading. No argument or "-" returns standard
## output and if a file handle is given as argument
## it is passed along.

### getOutputFile ######

sub getOutputFile {
  my $fname=shift;

  $fname="-" if (!defined $fname || $fname eq "");
  
  if (ref $fname) {
    return $fname;
  } else {
    my $newfile=new IO::File;
    my $ret=$newfile->open(">$fname");
    
    printf STDERR "cannot open file $fname for writing\n"
      if (!$ret);

    die "THE END" if (!$ret);

    return $newfile;
  }
}

## function: $handle = getAppendFile(file)
## returns an output file handle to append to a given file
## handle in a similar manner as <mark>getOutputFile</mark>.

### getAppendFile ######

sub getAppendFile {
  my $fname=shift;

  $fname="-" if (!defined $fname || $fname eq "");
  
  if (ref $fname) {
    return $fname;
  } else {
    my $newfile=new IO::File;
    my $ret=$newfile->open(">>$fname");
    
    printf STDERR "cannot open file $fname for writing\n"
      if (!$ret);

    die "THE END" if (!$ret);

    return $newfile;
  }
}

## function: $list = fragListFromOption(string)
## converts a fragment list command line argument
## into a list data structure

sub fragListFromOption {
  my $fragoption=shift;
  my $fraglist=();

  die "invalid list specification: >$fragoption<"
    if ($fragoption !~ /^[A-Za-z0-9:\.=_\-\+]+$/);

  my $fdef=undef;
  foreach my $f (split(/=/,$fragoption)) {
    my @fl=split(/_/,$f);
    my @l=split(/:/,$fl[0]);
    my $rec={};

    if ($#l>=2) {
      $rec->{segid}=shift @l;
    }

    my ($tchain,$tfrom)=($l[0]=~/([A-Za-z]*)([\-0-9]*)/);
    $rec->{chain}=$tchain if ($tchain ne "");
    $rec->{from}=$tfrom if ($tfrom ne "");
    if (defined $l[1]) {
      ($rec->{to}=$l[1])=~s/^[A-Za-z]+//;
    } else {
      $rec->{to}=$rec->{from};
    }
    if (defined $fl[1]) {
      $rec->{force}=$fl[1];
      $fdef=$fl[1];
    } else {
      $rec->{force}=$fdef;
    }
    push(@{$fraglist},$rec);
  }
  return $fraglist;
}

## function: $string = fragOptionFromList(list[,force])
## generates a command line option string from
## a fragment list data structure. A force constant may be
## given separately to be included in the output 
## if the force field is not defined in the input list.

sub fragOptionFromList {
  my $fraglist=shift;
  my $force=shift;

  my $lastforce;

  if (defined $fraglist) {
    my @list=();
    foreach my $f (@{$fraglist}) {
      my $str=($f->{from}==$f->{to})?"$f->{from}":"$f->{from}:$f->{to}";
      
      $str=(uc $f->{chain}).$str if (defined $f->{chain} && $f->{chain} ne "");

      my $fval=(defined $force)?
	$force:((defined $f->{force})?$f->{force}:undef);
      
      if (defined $fval) {
	if (!defined $lastforce || $fval ne $lastforce) {
	  $str.="_$fval";
	  $lastforce=$fval;
	}
      }
      $str=$f->{segid}.":".$str if (defined $f->{segid} && $f->{segid} ne "");
      push(@list,$str);
    }
    return join("=",@list); 
  } else {
    return "";
  }
}

## function: $list = gradForceList($list,$force) 
## generates a list with gradually increasing 
## force constants away from the list boundaries

sub gradForceList {
  my $list=shift;
  my $force=shift;

  my $nlist=();

  my @fac=(0.1,0.2,0.4,0.7);

  my $i;
  foreach my $l ( @{$list} ) {
    my $range=$l->{to}-$l->{from}+1;
    for ($i=0; $i<$range/2 && $i<=$#fac; $i++) {
      my $nrec={};
      $nrec->{from}=$l->{from}+$i;
      $nrec->{to}=$nrec->{from};
      $nrec->{chain}=$l->{chain};
      $nrec->{force}=$force*$fac[$i];
      push (@{$nlist},$nrec);
    }

    my $flat=$range-2*($#fac+1);
    if ($flat>0) {
      my $nrec={};
      $nrec->{from}=$l->{from}+$i;
      $nrec->{to}=$nrec->{from}+$flat-1;
      $nrec->{chain}=$l->{chain};
      $nrec->{force}=$force;
      push (@{$nlist},$nrec);
      $i+=$flat;
    }
    
    for (; $i<$range; $i++) {
      my $nrec={};
      $nrec->{from}=$l->{from}+$i;
      $nrec->{to}=$nrec->{from};
      $nrec->{chain}=$l->{chain};
      $nrec->{force}=$force*$fac[$range-$i-1];
      push (@{$nlist},$nrec);
    }
  }

  return $nlist;
}

## function: $list = fragListFromArray(arr[,chain])
## generates a fragment list data structure
## from an array of residues

sub fragListFromArray {
  my $arr=shift;
  my $chain=shift;

  my $list=();
  
  my $trec={};
  $trec->{chain}=$chain;
  $trec->{from}=$arr->[0];
  for (my $i=1; $i<=$#{$arr}; $i++) {
    if ($arr->[$i]-$arr->[$i-1]>1) {
      $trec->{to}=$arr->[$i-1];
      push(@{$list},$trec);
      $trec={};
      $trec->{chain}=$chain;
      $trec->{from}=$arr->[$i];
    }
  }
  $trec->{to}=$arr->[$#{$arr}];
  push(@{$list},$trec);

  return $list;
}

## function: $list = readFormat(handle,reclen,numrec)
## reads a formatted line with fields of a specific record length

sub readFormat {
  my $handle=shift;
  my $reclen=shift;
  my $numrec=shift;

  my @arr=();

 FORMATINP:
  while (<$handle>) {
    for (my $i=0; $i+$reclen<=length($_); $i+=$reclen) {
      push(@arr,substr($_,$i,$reclen));
      last FORMATINP if ($#arr+1>=$numrec);
    }      
  }
  printf STDERR "read only %d elements (%d requested)\n",$#arr+1,$numrec
    if ($#arr+1<$numrec);

  return @arr;
}

## function: $buffer = readData(handle,pattern)
## reads data from the given handle into a buffer until 
## a pattern is found.

sub readData {
  my $handle=shift;
  my $pattern=shift;

  $pattern="" if (!defined $pattern);

  my $buffer="";
  while (<$handle>) {
    return $buffer 
      if ($pattern ne "" && $_=~/$pattern/);
    $buffer.=$_;
  }
  return $buffer;
}

## function: $status = validGzip(filename)
## checks whether a file is a valid compressed gzip file.

sub validGzip {
  my $f=shift;

  my $ret=`gzip -t $f`;

  return ($ret!~/not/);
}

## function: $status = checkFile(filename)
## checks whether a file is available and readable

sub checkFile {
  my $f=shift;
  my $dontcheckzip=shift;
  $dontcheckzip=0 if (!defined $dontcheckzip);

  return ($f eq "-" || (-r $f && !-z $f) || (-r "$f.gz" && ($dontcheckzip || validGzip("$f.gz"))));
}

## function: remove(filename)
## removes a file if it exists

sub remove {
  my $f=shift;
  
  unlink $f if (-r $f);
}

## function: compress(filename)
## compresses a file with gzip if it
## exists

sub compress {
  my $f=shift;

  system "gzip -f $f" if (-r $f);
}

## function: $value = dihedral(coor1,coor2,coor3,coor4)
## calculates the dihedral value for four points given
## by the <mark>coor*</mark> arguments. 

sub dihedral {
  my $c1=shift;
  my $c2=shift;
  my $c3=shift;
  my $c4=shift;

  my $dx12=$c1->{xcoor}-$c2->{xcoor};
  my $dy12=$c1->{ycoor}-$c2->{ycoor};
  my $dz12=$c1->{zcoor}-$c2->{zcoor};

  my $dx23=$c2->{xcoor}-$c3->{xcoor};
  my $dy23=$c2->{ycoor}-$c3->{ycoor};
  my $dz23=$c2->{zcoor}-$c3->{zcoor};

  my $dx43=$c4->{xcoor}-$c3->{xcoor};
  my $dy43=$c4->{ycoor}-$c3->{ycoor};
  my $dz43=$c4->{zcoor}-$c3->{zcoor};

  my $px1=$dy12*$dz23-$dy23*$dz12;
  my $py1=$dz12*$dx23-$dz23*$dx12;
  my $pz1=$dx12*$dy23-$dx23*$dy12;
  
  my $np1=sqrt($px1*$px1+$py1*$py1+$pz1*$pz1);

  if ($np1<=0.00001) {
    printf STDERR "np1 zero: %s:%s:%s:%d:%s %s:%s:%s:%d:%s %s:%s:%s:%d:%s %s:%s:%s:%d:%s\n",
     $c1->{seg},$c1->{chain},$c1->{resname},$c1->{resnum},$c1->{atomname},
     $c2->{seg},$c2->{chain},$c2->{resname},$c2->{resnum},$c2->{atomname},
     $c3->{seg},$c3->{chain},$c3->{resname},$c3->{resnum},$c3->{atomname},
     $c4->{seg},$c4->{chain},$c4->{resname},$c4->{resnum},$c4->{atomname};
     return 0.0;
  }

  $px1/=$np1;
  $py1/=$np1;
  $pz1/=$np1;

  my $px2=$dy43*$dz23-$dy23*$dz43;
  my $py2=$dz43*$dx23-$dz23*$dx43;
  my $pz2=$dx43*$dy23-$dx23*$dy43;
  my $np2=sqrt($px2*$px2+$py2*$py2+$pz2*$pz2);
 
  if ($np2<=0.00001) {
    printf STDERR "np2 zero: %s:%s:%s:%d:%s %s:%s:%s:%d:%s %s:%s:%s:%d:%s %s:%s:%s:%d:%s\n",
     $c1->{seg},$c1->{chain},$c1->{resname},$c1->{resnum},$c1->{atomname},
     $c2->{seg},$c2->{chain},$c2->{resname},$c2->{resnum},$c2->{atomname},
     $c3->{seg},$c3->{chain},$c3->{resname},$c3->{resnum},$c3->{atomname},
     $c4->{seg},$c4->{chain},$c4->{resname},$c4->{resnum},$c4->{atomname};
     return 0.0;
  }

  $px2/=$np2;
  $py2/=$np2;
  $pz2/=$np2;
  
  my $dp12=$px1*$px2+$py1*$py2+$pz1*$pz2;

  my $ts=1.0-$dp12*$dp12;
  $ts=($ts<0.0)?0.0:sqrt($ts);
  my $angle=$pi/2.0-atan2($dp12,$ts);

  my $px3=$py1*$pz2-$py2*$pz1;
  my $py3=$pz1*$px2-$pz2*$px1;
  my $pz3=$px1*$py2-$px2*$py1;

  my $dp233=$px3*$dx23+$py3*$dy23+$pz3*$dz23;
  
  if ($dp233>0.0) {
    $angle=-$angle;
  }

  return $angle/$pi*180.0;
}

## function: $value = distance(coor1,coor2)
## calculates the distance between two atoms
## by the <mark>coor*</mark> arguments. 

sub distance {
  my $c1=shift;
  my $c2=shift;

  my $dx=$c1->{xcoor}-$c2->{xcoor};
  my $dy=$c1->{ycoor}-$c2->{ycoor};
  my $dz=$c1->{zcoor}-$c2->{zcoor};

  my $d=sqrt($dx*$dx+$dy*$dy+$dz*$dz);
  return $d;
}

## function: $value = angle(coor1,coor2,coor3)
## calculates the distance between two atoms
## by the <mark>coor*</mark> arguments. 

sub angle {
  my $c1=shift;
  my $c2=shift;
  my $c3=shift;

  my $dx12=$c1->{xcoor}-$c2->{xcoor};
  my $dy12=$c1->{ycoor}-$c2->{ycoor};
  my $dz12=$c1->{zcoor}-$c2->{zcoor};

  my $dx32=$c3->{xcoor}-$c2->{xcoor};
  my $dy32=$c3->{ycoor}-$c2->{ycoor};
  my $dz32=$c3->{zcoor}-$c2->{zcoor};

  my $d12=sqrt($dx12*$dx12+$dy12*$dy12+$dz12*$dz12);
  my $d32=sqrt($dx32*$dx32+$dy32*$dy32+$dz32*$dz32);
  
  my $d123=$dx12*$dx32+$dy12*$dy32+$dz12*$dz32;
  $d123/=($d12*$d32);
  my $ts=1.0-$d123*$d123;
  $ts=($ts<0.0)?0.0:sqrt($ts);
  my $angle=$pi/2.0-atan2($d123,$ts);

  return $angle/$pi*180.0;
}

## function: $dir = dataDir(index)
## returns the ensemble data directory for a given index

sub dataDir {
  my $n=shift;
  my $a=int($n/100);
  my $b=$n%100;
  return "$a/$b";
}

## function: ($avg,stddev,$n) = average(list)
## calculates average and standard deviation of
## values given in the list argument

sub average {
  my $list=shift;

  my $sum=0.0;
  my $sumsq=0.0;
  my $n=0;

  my $stddev;

  die "list not defined"
    if (!defined $list || $#{$list}<0);

  return ($list->[0]->{val},0.0,1)
    if ($#{$list}==0);
  
  foreach my $s ( @{$list} ) {
    my $val=$s->{val};
    $sum+=$val;
    $sumsq+=$val*$val;
    $n++;
  }
  
  my $avg=($n==0)?0.0:$sum/$n;
  if ($sumsq-$n*($avg*$avg) < 0){$stddev=0.0;}  #Necessary to prevent tolerance errors
  else {$stddev=($n==0)?0.0:sqrt(($sumsq-$n*($avg*$avg))/($n-1));}
  
  return ($avg,$stddev,$n);
}

## function: $list = limRange(list,min,max)
## generates a sublist from the list given
## as argument with values that lie within
## the range given by <mark>min</mark> and
## <mark>max</mark>.

sub limRange {
  my $list=shift;
  my $min=shift;
  my $max=shift;

  my $outlist=();

  foreach my $s ( @{$list} ) {
    my $val=$s->{val};
    push(@{$outlist},$s)
      if ((!defined $min || $val>=$min) && 
	  (!defined $max || $val<=$max));
  }
  
  return $outlist;
}

## function: ($list,$lastlimit,$average,$stddev,$n) = 
## function:  limCore(list[,parameters])
## generates a sublist by excluding data points that
## are further away from the mean than a multiple of the standard
## deviation. Parameters in hash-style key=>value pair form can be given 
## to select whether all values below (<mark>lower</mark>)
## and above (<mark>upper</mark>) are cropped, to change
## the multiple <mark>mult</mark> (default: 3) and to change a 
## tolerance limit <mark>tol</mark> for terminating the iterative 
## solution.<BR>
## Return values are the reduced list, the last limit value used,
## and the average, standard deviation and number of points for the
## reduced list.

sub limCore {
  my $list=shift;

  die "list not defined"
    if (!defined $list || $#{$list}<0);

  return ($list,0.0,$list->[0]->{val},0.0,1)
    if ($#{$list}==0);

  my %par=@_;

  my $lower=(defined $par{lower})?$par{lower}:1;
  my $upper=(defined $par{upper})?$par{upper}:1;
  my $mult=(defined $par{mult})?$par{mult}:3.0;
  my $tol=(defined $par{tol})?$par{tol}:0.05;

  my $avg;
  my $stddev;
  my $limit=1.0E20;
  my $lastlimit;

  my $n; 

  do {
    $lastlimit=$limit;

    $n=0;
    my $sum=0.0;
    my $sumsq=0.0;

    foreach my $vl ( @{$list} ) {
      my $v=$vl->{val};
      if ((!$upper || $v<($avg+$limit)) && 
	  (!$lower || $v>($avg-$limit))) {
	$sum+=$v;
	$sumsq+=$v*$v;
	$n++;
      }
    }
    $avg=($n==0)?0.0:$sum/$n;
    if (($sumsq-$n*($avg*$avg))<0) {$stddev=0.0;} #Necessary because of tolerance errors
    else {$stddev=($n==0)?0.0:sqrt(($sumsq-$n*($avg*$avg))/($n-1));}
    $limit=$mult*$stddev;
  } while($lastlimit-$limit>$tol && $n>0);

  my $outlist=&GenUtil::limRange($list,
     ($lower)?($avg-$lastlimit):undef,($upper)?($avg+$lastlimit):undef);
  
  return ($outlist,$lastlimit,$avg,$stddev,$n);
}

## function: $corr = spRank($xlist,$ylist)
## calculates a Spearman rank coefficient for x and y values 
## given as list arguments. 

sub spRank {
  my $x=shift;
  my $y=shift;

  my $spearbin=&GenUtil::findExecutable("spear");
  die "cannot find spear executable"
    if (!defined $spearbin);

  local (*READ,*WRITE);
  my $pid=open2(*READ,*WRITE,"$spearbin");

  for (my $i=0; $i<=$#{$x}; $i++) {
    printf WRITE "%f %f\n",$x->[$i]->{val},$y->[$i]->{val};
  }
  close WRITE;
  my $ret=<READ>;
  close READ;
  waitpid($pid,0);

  $ret=~s/^ +//;
  my @fret=split(/ +/,$ret);
  
  return $fret[0];
}

## function: ($m,$n,$corr) = linearFit($xlist,$ylist)
## calculates a linear fit for x and y values given as
## list arguments. The slope, intercept and correlation
## coefficient are returned.

sub linearFit {
  my $x=shift;
  my $y=shift;

  my $nval=$#{$x}+1;
  $nval=$#{$y}+1 if ($#{$y}<$#{$x});

  my $m=0.0;
  my $n=0.0;

  my $sx=0.0;
  my $sy=0.0;
  for (my $i=0; $i<$nval; $i++) {
    $sx+=$x->[$i]->{val};
    $sy+=$y->[$i]->{val};
  }
  my $ss=$nval;
  my $sxoss=$sx/$ss;
  
  my $st2=0.0;
  for (my $i=0; $i<$nval; $i++) {
    my $t=$x->[$i]->{val}-$sxoss;
    $st2+=$t*$t;
    $m+=$t*$y->[$i]->{val};
  }
  $m/=$st2;
  $n=($sy-$sx*($m))/$ss;
  my $siga=sqrt((1.0+$sx*$sx/($ss*$st2))/$ss);
  my $sigb=sqrt(1.0/$st2);

  my $chi2=0.0;
  my $vary=0.0;
  for (my $i=0; $i<$nval; $i++) {
    my $tc=($y->[$i]->{val}-$n-$m*$x->[$i]->{val});
    $chi2+=$tc*$tc;
    $tc=($y->[$i]->{val}-$sy/$ss);
    $vary+=$tc*$tc;
  }
      
  my $c=sqrt(1.0-($chi2/$vary));

  return ($m,$n,$c);
}

## function: $status = dataAvailable(handle,timeout)
## checks whether input data is available at the given
## handle.

sub dataAvailable {
  my $handle=shift;
  my $timeout=shift;

  $timeout=0.5 if ($timeout <= 0);
  
  my ($nfound,$rmask,$nread,$line);

  $rmask="";
  vec($rmask,fileno($handle),1)=1;
  ($nfound,$rmask)=select($rmask,undef,undef,$timeout);
    
  return $nfound;
}

## function: ret = safeToSend(filename)
## determines whether a file is safe to be sent 
## to a client from a replica exchange or ensemble server.

sub safeToSend {
  my $fname=shift;

  return 0 if (!checkFile($fname));

  $fname="$fname.gz" if (!-r "$fname" && -r "$fname.gz");

  my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
      $atime,$mtime,$ctime,$blksize,$blocks)=stat($fname);
  
#  return (($mode & 00004) && $> == $uid);
  return ($>==$uid);
}

## function: str = HTMLFont

sub HTMLFont {
  my $size=shift;
  my $color=shift;
  return "<FONT FACE=courier COLOR=$color SIZE=$size>";
}

## function: $port = getServerPort(base)
## finds an available TCP/IP port starting from <mark>base</mark>

sub getServerPort {
  my $baseaddress=shift;
  $baseaddress=4000 if (!defined $baseaddress);

  my $port=$ENV{PSPORT};
  if (!defined $port) {
    for ($port=$baseaddress; &alive("localhost",$port) && $port<$baseaddress+100; $port++) 
      {}
  }

  return $port;
}

## function: $socket = getServerSocket(port)
## obtains a server socket handle for incoming
## connections at the given port address

sub getServerSocket {
  my $port=shift;
  my $server=IO::Socket::INET->new(Proto     => 'tcp',
				   LocalPort => $port,
				   Listen    => SOMAXCONN,
				   Reuse     => 1);
  die "cannot start job server" if (!$server);
  return $server;
}

## function: $num = getServerID()
## generates a random server ID

sub getServerID {
  return int(rand(1000000));
}
  
## function: $status = alive(host,port)
## checks whether a server is alive at the
## given host and port.

sub alive {
  my $hostname=shift;
  my $port=shift;
  
  my $tst=&connectToServer($hostname,$port);
  return 0 if (!defined $tst);

  my $ret=&sendCommand($tst,"TEST");
  return 0 if (!defined $ret || $ret !~ /alive/);

  return 1;
}

## function: $socket = connectToServer(host,port)
## obtains a client socket to connect to the given
## host and port address.
  
sub connectToServer {
  my $hostname=shift;
  my $port=shift;

  my $remoteSocket=IO::Socket::INET->
    new(Proto    => 'tcp',
	PeerAddr => $hostname,
	PeerPort => $port);


  return $remoteSocket;
}

## function: waitForData(socket)
## blocks until data is available on socket

sub waitForData {
  my $socket=shift;
  
  my $mask='';
  vec($mask,fileno($socket),1)=1;
  my $ret=select($mask,undef,undef,undef);
}

## function: readyToWrite(socket)
## blocks until socket is ready for writing

sub readyToWrite {
  my $socket=shift;
  
  my $mask='';
  vec($mask,fileno($socket),1)=1;
  my $ret=select(undef,$mask,undef,undef);
}

## function: $data = readFromSocket(socket,nbytes)
## reads data from a socket

sub readFromSocket {
  my $socket=shift;
  my $expected=shift;
  my $toread=shift;

  my $buffer="";

  do {
    &waitForData($socket);
    my $tbuf='';
    my $nread=(defined $toread)?$toread-length($buffer):100000;
    my $ret=recv $socket,$tbuf,$nread,0;
    return undef if (!defined $ret || length($tbuf)<=0);
    $buffer.=$tbuf;
  } while ( (defined $toread && length($buffer)<$toread) ||
	    (defined $expected && $buffer !~ /$expected/) );

  return $buffer;
}

## function: $data = recvFixedData(socket)
## reads next data item from socket

sub recvFixedData {
  my $socket=shift;
  my $datalength=&readFromSocket($socket,undef,8);
  return undef if (!defined $datalength || $datalength<=0);
  return &readFromSocket($socket,undef,$datalength);
}

## function: $ret = writeToSocket(socket,data)
## writes data to socket

sub writeToSocket {
  my $socket=shift;
  my $data=shift;
  my $towrite=length($data);
  my $havewritten=0;
  
  do {
    &readyToWrite($socket);
    my $nwrite=send $socket,substr($data,$havewritten),0;
    return 0 if (!defined $nwrite || $nwrite<=0);
    $havewritten+=$nwrite;
  } while ($havewritten<$towrite);
  return 1;
}

## function: $ret = sendFixedData(socket,data)
## writes data item to socket

sub sendFixedData {
  my $socket=shift;
  my $data=shift;

  die "data item too large"
    if (length($data)>99999999);

  return &writeToSocket($socket,sprintf("%8d",length($data)).$data);
}

## function: $ret = sendCommand(socket,command[,keepOpen])
## sends a command to a server through the previously opened 
## socket handle. It returns the answer from the server.
## The socket connection is closed unless the <mark>keepOpen</mark>
## flag is set.

sub sendCommand {
  my $remoteSocket=shift;
  my $cmd=shift;
  my $keepOpen=shift;

  return undef if (!defined $remoteSocket);

  fcntl($remoteSocket,F_SETFL,O_NONBLOCK);

  &writeToSocket($remoteSocket,"$cmd\n") || 
    die "ERROR: cannot send >$cmd< to server";

  my $buffer=&readFromSocket($remoteSocket,"[\n\r]\$");
  
  close $remoteSocket unless (defined $keepOpen && $keepOpen && defined $buffer);

  if (!defined $buffer) {
#    print STDERR "ERROR: lost socket connection while sending >$cmd<\n";
    return undef;
  } else {
    chomp $buffer;
    if ($buffer =~/ERROR/) {
      print STDERR $buffer," while sending >$cmd<\n";
      return undef;
    } else {
      return $buffer;
    }
  }
}

## function: $list = readHostFile(file)
## reads a remote host file and returns
## the data as a list.

sub readHostFile {
  my $fname=shift;
  my $list=();

  my %haveHost;

  my $inp=&getInputFile($fname);

  while (<$inp>) {
    chomp;
    s/^ +//;
    my @f=split(/ +/);

    if (defined $haveHost{$f[0]} && !defined $f[3]) {
      $haveHost{$f[0]}->{cpus}+=$f[1];
    } else {
      my $rec={};
      $rec->{name}=$f[0];
      if (defined $f[1]) {
	$rec->{cpus}=($f[1]>0)?$f[1]:1;
	$rec->{localdir}=(defined $f[2])?$f[2]:undef;
        $rec->{gpu}=(defined $f[3])?$f[3]:undef;
      } else {
	$rec->{cpus}=1;
      }
      push(@{$list},$rec);
      $haveHost{$f[0]}=$list->[$#{$list}];
    }
  }

  undef $inp;

  return $list;
}

## function: $maxavail = remoteCPUs(hostrec)
## returns the maximum number of CPUs available on a remote machine

sub remoteCPUs {
  my $hrec=shift;
  return $hrec->{cpus};
}

## function: ($pid,$cpusleft) = 
## function:   submitRemote(hostrec,serverhost,serverport,serverid,
## function:                cpus,directory,command,options)
## submits a remote client job on the host given by <mark>hostrec</mark>.
## On the remote machine the given command is executed in the given 
## directory with the given options. Additional options are set to
## connect to the job server using the <mark>serverhost</mark>, 
## <mark>serverport</mark>, and <mark>serverid</mark> variables
## and request up to the given number of CPUs.<BR>
## Return values are the process ID of the subprocess running <mark>rsh</mark>
## and the number of CPUs left. 

sub submitRemote {
  my $hrec=shift;
  my $srec=shift;
  my $cpus=shift;
  my $cmd=shift;
  my $options=shift;
  my $subdir=shift;
  my $remsubdir=shift;
  my $logfile=shift;
  my $workdir=shift;

  $subdir=0 if (!defined $subdir);
  $remsubdir=1 if (!defined $remsubdir);

  my $host=$hrec->{name};
  my $scpus=($cpus>$hrec->{cpus})?$hrec->{cpus}:$cpus;

  my $dir=$hrec->{localdir};
  $dir="." if (!defined $dir);

  my $rsh=$ENV{'REMOTESHELL'};
  if (!defined $rsh || $rsh eq "") {
    $rsh="ssh";
  }

  my $pid=fork();
  if (!$pid) {
    my $rdir;
    my $sdir;
    if ($subdir) {
      $sdir=$hrec->{name}."-".sprintf("%d",int(rand(1000000)));
      system "$rsh -n $host '( mkdir -p $dir; cd $dir; mkdir -p $sdir )'";
      $rdir="$dir/$sdir";
    } else {
      $rdir=$dir;
    }

    system "$rsh -n $host '( cd $rdir; echo $srec->{name}:$srec->{port}:$srec->{id} > server-$host-$$.info; /bin/chmod 600 server-$host-$$.info )'";

#printf STDERR "running: $cmd -cpus $scpus -rserv server-$host-$$.info $options\n";

    system "$rsh -n $host '( cd $rdir; $cmd -cpus $scpus -rserv server-$host-$$.info $options )'";
    system "$rsh -n $host '( cd $rdir; /bin/rm server-$host-$$.info )'";

    if ($subdir && defined $logfile) {
      system "cd $workdir; $rsh $host '( cd $rdir/$workdir; tar cf - */$logfile )' | tar xf -";
    }

    if ($subdir && $remsubdir) {
      system "$rsh -n $host '( cd $dir; /bin/rm -rf $sdir )'";
    } 
    exit 0;
  }

  return ($pid,$cpus-$scpus);
}

## function: readCustomFile(custom, file)
## reads custom commands from a custom file

sub readCustomFile {
  my $custom=shift;
  my $fname=shift;

  return undef if (!checkFile($fname));

  my $inp=getInputFile($fname);
  
  my $key=undef;
  while (<$inp>) {
    if (/\#CUSTOM +(.+)\n$/) {
      $key=$1;
      $custom->{$key}="";
    } elsif (defined $key) {
      $custom->{$key}.=$_;
    }
  }

  close $inp;
}

## function: writeCustomFile(custom,file)
## writes custom commands to a custom file

sub writeCustomFile {
  my $custom=shift;
  my $fname=shift;

  return if (!defined $custom);
  
  my $out=getOutputFile($fname);

  foreach my $c ( keys %{$custom} ) {
    print $out "#CUSTOM $c\n",$custom->{$c};
  }

  close $out;
}

## function: writeArchiveFile(file,data,inx)
## writes data set to archive file

sub writeArchiveFile {
  my $file=shift;
  my $data=shift;
  my $inx=shift;
  my $ftype=shift;

  my $len=length($data);

  return if ($len==0);
  
  my $offset;

  my $header;
  my $dlen;
  my $dcnt;
  my $dft;

  sysopen(DAD,$file,O_RDWR|O_CREAT);
  my $flen=sysseek DAD,0,2;
  if ($flen>0) {
    sysseek(DAD,0,0);
    
    $header="";
    $offset=0;
    do {
      my $nread=sysread(DAD,$header,40-$offset,$offset);
      die "cannot read from archive file" unless (defined $nread);
      $offset+=$nread;
      printf STDERR "read header incomplete" if ($offset<40);
    } while ($offset<40);

#    sysread(DAD,$header,40);

    $dlen=substr($header,0,10)+0;
    $dcnt=substr($header,10,10)+0;
    $dft=substr($header,20,10)+0;

    die "data length mismatch: header says >$dlen<, data length is >$len<" if ($dlen != $len);
    die "file type mismatch" if (defined $ftype && $ftype != $dft);

    $dcnt++ if ($inx>$dcnt);

    sysseek(DAD,0,0);
  } else {
    $dlen=$len;
    $dcnt=1;
    $dft=$ftype;
  }
  
  $header=sprintf("%10d%10d%10d%10d",$dlen,$dcnt,$dft,0);

  $offset=0;
  do {
    my $nwrite=syswrite(DAD,$header,40-$offset,$offset);
    die "cannot write to archive file" unless (defined $nwrite);
    $offset+=$nwrite;
      printf STDERR "write header incomplete" if ($offset<40);
  } while ($offset<40);

#  syswrite(DAD,$header,40);

  if ($inx<=$dcnt) {
    sysseek(DAD,$dlen*($inx-1),1);
  } else {
    sysseek(DAD,0,2);
  }

  $offset=0;
  do {
    my $nwrite=syswrite(DAD,$data,length($data)-$offset,$offset);
    die "cannot write to archive file" unless (defined $nwrite);
    $offset+=$nwrite;
      printf STDERR "write data incomplete" if ($offset<length($data));
  } while ($offset<length($data));

#  syswrite DAD,$data,length($data);

  close DAD;
}

## function: archiveFile(arfile,file,inx)
## copies a given file into an archive file at the given index

sub archiveFile {
  my $arfile=shift;
  my $file=shift;
  my $inx=shift;
  
  my $tfile=getInputFile($file);
  my $len=sysseek($tfile,0,2);
  sysseek($tfile,0,0);
  my $data="";
  sysread($tfile,$data,$len+0);
  writeArchiveFile($arfile,$data,$inx,0);
}

## function: $data=readArchiveFile(file,inx)
## reads data set from archive file

sub readArchiveFile {
  my $file=shift;
  my $inx=shift;

  my $header;
  my $dlen;
  my $dcnt;
  my $dft;

  sysopen(DAD,$file,O_RDONLY) || die "cannot open archive file $file";
  sysread(DAD,$header,40);
  $dlen=substr($header,0,10)+0;
  $dcnt=substr($header,10,10)+0;

  die "cannot find data record $inx in archive file $file" if ($inx>$dcnt);

  sysseek(DAD,$dlen*($inx-1),1);

  my $data;
  sysread DAD,$data,$dlen;

  close DAD;

  return $data;
}

## function: ($length,$count,$filetype,$aux)=readArchiveHeader(file)
## reads header from archive file

sub readArchiveHeader {
  my $file=shift;
  my $inx=shift;

  my $header;
  my $dlen;
  my $dcnt;
  my $dft;
  my $daux;

  sysopen(DAD,$file,O_RDONLY) || die "cannot open archive file $file";
  sysread(DAD,$header,40);
  $dlen=substr($header,0,10)+0;
  $dcnt=substr($header,10,10)+0;
  $dft=substr($header,20,10)+0;
  $daux=substr($header,30,10)+0;

  close DAD;

  return ($dlen,$dcnt,$dft,$daux);
}

## function: parsePar(hash,string) 
## parses the given parameter string and stores
## the values in the hash data structure

sub parsePar {
  my $hdat=shift;
  my $str=shift;
  
  foreach my $p ( split(/,/,$str) ) {
    my ($key,$val)=split(/=/,$p,2);
    if (defined $val) {
      $hdat->{$key}=$val;
    } else {
      if ($key!~/^noe/ && $key=~/^no(.+)$/) {
	$hdat->{$1}=0;
      } else {
	$hdat->{$key}=1;
      }
    }
  }
}


## function: zPad(num,len) 
## pads the input number with leading zeros
## to return a string of the desired length
sub zPad {
    my $num=shift;
    my $newlen=shift;

    my $origlen=length($num);
    for (my $i=0; $i<($newlen-$origlen); $i++) {
	$num="0".$num;
    }
    return $num;
}

## function: readFortran(handle)

sub readFortran {
  my $handle=shift;

  my $dat;
  my $tdat;

  my $len;
  if ($GenUtil::fortranheadersize<0) {
    read($handle,$tdat,8) || die "cannot read data";
    if (substr($tdat,4,1) eq "C") {
       $GenUtil::fortranheadersize=4;
       $len=unpack("L",substr($tdat,0,4));
       read($handle,$dat,$len-4) || die "cannot read data";
       $dat=substr($tdat,4,4).$dat;
    } else {
       $GenUtil::fortranheadersize=8;
       $len=unpack("L",$tdat);
       read($handle,$dat,$len) || die "cannot read data";
    } 
  } else {
    read($handle,$tdat,$GenUtil::fortranheadersize) || die "cannot read data";
    $len=unpack("L",$tdat);
    read($handle,$dat,$len) || die "cannot read data";
  }
  read($handle,$tdat,$GenUtil::fortranheadersize) || die "cannot read data";

  return ($dat,$len);
}

## function: acos(value) 

sub acos {
  my $val=shift;
  
  return $pi/2.0-atan2($val,sqrt(1.0-$val*$val));
}

sub asin {
  my $val=shift;
  
  return atan2($val,sqrt(1.0-$val*$val));
}

sub nint {
  my $x = $_[0];
  my $n = int($x);
  if ( $x > 0 ) {
    if ( $x-$n > 0.5) {
      return $n+1;
    }
    else {
      return $n;
    }
  }
  else {
    if ( $n-$x > 0.5) {
      return $n-1;
    }
    else {
      return $n;
    }
  }
}
 

1;
