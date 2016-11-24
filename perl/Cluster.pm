# Cluster package
# generate/manage clusters
#
# http://mmtsb.scripps.edu/doc/Cluster.pm.html
# 2000, Michael Feig, Brooks group, TSRI

package Cluster;

require 5.004;

use strict;

use FileHandle;
use IPC::Open2;

use GenUtil;
use Molecule;
use Analyze;

## data: par -> { selmode, lsqfit, fraglist, fitfraglist,
## data:          filetype, clustermode, contmaxdist,
## data:          maxnum, minsize, method,
## data:          radius, maxerr, iterate, mixfactor }
## clustering paramters

## data: topcl   -> { tag, level, 
## data:              element[] -> { name, inx, dist },
## data:              subcl[] -> { tag, level, element, subcl } }
## cluster hierarchy data structure

## data: allcl[] -> { tag, level, element, subcl }  
## list of all clusters 

## constructor: new([parameters])
## creates a new Cluster object. Clustering parameters
## may given as hash-style key=>value pair arguments

sub new {
  my $self={};

  bless $self;

  $self->{par}->{selmode}="heavy";
  $self->{par}->{lsqfit}=1;
  $self->{par}->{fraglist}=undef;
  $self->{par}->{fitfraglist}=undef;
  $self->{par}->{filetype}="pdb";
  $self->{par}->{clustermode}="rmsd";
  $self->{par}->{contmaxdist}=12.0;
  $self->{par}->{maxnum}=4;
  $self->{par}->{minsize}=100;
  $self->{par}->{radius}=2.5;
  $self->{par}->{maxerr}=0.01;
  $self->{par}->{iterate}=1;
  $self->{par}->{mixfactor}=0.3;
  $self->{par}->{method}="jclust";

  $self->setPar(@_);
  $self->{topcl}->{tag}="t";
  $self->{topcl}->{level}=0;

  $self->{allcl}=();
  $self->{_getall}=0;
 
  return $self;
}

## method: setPar(parameters)
## sets clustering parameters in hash-style key=>value
## format

sub setPar {
  my $self=shift;
  my %par=@_;

  foreach my $p ( keys %par ) {
    $self->{par}->{$p}=$par{$p}
       if (defined $par{$p} && $par{$p} ne "");
  }
}

## method: setFileList(filelist)
## resets the top level cluster from a list of files.
## This method needs to be called before clustering
## can be done.

sub setFileList {
  my $self=shift;
  my $filelist=shift;
  
  undef $self->{topcl};
  
  $self->{topcl}={};
  $self->{topcl}->{element}=();
  $self->{topcl}->{tag}="t";
  $self->{topcl}->{level}=0;

  my $inx=1;
  foreach my $f ( @{$filelist} ) {
    my $elrec={};
    $elrec->{name}=(ref $f && defined $f->{name})?$f->{name}:$f;
    $elrec->{inx}=(ref $f && defined $f->{inx})?$f->{inx}:$inx++;
    push(@{$self->{topcl}->{element}},$elrec);
  }

  $self->{_getall}=1;
}  

## method: runCluster([cluster][,method])
## runs the jclust clustering program for the set of files in
## the top level cluster or in <mark>cluster</mark> if
## given as argument. 

sub runCluster {
  my $self=shift;
  my $cluster=shift;

  $cluster=$self->{topcl} if (!defined $cluster);

  &GenUtil::log("Cluster::runCluster",$self->{par}->{method});

  my $clustexe;
  my $ucmethod=uc $self->{par}->{method};
  if (defined $ENV{$ucmethod} && -x $ENV{$ucmethod}) {
    $clustexe=$ENV{$ucmethod};
  } else {
    $clustexe=&GenUtil::findExecutable($self->{par}->{method});
  }

  die "cannot find clustering executable"
    if (!defined $clustexe);

  my $options="-$self->{par}->{filetype}";

  $options.=" -centroid";

 
  if ($self->{par}->{method} eq "jclust") {
    $options.=" -max $self->{par}->{maxnum}";
  } elsif ($self->{par}->{method} eq "kclust") {
    $options.=" -cdist";
    $options.=" -radius $self->{par}->{radius}";
    $options.=" -iterate -maxerr $self->{par}->{maxerr}"
      if ($self->{par}->{iterate});
  }

  if ($self->{par}->{clustermode} eq "rmsd") {
    $options.=" -mode rmsd";
    $options.=" -$self->{par}->{selmode}";
    $options.=" -lsqfit" if ($self->{par}->{lsqfit});
  } elsif ($self->{par}->{clustermode} eq "contact") {
    die "cannot use kclust for contact map clustering"
      if ($self->{par}->{method} eq "kclust");
    $options.=" -mode contact";
    $options.=" -maxdist $self->{par}->{contmaxdist}" 
      if (defined $self->{par}->{contmaxdist});
  } else {
    if ($self->{par}->{clustermode} eq "mix") {
      die "cannot use jclust for cluster mode $self->{par}->{clustermode}"
        if ($self->{par}->{method} eq "jclust");
      $options.=" -mode mix $self->{par}->{mixfactor}";
      $options.=" -$self->{par}->{selmode}";
      $options.=" -lsqfit" if ($self->{par}->{lsqfit});
    } else {
      $options.=" -mode $self->{par}->{clustermode}";
    }
  }
    
  if (defined $self->{par}->{fraglist}) {
    $options.=" -l $self->{par}->{fraglist}";
    if ($self->{par}->{lsqfit} && 
	($self->{par}->{clustermode} eq "rmsd" || 
	 $self->{par}->{clustermode} eq "mix")) {
      if (defined $self->{par}->{fitfraglist}) {
	$options.=" -fit $self->{par}->{fitfraglist}";
      } else {
	$options.=" -fitxl";
      }
    }
  } 

  &GenUtil::log("Cluster::runCluster","running $clustexe $options");

  printf STDERR "$clustexe $options\n";
  

  local (*READ,*WRITE);
  my $pid=open2(*READ,*WRITE,"$clustexe $options");
  
  foreach my $n ( @{$cluster->{element}} ) {
    if (!&GenUtil::checkFile($n->{name},1)) {
      printf STDERR "cannot read file $n->{name} in Cluster::runCluster\n";
    } else {
      print WRITE $n->{name},"\n";
#      print STDERR $n->{name},"\n";
    }
  }
  close WRITE;

  $cluster->{subcl}=();

  my $clrec;
  while (<READ>) {
    chomp;
    if (/^\#Cluster +([0-9]+)/) {
      $clrec={};
      $clrec->{tag}="$cluster->{tag}.$1";
      $clrec->{element}=();
      $clrec->{level}=$cluster->{level}+1;
    } elsif (/^([0-9]+) +([\S]+)\s+([\S]+)/) {
      my $el={};
      $el->{inx}=$cluster->{element}->[($1)-1]->{inx};
      $el->{name}=$2;
      $el->{dist}=$3;
      push(@{$clrec->{element}},$el);
    } elsif (/^([0-9]+) +([\S]+)/) {
      my $el={};
      $el->{inx}=$cluster->{element}->[($1)-1]->{inx};
      $el->{name}=$2;
      push(@{$clrec->{element}},$el);
    } elsif (/^\#Centroid +([0-9]+)/) {
      if ($self->{par}->{clustermode} eq "rmsd" ||
	  $self->{par}->{clustermode} eq "mix") {
	$clrec->{centroid}=Molecule::new(\*READ);
      } elsif ($self->{par}->{clustermode} eq "contact") {
	$clrec->{centroid}=&GenUtil::readData(\*READ,"^\#End");
      } else {
	$clrec->{centroid}=&GenUtil::readData(\*READ,"^\#End");
      }
      push(@{$cluster->{subcl}},$clrec);
    } elsif (/^\#End/) {
      push(@{$cluster->{subcl}},$clrec);
    }
  }
  close READ;
  waitpid($pid,0);

  $self->{_getall}=1;
}

## method: clusterHierarchy([level[,cluster]])
## generates a cluster hierarchy beginning at the
## top level cluster or at a cluster given as argument
## by calling <mark>runCluster</mark> recursively.
## The level argument sets the maximum number of cluster
## levels that may be performed recursively

sub clusterHierarchy {
  my $self=shift;
  my $level=shift;
  my $cluster=shift;
  my $subcluster=shift;

  return if (defined $level && $level<=0);

  $level=1 if (!defined $level || $self->{par}->{method} eq "kclust");

  $cluster=$self->{topcl} if (!defined $cluster);

  my $elements=$#{$cluster->{element}}+1;
  &GenUtil::log("Cluster::clusterHierarchy","level $level, elements: $elements");
  if ($elements>=$self->{par}->{minsize} || !defined $subcluster) {
    $self->runCluster($cluster);

    my $nsub=$#{$cluster->{subcl}}+1;
    &GenUtil::log("Cluster::cluster","$nsub sub clusters");
    for (my $i=0; $i<=$#{$cluster->{subcl}}; $i++) {
      $self->clusterHierarchy($level-1,$cluster->{subcl}->[$i],1);
    }
  }
}

## method: top()
## returns the top level cluster

sub top {
  my $self=shift;
  return $self->{topcl};
}

## method: readFile(file)
## reads cluster information from a file

sub readFile {
  my $self=shift;
  my $fin=&GenUtil::getInputFile(shift);

  undef $self->{topcl};
  $self->{topcl}=&_readCluster($fin,0)
    if (defined $fin);

  undef $fin;

  $self->{_getall}=1;
}

sub _readCluster {
  my $fin=shift;
  my $level=shift;

  my $c={};
  while (<$fin>) {
    if (/^\@cluster/) {
      my @f=split(/ +/);

      $c->{tag}=$f[1];
      $c->{level}=$level;
      my ($nelem,$nsub)=($f[3],$f[5]);

      $c->{element}=();
      for (my $i=0; $i<$nelem; $i++) {
	my $line=<$fin>;
	if (defined $line) {
	  my $trec={};
	  my @f=split(/\s+/,$line);
	  $trec->{inx}=$f[0];
	  $trec->{name}=$f[1];
	  $trec->{dist}=$f[2] if (defined $f[2] && $f[2] ne "");
	  push (@{$c->{element}},$trec);
	}
      }
      $c->{subcl}=();
      for (my $i=0; $i<$nsub; $i++) {
	my $sc=&_readCluster($fin,$c->{level}+1);
	push (@{$c->{subcl}},$sc);
      }
      return $c;
    }
  }
}

## method: writeFile(file[,centflag,centoutname])
## writes the current cluster hierarchy to a file.
## Centroids are generated if <mark>centflag</mark>
## is set. The name used in the file template 
## may be given through <mark>centoutname</mark>.

sub writeFile {
  my $self=shift;
  my $fout=&GenUtil::getOutputFile(shift);
  my $centflag=shift;
  my $centoutname=shift;

  printf $fout "# cluster file\n";
  my $now=localtime();
  printf $fout "# automatically generated on: %s\n",$now;
  printf $fout "# mode: %s, filetype: %s, lsqfit: %d, selmode: %s\n",
  $self->{par}->{clustermode},$self->{par}->{filetype},$self->{par}->{lsqfit},
  $self->{par}->{selmode};
  
  &_writeCluster($fout,$self->{topcl},$self->{par}->{clustermode},
		 $centflag,$centoutname);

  undef $fout;
}  

sub _writeCluster {
  my $fout=shift;
  my $c=shift;
  my $clmode=shift;
  my $centflag=shift;
  my $centoutname=shift;

  printf $fout "\@cluster %s has %d elements, %d subclusters\n",
  $c->{tag},$#{$c->{element}}+1,$#{$c->{subcl}}+1;

  return if ($#{$c->{element}}<0);

  for (my $i=0; $i<=$#{$c->{element}}; $i++) {
    my $ce=$c->{element}->[$i];
    if (defined $ce->{dist}) {
      printf $fout "%d %s %f\n",$ce->{inx},$ce->{name},$ce->{dist};
    } else {
      printf $fout "%d %s\n",$ce->{inx},$ce->{name};
    }
  }

  my $havecentflag=(defined $centflag && $centflag);
  my $havecentroid=(defined $c->{centroid});
  &GenUtil::log("_writeCluster","centflag: $havecentflag, have centroid: $havecentroid"); 

  if (defined $centflag && $centflag && defined $c->{centroid}) {
    if ($clmode eq "rmsd") {
      my $filename="$centoutname-$c->{tag}.pdb";
      $c->{centroid}->writePDB($filename,ssbond=>0);
    } elsif ($clmode eq "contact") {
      my $filename="$centoutname-$c->{tag}.rgb";
      open OUTCENT,">$filename";
      syswrite OUTCENT, $c->{centroid}, length($c->{centroid});
      close OUTCENT;
    }
  }
  
  for (my $i=0; $i<=$#{$c->{subcl}}; $i++) {
    &_writeCluster($fout,$c->{subcl}->[$i],$clmode,$centflag,$centoutname);
  }
}


## method: writeCentroids(centoutname)
## writes the current cluster hierarchy to a file.
## Centroids are generated if <mark>centflag</mark>
## is set. The name used in the file template 
## may be given through <mark>centoutname</mark>.

sub writeCentroids {
  my $self=shift;
  my $centoutname=shift;

  &_writeCentroids($self->{topcl},$self->{par}->{clustermode},$centoutname);
}  

sub _writeCentroids {
  my $c=shift;
  my $clmode=shift;
  my $centoutname=shift;

  return if ($#{$c->{element}}<0);

  if (defined $c->{centroid}) {
    if ($clmode eq "rmsd") {
      my $filename="$centoutname-$c->{tag}.pdb";
      $c->{centroid}->writePDB($filename,ssbond=>0);
    } elsif ($clmode eq "contact") {
      my $filename="$centoutname-$c->{tag}.rgb";
      open OUTCENT,">$filename";
      syswrite OUTCENT, $c->{centroid}, length($c->{centroid});
      close OUTCENT;
    }
  }
  
  for (my $i=0; $i<=$#{$c->{subcl}}; $i++) {
    &_writeCentroids($c->{subcl}->[$i],$clmode,$centoutname);
  }
}


sub _addClusterRef {
  my $self=shift;
  my $cl=shift;

  push(@{$self->{allcl}},$cl);
  foreach my $c ( @{$cl->{subcl}} ) {
    $self->_addClusterRef($c);
  }
}

## method: $list = allClusters()
## return a list of all clusters

sub allClusters {
  my $self=shift;
  
  return $self->{allcl} 
    if (!$self->{_getall});

  $self->{allcl}=();
  $self->_addClusterRef($self->{topcl});

  $self->{_getall}=0;
  return $self->{allcl};
}

## method: $list = clusterList(level)
## return all clusters found at the given level

sub clusterList {
  my $self=shift;
  my $level=shift;

  my $list=();

  foreach my $c ( @{$self->allClusters()} ) { 
    push (@{$list},$c)
      if ( ($#{$c->{subcl}}<0 && $level<0) ||
	   ($c->{level}==$level) );
  }

  return $list;
}

## method: $list = fileList(level)
## returns the list of files at a given level.

sub fileList {
  my $self=shift;
  my $level=shift;

  foreach my $c ( @{$self->allClusters()} ) { 
    if ($level eq $c->{tag}) {
      return $c->{element};
    }
  }
  return undef;
}

1;
