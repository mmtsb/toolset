#!/usr/bin/env perl

# runs BLAST, PSI-BLAST, or PDB-BLAST
#
# 2006, Michael Feig, Michigan State University

sub usage {
  printf STDERR "usage:   psiblast.pl [options] [seqFile]\n";
  printf STDERR "options:             [-nr] [-pdb] [-db name]\n";
  printf STDERR "                     [-eshow value] [-ekeep value]\n";
  printf STDERR "                     [-n alignments] [-i iterations]\n";
  printf STDERR "                     [-gap initial extension]\n";
  printf STDERR "                     [-matrix BLOSUM45|BLOSUM50|BLOSUM62|BLOSUM80|PAM250|PAM70|PAM30]\n";
  printf STDERR "                     [-verbose] [-log file]\n";
  printf STDERR "                     [-readlog file] [-no value]\n";
  printf STDERR "                     [-csblast]\n";
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

my $seqname;

my $db="nr";
my $eshow=10;
my $ekeep=0.001;
my $nalign=250;
my $niter=3;

my $gapinitial;
my $gapextend;

my $matrix="BLOSUM62";

my $logfile;

my $verbose=0;

my $readlog;
my $no;

my $csblast=0;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-nr") {
    shift @ARGV;
    $db="nr";
  } elsif ($ARGV[0] eq "-pdb") {
    shift @ARGV;
    $db="pdb";
  } elsif ($ARGV[0] eq "-db") {
    shift @ARGV;
    $db=shift @ARGV;
  } elsif ($ARGV[0] eq "-matrix") {
    shift @ARGV;
    $matrix=shift @ARGV;
  } elsif ($ARGV[0] eq "-eshow") {
    shift @ARGV;
    $eshow=shift @ARGV;
  } elsif ($ARGV[0] eq "-ekeep") {
    shift @ARGV;
    $ekeep=shift @ARGV;
  } elsif ($ARGV[0] eq "-gap") {
    shift @ARGV;
    $gapinitial=shift @ARGV;
    $gapextend=shift @ARGV;
  } elsif ($ARGV[0] eq "-no") {
    shift @ARGV;
    $no=shift @ARGV;
  } elsif ($ARGV[0] eq "-readlog") {
    shift @ARGV;
    $readlog=shift @ARGV;
  } elsif ($ARGV[0] eq "-n") {
    shift @ARGV;
    $nalign=shift @ARGV;
  } elsif ($ARGV[0] eq "-i") {
    shift @ARGV;
    $niter=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    $logfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-verbose") {
    shift @ARGV;
    $verbose=1;
  } elsif ($ARGV[0] eq "-csblast") {
    shift @ARGV;
    $csblast=1;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $seqname = shift @ARGV;
    $done=1;
  }
}

my $logout;

my $inp=&GenUtil::getInputFile($seqname);

open OUT,">$$.seq";

my $seqstr="";
my $tag;
my $printedheader;
while (<$inp>) {
  chomp;
  if (/^>(.+)$/ && !defined $tag) {
    $tag=$1;
  } elsif (/^>/ || /^\#/) {
  } else {
    if (!defined $printedheader) {
      if (defined $tag) {
	printf OUT "> %s\n",$tag;
      } else {
	printf OUT "> unknown\n";
      }
      $printedheader=1;
    }
    
    s/[^A-Z]//g;
    
    printf OUT "%s\n",$_;
    
    $seqstr.=$_;
  }
}
close OUT;
  
undef $inp;

if (!defined $readlog) {  
  my $blastpgp=&GenUtil::findExecutable("blastpgp");
  die "cannot find blastpgp" if (!defined $blastpgp); 
  my $gapstr="";
  $gapstr=sprintf("-G %s -E %s",$gapinitial,$gapextend) 
    if (defined $gapinitial && defined $gapextend);

  my $csblastexec;
  my $blastpath;
  my $cslib;

  if ($csblast) {
    $csblastexec=&GenUtil::findExecutable("csblast");
    die "cannot find csblast executable" if (!defined $csblastexec);
    $blastpath=$blastpgp;
    $blastpath=~s/\/blastpgp//;
    if (exists $ENV{CSBLASTLIB}) {
      $cslib=$ENV{CSBLASTLIB};
    } elsif (-r "/apps/csblast/K4000.lib") {
      $cslib="/apps/csblast/K4000.lib";
    } elsif (-r "K4000.lib") {
      $cslib="K4000.lib";
    }
    die "cannot find K4000 library" if (!defined $cslib || !-r $cslib);
  }


  my $cmd;
  if ($db eq "pdb") {
    if ($csblast) {
      $cmd=sprintf("%s -M %s %s -i %s -d nr -j %d -b %d -e %f -h %f -C $$.chk --blast-path %s -D %s 2>&1 > $$.blastoutnr\n",
		   $csblastexec,$matrix,$gapstr,"$$.seq",$niter,$nalign,$eshow,$ekeep,$blastpath,$cslib);
      system($cmd);
      $cmd=sprintf("%s -M %s %s -i %s -d pdbaa -j 1 -b %d -e %f -R $$.chk --blast-path %s -D %s 2>&1 > $$.blastout",
		   $csblastexec,$matrix,$gapstr,"$$.seq",$nalign,$eshow,$blastpath,$cslib);
    } else {
      $cmd=sprintf("%s -M %s %s -i %s -d nr -j %d -b %d -e %f -h %f -C $$.chk 2>&1 > $$.blastoutnr\n",$blastpgp,$matrix,$gapstr,"$$.seq",$niter,$nalign,$eshow,$ekeep);
      system($cmd);
      $cmd=sprintf("%s -M %s %s -i %s -d pdbaa -j 1 -b %d -e %f -R $$.chk 2>&1 > $$.blastout",$blastpgp,$matrix,$gapstr,"$$.seq",$nalign,$eshow,(defined $logfile)?$logfile:"$$.blastout");
    }
  } else {
    if ($csblast) {
      $cmd=sprintf("%s -M %s %s -i %s -d nr -j %d -b %d -e %f -h %f --blast-path %s -D %s 2>&1 > $$.blastout\n",
		   $csblastexec,$matrix,$gapstr,"$$.seq",$niter,$nalign,$eshow,$ekeep,$blastpath,$cslib);
    } else {
      $cmd=sprintf("%s -M %s %s -i %s -d %s -j %d -b %d -e %f 2>&1 > $$.blastout",$blastpgp,$matrix,$gapstr,"$$.seq",$db,$niter,$nalign,$eshow);
    }
  }

  printf STDERR "$cmd\n" if ($verbose);

  system($cmd);
  
  if (defined $logfile) {
    $logout=&GenUtil::getOutputFile($logfile);
    if ($db eq "pdb") {
      open INP,"$$.blastoutnr";
      while (<INP>) {
	print $logout $_;
      }
      close INP;
    }
  }

  $readlog="$$.blastout";
} 

my $data;
my $trec;
my $skip;
open INP,$readlog;
while(<INP>) {
  print $logout $_ if (defined $logout);

  chomp;
  if (/^>([^\|]+)\|(.+)/) {
    if (defined $trec && exists $trec->{id}) {
      push(@{$data},$trec);
    }
    $trec={};
    $skip=0;
    $trec->{type}=$1;
    my $rest=$2;

    if ($trec->{type} eq "pdb" && $rest=~/([^\|]+)\|([A-Z]*) (.+)/) {
      $trec->{id}=$1;
      $trec->{chain}=$2;
      $trec->{comment}=$3;
    } elsif ($rest=~/([^ ]+)\s+(.*)/) {
      $trec->{id}=$1;
      $trec->{comment}=$2;
    } else {
      printf STDERR "entry ::%s:: ignored\n",$_;
    }
  } elsif (/Expect\s+=\s+([^,]+),/) {
    if (exists $trec->{evalue}) {
      $skip=1;
    } else {
      $trec->{evalue}=$1;
    }
  } elsif (!$skip && /Identities = [0-9]+\/[0-9]+ \(([0-9]+).\)/) {
    $trec->{seqident}=$1;
    $trec->{query}=();
    $trec->{result}=();
  } elsif (!$skip && /Query: ([0-9]+) +([^ ]+)/) {
    my $srec={};
    $srec->{inx}=$1;
    $srec->{seq}=$2;
    push(@{$trec->{query}},$srec);
  } elsif (!$skip && /Sbjct: ([0-9]+) +([^ ]+)/) {
    my $srec={};
    $srec->{inx}=$1;
    $srec->{seq}=$2;
    push(@{$trec->{result}},$srec);
  } elsif (/round/ || /significant alignment/) {
    $data=();
    undef $trec;
  }
}

push(@{$data},$trec);

close INP;

undef $logout if (defined $logout);

my $seqlen=length($seqstr);

my $empty=$seqstr;
for (my $i=0; $i<$seqlen; $i++) {
  substr($empty,$i,1)=".";
}

foreach my $d (@{$data}) {
  my $query="";
  my $result="";
  
  my $qinx=$d->{query}->[0]->{inx};
  my $rinx=$d->{result}->[0]->{inx};

  for (my $na=0; $na<=$#{$d->{query}}; $na++) {
    my $aq=$d->{query}->[$na];
    my $ar=$d->{result}->[$na];

    $query.=$aq->{seq};
    $result.=$ar->{seq};
  } 

  my $alen=length($result);

  my $de=substr($empty,0,$qinx-1);
  my $deq=substr($seqstr,0,$qinx-1);
  my $nm=0;

  my $qoff=$qinx;
  for (my $ia=0; $ia<$alen; $ia++) {
	substr($de,$ia+$qoff-1,1)=substr($result,$ia,1);
        $deq.=substr($query,$ia,1);
	$nm++ if (substr($query,$ia,1) ne "-");
  }
  $deq.=substr($seqstr,$qinx+$nm-1);
  $de.=substr($empty,$qinx+$nm-1);

  $d->{alignment}=$de;
  $d->{querymatch}=$deq;
}

my $nn=1;
foreach my $d (@{$data}) {
  my $name=$d->{id};
  if (defined $d->{chain} && $d->{chain} ne "") {
    $name.="_".$d->{chain};
  }

  if (!defined $no || $no==$nn++) {
    printf ">TARGET \n";
    printf "%s\n",$d->{querymatch};
    printf ">%s ::%s::%s::%s\n",$name,$d->{evalue},$d->{seqident},$d->{comment};
    printf "%s\n",$d->{alignment};
  }
}

&GenUtil::remove("$$.chk") if (-r "$$.chk");
&GenUtil::remove("$$.seq") if (-r "$$.seq");
&GenUtil::remove("$$.blastout") if (-r "$$.blastout");
&GenUtil::remove("$$.blastoutnr") if (-r "$$.blastoutnr");

