#!/usr/bin/env perl

# generates MONSSTER sequence file
#
# http://mmtsb.scripps.edu/doc/genseq.pl.html
# 2000, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:   genseq.pl [options] [-pdb | -monsster | -one] [file]\n";
  printf STDERR "options: [-out [one][sec] | monsster]\n";
  printf STDERR "         [-sel from:to]\n";
  printf STDERR "         [-s inx:sequence[=inx:sequence]]\n";
  printf STDERR "         [-2ndpred file[:file...]\n";
  printf STDERR "         [-2ndone file]\n";
  printf STDERR "         [-dssp] [-dsspfull] [-dsspnum]\n";
  printf STDERR "         [-fill]\n";
  exit 1;
}

require 5.004;

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;

use GenUtil;
use Molecule;
use Sequence;

my $inpfile;
my $input="pdb";
my $output="monsster";
my @predfiles=();
my $secondary="";
my $slist=();
my $secfile;
my $fraglist;
my $fill;
my $header="sp";
my $fasta=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-pdb") {
    shift @ARGV;
    $input="pdb";
  } elsif ($ARGV[0] =~ /^-(monsster|one|pdb)$/) {
    ($input=shift @ARGV)=~s/^-//;
  } elsif ($ARGV[0] eq "-fasta") {
    shift @ARGV;
    $fasta=1;
    $output="one";
    $input="one";
  } elsif ($ARGV[0] eq "-out") {
    shift @ARGV;
    $output=shift @ARGV;
  } elsif ($ARGV[0] eq "-sel") {
    shift @ARGV;
    $fraglist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-dssp") {
    shift @ARGV;
    $secondary="dssp";
  } elsif ($ARGV[0] eq "-dsspfull") {
    shift @ARGV;
    $secondary="dsspfull";
  } elsif ($ARGV[0] eq "-dsspnum") {
    shift @ARGV;
    $secondary="dsspnum";
  } elsif ($ARGV[0] eq "-fill") {
    shift @ARGV;
    $fill=1;
  } elsif ($ARGV[0] eq "-2ndone") {
    shift @ARGV;
    $secfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-2ndpred") {
    shift @ARGV;
    $secondary="files";
    @predfiles=split(/:/,shift @ARGV);
  } elsif ($ARGV[0] eq "-s") {
    shift @ARGV;
    foreach my $t ( split(/=/,shift @ARGV) ) {
      my ($tinx,$tseq)=split(/:/,$t);
      my $srec={};
      $srec->{inx}=$tinx;
      $srec->{seq}=$tseq;
      push(@{$slist},$srec);    
    }
  } elsif ($ARGV[0] =~ /^-.*/) { 
    printf STDERR "invalid option %s\n",shift @ARGV;
    &usage();
  } else {
    $inpfile=shift @ARGV;
  }
}

my $seq;
if ($input eq "pdb") {
  my $mol=Molecule::new();
  $mol->readPDB($inpfile);
  $mol->selectChain("");
  $seq=Sequence::new($mol,$slist,$fill);
  $seq->secFromDSSP($mol)
    if ($secondary eq "dssp" && $#{$slist}<0);
  $seq->secFromDSSP($mol,1)
    if ($secondary eq "dsspfull" && $#{$slist}<0);
  $seq->secFromDSSP($mol,1,1)
    if ($secondary eq "dsspnum" && $#{$slist}<0);
} elsif ($input eq "monsster") {
  $seq=Sequence::new();
  $seq->readMONSSTER($inpfile);
} elsif ($input eq "one") {
  my $inp=&GenUtil::getInputFile($inpfile);
  my $seqstr="";
  while(<$inp>) { 
    chomp $_; 
    if (/^\>(.*)/) {
      $header=$1;
    } else {
      $seqstr.=$_;
    }
  }
  $seqstr=~s/ +//g;
  $seq=Sequence::new($seqstr);
} else {
  die "invalid input mode $input";
}

if (defined $secfile) {
  $seq->secFromOneFile($secfile);
} elsif ($secondary eq "files") {
  $seq->secFromPredict(@predfiles);
}

if (defined $fraglist) {
  $seq->setValidResidues($fraglist);
}

if ($output eq "monsster") {
  $seq->writeMONSSTER(\*STDOUT,1);
} elsif ($output =~ /^one/) {
  printf STDOUT ">%s\n",$header if ($fasta);
  printf STDOUT "%s\n",$seq->abbrevSeq(1);
  if ($output eq "onesec") {
    printf STDOUT "%s\n",$seq->abbrevSec(1);
  }
} elsif ($output =~ /sec/) {
  printf STDOUT "%s\n",$seq->abbrevSec(1);
} else {
  die "unknown output mode\n";
}

