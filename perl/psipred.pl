#!/usr/bin/env perl

# runs PSIPRED
#
# 2006, Michael Feig, Michigan State University

sub usage {
  printf STDERR "usage:   psipred.pl [options] [seqFile]\n";
  printf STDERR "options:            [-long]\n";
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

my $seqname;
my $long=0;

my $done=0;
while ($#ARGV>=0 && !$done) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-long") {
    shift @ARGV;
    $long=1;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "invalid option\n";
    &usage();
  } else {
    $seqname = shift @ARGV;
    $done=1;
  }
}

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

my $blastpgp=&GenUtil::findExecutable("blastpgp");
my $makemat=&GenUtil::findExecutable("makemat");

die "cannot find blastpgp or makemat" if (!defined $blastpgp || !defined $makemat);

system("$blastpgp -b 0 -j 3 -h 0.001 -d nr -i $$.seq -C $$.chk 2>&1 > /dev/null");
system("echo $$.chk > $$.pn");
system("echo $$.seq  > $$.sn");
system("$makemat -P $$");

my $psipreddir=$ENV{'PSIPREDDIR'};
if (!defined $psipreddir || $psipreddir eq "" || !-d $psipreddir) {
  $psipreddir="/apps/psipred";
  $psipreddir="/usr/local/psipred" if (!defined $psipreddir || !-d $psipreddir);
} 

die "cannot find PSIPRED directory" if (!defined $psipreddir || !-d $psipreddir);

system("$psipreddir/bin/psipred $$.mtx $psipreddir/data/weights.dat $psipreddir/data/weights.dat2 $psipreddir/data/weights.dat3 $psipreddir/data/weights.dat4 > $$.ss");

system("$psipreddir/bin/psipass2 $psipreddir/data/weights_p2.dat 1 1.0 1.3 $$.ss2 $$.ss > $$.horiz");

open INP,"$$.ss2";
while (<INP>) {
  if (/[0-9]+\s([A-Z])\s([A-Z])\s+[0-9\.]+\s+[0-9\.]+\s+[0-9\.]+/) {
    if ($long) {
      print;
    } else {
      printf "%s",$2;    
    }
  } 
}
close INP;

printf "\n" unless ($long);

system("rm -f $$.*");
