#!/usr/bin/env perl
#
# run single hybrid lattice/MD simulations with CHARMM
# http://mmtsb.scripps.edu/doc/hlamc.pl.html
# 2002, Michael Feig, Brooks group, TSRI
#

sub usage {
  printf STDERR "usage:    hlamc.pl [options] PDBfile\n";
  printf STDERR "options:  [-n runs]\n";
  printf STDERR "          [-seq file]\n";
  printf STDERR "          [-par mctemp=value,mode=coupled|free,\n";
  printf STDERR "                [no]keepchain,[no]eval,\n";
  printf STDERR "                caforce=value,hmcmforce=value]\n";
  printf STDERR "          [-latpar ncycle=value,icycle=value,temp=value]\n";
  printf STDERR "          [-aapar CHARMMparams]\n";
  printf STDERR "          [-evalpar CHARMMparams]\n";
  printf STDERR "          [-dir name]\n";
  printf STDERR "          [-enstag name]\n";
  printf STDERR "          [-log file]\n";
  printf STDERR "          [-evallog file] [-aalog file]\n";
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
use CHARMM;
use MONSSTER;
use SICHO;
use Sequence;
use Ensemble;

my %aapar;
my %defaapar = ( 
 sdsteps    =>  20,
 minsteps   =>  100,
 dynsteps   =>  0,
 dyntemp    =>  298,
 dynoutfrq  =>  20
); 

my %evalpar;
my %defevalpar= (); 

my %latpar;
my %deflatpar = (
 icycle  =>  20,
 temp    =>  1.0,
 ncycle  =>  1
);

my %genpar;
my %defgenpar = ( 
  mctemp  => 30.0,
  caforce   => 0.5,
  hmcmforce => 0.2,
  mcruns    => 10,
  mode      => "free",
  keepchain => 1,		 
  eval      => 1
);

my $dir=".";

my $pdbfile;

my $evallog;
my $aalog;

my $enstag="hla";

my $seqfile;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-aalog") {
    shift @ARGV;
    $aalog=shift @ARGV;
  } elsif ($ARGV[0] eq "-evallog") {
    shift @ARGV;
    $evallog=shift @ARGV;
  } elsif ($ARGV[0] eq "-log") {
    shift @ARGV;
    &GenUtil::setLogFile(shift @ARGV);
  } elsif ($ARGV[0] eq "-enstag") {
    shift @ARGV;
    $enstag=shift @ARGV;
  } elsif ($ARGV[0] eq "-latpar") {
    shift @ARGV;
    &GenUtil::parsePar(\%latpar,shift @ARGV);
  } elsif ($ARGV[0] eq "-aapar") {
    shift @ARGV;
    &GenUtil::parsePar(\%aapar,shift @ARGV);
  } elsif ($ARGV[0] eq "-evalpar") {
    shift @ARGV;
    &GenUtil::parsePar(\%evalpar,shift @ARGV);
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%genpar,shift @ARGV);
  } elsif ($ARGV[0] eq "-seq") {
    shift @ARGV;
    $seqfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-dir") {
    shift @ARGV;
    $dir=shift @ARGV;
  } elsif ($ARGV[0] eq "-n") {
    shift @ARGV;
    $genpar{mcruns}=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "unknown option %s\n",shift @ARGV;
    &usage();
  } else {
    $pdbfile=shift @ARGV;
  }
}

&GenUtil::makeDir("$dir");

my $latens=Ensemble->new("lat",$dir,"lat.cfg");

my $ens=Ensemble->new($enstag,$dir);
my $at=$ens->{par}->{runs}+1;

$ens->set(seq=>$seqfile);

die "need sequence file"
  if (!defined $ens->getSeq());

$ens->setOption(%genpar);
$ens->setDefOption(%defgenpar);

$ens->setPar(%aapar);
$ens->setDefPar(%defaapar);

$latens->setPar(%latpar);
$latens->setDefPar(%deflatpar);

$ens->save();
$latens->saveOptions();

my $mcharmm=&CHARMM::new((defined $aalog)?"$dir/$aalog":undef);
$mcharmm->loadParameters(%{$ens->getPar()});
$mcharmm->setupFromPDB($pdbfile);

my $echarmm;

if ($ens->{opt}-{"eval"} != 0) {
  my $evalens=Ensemble->new("eval",$dir,"eval.cfg");
  $evalens->setPar(%evalpar);
  $evalens->setDefPar(%defevalpar);
  $evalens->saveOptions();

  $echarmm=&CHARMM::new((defined $evallog)?"$dir/$evallog":undef);
  $echarmm->loadParameters(%{$evalens->getPar()});
  $echarmm->setupFromPDB($pdbfile);
  $echarmm->setupEnergy();
}

my $first=1;
my $outmol=Molecule::new();

my $monsster=MONSSTER::new($ens->getSeq());
$monsster->setDirectory($dir);
$monsster->setParameters(%{$latens->getPar()});
$monsster->setParameters(tsteps=>1);

my $inpchain=SICHO::new();
my $initmol=Molecule::new($pdbfile);
$initmol->center();
$inpchain->genMONSSTERFromAllAtom($initmol);

my $cons=();

if ($ens->{opt}->{mode} eq "coupled" && $ens->{opt}->{hmcmforce}>0.0) {
  my $c={};
  $c->{type}="hmcm";
  $c->{list}=&GenUtil::fragListFromOption(sprintf("-9999:9999_%1.1f",$ens->{opt}->{hmcmforce}));
  $c->{reffile}="$dir/$MONSSTER::fileName{finalchain}";
  push(@{$cons},$c);
} elsif ($ens->{opt}->{mode} eq "free" && $ens->{opt}->{caforce}>0.0) {
  my $c={};
  $c->{sel}="ca";
  $c->{type}="self";
  $c->{list}=&GenUtil::fragListFromOption(sprintf("-9999:9999_%1.1f",$ens->{opt}->{caforce}));
  $c->{exclmode}=0;
  push(@{$cons},$c);
}

my $lastener=0.0;

for (my $job=1; $job<=$ens->{opt}->{mcruns}; $job++) {
  $monsster->run($inpchain);

  $mcharmm->clearEnergy() if (!$first);
  $mcharmm->clearRestraints() if (!$first);
  $mcharmm->clearShake() if (!$first);

  my $outchain;
  if ($ens->{opt}->{mode} eq "free") {
    $outchain=SICHO::new();
    $outchain->readChain("$dir/$MONSSTER::fileName{finalchain}");
  
    my $rebmol=Molecule::new();
    $rebmol->rebuildFromSICHO($ens->getSeq(),$outchain);
    $rebmol->writePDB("$dir/rebuilt.pdb",translate=>"CHARMM22",ssbond=>0);

    $mcharmm->initCoordinates();
    $mcharmm->readFromPDB("$dir/rebuilt.pdb");
  }

  $mcharmm->setupEnergy();

  if ($ens->{opt}->{mode} eq "free") {
    $mcharmm->setupRestraints(1.0,$cons);
  } else {
    $mcharmm->hmcmRestraint($cons->[0]);
  }    

  $mcharmm->minimizeSD() if ($mcharmm->{par}->{sdsteps}>0);
  $mcharmm->minimize() if ($mcharmm->{par}->{minsteps}>0);

  my $me=$mcharmm->_processEneOutput() if ($mcharmm->{par}->{minsteps}+$mcharmm->{par}->{sdsteps}>0);
    
  if ($mcharmm->{par}->{dynsteps}>0 && 
      (!defined $me ||
       ($me->[$#{$me}]->{total}<0.0 && $mcharmm->{_lastOutput}!~/far from minimum/))) {
    $mcharmm->shake();
    $mcharmm->runDynamics(undef,undef,undef,undef,dyntwin=>0.02*$mcharmm->{par}->{dyntemp});
  }
    
  my $chmoutpdb="$dir/chmout.pdb";
  $mcharmm->writePDB($chmoutpdb);
  
  my $etrial;
  if ($ens->{opt}->{"eval"}!=0) {
    $echarmm->clearEnergy();
    $echarmm->initCoordinates();
    $echarmm->readFromPDB("$dir/chmout.pdb");
    $echarmm->setupEnergy();
    $etrial=$echarmm->getEnergy()->{total}+0.0;
  } else {
    my $e=$mcharmm->getEnergy();
    $etrial=$e->{total}-$e->{constr};
  }

  if ($first || $etrial<$lastener ||
      rand()<exp(-($etrial-$lastener)/$ens->{opt}->{mctemp})) {
    $lastener=$etrial;
    $outmol->readPDB($chmoutpdb,
     translate=>&CHARMM::getConvType($mcharmm->{par}->{param}),chainfromseg=>1);
    &GenUtil::remove($chmoutpdb);

    if ($ens->{opt}->{keepchain}) {
      $inpchain=$outchain;
    } else {
      $inpchain=SICHO::new();
      $outmol->center();
      $inpchain->genMONSSTERFromAllAtom($outmol);
      if ($#{$inpchain->{sidechain}}<0) {
	my $tf=sprintf "tmp-%d.pdb",$job;
	printf STDERR "no chain! writing out $tf\n";
	$outmol->setSSBonds($mcharmm->{molecule}->getSSBonds());
	$outmol->writePDB($tf);
      }
    }

    &GenUtil::log("hlamc",sprintf("mcRun %d etrial: %f accepted",$job,$etrial));
    $ens->checkinPDB($at++,$outmol,$etrial);
  } else {
    &GenUtil::log("hlamc",sprintf("mcRun %d etrial: %f rejected",$job,$etrial));
  }
  
  $first=0;
}

$ens->save();
 
$mcharmm->finish();
$echarmm->finish()
  if ($ens->{opt}->{"eval"}!=0);
