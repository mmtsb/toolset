#!/usr/bin/env perl
#
# start and equilibrate explicit solvent MD simulations with CHARMM
# 2006, Michael Feig, Michigan State University
#

sub usage {
  printf STDERR "usage:    equiCHARMM.pl [options] PDBfile\n";
  printf STDERR "options:  [-prefix filename] [-logs] [-trajs] [-cmds]\n";
  printf STDERR "          [-par CHARMMparams]\n";
  printf STDERR "          [-cons [ca|cb|cab|heavy] refpdb|self min:max[_force][=...]]\n";
  printf STDERR "          [-fixsolute] [-densitysteps value]\n";
  printf STDERR "          [-equi temp:steps[=temp:steps]]\n";
  printf STDERR "          [-cubicbox] [-octabox] [-cutoff value]\n";
  printf STDERR "          [-neutralize] [-addions concentration]\n";
  printf STDERR "          [-[no]center]\n";
  printf STDERR "          [-custom file]\n";
  printf STDERR "          [-verbose]\n";
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

my %par = ( 
 shake      =>  1,
 dyntemp    =>  298,
 dynens     =>  "NVT",
 dyneqfrq   =>  100,
 dynoutfrq  =>  50,
 sdsteps    =>  50,
 minsteps   =>  500,
 param      =>  36,
 cmap       =>  1
); 

my $cons=();

my $pdbfile;
my $prefix="md";

my $logs;
my $trajs=0;
my $cmds;

my $fixsolute=0;

my $customfile;

my $shape;
my $cutoff=9;
my $center=1;

my $neutralize=0;

my $verbose=0;

my $addionsconc=0;

my @equisteps=();

my $densitysteps=5000;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-prefix") {
    shift @ARGV;
    $prefix=shift @ARGV;
  } elsif ($ARGV[0] eq "-logs") {
    shift @ARGV;
    $logs=1;
  } elsif ($ARGV[0] eq "-cmds") {
    shift @ARGV;
    $cmds=1;
  } elsif ($ARGV[0] eq "-trajs") {
    shift @ARGV;
    $trajs=1;
  } elsif ($ARGV[0] eq "-custom") {
    shift @ARGV;
    $customfile=shift @ARGV;
  } elsif ($ARGV[0] eq "-center") {
    shift @ARGV;
    $center=1;
  } elsif ($ARGV[0] eq "-nocenter") {
    shift @ARGV;
    $center=0;
  } elsif ($ARGV[0] eq "-par") {
    shift @ARGV;
    &GenUtil::parsePar(\%par,shift @ARGV);
  } elsif ($ARGV[0] eq "-cons") {
    shift @ARGV;
    my $c={};
    if ($ARGV[0] =~ /^(ca|cb|cab|cabp|cap|heavy)$/) {
      $c->{sel}=shift @ARGV;
    } else {
      $c->{sel}="heavy";
    }
    if ($ARGV[0] eq "self") {
      shift @ARGV;
      $c->{type}="self";
    } else {
      $c->{type}="ref";
      $c->{reffile}=shift @ARGV;
    }

    $c->{list}=&GenUtil::fragListFromOption(shift @ARGV);
    $c->{exclmode}=0;
    push(@{$cons},$c);
  } elsif ($ARGV[0] eq "-verbose") {
    shift @ARGV;
    $verbose=1;
  } elsif ($ARGV[0] eq "-fixsolute") {
    shift @ARGV;
    $fixsolute=1;
  } elsif ($ARGV[0] eq "-equi") {
    shift @ARGV;
    foreach my $e ( split(/=/,shift @ARGV) ) {
      my ($t,$st)=split(/:/,$e);
      push(@equisteps,{temp=>$t, steps=>$st});
    }
  } elsif ($ARGV[0] eq "-densitysteps") {
    shift @ARGV;
    $densitysteps=shift @ARGV;
  } elsif ($ARGV[0] eq "-cubicbox") {
    shift @ARGV;
    $shape="cubic";
  } elsif ($ARGV[0] eq "-octabox") {
    shift @ARGV;
    $shape="octahedron";
  } elsif ($ARGV[0] eq "-neutralize") {
    shift @ARGV;
    $neutralize=1;
  } elsif ($ARGV[0] eq "-addions") {
    shift @ARGV;
    $addionsconc=shift @ARGV;
  } elsif ($ARGV[0] eq "-cutoff") {
    shift @ARGV;
    $cutoff=shift @ARGV;
  } elsif ($ARGV[0] =~ /^-/) {
    printf STDERR "unknown option %s\n",shift @ARGV;
    &usage();
  } else {
    $pdbfile=shift @ARGV;
  }
}

if ($#equisteps<0) {
  push(@equisteps,{temp=>50, steps=>1000});
  push(@equisteps,{temp=>100, steps=>1000});
  push(@equisteps,{temp=>150, steps=>1000});
  push(@equisteps,{temp=>200, steps=>2000});
  push(@equisteps,{temp=>250, steps=>2000});
  push(@equisteps,{temp=>275, steps=>2000});
  push(@equisteps,{temp=>300, steps=>2000});
}

my $mol=Molecule::new();
$mol->readPDB($pdbfile,ignoreseg=>0);
$mol->translate("CHARMM22");

my $nmol=$mol->completeWater();
$nmol->completeResidue();
$nmol->fixCOO() if (!$par{blocked});
$nmol->translate("CHARMM22");
$mol=$nmol;

my $ss=$mol->getSSBonds();

my $charmm=&CHARMM::new();
$charmm->loadParameters(%par);
$mol->fixHistidine($charmm->{par}->{hsd},$charmm->{par}->{hse},$charmm->{par}->{hsp});
$mol->changeResName($charmm->{par}->{resmod});
$mol->generateSegNames();

$charmm->setupFromMolecule($mol);
my $chmoutpdb="t$$.out.pdb";
$charmm->writePDB($chmoutpdb);
$mol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($par{param}),chainfromseg=>1);
$mol->setSSBonds($ss);

&GenUtil::remove($chmoutpdb);

my $chg=0;
if ($neutralize) {
  $charmm->setupEnergy();
  if (defined $customfile) {
    foreach my $c ( split(/:/,$customfile))  {
      if (&GenUtil::checkFile($c)) {
	my $custom=&GenUtil::readData(&GenUtil::getInputFile($c));
	$charmm->stream($custom);
      }
    }
  }
  $chg=$charmm->getTotalCharge();
  printf STDERR "charge of solute: %s\n",$chg if ($verbose);
}
$charmm->finish();

$mol->center() if ($center);

my $err=$mol->solvate($cutoff,$shape,undef,undef,(!$center)?0:undef);

my $boxa;
my $boxb;
my $boxc;

if ($err=~/box size: ([0-9\.]+) x ([0-9\.]+) x ([0-9\.]+)/) {
  $boxa=$1;
  $boxb=$2;
  $boxc=$3;
} else {
  die "could not solvate molecule\n".$err;
}

my $naions=0;
my $clions=0;

if ($neutralize) {
  if ($chg<0) {
    $naions=int(-$chg+0.5);
  } elsif ($chg>0) {
    $clions=int($chg+0.5);
  }
}

my $nwat=0;
foreach my $c ( @{$mol->{chain}} ) {
  foreach my $r ( @{$c->{res}} ) {
    if ($r->{name} eq "TIP3" || $r->{name} eq "HOH") {
      $nwat++;
    }
  }
}

my $watdens=$nwat/($boxa*$boxb*$boxc);
my $watndens=0.0333692;
my $watfrac=$watdens/$watndens;

if ($addionsconc>0) {
# volume*Avogadro's number*1000
  my $numions=$watfrac*$addionsconc*$boxa*$boxb*$boxc*6.022141E-4;
  my $ionpairs=int($numions+0.5);
  $naions+=$ionpairs;
  $clions+=$ionpairs;
}

if ($naions>0 || $clions>0) {
  my $fac=$watfrac*$boxa*$boxb*$boxc*6.022141E-4;
  printf STDERR "adding ions: %d (%1.2fM) Na+, %d (%1.2fM) Cl-\n",
    $naions,$naions/$fac,$clions,$clions/$fac if ($verbose);

  my @ions=();
  if ($naions>0) {
    my $trec={};
    $trec->{name}="SOD";
    $trec->{num}=$naions;
    push(@ions,$trec);
  } 
  if ($clions>0) {
    my $trec={};
    $trec->{name}="CLA";
    $trec->{num}=$clions;
    push(@ions,$trec);
  } 

  $mol->replaceIons(\@ions);
  $mol=$mol->clone(1);
}

$mol->writePDB(sprintf("%s.solvated.pdb",$prefix),translate=>"CHARMM22");

my $charmm=&CHARMM::new($logs?sprintf("%s.heat.log",$prefix):undef,$cmds?sprintf("%s.heat.cmd",$prefix):undef);

if ($shape eq "cubic") {
  $par{boxshape}="cubic";
  $par{boxsize}=$boxa;
} elsif ($shape eq "octahedron") {
  $par{boxshape}="octahedron";
  $par{boxx}=$boxa;
  $par{boxy}=$boxb;
  $par{boxz}=$boxc;
} else {
  $par{boxshape}="ortho";
  $par{boxx}=$boxa;
  $par{boxy}=$boxb;
  $par{boxz}=$boxc;
}

$charmm->loadParameters(%par);

$mol->fixHistidine($charmm->{par}->{hsd},$charmm->{par}->{hse},$charmm->{par}->{hsp});
$mol->changeResName($charmm->{par}->{resmod});
$mol->generateSegNames();

$charmm->setupFromMolecule($mol);

$charmm->setupEnergy();

if ($fixsolute) {
  $charmm->_sendCommand("cons fix select ( .not. (resname TIP3 .or. resname CLA .or. resname SOD) ) end");
} else {  
  $charmm->_sendCommand("cons harm force 5.0 mass select ( .not. (resname TIP3 .or. resname CLA .or. resname SOD) ) end");
}

if (defined $customfile) {
  foreach my $c ( split(/:/,$customfile))  {
    if (&GenUtil::checkFile($c)) {
      my $custom=&GenUtil::readData(&GenUtil::getInputFile($c));
      $charmm->stream($custom);
    }
  }
}

if ($charmm->{par}->{sdsteps}>0) {
  $charmm->minimizeSD();
}

if ($charmm->{par}->{minsteps}>0) {
  $charmm->minimize();
}

my $chmoutpdb=lc "pdb$$-out";
$charmm->writePDB($chmoutpdb);

my $outmol=Molecule::new();
$outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),
		 chainfromseg=>1);

$outmol->setSSBonds($charmm->{molecule}->getSSBonds());
$outmol->writePDB(sprintf("%s.min.pdb",$prefix),translate=>"CHARMM22");
&GenUtil::remove($chmoutpdb);

printf STDERR "minimization completed\n" if ($verbose);

if (!$fixsolute) {
  $charmm->_sendCommand("cons harm clear");
  $charmm->setupRestraints(1.0,$cons)
  if ($#{$cons}>=0);
} 

$charmm->shake();

my $finaltemp=$charmm->{par}->{dyntemp};

my $nrun=0;
foreach my $e ( @equisteps ) {
  $nrun++;

  $charmm->runDynamics(undef,undef,
  $trajs?sprintf("%s.heat.%d.dcd",$prefix,$e->{temp}):undef,undef,
  dyntemp=>$e->{temp},dyntwin=>0.02*$e->{temp},dynsteps=>$e->{steps});
  
  $charmm->writePDB($chmoutpdb);
  my $outmol=Molecule::new();
  $outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),
		   chainfromseg=>1);
 
  $outmol->setSSBonds($charmm->{molecule}->getSSBonds());
  $outmol->writePDB(sprintf("%s.heat.%d.pdb",$prefix,$e->{temp}),translate=>"CHARMM22");
  $mol=$outmol;
  &GenUtil::remove($chmoutpdb);
  
  printf STDERR "warm-up run %d (%d steps @ %fK) completed\n",$nrun,$e->{steps},$e->{temp} if ($verbose);
}
$charmm->finish();

if ($densitysteps>0) {
 for ($nrun=1; $nrun<=3; $nrun++) {
  my $charmm=&CHARMM::new($logs?sprintf("%s.equi.%d.log",$prefix,$nrun):undef,
                          $cmds?sprintf("%s.equi.%d.cmd",$prefix,$nrun):undef);

  if ($shape eq "cubic") {
    $par{boxshape}="cubic";
    $par{boxsize}=$boxa;
  } elsif ($shape eq "octahedron") {
    $par{boxshape}="octahedron";
    $par{boxx}=$boxa;
    $par{boxy}=$boxb;
    $par{boxz}=$boxc;
  } else {
    $par{boxshape}="ortho";
    $par{boxx}=$boxa;
    $par{boxy}=$boxb;
    $par{boxz}=$boxc;
  }
    
  $charmm->loadParameters(%par);

  $mol->fixHistidine($charmm->{par}->{hsd},$charmm->{par}->{hse},$charmm->{par}->{hsp});
  $mol->changeResName($charmm->{par}->{resmod});
  $mol->generateSegNames();

  $charmm->setupFromMolecule($mol);
  $charmm->setupEnergy();

  if ($fixsolute) {
    $charmm->_sendCommand("cons fix select ( .not. (resname TIP3 .or. resname CLA .or. resname SOD) ) end");
  } else {
    $charmm->setupRestraints(1.0,$cons)
     if ($#{$cons}>=0);
  }

  if (defined $customfile) {
    foreach my $c ( split(/:/,$customfile))  {
      if (&GenUtil::checkFile($c)) {
	my $custom=&GenUtil::readData(&GenUtil::getInputFile($c));
	$charmm->stream($custom);
      }
    }
  }

  $charmm->shake();

  $charmm->minimizeSD(sdsteps=>20);

  my $restart=sprintf("%s.equi.%d.restart",$prefix,$nrun);
  my $traj=sprintf("%s.equi.%d.dcd",$prefix,$nrun);
  $charmm->runDynamics(undef,$restart,$traj,undef,
		       dyntwin=>0.02*$charmm->{par}->{dyntemp},dynsteps=>$densitysteps);
  
  $charmm->writePDB($chmoutpdb);
  my $outmol=Molecule::new();
  $outmol->readPDB($chmoutpdb,translate=>&CHARMM::getConvType($charmm->{par}->{param}),
		   chainfromseg=>1);
  
  $outmol->setSSBonds($charmm->{molecule}->getSSBonds());
  $outmol->writePDB(sprintf("%s.equi.%d.pdb",$prefix,$nrun),translate=>"CHARMM22");
  &GenUtil::remove($chmoutpdb);

  $mol=$outmol;

  printf STDERR "equilibration run %d (%d steps @ %fK) completed\n",$nrun,$densitysteps,$finaltemp if ($verbose);

  my $maxboxsize=-999;
  $maxboxsize=$boxa if ($boxa>$maxboxsize);
  $maxboxsize=$boxb if ($boxb>$maxboxsize);
  $maxboxsize=$boxc if ($boxc>$maxboxsize);

  $charmm->{par}->{rdfrmax}=$maxboxsize/2.0;
  my ($nframes,$profile)=$charmm->analyzeRadialDistribution($traj,cofm=>1,ndens=>1);

  my $avgw=0.0;
  my $navgw=0;
  foreach my $p ( @{$profile} ) {
    if ($p->{x}>($maxboxsize/2.0-5)) {
       $avgw+=$p->{val};
       $navgw++;
    }
  }
  die "unexpected radial distribution result" if ($navgw<=0);

  $charmm->finish();
 
  $avgw/=$navgw;

  my $oa=$boxa;
  my $ob=$boxb; 
  my $oc=$boxc; 

  my $scale=($avgw/$watndens)**(1/3);

  $boxa*=$scale;
  $boxb*=$scale;
  $boxc*=$scale;

  if ($boxa>$oa+0.2) {
    $boxa=$oa+0.2;
  } elsif ($boxa<$oa-0.2) {
    $boxa=$oa-0.2;
  }
  if ($boxb>$ob+0.2) {
    $boxb=$ob+0.2;
  } elsif ($boxb<$ob-0.2) {
    $boxb=$ob-0.2;
  }
  if ($boxc>$oc+0.2) {
    $boxc=$oc+0.2;
  } elsif ($boxc<$oc-0.2) {
    $boxc=$oc-0.2;
  }

  printf STDERR "avg. density: %f; changed box size %f x %f x %f ---> %f x %f x %f\n",
    $avgw,$oa,$ob,$oc,$boxa,$boxb,$boxc if ($verbose);
 }
}

