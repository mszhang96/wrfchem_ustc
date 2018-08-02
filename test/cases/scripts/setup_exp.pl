#!/usr/bin/perl -w
#########################################################################
# This script setup a series of WRF runs, doing cold starts on the
# meteorology and hot starts on the chemistry 
#
# Chun Zhao  ; 5-Aug-2016 
#
# History:
#   29-Sep-2005  Working copy for WRF v2.1.
#   10-Nov-2006  Added option to cycle clouds.
#   22-Apr-2008  Added option to choose either PBS or SLURM queueing
#   25-Aug-2009  Modified for WRF-CHEM v3.1.1
#   5-Aug-2016   Modified for setting the WRF-CHEM run cases automatically 
#########################################################################

#######################
#  USER DEFINED INFO  #
#######################
# Queing system to use...
#     1 = SLURM (currently used only for wrf.exe and not real.exe so 
#                set $queuereal=0)
#$machine       = "yellowstone"; #yellowstone,eastwind,tianhe
$machine       = "cheyenne";
$type          = "wrfchem" ; 
$boundarytype  = "ms37_to_ms34" ;
#$boundarytype  = "moz_to_ms34" ;
$version       = "v3.5.1_20161010" ; 
$create_case   = 0 ; 

$mode          = "forecast" ; #forecast or continue 

# Case directory path...
$region        = "hemisphere" ;
$resolution    = "1deg" ;
$grids         = "360x145" ;
$casedir       = "$region/$resolution/$grids/$version/2012/ms37_std";

# wrf.exe queueing script name... 
$queuewrf      = "wrf_batch.$machine";
$queuereal     = "real_batch.$machine";
$queueupdatebc   = "updatebc_batch.$machine";

#output stream list
$iofields_filename  = "streams_io_minfields.txt";

# Script directory path...
$scriptdir = "/glade/p/work/chunzhao/WRFCHEM/scripts" ; 

# Model directory path...
$modeldir       = "/glade/p/work/chunzhao/WRFCHEM/model/$version";

# Experiment directory head...
$expdir        = "/glade/scratch/chunzhao/WRFCHEM/run_cases/exps/$casedir";

# Output log directory...
$outputdir     = "$expdir/output";

# MET directory 
$metdir       = "/glade/scratch/chunzhao/WRFCHEM/input/EChina/12km/120x100/ECMWF/2014" ;

# Emission directory 
$emisdir       = "/glade/scratch/chunzhao/WRFCHEM/run_cases/emissions/EChina/12km/120x100" ;

# Archive directory (for data merge to HPSS) 
$transdir     = "/glade/scratch/chunzhao/WRFCHEM/output/$casedir";
system("mkdir -p '$transdir/wrfout'");
system("mkdir -p '$transdir/wrfrst'");
system("mkdir -p '$transdir/wrfdetail'");

# Boundary data directory for regional simulation 
$boundarycode  = "/glade/p/work/chunzhao/WRFCHEM/boundary/wrfchem/codes/code4" ;
$boundarydir   = "/glade/p/work/chunzhao/WRFCHEM/boundary/wrfchem/cases/$casedir" ;
$parentdir   = "/glade/scratch/chunzhao/WRFCHEM/output/hemisphere/1deg/360x145/ms37_old/2014/wrfout";
$boundarycolon = 0 ; #parent files with colon or not, default is 1

# Run start date (use a 4 digit year)...
$runsyr        = 2014;
$runsmn        = 1;
$runsdy        = 5;
$runshr        = 0;
# Run end date (use a 4 digit year)...
$runeyr        = 2014;
$runemn        = 2;
$runedy        = 5;
$runehr        = 0;
# Blocking options for each cycle of the run...
#     block unit: 0 = Hours
#     block length is the number of units per block
$blkunit       = 0; #do not change this for now, it will require lots of other changes
$blkday        = 5; # in days 
$blklen        = $blkday*24; # in hours 
# Number of hours into current cycle to start next cycle...
# This allows for overlap at the end of each cycle while spinup occurs
# for the next cycle. Make sure there is history output at this time
# and the path to it is set in $cyclesoilpath if cyclesoil (and also
# for cyclecloud) is turned on.
if ($mode eq "forecast") { 
$nextcyclehr   = $blklen-1*24 ; # days to restart 
$restart   = 2 ; # 0:restart since second cycle, 1:restart since first cycle, 2:never restart 
$firstforecast = 1 ; #0: not the first forecast cycle, 1: this is the first forecast cycle
} elsif ($mode eq "continue") {
$nextcyclehr   = $blklen ; # days to restart 
$restart   = 0 ; # 0:restart since second cycle, 1:restart since first cycle, 2:never restart 
}
# Number of domains to run...
$max_dom       = 1;
# Parts of the model to run (0=No, 1=Yes, 2=Yes, but not 1st time through loop)
$runReal       = 1;
$updatePFT     = 0; #(0=No PFT input file, 1=update PFT from the input file)
$pftfile       = "surfdata_0100x0120_chunchina12km.NCARpft_c03172017.nc" ;
$updatebdy     = 1; #(0=No chemical boundary, 1=wrfchem boundary, 2=MOZART boundary
$updateic      = 1; #(0=No chemical initial condition, 1=Global WRFChem chemical initial condition, 2=MOZART chemical initial condition 
$runEmis       = 1;
$hourlyAnthro  = 1;
$monthlyFire   = 0;
$weekFire      = 0;
$hourlyFire    = 1;
$runWRF        = 1;
$runMoveOutput = 1;

###################
#  PROGRAM START  #
###################
use File::Copy;
use FileHandle;
use POSIX;

###  Unbuffer the output so that the log is always up to date
STDOUT->autoflush(1);
STDERR->autoflush(1);

####Create case
if (-d $expdir) {
# exp directory exists
$create_case = 0 ;
}
if ($create_case eq 1) {

#mkdir -p "$expdir" or die "STOP: Could not create $expdir\n";
system("mkdir -p '$expdir'");
system("mkdir -p '$outputdir'");

copy "namelist/namelist.input.$type", "$expdir/namelist.input"
        or die "ERROR: Failed copying namelist.input: $?";
copy "namelist/$iofields_filename", "$expdir/streams_io.txt"
        or die "ERROR: Failed copying $iofields_filename: $?";
copy "machine/$queuereal", "$expdir/$queuereal"
        or die "ERROR: Failed copying $queuereal: $?";
copy "machine/$queuewrf", "$expdir/$queuewrf"
        or die "ERROR: Failed copying $queuewrf: $?";

### link input files 
&linkfiles or die "linkfiles crashed\n";

system("ln -sf $metdir $expdir/metdata");
exit 0;
} 

###  Figure out the number of time increments in the given run period
my $rstime = POSIX::mktime(0,0,$runshr,$runsdy,$runsmn-1,$runsyr-1900,0,0,-1);
my $retime = POSIX::mktime(0,0,$runehr,$runedy,$runemn-1,$runeyr-1900,0,0,-1);

if ( $retime <= $rstime ) {die "Ill configured dates were selected.\n"}
my $runlen = POSIX::ceil( ($retime-$rstime)/3600 ); #number of hr increments
print "number of hour increments in the run = $runlen\n";

###  Start looping through each block until the run is over
my ($blksyr,$blksmn,$blksdy,$blkshr,
    ,$blkeyr,$blkemn,$blkedy,$blkehr,$betime);
  $blksyr  = $runsyr;
  $blksmn  = $runsmn;
  $blksdy  = $runsdy;
  $blkshr  = $runshr;

$betime  = &IncBlk($blkshr,$blksdy,$blksmn,$blksyr);  # set the end time of this block  -czhao

my $theend  = 0;
while ( not $theend ) {

  printf "\nStarting block %2.2d%2.2d%2.2d%2.2d-%2.2d%2.2d%2.2d%2.2d\n"
    ,&yr2($blksyr),$blksmn,$blksdy,$blkshr
    ,&yr2($blkeyr),$blkemn,$blkedy,$blkehr;

 ### chdir "$expdir" or die "Could not cd to $expdir: $!\n";

###  link the initial and restart files 
  &linkinit or die "linkinit crashed\n";
  if ($restart eq 1 ) {
  $restartvalue = ".true."
  } else {
  $restartvalue = ".false."
  }

###  Do emission update 
  if ( $runEmis==1 ) {
    &updateEmis or die "runEmis crashed\n";
  } 

###  Do real.exe
  if ( $runReal==1 ) {
    &Realexe or die "Real crashed\n";
  }

###  Do chemical boundary 
  if ($updatebdy ge 1 ) {
  $bcvalue = ".true."
  } else {
  $bcvalue = ".false."
  }
  if ( $updatebdy ge 1 or $updateic ge 1) {
    &updatebc or die "updatebc crashed\n";
  } 

###  Update PFT
  $inputpftvalue = ".false."  ;
  if ($updatePFT==1) {
    copy "namelist/add_pft.ncl", "$expdir/add_pft.ncl"
        or die "ERROR: Failed copying add_pft.ncl: $?";
    chdir "$expdir" or die "Could not cd to $expdir: $!\n";
    system("ln -sf $emisdir/../$pftfile ./surfdata.nc");
    system("ncl add_pft.ncl >& add_pft.log");
    chdir "$scriptdir";
    $inputpftvalue = ".true."  ;
  }

###  Do WRF
  if ( $runWRF==1 ) {
    &WRF or die "WRF crashed\n";
  }

###  Do MoveOutput
  if ( $runMoveOutput==1 ) {
    &MoveOutput or die "MoveOutput crashed\n";
  }

  printf "\nCompleted block %2.2d%2.2d%2.2d%2.2d-%2.2d%2.2d%2.2d%2.2d\n"
    ,&yr2($blksyr),$blksmn,$blksdy,$blkshr
	,&yr2($blkeyr),$blkemn,$blkedy,$blkehr;

###  Increment the times
  if ( $betime == $retime ) { $theend = 1; }
  ($blkshr,$blksdy,$blksmn,$blksyr) = &IncDate($nextcyclehr,$blkshr,$blksdy,$blksmn,$blksyr); # set the start date of the following block
  $betime = &IncBlk($blkshr,$blksdy,$blksmn,$blksyr); #the end of the following run

  if ($restart ne 2 ) {
  $restart = 1 ;
  } else { 
  $restart = 2 ;
  }

  if ($mode eq "forecast") {
  $firstforecast = 0 ;
  }

}  #$theend
print "run completed\n";

#########################################################################
#  Sets the next block calendar end time based on
#  the time passed plus the block increment.
#########################################################################
sub IncBlk {
  my ($hr,$dy,$mn,$yr) = @_;
  my $newtime;

 SWITCH: {
# Hour blocks
    if ( $blkunit == 0 ) { 
      $newtime =  POSIX::mktime(0,0,$hr+$blklen,$dy,$mn-1,$yr-1900,0,0,-1);
      last SWITCH;
    }
  }

  if ( $newtime > $retime ) { $newtime = $retime; }
  $blkeyr = (localtime($newtime))[5] + 1900;
  $blkemn = (localtime($newtime))[4] + 1;
  $blkedy = (localtime($newtime))[3];
  $blkehr = (localtime($newtime))[2];
  return $newtime;
}


#########################################################################
#  Returns the date incremented by the number of hours requested
#########################################################################
sub IncDate {
  my ($hinc,$ihr,$idy,$imn,$iyr) = @_;
  my ($newtime,$ohr,$ody,$omn,$oyr);

  $newtime =  POSIX::mktime(0,0,$ihr+$hinc,$idy,$imn-1,$iyr-1900,0,0,-1);

  $oyr = (localtime($newtime))[5] + 1900;
  $omn = (localtime($newtime))[4] + 1;
  $ody = (localtime($newtime))[3];
  $ohr = (localtime($newtime))[2];

  return ($ohr,$ody,$omn,$oyr);
}

#########################################################################
# Moves WRF output, restart, and rsl files and then copies the
# namelist.input file to a subdirectory of $outputdir that is
# created and named based on the block start time.
#########################################################################
sub MoveOutput {
  my ($currentdir, $dirname);

  printf "\nMoveOutput for block %2.2d%2.2d%2.2d%2.2d-%2.2d%2.2d%2.2d%2.2d\n"
    ,&yr2($blksyr),$blksmn,$blksdy,$blkshr
    ,&yr2($blkeyr),$blkemn,$blkedy,$blkehr;

# $currentdir = `pwd`;
  chdir "$expdir" or die "Could not cd to $expdir: $!\n";

  $dirname = sprintf("$outputdir/wrf%4.4d%2.2d%2.2d%2.2d",
					 $blksyr,$blksmn,$blksdy,$blkshr);
  mkdir "$dirname" or die "STOP: Could not create $dirname\n";

  system("mv rsl.* $dirname") == 0
	or die "ERROR: Failed moving rsl files: $?";

  if ($mode eq "forecast") {
  $firstfile = sprintf("wrfout_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00",
                       $blksyr,$blksmn,$blksdy,$blkshr);
  system("rm -f $firstfile");
  }
  
  $lastfile = sprintf("wrfout_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00", 
                       $blkeyr,$blkemn,$blkedy,$blkehr);
  system("mv $lastfile $dirname"); # for next simulation
  system("mv wrfout_d* $transdir/wrfout") == 0
	or die "ERROR: Failed moving wrfout files: $?";

  if($mode eq "forecast") {
  $endtime =  POSIX::mktime(0,0,$blkehr,$blkedy,$blkemn-1,$blkeyr-1900,0,0,-1);
  $endtime = $endtime - (1*24*3600) ;
  $endtime = $endtime + 3600;  # one hour uncertainty
  $endyr = (localtime($endtime))[5] + 1900;
  $endmn = (localtime($endtime))[4] + 1;
  $enddy = (localtime($endtime))[3];

  $lastfile1 = sprintf("wrfout_d0?_%4.4d-%2.2d-%2.2d_00:00:00",
                       $endyr,$endmn,$enddy);
  system("ln -s $transdir/wrfout/$lastfile1 $dirname/") ;
  }

  $firstfile = sprintf("wrfrst_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00", 
                       $blksyr,$blksmn,$blksdy,$blkshr);
  $lastfile = sprintf("wrfrst_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00", 
                       $blkeyr,$blkemn,$blkedy,$blkehr);
  system("rm -f $firstfile");
  system("mv wrfrst_d* $transdir/wrfrst"); #== 0;
  system("ln -s $transdir/wrfrst/$lastfile $dirname/");

  $lastfile1 = sprintf("wrf1hr_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00", 
                       $blkeyr,$blkemn,$blkedy,$blkehr);
  $lastfile2 = sprintf("wrf3hr_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00", 
                       $blkeyr,$blkemn,$blkedy,$blkehr);
  $lastfile3 = sprintf("wrfchemistry_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00", 
                       $blkeyr,$blkemn,$blkedy,$blkehr);
  system("rm $lastfile1 $lastfile2 $lastfile3");

  if ($mode eq "forecast") {
  $firstfile1 = sprintf("wrf1hr_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00",
                       $blksyr,$blksmn,$blksdy,$blkshr);
  $firstfile2 = sprintf("wrf3hr_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00",
                       $blksyr,$blksmn,$blksdy,$blkshr);
  $firstfile3 = sprintf("wrfchemistry_d0?_%4.4d-%2.2d-%2.2d_%2.2d:00:00",
                       $blksyr,$blksmn,$blksdy,$blkshr);
  system("rm $firstfile1 $firstfile2 $firstfile3");
  }

  system("mv wrf?hr_d* $transdir/wrfdetail"); #== 0;
  system("mv wrfchemistry* $transdir/wrfdetail"); #== 0;

  system("rm -f auxhist*");

  copy "namelist.input", "$dirname/namelist.input"
	or die "ERROR: Failed copying namelist.input: $?";

  chdir "$scriptdir";
  return 1;
}


#########################################################################
# Submits the SLURM script via qsub and waits until the job completes
# before returning control to the calling routine. Checks of the queue
# are made every $sec seconds.
#
# wig; 22-Apr-2008
#########################################################################
sub QsubAndWaitSLURM {
  my ($cmd,$sec) = @_;
  my ($foundit, $id, $result);

  if ($machine eq 'yellowstone') {
  $result = `bsub < $cmd 2>&1`;
  $id0 = (split(/ /,$result))[1]; #get job ID number
  $id1= (split(/</,$id0))[1]; #get job ID number
  $id = (split(/>/,$id1))[0]; #get job ID number
  print "Waiting for Bsub job $id to complete ($cmd).";
  } elsif ($machine eq 'cheyenne') {
  $result = `qsub < $cmd 2>&1`;
  chomp $result ;
  $idx = index($result,'.') ;
  $id = substr $result, 0, $idx; #get job ID number
  print "Waiting for QBS job $id to complete ($cmd).";
  } elsif ($machine eq 'aemol') {
  $result = `qsub < $cmd 2>&1`;
  chomp $result ;
  $idx = index($result,'.') ;
  $id = substr $result, 0, $idx; #get job ID number
  print "Waiting for QBS job $id to complete ($cmd).";
  } elsif ($machine eq 'tianhe') {
  $result = `sbatch $cmd 2>&1`;
  chomp $result;
  $id = (split(/ /,$result))[-1]; #get job ID number
  print "Waiting for Sbatch job $id to complete ($cmd).";
  } else {
  $result = `sbatch $cmd 2>&1`;
  chomp $result;
  $id = (split(/ /,$result))[-1]; #get job ID number
  print "Waiting for Sbatch job $id to complete ($cmd).";
  }

  while (1) {
     #-------------------------------------------------------------
     # Check instruction file for changes. Any lines beginning
     # with "#" are ignored.
     if ( -e "AutowrfInstruct" ) {
     open(FH,"< AutowrfInstruct");
     my $stopcmd = 0;
     foreach (<FH>) {
     if ( not /^\#/ ) {
         if ( /stop|die|end|quit/i ) {
              print "Instructed to stop.\n";
              $stopcmd = 1;
              last;
              }
       }
     }
     if ( $stopcmd ) { last; }
     }
     #-------------------------------------------------------------
	sleep $sec;
        if ($machine eq 'eastwind') {
	$foundit = `squeue -a | grep $id`;
        $status = (split(/ /,$foundit))[-2];
        }
        if ($machine eq 'edison') {
	$foundit = `squeue -a | grep $id`;
        $status = (split(/ /,$foundit))[-9];
        }
        if ($machine eq 'yellowstone') {
        $foundit = `bjobs | grep $id`;
        $status = (split(/ /,$foundit))[-2];
        }
        if ($machine eq 'cheyenne') {
        $foundit = `qstat -a | grep $id`;
        $status = (split(/ /,$foundit))[-1];
        }
        if ($machine eq 'aemol') {
        $foundit = `qstat -a | grep $id`;
        $status = (split(/ /,$foundit))[-1];
        }
        if ($machine eq 'tianhe') {
        $foundit = `yhqueue -a | grep $id`;
        $status = (split(/ /,$foundit))[-2];
        }
        if (not $foundit) {last;}
        elsif ($status eq 'C') {last;}
	print ".";
  }
  print " Complete. Moving on...\n";
}

#########################################################################
# Update emission 
#########################################################################
sub updateEmis {

  my ($currentdir, $success);
  my ($semistime, $eemistime, $emistime, $id, $emisyr, $emismn, $emisdy, $inname, $outname);

  $currentdir = `pwd`;
  chdir "$expdir" or die "Could not cd to $expdir: $!\n";

  printf "\nupdateEmis for block %2.2d%2.2d%2.2d%2.2d-%2.2d%2.2d%2.2d%2.2d\n"
    ,&yr2($blksyr),$blksmn,$blksdy,$blkshr
        ,&yr2($blkeyr),$blkemn,$blkedy,$blkehr;

  system("rm -f wrfchemi_*");
  system("rm -f wrffirechemi_*");
  system("rm -f wrfchemi_gocart_bg_d0?");
  
  for ($i=1; $i<=$max_dom; $i++) {

  # Gocart dust 
     $inname  = sprintf("$emisdir/wrfchemi_gocart_bg_d%2.2d", $i);
     $outname  = sprintf("wrfchemi_gocart_bg_d%2.2d", $i);
     system("ln -s $inname $outname");

  # anthropogenic
  if ($hourlyAnthro == 1) {
  $semistime = POSIX::mktime(0,0,$blkshr,$blksdy,$blksmn-1,$blksyr-1900,0,0,-1);
  $eemistime = POSIX::mktime(0,0,$blkehr,$blkedy,$blkemn-1,$blkeyr-1900,0,0,-1);
  $emistime = $semistime ;
  while ($emistime <= $eemistime) {
  $emisyr = (localtime($emistime))[5] + 1900;
  $emismn = (localtime($emistime))[4] + 1;
  $emisdy = (localtime($emistime))[3];
  $emishr = (localtime($emistime))[2];
  printf "\nCreate anthropogenic emission file for the day of %4.4d%2.2d%2.2d%2.2d\n"
    ,$emisyr,$emismn,$emisdy,$emishr;
  $inname  = sprintf("$emisdir/wrfchemi_d%2.2d_%4.4d-%2.2d-%2.2d_*",$i, $emisyr,$emismn,$emisdy);
  system("ln -s $inname ./");
  $emistime = $emistime + 24*3600 ; # one day increment as emission created in daily files
  }
  } else {
  $semistime = POSIX::mktime(0,0,$blkshr,$blksdy,$blksmn-1,$blksyr-1900,0,0,-1);
  $eemistime = POSIX::mktime(0,0,$blkehr,$blkedy,$blkemn-1,$blkeyr-1900,0,0,-1);
  $emistime = $semistime ;
  while ($emistime <= $eemistime) {
  $emisyr = (localtime($emistime))[5] + 1900;
  $emismn = (localtime($emistime))[4] + 1;
  $emisdy = (localtime($emistime))[3];
  printf "\nCreate anthropogenic emission file for the day of %4.4d%2.2d%2.2d\n"
    ,$emisyr,$emismn,$emisdy;
  $inname  = sprintf("$emisdir/wrfchemi_d%2.2d_%4.4d%2.2d", $i,$emisyr,$emismn);
  #$inname  = sprintf("$emisdir/wrfchemi_d01");
  $outname = sprintf("wrfchemi_d%2.2d_%4.4d-%2.2d-%2.2d_00:00:00",$i,$emisyr,$emismn,$emisdy);
  system("ln -s $inname $outname");
  $emistime = $emistime + 24*3600 ; # one day increment as emission created in daily files
  }
  } # hourlyAnthro

  # fire 
  if ($monthlyFire == 1 ) {
  chdir "$emisdir" or die "Could not cd to $emisdir: $!\n";
  system("ls wrffirechemi* >& fireemis.list");
  $semistime = POSIX::mktime(0,0,$blkshr,$blksdy,$blksmn-1,$blksyr-1900,0,0,-1);
  $eemistime = POSIX::mktime(0,0,$blkehr,$blkedy,$blkemn-1,$blkeyr-1900,0,0,-1);
  $emistime = $semistime;
  while ($emistime <= $eemistime) {
  $emisyr = (localtime($emistime))[5] + 1900;
  $emismn = (localtime($emistime))[4] + 1;
  $emisdy = (localtime($emistime))[3];

  $emistime0 = $emistime ;
    while (1) {
     chdir "$emisdir" or die "Could not cd to $emisdir: $!\n";
     $emisyr0 = (localtime($emistime0))[5] + 1900;
     $emismn0 = (localtime($emistime0))[4] + 1;
     $emisdy0 = (localtime($emistime0))[3];
     $emisfile = sprintf("wrffirechemi_d%2.2d_%4.4d%2.2d", $i,$emisyr0,$emismn0);
     $foundit = `grep "$emisfile" fireemis.list`;
      if ($foundit) {
       chdir "$expdir" or die "Could not cd to $expdir: $!\n";
       printf "\nCreate fire emission file for the day of %4.4d%2.2d%2.2d\n"
              ,$emisyr,$emismn,$emisdy;
       $inname  = sprintf("$emisdir/$emisfile") ;
       $outname = sprintf("wrffirechemi_d%2.2d_%4.4d-%2.2d-%2.2d_00:00:00",$i,$emisyr,$emismn,$emisdy);
       system("ln -s $inname $outname");
       last ;
      }
      else {
        $emistime0 = $emistime0 - 1*24*3600 ;
        #$emisyr0 = (localtime($emistime0))[5] + 1900;
        #$emismn0 = (localtime($emistime0))[4] + 1;
        #$emisdy0 = (localtime($emistime0))[3];
        if ($emistime0 < $emistime - 10*24*3600) {    # back searching for 10 days 
          chdir "$expdir" or die "Could not cd to $expdir: $!\n";
          die "Could not find any on-date fire emission data\n";
         }
      }
    } # while (1)

  $emistime = $emistime + 24*3600 ; # one day increment as emission created in daily files
  } #while emistime
  # this is only for reading by real.exe, not used in simulation
  $inname  = sprintf("wrffirechemi_d%2.2d_%4.4d-%2.2d-%2.2d_00:00:00",$i,$blksyr,$blksmn,$blksdy) ;
  $outname  = sprintf("wrffirechemi_d%2.2d",$i) ;
  system("ln -s $inname $outname");
  } #monthlyFire

  if ($weekFire == 1 ) {
  chdir "$emisdir" or die "Could not cd to $emisdir: $!\n";
  system("ls wrffirechemi* >& fireemis.list");
  $semistime = POSIX::mktime(0,0,$blkshr,$blksdy,$blksmn-1,$blksyr-1900,0,0,-1);
  $eemistime = POSIX::mktime(0,0,$blkehr,$blkedy,$blkemn-1,$blkeyr-1900,0,0,-1);
  $emistime = $semistime;
  while ($emistime <= $eemistime) {
  $emisyr = (localtime($emistime))[5] + 1900;
  $emismn = (localtime($emistime))[4] + 1;
  $emisdy = (localtime($emistime))[3];

  $emistime0 = $emistime ;
    while (1) {
     chdir "$emisdir" or die "Could not cd to $emisdir: $!\n";
     $emisyr0 = (localtime($emistime0))[5] + 1900;
     $emismn0 = (localtime($emistime0))[4] + 1;
     $emisdy0 = (localtime($emistime0))[3];
     $emisfile = sprintf("wrffirechemi_d%2.2d_%4.4d%2.2d%2.2d",$i, $emisyr0,$emismn0,$emisdy0);
     $foundit = `grep "$emisfile" fireemis.list`;
      if ($foundit) {
       chdir "$expdir" or die "Could not cd to $expdir: $!\n";
       printf "\nCreate fire emission file for the day of %4.4d%2.2d%2.2d\n"
              ,$emisyr,$emismn,$emisdy;
       $inname  = sprintf("$emisdir/$emisfile") ;
       $outname = sprintf("wrffirechemi_d%2.2d_%4.4d-%2.2d-%2.2d_00:00:00",$i,$emisyr,$emismn,$emisdy);
       system("ln -s $inname $outname");
       last ;
      }
      else {
        $emistime0 = $emistime0 - 1*24*3600 ;
        #$emisyr0 = (localtime($emistime0))[5] + 1900;
        #$emismn0 = (localtime($emistime0))[4] + 1;
        #$emisdy0 = (localtime($emistime0))[3];
        if ($emistime0 < $emistime - 10*24*3600) {    # back searching for 10 days 
          chdir "$expdir" or die "Could not cd to $expdir: $!\n";
          die "Could not find any on-date fire emission data\n";
         }
      }
    } # while (1)
  $emistime = $emistime + 24*3600 ; # one day increment as emission created in daily files
  } #while emistime
  # this is only for reading by real.exe, not used in simulation
  $inname  = sprintf("wrffirechemi_d%2.2d_%4.4d-%2.2d-%2.2d_00:00:00",$i,$blksyr,$blksmn,$blksdy) ;
  $outname  = sprintf("wrffirechemi_d%2.2d",$i) ;
  system("ln -s $inname $outname");
  } #weekFire

  if ($hourlyFire == 1 ) {
  chdir "$emisdir" or die "Could not cd to $emisdir: $!\n";
  $semistime = POSIX::mktime(0,0,$blkshr,$blksdy,$blksmn-1,$blksyr-1900,0,0,-1);
  $eemistime = POSIX::mktime(0,0,$blkehr,$blkedy,$blkemn-1,$blkeyr-1900,0,0,-1);
  $eemistime = $eemistime + 86400;  #one more day to link
  $emistime = $semistime ;
  while ($emistime <= $eemistime) {
   $emisyr = (localtime($emistime))[5] + 1900;
   $emismn = (localtime($emistime))[4] + 1;
   $emisdy = (localtime($emistime))[3];

   $emistime0 = $emistime ;
  #  chdir "$emisdir" or die "Could not cd to $emisdir: $!\n";
     $emisyr0 = (localtime($emistime0))[5] + 1900;
     $emismn0 = (localtime($emistime0))[4] + 1;
     $emisdy0 = (localtime($emistime0))[3];
     $emisfile = sprintf("wrffirechemi_d%2.2d_%4.4d-%2.2d-%2.2d_*", $i,$emisyr0,$emismn0,$emisdy0);
     chdir "$expdir" or die "Could not cd to $expdir: $!\n";
     printf "\nCreate fire emission file for the day of %4.4d%2.2d%2.2d\n"
              ,$emisyr,$emismn,$emisdy;
     system("ln -s $emisdir/$emisfile ./");
   $emistime = $emistime + 24*3600 ; # one day increment 
  } #while emistime
  # this is only for reading by real.exe, not used in simulation
  $inname  = sprintf("wrffirechemi_d%2.2d_%4.4d-%2.2d-%2.2d_00:00:00",$i,$blksyr,$blksmn,$blksdy) ;
  $outname  = sprintf("wrffirechemi_d%2.2d",$i) ;
  system("ln -s $inname $outname");
  } #hourlyFire


  } #max_dom

  chdir "$currentdir";

  $success = 1;
  return $success;

}

#########################################################################
# Update chemistry boundary 
# Currently, the default namelist includes using global data for initial 
# condition as well. This requires 'chem_in_opt=1'. For restart run, leaving
# this on should not matter.
#########################################################################
sub updatebc {

  my ($currentdir, $success);
  my ($sbdtime, $ebdtime, $bdtime, $id, $bdyr, $bdmn, $bddy, $inname, $outname);

  printf "\nupdatebc for block %2.2d%2.2d%2.2d%2.2d-%2.2d%2.2d%2.2d%2.2d\n"
    ,&yr2($blksyr),$blksmn,$blksdy,$blkshr
        ,&yr2($blkeyr),$blkemn,$blkedy,$blkehr;

  $currentdir = `pwd`;

  #Create boundary processing directory
  if (-d $boundarydir) {

   printf "\nUse the existing boundary processing case\n"

  } else {

  chdir "$scriptdir" or die "Could not cd to $scriptdir: $!\n";
  system("mkdir -p '$boundarydir'/data/parent");
  system("mkdir -p '$boundarydir'/data/nest");
  if ($updatebdy eq 1 or $updateic == 1 ) { 
  system("ln -sf $boundarycode/wrfchembc $boundarydir/wrfchembc");
  }
  if ($updatebdy eq 2 or $updateic == 2 ) { 
  system("ln -sf $boundarycode/mozbc $boundarydir/wrfchembc");
  }
  copy "machine/$queueupdatebc", "$boundarydir/$queueupdatebc"
        or die "ERROR: Failed copying $queueupdatebc: $?";
  if ($updateic eq 0) {
  copy "namelist/namelist.input.updatebc.template.$boundarytype", "$boundarydir/namelist.input_template"
         or die "ERROR: Failed copying boundary namelist.input: $?";
  } else {
   if ($mode eq "forecast") {
   copy "namelist/namelist.input.updatebcic.template.$boundarytype", "$boundarydir/namelist.input_template"
         or die "ERROR: Failed copying boundary namelist.input: $?";
   $updateic = 0 ;
   } elsif ($mode eq "continue") {
   copy "namelist/namelist.input.updatebcic.template.$boundarytype", "$boundarydir/namelist.input_template"
         or die "ERROR: Failed copying boundary namelist.input: $?";
   }
  }
 
  } ##boundary case directory exist or not

  if ( $updatebdy==1 or $updateic == 1 ) { #chemical boundary from WRF-Chem 
  chdir "$boundarydir/data/parent" or die "Could not cd to $boundarydir/data/parent: $!\n";
  $sbdtime = POSIX::mktime(0,0,$blkshr,$blksdy,$blksmn-1,$blksyr-1900,0,0,-1);
  $ebdtime = POSIX::mktime(0,0,$blkehr,$blkedy,$blkemn-1,$blkeyr-1900,0,0,-1);
  $ebdtime = $ebdtime+24*3600 ; # one day after, disabled in this case
  system("rm -f wrfout_d01*");
  $bdtime = $sbdtime;
  while ($bdtime <= $ebdtime) {
   $bdyr = (localtime($bdtime))[5] + 1900;
   $bdyr2005 = 2005 ;
   $bdmn = (localtime($bdtime))[4] + 1;
   $bddy = (localtime($bdtime))[3];
   printf "\nCreate parent file for the day of %4.4d%2.2d%2.2d\n"
    ,$bdyr,$bdmn,$bddy;
   if ($boundarycolon eq 1) {
   $inname  = sprintf("wrfout_d01_%4.4d-%2.2d-%2.2d_00:00:00", $bdyr,$bdmn,$bddy);
   } else {
   $inname  = sprintf("wrfout_d01_%4.4d-%2.2d-%2.2d_00.00.00", $bdyr,$bdmn,$bddy);
   }
   $outname  = sprintf("wrfout_d01_%4.4d-%2.2d-%2.2d_00:00:00", $bdyr,$bdmn,$bddy);
   system("ln -s $parentdir/$inname ./$outname"); # need to change the variable "Times" within each file
   $bdtime = $bdtime + 24*3600 ; # one day increment 
  }

  chdir "$boundarydir/data/nest" or die "Could not cd to $boundarydir/data/nest: $!\n";
  system("rm -f met_em.d0?*");
  $bdtime = $sbdtime;
  $sbdyr = (localtime($bdtime))[5] + 1900;
  while ($bdtime <= $ebdtime) {
   $bdyr = (localtime($bdtime))[5] + 1900;
   $bdmn = (localtime($bdtime))[4] + 1;
   $bddy = (localtime($bdtime))[3];
   printf "\nCreate nest met file for the day of %4.4d%2.2d%2.2d\n"
    ,$bdyr,$bdmn,$bddy;
   $inname  = sprintf("met_em.d0?.%4.4d-%2.2d-%2.2d_*", $bdyr,$bdmn,$bddy);
   system("ln -s $expdir/metdata/$inname ./");
   $bdtime = $bdtime + 24*3600 ; # one day increment 
  }
  system("rm -f wrfbdy_d0?");
  system("rm -f wrfinput_d0?");
  system("ln -s $expdir/wrfbdy_d0? ./");
  system("ln -s $expdir/wrfinput_d0? ./");

  chdir "$boundarydir" or die "Could not cd to $boundarydir: $!\n";
 

 for ($i=1; $i<=$max_dom; $i++) {
 
# print "\nChanging time in parent file\n";
# system("ncl changetime.ncl");  # change "Times" within each file
# print "\nUpdating boundary\n";
# system("./wrfchembc < namelist.input > updatebc.out");

  #modify namelist.inupt 
  if ($i eq 1) {
  copy "namelist.input_template", "namelist.input";
  open(FI,"< namelist.input_template")
        or die "Problem opening namelist.input_template, working on $i: $!\n";
  open(FO,"> namelist.input")
        or die "Problem opening namelist.input, working on $i: $!\n";
  foreach (<FI>) {
        if (/simulation_year/) {
          $newline = "simulation_year       = $sbdyr,    \n";
        }
        elsif (/domain/) {
          $newline = "domain      = $i,    \n";
        } 
        elsif (/do_bc/)  {
          $newline = "do_bc       =  $bcvalue, \n";
        }
        else { $newline = $_; }
        print FO "$newline";
  }
  close FI;
  close FO;
  } else {
  copy "namelist.input_template", "namelist.input";
  open(FI,"< namelist.input_template")
        or die "Problem opening namelist.input_template, working on $i: $!\n";
  open(FO,"> namelist.input")
        or die "Problem opening namelist.input, working on $i: $!\n";
  foreach (<FI>) {
        if (/simulation_year/) {
          $newline = "simulation_year       = $sbdyr,    \n";
        }
        elsif (/domain/) {
          $newline = "domain      = $i,    \n";
        }
        elsif (/do_bc/)  {
          $newline = "do_bc       =    .false., \n";
        }
        else { $newline = $_; }
        print FO "$newline";
  }
  } #if $i
  printf "\nUpdating boundary for domain %2.2d\n", $i;


   # Submit the script to the queue and check the results.
   &QsubAndWaitSLURM($queueupdatebc,300);
 } #$i, max_dom
 
  } elsif ( $updatebdy==2 or $updateic == 2 ) { #chemical boundary from MOZART
  chdir "$boundarydir/data/parent" or die "Could not cd to $boundarydir/data/parent: $!\n";
  system("rm -f mozart00??.nc");
  $sbdtime = POSIX::mktime(0,0,$blkshr,$blksdy,$blksmn-1,$blksyr-1900,0,0,-1);
  $ebdtime = POSIX::mktime(0,0,$blkehr,$blkedy,$blkemn-1,$blkeyr-1900,0,0,-1);
  $bdtime = $sbdtime - 24*3600 ; # one day ahead 
  $bdyr0 = 0 ;
  $bdmn0 = 0 ;
  $id = 1 ;
  while ($bdtime <= $ebdtime) {
   $bdyr = (localtime($bdtime))[5] + 1900;
   $bdmn = (localtime($bdtime))[4] + 1;
   $bddy = (localtime($bdtime))[3];
   printf "\nCreate boundary file for the month of %4.4d%2.2d\n"
    ,$bdyr,$bdmn;
   $inname  = sprintf("mozart4geos5_%4.4d%2.2d.nc", $bdyr,$bdmn);
   if ( $id < 10) {
   $outname = sprintf("mozart000%1.1d.nc",$id);
   }
   else {$outname = sprintf("mozart00%2.2d.nc",$id);}
   if ($bdyr ne $bdyr0 or $bdmn ne $bdmn0) {
   system("ln -s $parentdir/$inname $outname");
   $id = $id + 1;
   $bdyr0 = $bdyr ; 
   $bdmn0 = $bdmn ;
   }

   $bdtime = $bdtime + 3600*24 ; # one day increment 
  }

  chdir "$boundarydir/data/nest" or die "Could not cd to $boundarydir/data/nest: $!\n";
  system("rm -f met_em.d0?.*");
  $bdtime = $sbdtime;
  $sbdyr = (localtime($bdtime))[5] + 1900;
  while ($bdtime <= $ebdtime) {
   $bdyr = (localtime($bdtime))[5] + 1900;
   $bdmn = (localtime($bdtime))[4] + 1;
   $bddy = (localtime($bdtime))[3];
   printf "\nCreate nest met file for the day of %4.4d%2.2d%2.2d\n"
    ,$bdyr,$bdmn,$bddy;
   $inname  = sprintf("met_em.d0?.%4.4d-%2.2d-%2.2d_*", $bdyr,$bdmn,$bddy);
   system("ln -s $expdir/metdata/$inname ./");
   $bdtime = $bdtime + 24*3600 ; # one day increment 
  }
  system("rm -f wrfbdy_d0?");
  system("rm -f wrfinput_d0?");
  system("ln -s $expdir/wrfbdy_d0? ./");
  system("ln -s $expdir/wrfinput_d0? ./");

  chdir "$boundarydir" or die "Could not cd to $boundarydir: $!\n";

 for ($i=1; $i<=$max_dom; $i++) {

  if ($i eq 1) {
  #modify namelist.inupt 
  copy "namelist.input_template", "namelist.input";
  open(FI,"< namelist.input_template")
        or die "Problem opening namelist.input_template, working on $i: $!\n";
  open(FO,"> namelist.input")
        or die "Problem opening namelist.input, working on $i: $!\n";
  foreach (<FI>) {
        if (/simulation_year/) {
          $newline = "simulation_year       = $sbdyr,    \n";
        }
        elsif (/domain/) {
          $newline = "domain      = $max_dom,    \n";
        }
        elsif (/do_bc/)  {
          $newline = "do_bc       =  $bcvalue, \n";
        }
        else { $newline = $_; }
        print FO "$newline";
  }
  close FI;
  close FO;
  } else {
  copy "namelist.input_template", "namelist.input";
  open(FI,"< namelist.input_template")
        or die "Problem opening namelist.input_template, working on $i: $!\n";
  open(FO,"> namelist.input")
        or die "Problem opening namelist.input, working on $i: $!\n";
  foreach (<FI>) {
        if (/simulation_year/) {
          $newline = "simulation_year       = $sbdyr,    \n";
        }
        elsif (/domain/) {
          $newline = "domain      = $max_dom,    \n";
        }
        elsif (/do_bc/)  {
          $newline = "do_bc       =    .false., \n";
        }
        else { $newline = $_; }
        print FO "$newline";
  }
  close FI;
  close FO;
  } #if $i,

  printf "\nUpdating boundary for domain %2.2d\n", $i;

   # Submit the script to the queue and check the results.
# # system("mozbc < namelist.input > updatebc.out");
  # Submit the script to the queue and check the results.
   &QsubAndWaitSLURM($queueupdatebc,300);
  } #$i, max_dom

  } #if $updatebdy

  open (FH,"<updatebc.out")
        or die "STOP: did not produce an updatebc.out file.\n";
  $success = 0;
  foreach (<FH>) {
        if (/bc_wrfchem completed successfully/) {
          $success = 1;
          if ($updatebdy eq 1 or $updateic == 1 ){
          system("rm $boundarydir/data/parent/wrfout_d01*") ;
          }
          if ($updatebdy eq 2 or $updateic == 2 ){
          system("rm $boundarydir/data/parent/mozart*") ;
          }
          last;
        }
  }
  close FH;

  if ($success == 0) {system "rm -f $expdir/wrfbdy_d01";}
  if ($success == 0) { die "updatebc crashed\n";}

#  chdir "$expdir" or die "Could not cd to $expdir: $!\n";
#  chdir "$currentdir";

  chdir "$scriptdir" or die "Could not cd to $scriptdir: $!\n";

  return $success;
}

#########################################################################
# Runs real.exe
#########################################################################
sub Realexe {
  my ($currentdir, $i, $j, $k, $name, $newline, $pre, $post, $sifile, $siprd);
  my (@files);

  printf "\nRealexe for block %2.2d%2.2d%2.2d%2.2d-%2.2d%2.2d%2.2d%2.2d\n"
    ,&yr2($blksyr),$blksmn,$blksdy,$blkshr
    ,&yr2($blkeyr),$blkemn,$blkedy,$blkehr;

  $currentdir = `pwd`;
  chdir "$expdir" or die "Could not cd to $expdir: $!\n";

  &SetNamelistDates;

  #modify namelist for processors & restart
  copy "namelist.input", "namelist.input.tmp";
  open(FI,"< namelist.input.tmp")
        or die "Problem opening namelist.input.tmp, working on $i: $!\n";
  open(FO,"> namelist.input")
        or die "Problem opening namelist.input, working on $i: $!\n";
  foreach (<FI>) {
        if (/nproc_x/) {
          $newline = "nproc_x                             = 6,    \n";
        }
        elsif (/nproc_y/) {
          $newline = "nproc_y                             = 4,    \n";
        }
        elsif (/io_form_restart/) {
          $newline = "io_form_res                     = 2, \n";
        }
        elsif (/restart /) {
          $newline = "restart                             = $restartvalue, \n";
        }
        else { $newline = $_; }
        print FO "$newline";
  }
  close FI;
  close FO;

  copy "namelist.input", "namelist.input.tmp";
  open(FI,"< namelist.input.tmp")
        or die "Problem opening namelist.input.tmp, working on $i: $!\n";
  open(FO,"> namelist.input")
        or die "Problem opening namelist.input, working on $i: $!\n";
  foreach (<FI>) {
        if (/io_form_res/) {
          $newline = "io_form_restart                     = 2, \n";
        }
        else { $newline = $_; }
        print FO "$newline";
  }
  close FI;
  close FO;

  if ($type eq "wrfchem") {
  my $auxinput12 = 0;
  copy "namelist.input", "namelist.input.tmp";
  open FI, "<namelist.input.tmp" or die $!;
  my @lines = <FI>;
  close FI or die $!;
  my $idx = 0;
  do {
    if($lines[$idx] =~ /io_form_auxinput12/) {
       $auxinput12 = 1; 
    }
    $idx++;
  } until($idx >= @lines);

  if ($mode eq "forecast" and $updateic eq 0 and $auxinput12 eq 0 and $firstforecast ne 1) {
  copy "namelist.input", "namelist.input.tmp";
  open FI, "<namelist.input.tmp" or die $!;
  my @lines = <FI>;
  close FI or die $!;
  my $idx = 0;
  do {
    if($lines[$idx] =~ /debug_level/) {
        splice @lines, $idx, 0, "auxinput12_inname  = 'wrf_chem_input',\n";
        $idx++;
    } elsif ($lines[$idx] =~ /io_form_auxinput10/) {
        splice @lines, $idx, 0, "io_form_auxinput12  = 2,\n";
        $idx++;
    }
    $idx++;
  } until($idx >= @lines);
  open FO, ">namelist.input" or die $!;
  print FO join("",@lines);
  close FO;
  }
  }

  system "rm -f rsl.*";
      &QsubAndWaitSLURM($queuereal,30);
  open (FH,"<rsl.error.0000")
	or die "STOP: real.exe did not produce an rsl.error.0000 file.\n";
  $i = 0;
  foreach (<FH>) {
	if (/SUCCESS COMPLETE REAL_EM/) {
	  $i = 1;
	  last;
	}
  }
  close FH;

  chdir "$scriptdir";
  return $i;
}

#########################################################################
# Sets the start and end dates in namelist.input given the block start
# and end dates. Note that this is currently hard coded for a maximum of
# 3 domains. The format strings need to be extended if more domains are
# going to be used.
#########################################################################
sub SetNamelistDates {
  my $newline;

  if ($max_dom > 3) {
	die "STOP: The format strings for the dates in SetNamelistDates need to be extended to handle $max_dom domains.\n";
  }

  copy "namelist.input", "namelist.input.tmp";
  open(FI,"< namelist.input.tmp")
	or die "Problem opening namelist.input.tmp, working on $i: $!\n";
  open(FO,"> namelist.input")
	or die "Problem opening namelist.input, working on $i: $!\n";
  foreach (<FI>) {
        if (/run_days/) {
          $newline = sprintf("run_days                            = %2.2d,\n",
                                                 $blkday);
        }
	elsif (/start_year/) { 
	  $newline = sprintf("start_year                          = %4.4d, %4.4d, %4.4d,\n",
						 $blksyr, $blksyr, $blksyr);
	}
	elsif (/start_month/) { 
	  $newline = sprintf("start_month                         = %2.2d,   %2.2d,   %2.2d,\n",
						 $blksmn, $blksmn, $blksmn);
	}
	elsif (/start_day/) { 
	  $newline = sprintf("start_day                           = %2.2d,   %2.2d,   %2.2d,\n",
						 $blksdy, $blksdy, $blksdy);
	}
	elsif (/start_hour/) { 
	  $newline = sprintf("start_hour                          = %2.2d,   %2.2d,   %2.2d,\n",
						 $blkshr, $blkshr, $blkshr);
	}
	elsif (/end_year/) { 
	  $newline = sprintf("end_year                            = %4.4d, %4.4d, %4.4d,\n",
						 $blkeyr, $blkeyr, $blkeyr);
	}
	elsif (/end_month/) { 
	  $newline = sprintf("end_month                           = %2.2d,   %2.2d,   %2.2d,\n",
						 $blkemn, $blkemn, $blkemn);
	}
	elsif (/end_day/) { 
	  $newline = sprintf("end_day                             = %2.2d,   %2.2d,   %2.2d,\n",
						 $blkedy, $blkedy, $blkedy);
	}
	elsif (/end_hour/) { 
	  $newline = sprintf("end_hour                            = %2.2d,   %2.2d,   %2.2d,\n",
						 $blkehr, $blkehr, $blkehr);
	}
	  else { $newline = $_; }
	  print FO "$newline";
	}
	close FI;
	close FO;
}


#########################################################################
# Runs WRF
#########################################################################
sub WRF {
  my ($currentdir, $foundit, $newline, $sucess);

  printf "\nWRF for block %2.2d%2.2d%2.2d%2.2d-%2.2d%2.2d%2.2d%2.2d\n"
    ,&yr2($blksyr),$blksmn,$blksdy,$blkshr
	,&yr2($blkeyr),$blkemn,$blkedy,$blkehr;

  $currentdir = `pwd`;
  chdir "$expdir" or die "Could not cd to $expdir: $!\n";

### Setup the namelist. Note that the history interval is forced to
  # hourly in order to simplfy handling of the soil cycling. I do not
  # want to take the time to code up parsing times in netCDF files. If
  # more than 3 domains are used, the history line needs to be extended.
  &SetNamelistDates;
  copy "namelist.input", "namelist.input.tmp";
  open(FI,"< namelist.input.tmp")
	or die "Problem opening namelist.input.tmp, working on $i: $!\n";
  open(FO,"> namelist.input")
	or die "Problem opening namelist.input, working on $i: $!\n";
  foreach (<FI>) {
	if (/max_dom/) {
	  $newline = "max_dom                             = $max_dom,\n";
	}
#       elsif (/history_interval/) {
#         $newline = "history_interval                    = 1440, 60, 60,\n";
#       }
        elsif (/input_pft/) {
          $newline = "input_pft                           = $inputpftvalue, \n";
        }
        elsif (/nproc_x/) {
          $newline = "nproc_x                             = 32,    \n";
        }
        elsif (/nproc_y/) {
          $newline = "nproc_y                             = 15,    \n";
        }
        else { $newline = $_; }
	print FO "$newline";
  }
  close FI;
  close FO;

### Submit the WRF script to the queue and check the results.
    &QsubAndWaitSLURM($queuewrf,300);

  if (not -e "rsl.error.0000") {
	die "STOP: wrf.exe did not produce an rsl.error.0000 file.\n";
  }
  $foundit = `grep "SUCCESS COMPLETE WRF" rsl.error.0000`;
  if ($foundit) {
	$sucess = 1;
  }
  else {
	$sucess = 0;
  }

  chdir "$scriptdir";
  return $sucess;
}

#########################################################################
# Convert a 4 digit date to 2 digits
#########################################################################
sub yr2 {
  if ( $_[0] < 2000 ) {
    return $_[0]-1900;
  }
  else {
    return $_[0]-2000;
  }
}

#########################################################################
# link the inital files 
#########################################################################
sub linkinit {

  my ($success);

  printf "\nlink initial files for case $expdir\n" ;

  chdir "$expdir" or die "Could not cd to $expdir: $!\n";

  ### link the initial files 
  #system("rm -f wrf_chem_input*");
  # $initime0 =  POSIX::mktime(0,0,$blkshr,$blksdy,$blksmn-1,$blksyr-1900,0,0,-1);
  $initime =  POSIX::mktime(0,0,$blkshr,$blksdy,$blksmn-1,$blksyr-1900,0,0,-1);
  $initime = $initime - (3600*$nextcyclehr) ;
  # before 2007, time always miss the 2am in April 2 for unknown reason 
  # check if the time is qualified for one hour back
  #$marktime =  POSIX::mktime(0,0,2,2,4-1,$blksyr-1900,0,0,-1); # April 2 2:00
  #if ($blksyr < 2007 && $initime0 > $marktime && $initime < $marktime ) 
  #{$initime = $initime + 3600;
  $initime = $initime + 3600;
  $inityr = (localtime($initime))[5] + 1900;
  $initmn = (localtime($initime))[4] + 1;
  $initdy = (localtime($initime))[3];

  $initdir = sprintf("$outputdir/wrf%4.4d%2.2d%2.2d%2.2d",$inityr,$initmn,$initdy,$blkshr ) ;
  $inname = sprintf("$initdir/wrfout_d01_%4.4d-%2.2d-%2.2d_%2.2d:00:00", $blksyr,$blksmn,$blksdy,$blkshr);
  $outname = sprintf("wrf_chem_input_d01");
  if ($restart eq 2) { 
    system("rm -f wrf_chem_input*");
    system("ln -s $inname $outname"); 
  }

  ### link restart files 
  system("rm -f wrfrst*");
  $initdir = sprintf("$outputdir/wrf%4.4d%2.2d%2.2d%2.2d",$inityr,$initmn,$initdy,$blkshr ) ;
  $inname = sprintf("$initdir/wrfrst_d01_%4.4d-%2.2d-%2.2d_%2.2d:00:00", $blksyr,$blksmn,$blksdy,$blkshr);
  $outname = sprintf("wrfrst_d01_%4.4d-%2.2d-%2.2d_%2.2d:00:00", $blksyr,$blksmn,$blksdy,$blkshr);
  if ($restart eq 1) { system("ln -s $inname $outname"); } 

  chdir "$scriptdir" or die "Could not cd to $scriptdir: $!\n";

  $success = 1;
  return $success;
}


#########################################################################
# link the run files 
#########################################################################
sub linkfiles {

  my ($success);

  printf "\nlink files for case $expdir \n" ;

  chdir "$expdir" or die "Could not cd to $expdir: $!\n";
  
  $modelmaindir = "$modeldir/main";
  $modelrundir  = "$modeldir/run";

  system("rm -f wrf.exe ");
  system("ln -s $modelmaindir/wrf.exe ./");
  system("rm -f real.exe");
  system("ln -s $modelmaindir/real.exe ./");
  system("rm -f tc.exe");
  system("ln -s $modelmaindir/tc.exe ./");
  system("rm -f ndown.exe");
  system("ln -s $modelmaindir/ndown.exe ./");

  system("rm -f bulk*");
  system("ln -s $modelrundir/bulk* ./");
  system("rm -f capacity.asc");
  system("ln -s $modelrundir/capacity.asc ./");
  system("rm -f coeff*");
  system("ln -s $modelrundir/coeff* ./");
system("rm -f constants.asc");
system("ln -s $modelrundir/constants.asc ./");
system("rm -f GENPARM.TBL");
system("ln -s $modelrundir/GENPARM.TBL ./");
system("rm -f kernels*");
system("ln -s $modelrundir/kernels* ./");
system("rm -f LANDUSE.TBL");
system("ln -s $modelrundir/LANDUSE.TBL ./");
system("rm -f masses.asc");
system("ln -s $modelrundir/masses.asc ./");
system("rm -f MPTABLE.TBL");
system("ln -s $modelrundir/MPTABLE.TBL ./");
system("rm -f SOILPARM.TBL");
system("ln -s $modelrundir/SOILPARM.TBL ./");
system("rm -f termvels.asc");
system("ln -s $modelrundir/termvels.asc ./");
system("rm -f URBPARM.TBL");
system("ln -s $modelrundir/URBPARM.TBL ./");
system("rm -f VEGPARM.TBL");
system("ln -s $modelrundir/VEGPARM.TBL ./");
system("rm -f CAM_ABS_DATA");
system("ln -s $modelrundir/CAM_ABS_DATA ./");
system("rm -f CAM_AEROPT_DATA");
system("ln -s $modelrundir/CAM_AEROPT_DATA ./");
system("rm -f ETAMPNEW_DATA");
system("ln -s $modelrundir/ETAMPNEW_DATA ./");
system("rm -f grib2map.tbl");
system("ln -s $modelrundir/grib2map.tbl ./");
system("rm -f gribmap.txt");
system("ln -s $modelrundir/gribmap.txt ./");
system("rm -f ozone.formatted");
system("ln -s $modelrundir/ozone.formatted ./");
system("rm -f ozone_lat.formatted");
system("ln -s $modelrundir/ozone_lat.formatted ./");
system("rm -f ozone_plev.formatted");
system("ln -s $modelrundir/ozone_plev.formatted ./");
system("rm -f README.namelist");
system("ln -s $modelrundir/README.namelist ./");
system("rm -f RRTM_DATA");
system("ln -s $modelrundir/RRTM_DATA ./");
system("rm -f RRTMG_LW_DATA");
system("ln -s $modelrundir/RRTMG_LW_DATA ./");
system("rm -f RRTMG_SW_DATA");
system("ln -s $modelrundir/RRTMG_SW_DATA ./");
system("rm -f tr49t67");
system("ln -s $modelrundir/tr49t67 ./");
system("rm -f tr49t85");
system("ln -s $modelrundir/tr49t85 ./");
system("rm -f tr67t85");
system("ln -s $modelrundir/tr67t85 ./");
system("rm -f *.TBL");
system("ln -s $modelrundir/*.TBL ./");
system("rm -f CLM*");
system("ln -s $modelrundir/CLM* ./");
system("rm -f CAMtr_volume_mixing_ratio*");
system("ln -s $modelrundir/CAMtr_volume_mixing_ratio* ./");
system("rm -f CCN_ACTIVATE.BIN");
system("ln -s $modelrundir/CCN_ACTIVATE.BIN ./");
system("rm -f aerosol*formatted");
system("ln -s $modelrundir/aerosol*formatted ./");
system("rm -f ETAMPNEW_DATA.expanded_rain");
system("ln -s $modelrundir/ETAMPNEW_DATA.expanded_rain ./");
system("rm -f megan21_emis_factors_c20130304.nc");
system("ln -s $modelrundir/megan21_emis_factors_c20130304.nc ./");

  chdir "$scriptdir" or die "Could not cd to $scriptdir: $!\n";

  $success = 1;
  return $success;

}
