#!/usr/bin/perl -w
#########################################################################
# This script runs a series of WRF runs, doing cold starts on the
# meteorology and hot starts on the chemistry 
#
# William.Gustafson@pnl.gov; 29-Sep-2005
# Chun Zhao  ; 25-Aug-2009 
#
# History:
#   29-Sep-2005  Working copy for WRF v2.1.
#   10-Nov-2006  Added option to cycle clouds.
#   22-Apr-2008  Added option to choose either PBS or SLURM queueing
#   25-Aug-2009  Modified for WRF-CHEM v3.1.1
#########################################################################


#######################
#  USER DEFINED INFO  #
#######################
$version       = "v3.5.1_20180208" ;
$mode          = "forecast" ; #forecast or continue 


# Case directory path...
$casedir       = "EChina/12km/120x100/$version/2015/APR/ms34_clm.mynn_wrficbdy_2015emis_adjust";

# Experiment directory head...
$expdir        = "/home/chunzhao/Work/WRFCHEM/run_cases/exps/$casedir";

# Output log directory...
$outputdir     = "$expdir/output";

# Archive directory (for data merge to HPSS) 
$transdir     = "/home/chunzhao/Work/WRFCHEM/output/$casedir";

# Run start date (use a 4 digit year)...
$runsyr        = 2015;
$runsmn        = 4;
$runsdy        = 11;
$runshr        = 0;
# Run end date (use a 4 digit year)...
$runeyr        = 2015;
$runemn        = 4;
$runedy        = 16;
$runehr        = 0;
# Blocking options for each cycle of the run...
#     block unit: 0 = Hours
#     block length is the number of units per block
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
} elsif ($mode eq "continue") {
$nextcyclehr   = $blklen ; # days to restart 
$restart   = 0 ; # 0:restart since second cycle, 1:restart since first cycle, 2:never restart 
}

# Number of domains to run...
$max_dom       = 1;
# Parts of the model to run (0=No, 1=Yes, 2=Yes, but not 1st time through loop)
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

  $blkeyr  = $runeyr;
  $blkemn  = $runemn;
  $blkedy  = $runedy;
  $blkehr  = $runehr;

  $inittime = &DecBlk($blkehr,$blkedy,$blkemn,$blkeyr); #the end of the current run

chdir "$expdir";

  printf "\nStarting block %2.2d%2.2d%2.2d%2.2d-%2.2d%2.2d%2.2d%2.2d\n"
    ,&yr2($blksyr),$blksmn,$blksdy,$blkshr
    ,&yr2($blkeyr),$blkemn,$blkedy,$blkehr;

###  Do MoveOutput
  if ( $runMoveOutput==1 ) {
    &MoveOutput or die "MoveOutput crashed\n";
  }

  printf "\nCompleted block %2.2d%2.2d%2.2d%2.2d-%2.2d%2.2d%2.2d%2.2d\n"
    ,&yr2($blksyr),$blksmn,$blksdy,$blkshr
	,&yr2($blkeyr),$blkemn,$blkedy,$blkehr;

print "run completed\n";

#########################################################################
#  Sets the next block calendar start time
#########################################################################
sub DecBlk {
  my ($hr,$dy,$mn,$yr) = @_;
  my $newtime;

  $newtime =  POSIX::mktime(0,0,$hr-$blklen,$dy,$mn-1,$yr-1900,0,0,-1);
  $lastyr = (localtime($newtime))[5] + 1900;
  $lastmn = (localtime($newtime))[4] + 1;
  $lastdy = (localtime($newtime))[3];
  $lasthr = (localtime($newtime))[2];
  return $newtime;
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

  $currentdir = `pwd`;
  chdir "$expdir" or die "Could not cd to $expdir: $!\n";

  $dirname = sprintf("$outputdir/wrf%4.4d%2.2d%2.2d%2.2d",
					 $lastyr,$lastmn,$lastdy,$blkshr);
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

#turn off die for restarts since they are not always saved
#	or die "ERROR: Failed moving wrfrst files: $?";
  copy "namelist.input", "$dirname/namelist.input"
	or die "ERROR: Failed copying namelist.input: $?";

  chdir "$currentdir";
  return 1;
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
