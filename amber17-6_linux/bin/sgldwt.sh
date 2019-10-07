#!/bin/csh
# This script calculate SGLD weight from sander simulation output
# csh
if ($#argv != 1 &&$#argv != 2 &&$#argv != 3 )then
  echo " Usage: sgldwt.sh sander_output [treflf [trefhf]] "
  exit
endif

set z = $argv[1]
set trlf = 0
set trhf = 0
if ($#argv >= 2 )then
set trlf = $argv[2]
endif
if ($#argv >= 3 )then
set trhf = $argv[3]
endif

#echo $trlf $trhf

awk -v sanderfile=$z -v trefl=$trlf -v trefh=$trhf  'BEGIN{t=-1;nt=0; \
while(getline <sanderfile){if($1=="NSTEP"){if(t!=$6){nt++;t=$6;}}} \
close(sanderfile);sta=nt/2;t=-1;nt=0;m=0;} \
{if($1=="NSTEP"){skip=0;if(t==$6)skip=1; \
else {nt++;t=$6;temp=$9;if(nt<sta)skip=1;}} if(skip)next; \
if($1=="SGLF"){tsglf=$5;trlf=$6;cflf=$7;elf=$8;wsg=$9;} \
if($1=="SGHF"){tsghf=$5;trhf=$6;cfhf=$7;ehf=$8;vsg=$9; \
avgtmp+=temp;avgelf+=elf;avgtlf+=tsglf;avgrlf+=trlf;avgflf+=cflf; \
avgthf+=tsghf;avgrhf+=trhf;avgfhf+=cfhf;avgvsg+=vsg;m++;}} \
END{avgtmp/=m;kt=0.001987*avgtmp; avgelf/=m;avgtlf/=m;avgrlf/=m; \
avgflf/=m;avgehf/=m;avgthf/=m;avgrhf/=m;avgfhf/=m;avgvsg/=m;  \
if(trefl>0){avgrlf=trefl;} if(trefh>0)avgrhf=trefh; t=-1; \
while(getline <sanderfile){if($1=="NSTEP"){dowt=0;if(t!=$6){nt++;t=$6;dowt=1;}} \
if(dowt){if($1=="SGLF"){elf=$8;} \
if($1=="SGHF"){ehf=$8;vsg=$9; \
wt=exp(((avgflf*avgrlf/avgtlf-1)*(elf-avgelf)+(avgfhf*avgrhf/avgthf-1)*(ehf-avgehf)+vsg-avgvsg)/kt); \
printf("%10.2f %14.6f\n",t,wt);}}}}'   $z



