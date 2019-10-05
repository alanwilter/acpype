#!/bin/csh
# This script print SGLD properties from sander simulation output
# csh
if ($#argv != 1 )then
  echo " Usage: sgldinfo.sh sander_output "
  exit
endif

set z = $argv[1]

awk -v sanderfile=$z  'BEGIN{t=-1;nt=0; \
while(getline <sanderfile){if($1=="Local"&&$2=="averaging")tlavg=$4;\
if($1=="NSTEP"){if(t!=$6){nt++;t=$6;}}} \
close(sanderfile);sta=nt/2;t=-1;nt=0;m=0;} \
{if($1=="NSTEP"){skip=0;if(t==$6)skip=1; \
else {nt++;t=$6;temp=$9;if(nt<sta)skip=1;}} if(skip)next; \
if($1=="SGLF"){sgft=$3;tempsg=$4;tsglf=$5;trlf=$6;cflf=$7;elf=$8;wsg=$9;} \
if($1=="SGHF"){sgff=$3;sgfd=$4;tsghf=$5;trhf=$6;cfhf=$7;ehf=$8;vsg=$9; \
avgtmp+=temp;avgsgft+=sgft;avgsgff+=sgff;avgsgfd+=sgfd; avgtsg+=tempsg; \
avgelf+=elf;avgtlf+=tsglf;avgrlf+=trlf;avgcflf+=cflf; avgwsg+=wsg;\
avgehf+=ehf;avgthf+=tsghf;avgrhf+=trhf;avgcfhf+=cfhf;avgvsg+=vsg;m++;}} \
END{if(m==0){printf(" This is not a SGMD/SGLD simulation: %-s\n",sanderfile);\
exit;} \
avgtmp/=m;avgsgft/=m;avgsgff/=m;avgsgfd/=m; avgtsg/=m; \
avgelf/=m;avgtlf/=m;avgrlf/=m;avgcflf/=m; avgwsg/=m;\
avgehf/=m;avgthf/=m;avgrhf/=m;avgcfhf/=m;avgvsg/=m; \
printf("  SGLD properties from SANDER output: %-40s \n",sanderfile); \
printf("   NPRINT         TIME       TEMP    TSGAVG   \n");    \
printf(" %8d   %10.2f   %8.2f  %8.4f \n",nt,t,avgtmp,tlavg); \
printf("  SGFT    TEMPSG    TEMPLF     TREFLF    FRICLF         EPOTLF     WTSG\n"); \
printf("%6.4f  %8.2f  %8.2f   %8.2f  %8.4f %14.2f %8.4f\n",avgsgft,avgtsg,avgtlf,avgrlf,avgcflf,avgelf,avgwsg); \
printf("  SGFF      SGFD    TEMPHF     TREFHF    FRICHF         EPOTHF    VIRSG\n"); \
printf("%6.4f  %8.2f  %8.2f   %8.2f  %8.4f %14.2f %8.4f\n",avgsgff,avgsgfd,avgthf,avgrhf,avgcfhf,avgehf,avgvsg); \
printf("  TEMPLF=  %8.2f \n",avgtlf); \
printf("  TEMPHF=  %8.2f \n",avgthf); \
}'   $z



