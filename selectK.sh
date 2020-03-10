#!/bin/bash

# This script is aimed to automate the fastPASE K selection procedure as described at https://github.com/inzilico/imputeqc
# Contributed by Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

#Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

#Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

#Initialize default values.
ncores=2
kmin=10
kmax=15
repeats=30
perc=0.05

#Help function
function HELP {
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo "${REV}-i${NORM}  --${BOLD}input${NORM} file name ${BOLD}without the .inp${NORM} extension."
  echo "${REV}-c${NORM}  --number of ${BOLD}CPU cores${NORM}. Default is $ncores."
  echo "${REV}-m${NORM}  --minimal number of ${BOLD}cluster${NORM} to test with ${BOLD}fastPHASE${NORM}. Default is $kmin."
  echo "${REV}-x${NORM}  --maximal number of ${BOLD}cluster${NORM} to test with ${BOLD}fastPHASE${NORM}. Default is $kmax."
  echo "${REV}-r${NORM}  --number of ${BOLD}replicates${NORM} to make with ${BOLD}make_test_files.R${NORM}. Default is $repeats."
  echo "${REV}-p${NORM}  --percent of ${BOLD}masked sites${NORM} in test files. Default is $perc."
  echo -e "${REV}-h${NORM}  --Displays this help message."\\n
  echo -e "Example: ${BOLD}./$SCRIPT -i test -c 2 -n 10 -x 15 -r 30 -p 0.05${NORM}"\\n
  echo -e "${BOLD}NOTE:${NORM} Change the hard-coded location of R scripts by replacing '~/Scripts/git/imputeqc/inst/extdata/' with your path in ./$SCRIPT."\\n
  exit 1
}

#Check if any argument is passed.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  HELP
fi

#Input arguments.
while getopts i:c:m:x:r:p:h: option
    do
        case $option in
            i) input=$OPTARG;;
            c) ncores=$OPTARG;;
            m) kmin=$OPTARG;;
            x) kmax=$OPTARG;;
            r) repeats=$OPTARG;;
            p ) perc=$OPTARG;;
            h) HELP
            exit 0;;
        esac
    done

#Check the input file argument.
if [ ! -f $input.inp ]; then
    echo "The input file is not found! See below how to specify the input file name."
    HELP
    exit
fi

# find make_test_files.R
filesR="~/Scripts/git/imputeqc/inst/extdata/"

# Generate the test files
Rscript $filesR/make_test_files.R -n $repeats $input.inp
Rscript $filesR/make_test_files.R -n $repeats -p $perc -o masked/$input\_masked $input.inp 

# run fastPHASE for the range of K.
for k in $(seq $kmin $kmax)
    do
        parallel --jobs $ncores ~/Programs/fastPHASE/fastPHASE -K$k -H-4 -n -Z -o$k $input/$input.m{}.inp masked/$input_masked.m{}.inp ::: $(seq 1 $repeats);
    done

# estimate quality of each run.
echo > $input\_estimateQuality.txt
for f in *_genotypes.out
    do
        Rscript inst/extdata/estimateQuality.R $input.inp $f masked/masks.RDS >> $input\_estimateQuality.txt 2>&1
    done

# extract the discordance score.
grep Discordance $input\_estimateQuality.txt > $input\_discordance.txt

echo -e "\nSelection of K is done."
echo "Check the output files" $input"_estimateQuality.txt and" $input"_discordance.txt"
exit 0
