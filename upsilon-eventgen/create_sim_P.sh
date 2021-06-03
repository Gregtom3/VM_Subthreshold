#!/bin/bash
#----------------------------------------#
#  1. Global Variables 
#     (all directories begin at '/' )
#----------------------------------------#

# The directory this file is stored
CURRENT_DIR="/sphenix/user/gregtom3/Duke/VM_Subthreshold/upsilon-eventgen"

# Where the bulk simulation data will be saved
# NOTE: Combined files stored elsewhere
RUNCARD_DIR="/sphenix/user/gregtom3/data/Fall2020/sartre_smear_simulations"

# The email you would like the condor jobs to be sent to
JOB_EMAIL="gregcondorjobs@gmail.com"

# The location of your Sartre 1.33 'examples' directory
# NOTE: Add 'sartreSMEAR.cpp' to this directory
#       The file can be copied from user 'gregtom3'
#       The file will need to be appended to CMakeLists.txt
#SARTRE_EXAMPLES_DIR="/sphenix/user/gregtom3/matousek_sartre_1.33/new_sartre/sartre/examples"

# Your default home directory
HOME_DIR="/phenix/u/gregtom3"

# Where you have 'libeicsmear.so' locally installed
# INSTALL_DIR="/sphenix/user/gregtom3/install/lib"


# Where you have the includes for 'libeicsmear.so'
# INCLUDE_DIR="/sphenix/user/gregtom3/install/include"


# Location of the 'eictree' directory in this repo
#EICTREE_DIR="/sphenix/user/gregtom3/ephenix-sbu/analysis/gregory_matousek/vlad_dvmp/eictree"

# Location of the 'smear' directory in this repo
#DET_DIR="/sphenix/user/gregtom3/ephenix-sbu/analysis/gregory_matousek/vlad_dvmp/smear"
#----------------------------------------#
#  2. Simulation Parameters
#----------------------------------------#
ELECTRON_ENERGY=10
HADRON_ENERGY=100

NEVENTS=1000
NBATCHES=5

# User gregtom3 edited Sartre's source code 
# to permit T restrictions

Q2MIN=0.00
Q2MAX=1.00

YMIN=0.01
YMAX=0.10

#----------------------------------------#
#  3. DVMP Types to Simulate
#----------------------------------------#
#DO_EP=1
#DO_ECA=1
#DO_EAU=1

#DO_JPSI=1
#DO_PHI=1
#DO_RHO=0
#DO_UPSILON=0

#DO_EE=1
#DO_MU=0
#DO_PION=0
#DO_KAON=1

#DO_BNONSAT=1
#DO_BSAT=1
#----------------------------------------#
#  4. DETECTOR SYSTEM FOR SMEARING
#----------------------------------------#
#DO_HANDBOOK=1
#DO_BEAST=0
#DO_EPHENIX=0
#DO_PERFECT=0
#----------------------------------------#
#  NO NEED TO EDIT BEYOND THIS
#  NO NEED TO EDIT BEYOND THIS
#  NO NEED TO EDIT BEYOND THIS
#----------------------------------------#

# Get total number of simulations
#-----------------------------------
declare -i NSIMS
NSIMS=${NEVENTS}/${NBATCHES}

# Generate simulation title
#-----------------------------------
TITLE="P-${NEVENTS}-${NBATCHES}_${ELECTRON_ENERGY}x${HADRON_ENERGY}_Y_${YMIN}_${YMAX}_Q2_${Q2MIN}_${Q2MAX}"

# Go back and create condor directory
#-----------------------------------
cd $CURRENT_DIR
OUTPUT_DIR="${CURRENT_DIR}/${TITLE}"
rm -rf $OUTPUT_DIR
mkdir $OUTPUT_DIR
ULTIMATE_COMBINE="${OUTPUT_DIR}/ultimateCombine.sh"
ULTIMATE_COMBINE2="${OUTPUT_DIR}/ultimateCombine2.sh"
CLEAR_FILES="${OUTPUT_DIR}/clearFiles.sh"
# Create sims
#-----------------------------------
######################################################
# Path to simulation runcard

# Create necessary condor filenames
SARTRE_BASE="${OUTPUT_DIR}/test"
SARTRE_SH="${SARTRE_BASE}.sh"
SARTRE_OUT="${SARTRE_BASE}.out"
SARTRE_LOG="${SARTRE_BASE}.log"
SARTRE_ERR="${SARTRE_BASE}.err"
SARTRE_JOB="${SARTRE_BASE}.job"

# Generate files
touch $SARTRE_SH
touch $SARTRE_OUT
touch $SARTRE_LOG
touch $SARTRE_ERR
touch $SARTRE_JOB

# Append the job file
echo "Executable = ${SARTRE_SH} " >> $SARTRE_JOB
echo "PeriodicHold = (NumJobStarts>=3 && JobStatus == 1)" >> $SARTRE_JOB
echo "Output = ${SARTRE_OUT} " >> $SARTRE_JOB
echo "Error = ${SARTRE_ERR} " >> $SARTRE_JOB
echo "Log = ${SARTRE_LOG} " >> $SARTRE_JOB
echo "Universe = vanilla " >> $SARTRE_JOB
echo "Priority = +0" >> $SARTRE_JOB
echo "Input = /dev/null" >> $SARTRE_JOB
echo "GetEnv = False" >> $SARTRE_JOB
echo "Initialdir = ${SARTRE_EXAMPLES_DIR}" >> $SARTRE_JOB
echo "+Experiment= \"sphenix\"" >> $SARTRE_JOB
echo "+Job_Type = \"cas\"" >> $SARTRE_JOB
echo "Notify_user = ${JOB_EMAIL}" >> $SARTRE_JOB

# One job file will simulate many batches
declare -i VAR
VAR=$NBATCHES-1
for m in $(seq 0 1 ${VAR})
do
    echo "Arguments = ${m}" >> $SARTRE_JOB
    echo "Queue" >> $SARTRE_JOB
done


# Append the sh file
echo "#!/bin/tcsh" >> $SARTRE_SH
echo "setenv HOME ${HOME_DIR}" >> $SARTRE_SH
echo "source /etc/csh.login" >> $SARTRE_SH
echo "foreach p (/etc/profile.d/*.csh)" >> $SARTRE_SH
echo "    source \$p" >> $SARTRE_SH
echo "end" >> $SARTRE_SH
echo "source /opt/sphenix/core/bin/sphenix_setup.csh" >> $SARTRE_SH
echo "setenv LD_LIBRARY_PATH \"${INSTALL_DIR}:$LD_LIBRARY_PATH\"" >> $SARTRE_SH
echo "setenv GSL_DIR /opt/sphenix/core/gsl" >> $SARTRE_SH
echo "cd ${CURRENT_DIR}" >> $SARTRE_SH
echo "root -l electro_eic_P_condor.C\\(${NEVENTS},${NBATCHES},\$1,${ELECTRON_ENERGY},${HADRON_ENERGY},${YMIN},${YMAX},${Q2MIN},${Q2MAX}\\)" >> $SARTRE_SH

condor_submit ${SARTRE_JOB}

cd ${OUTPUT_DIR}

COMBINE_SH="combine.sh"
touch $COMBINE_SH

echo "#!/bin/tcsh" >> $COMBINE_SH
echo "root -l ../macros/combineTree.C\\(${NBATCHES},\\\"${OUTPUT_DIR}\\\",\\\"${OUTPUT_DIR}\\\"\\)" >> $COMBINE_SH
echo "rm run*.root" >> $COMBINE_SH

chmod +x ./*
