#!/bin/bash
#----------------------------------------#
#  1. Global Variables 
#     (all directories begin at '/' )
#----------------------------------------#
# The directory this file is stored
CURRENT_DIR="/sphenix/user/gregtom3/Duke/VM_Subthreshold/psi2S-eventgen"
# The email you would like the condor jobs to be sent to
JOB_EMAIL="gregcondorjobs@gmail.com"
# Your default home directory
HOME_DIR="/phenix/u/gregtom3"
#----------------------------------------#
#  2. Simulation Parameters
#----------------------------------------#
NEVENTS_PHOTO=5000000
NEVENTS_ELECTRO=100000000
EBEAM=(12 15 17)
# Get total number of simulations
#-----------------------------------

# Generate simulation title
#-----------------------------------
NUCLEUS=("p" "D")
PHOTON=("photo" "electro")

for ebeam in "${EBEAM[@]}"
do
    for nucleus in "${NUCLEUS[@]}"
    do
	for photon in "${PHOTON[@]}"
	do	    
	    declare -i NEVENTS
	    if [ "${photon}" == "photo" ]; then
		NEVENTS=${NEVENTS_PHOTO}
	    else
		NEVENTS=${NEVENTS_ELECTRO}
	    fi
	    TITLE="${nucleus}-${photon}-${NEVENTS}-${ebeam}"
	    # Go back and create condor directory
	    #-----------------------------------
	    cd $CURRENT_DIR
	    OUTPUT_DIR="${CURRENT_DIR}/data/${TITLE}"
	    rm -rf $OUTPUT_DIR
	    mkdir $OUTPUT_DIR
	    # Create sims
	    #-----------------------------------
	    ######################################################
	    # Path to simulation runcard

	    # Create necessary condor filenames
	    SARTRE_BASE="${OUTPUT_DIR}/${TITLE}"
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

	    echo "Queue" >> $SARTRE_JOB

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
	    declare -i DEUT
	    if [ "${nucleus}" == "D" ]; then
		DEUT=1
	    else
		DEUT=0
	    fi

	    echo "root -l ${photon}_solid_study_psi2S.C\\(${EBEAM},${DEUT},${NEVENTS},\\\"${OUTPUT_DIR}\\\"\\)" >> $SARTRE_SH
	    cd $OUTPUT_DIR
	    chmod +x ./*
	    condor_submit ${SARTRE_JOB}
	    done
    done
done
