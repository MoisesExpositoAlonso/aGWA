#!/bin/bash

#  Reserve 8 CPUs for this job
#$ -pe parallel 10
#
# Request  RAM
#$ -l h_vmem=10G
#
#Request not to send it to the biggest node
#$ -l mem_total_lt=600G
#
#  Request it to run this long HH:MM:SS
#$ -l h_rt=24:00:00
#
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#
#  Run job from current working directory
#$ -cwd
##  Send email when the job begins, ends, aborts, or is suspended
#$ -m beas


time bash sh-ANALYSES_ANCESTRYGWA.sh label_input_genomes.tsv discrete trialdiscrete10 > discrete10.log

