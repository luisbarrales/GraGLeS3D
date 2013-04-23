#!/usr/bin/env zsh

### Job name
#BSUB -J Levelset

### File / path where STDOUT will be written, the %J is the job id
#BSUB -o Levelset.%J

### (OFF) Different file for STDERR, if not to be merged with STDOUT
#BSUB -e Levelset.e%J

### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 15 minutes you could also use this: 00:15
#BSUB -W 0:15					

### Request vitual memory you need for your job in MB
#BSUB -M 16000

### (OFF) Specify your mail address
#BSUB -u miessen@imm.rwth-aachen.de

### Send a mail when job is done
#BSUB -B 
#BSUB -N

### Use nodes exclusive
### BSUB -x

### Use the hpcwork with lustre
###BSUB -R "select[hpcwork]"

### Use esub for OpenMP
#BSUB -a openmp  
#BSUB -R span[hosts=1]

### Request the number of compute slots you want to use
###BSUB -n 2

###export VT_MAX_FLUSHES=0
###export VT_BUFFER_SIZE=1024


###Execute the application
Levelset_IMM
