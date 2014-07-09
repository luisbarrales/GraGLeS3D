#!/usr/bin/env zsh

### Job name
#BSUB -J LEVELSET

### File / path where STDOUT will be written , the %J is the job id
#BSUB -o LEVELSET.%J

### (OFF ) Different file for STDERR , if not to be merged with STDOUT
#BSUB -e LEVELSET.e%J



###  Request vitual memory  you  need  for your job  in  MB / 50000 = 9GB
### Memory limit for ALL nodes
#BSUB -M  42000

### Request the time you need for execution in minutes
### The format for the parameter is: [ hour :] minute ,
### that means for 80 minutes you could also use this : 1:20
#BSUB -W 4:15

### Specify your mail address
#BSUB -u miessen@imm.rwth-aachen.de

### Send a mail when job is done
#BSUB -N

#### compute nodes
#BSUB -n 64
### for shared memory jobs (OpenMP)
#BSUB -a "bcs openmp"

### Change to the work directory
mkdir mySim64_test_50000 
cd mySim642_test_50000


cp ../LevelSet_IMM ./
cp ../parameters.xml ./

### Execute your application
time LevelSet_IMM parameters.xml
