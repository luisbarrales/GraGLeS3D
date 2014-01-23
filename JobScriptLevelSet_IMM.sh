#!/usr/bin/env zsh

### Job name
#BSUB -J LEVELSET

### File / path where STDOUT will be written , the %J is the job id
#BSUB -o LEVELSET.%J

### (OFF ) Different file for STDERR , if not to be merged with STDOUT
#BSUB -e LEVELSET.e%J
#BSUB -M 2000
#BSUB -W 1:00

### Specify your mail address
#BSUB -u miessen@imm.rwth-aachen.de

### Send a mail when job is done
#BSUB -N

#### compute nodes
#BSUB -n 1


### Execute your application
LevelSet_IMM
