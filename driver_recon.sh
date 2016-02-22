#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -e output_erromp8
#PBS -o outputomp8

PURIFYPATH="/home/lwolz/myproj/imaging/purify/"

SRCPATH="build/src/"

PRGM="testrun_m31"

cd $PURIFYPATH


echo "NNODES   is "`uniq $PBS_NODEFILE | wc -l`         
echo "NPROCS   is "`cat  $PBS_NODEFILE | wc -l`         
echo "Using the following nodes:"      
cat $PBS_NODEFILE | uniq   

OUTPUTDIR=$PURIFYPATH"/data/test/chirp0_ppn8/"

mkdir $OUTPUTDIR

./$SRCPATH$PRGM $OUTPUTDIR > outputomp8.txt


