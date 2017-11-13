#!/bin/sh

NCORES=4
DEFFNM=md

TPR=${DEFFNM}.tpr
PDB=${DEFFNM}.pdb

MDRUN_OPTS=""

gmx mdrun -nt $NCORES -v -deffnm ${DEFFNM} -c ${PDB} -cpi -append $MDRUN_OPTS
