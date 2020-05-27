#!/bin/bash

#BSUB -q short-serial
#BSUB -n 1
#BSUB -cwd /home/users/pcos/job_output/ 
#BSUB -oo %J.out
#BSUB -eo %J.err
#BSUB -J obs_assessment_job 
#BSUB -W 04:00

source /home/users/jvegas/miniconda3/bin/activate /gws/smf/j04/primavera/envs/esmval_primavera

esmvaltool -c /home/users/pcos/config-files/config-user.yml /home/users/pcos/recipes/CRU.yml

esmvaltool -c /home/users/pcos/config-files/config-user.yml /home/users/pcos/recipes/GHCN.yml

esmvaltool -c /home/users/pcos/config-files/config-user.yml /home/users/pcos/recipes/HadCRUT4.yml

esmvaltool -c /home/users/pcos/config-files/config-user.yml /home/users/pcos/recipes/NCEP.yml

esmvaltool -c /home/users/pcos/config-files/config-user.yml /home/users/pcos/recipes/ERAINTERIM.yml
