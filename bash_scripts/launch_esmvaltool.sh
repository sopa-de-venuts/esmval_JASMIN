#!/bin/bash 
#BSUB -n 1
#BSUB -q high-mem
#BSUB -W 24:00
#BSUB –R "select[maxmem > 128000] rusage[mem=150000]"
#BSUB -o /home/users/jvegas/esmval_ohc.%J.out
#BSUB -e /home/users/jvegas/esmval_ohc.%J.err 


source ~jvegas/miniconda3/bin/activate esmvaltool
esmvaltool -c $HOME/config-user.yml recipe_transports.yml --max-years 5 --max-datasets 5 --skip-nonexistent
