ls
cd CMIP6/
cd
top
exit
vim config-files/config-developer.yml 
cd /
ls
sci1
exit
ls
rsync -zaP jcos@bscearth001:/esarchive/scratch/jcos/esmvaltool/scripts .
pwd
rsync -zaP jcos@bscearth001:/esarchive/scratch/jcos/esmvaltool/scripts /home/users/pcos/
rsync -zaP jcos@bscearth001.int.bsc.es:/esarchive/scratch/jcos/esmvaltool/scripts /home/users/pcos/
ls
ls scripts/
ls
module list
vim .bashrc 
source .bashrc
vim .bashrc 
source .bashrc
sci
vim .bashrc 
source .bashrc
sci2
ls
cd ..
ls
cd jvegas/
ls
vim launch_esmvaltool.sh 
vim config-developer.yml 
ls miniconda3/
ls miniconda3/pkgs/
ls
vim launch_ohc_diag.sh 
vim launch_test.sh 
vim launch_blocking
vim launch_blocking.sh 
vim recipe_wiki.yml 
vim launch_esmvaltool.sh 
cd miniconda3/bin/
ls
cd activate 
ls
cd..
cd ..
ls
cd bin/
vim activate 
cd ..
cd 
cd /
ls
cd group_workspaces/
ls
cd ..
ls projects/
ls utils/
ls utils/Modules/3.2.3/
ls utils/Modules/3.2.3/modulefiles/
cd
cd ../jvegas
ls
cd
bsub -man
bsub --man
ls -lt ../jvegas/
top
module available
module available | grep esmv
module available | grep esmv*
module list
module add esmvaltool
module add jaspy
module list
module rm jaspy
ls
python
pyhton -version
pyhton --version
python --version
python3 --version
python2.7 --version
python3.8 -v
python2.8 -v
python2.7 -v
mkdir ex04
ls /badc
ls /badc/cmpi6
ls /badc/cmip6
ls /badc/cmip6/data/
ls /badc/cmip6/data/CMIP6/
ls -lt /badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/
ls -lt /badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/latest/tas_Amon_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_20
ls -lt /badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/latest/tas_Amon_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_20*
ls -lt /badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/latest/
ls -lt /badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/latest/tas_Amon_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_205001-205012.nc 
INPUT_FILE=/badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/latest/tas_Amon_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_205001-205012.nc
INPUT_FILE
pwd
OUTPUT_FILE=/home/users/pcos/ex04/out.nc
module load jaspy
$INPUT_FILE 
cdo $INPUT_FILE $OUTPUT_FILE 
cdo selname $INPUT_FILE $OUTPUT_FILE 
ncdump -h ex04/out.nc 
mv ex04/ exs/
ls
vim bash_scripts/extract_HRMIP_proj.sh
sh bash_scripts/extract_HRMIP_proj.sh 
vim bash_scripts/extract_HRMIP_proj.sh
vim bash_scripts/extract_HRMIP_proj.sh 2015
sh bash_scripts/extract_HRMIP_proj.sh 2015
vim bash_scripts/extract_HRMIP_proj.sh
sh bash_scripts/extract_HRMIP_proj.sh 2015
vim bash_scripts/extract_HRMIP_proj.sh
sh bash_scripts/extract_HRMIP_proj.sh 2015
vim bash_scripts/extract_HRMIP_proj.sh
sh bash_scripts/extract_HRMIP_proj.sh 2015
ls exs/
vim bash_scripts/extract_HRMIP_proj.sh
rm exs/out.nc 
rm exs/tas-tas_Amon_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_201501-201512.nc 
vim bash_scripts/submit-all.sh
echo $BASH_VERSION
vim bash_scripts/submit-all.sh
sh bash_scripts/submit-all.sh 
vim bash_scripts/submit-all.sh
sh bash_scripts/submit-all.sh 
vim bash_scripts/submit-all.sh
sh bash_scripts/submit-all.sh 
sh bash_scripts/extract_HRMIP_proj.sh 2015
ls exs/
vim bash_scripts/submit-all.sh
vim bash_scripts/demo.sh
bsub < bash_scripts/demo.sh 
bkill < bash_scripts/demo.sh 
bkill
bkill <1171396>
bjobs
bkill 1171396
bjobs
ls
vim 1171396.err 
vim 1171396.out
rm 1171396.*
cd /
ls
cd gws
ls
ls nopw/
cd
vim bash_scripts/esmvaltool_job.sh
module load idl
idl
cd /gws/smf/j04/primavera/
ls
ls envs/
cd ..
ls
cd ..
ls
cd ..
ls
cd ..
ls
cd gws/
cd j04
ls
cd smf/j04/primavera/conf-files/
ls
vim config-user.yml 
cp config-user.yml /home/users/pcos/config-files/
vim config-user.yml 
cd /work/scratch-nompiio/
ls
cd sloosvel/
ls
cd ..
cd jvegas/
ls
cd ..
mkdir pcos
ls
mv pcos jcos
ls
mv jcos pcos
ls
:q
cd
vim config-files/config-developer-jv.yml 
cd scripts/
ls
pwd
cd ../config-files/
ls
vim config-user-jv.yml 
rm config-user-jv.yml 
ls
vim config-developer.yml 
rm config-developer.yml 
mv config-developer-jv.yml config-developer.yml 
pwd
ls
vim config-user.yml 
vim /work/scratch-nompiio/pcos/recipe_JASMIN_20200521_114428/run/mean/mean/log.txt
vim .bashrc
exit
pwd
vim config-files/config-developer.yml 
vim recipes/recipe_JASMIN.yml 
sci2
top
vim config-files/config-user
vim config-files/config-user.yml 
vim config-files/config-developer.yml 
cd /gws/nopw/j04/
ls
ls rsg_reproc_lv1 
ls docile/
ls launchpad/
ls ukca_vol1/
cd 
cd /group_workspaces/jasmin4
ls
cd esmeval/
cd
vim recipes/recipe_JASMIN.yml 
vim config-files/config-developer.yml 
cd /badc/cmip6/data/CMIP6/CMIP/
find ./ FGOALS-*
ls CAS/CAS-ESM2-0/historical/r1i1p1f1/Amon/tas/gn/cd c
cd
vim config-files/config-developer.yml 
ls /badc/
vim config-files/config-developer.yml 
ls /badc/cmip6/data/CMIP6/CMIP/
vim config-files/config-user.yml 
vim config-files/config-developer.yml 
cd /badc/cmip6/data/CMIP6/CMIP/AWI/AWI-CM-1-1-MR/historical/r1i1p1f1/Amon/tas/gn/latest/
ls
cd
pwd
vim config-files/config-user.yml 
vim bash_scripts/execute_recipe.sh 
cd /gws/smf/j04/
ls
cd /gws/smf/j04/cd
pwd
cd
pwd
cd /badc/cmip6/data/CMIP6/CMIP/AWI/AWI-CM-1-1-MR/historical/r1i1p1f1/Amon/tas/gn/
cd latest/
ls
cd /badc/cmip6/data/CMIP6/CMIP
ls
find . -name GFDL-CM4
find . -name CIESM
find . -name MCM-UA
find . -name MCM-UA-1-0
vim config-files/config-developer.yml 
sci2
cd /badc/cmip6/data/CMIP6/
ls
ls CMIP
cd CMIP/IPSL
ls
cd IPSL-CM6A-LR/
ls
ls amip/
cd historical/
ls
cd ..
ls
cd ..
ls
find ./ ssp585
ls
ls CNRM-CERFACS/CNRM-CM6-1/ cd
cd
vim config-files/config-user
vim config-files/config-user.yml 
vim recipes/recipe_JASMIN.yml 
vim config-files/config-user.yml 
vim recipes/recipe_JASMIN.yml 
vim recipes/
vim bash_scripts/execute_recipe.sh 
vim recipes/recipe_JASMIN.yml 
vim config-files/config-user.yml 
vim recipes/recipe_JASMIN.yml 
vim config-files/config-user.yml 
vim config-files/config-developer.yml 
vim bash_scripts/execute_recipe.sh 
vim recipes/recipe_JASMIN.yml 
vim config-files/config-developer.yml 
vim config-files/config-user.yml 
vim config-files/config-developer.yml 
vim config-files/config-user.yml 
vim config-files/config-developer.yml 
vim config-files/config-user.yml 
vim recipes/recipe_JASMIN.yml 
vim config-files/config-user.yml 
vim config-files/config-developer.yml 
vim config-files/config-user.yml 
vim recipes/recipe_JASMIN.yml 
vim config-files/config-user.yml 
vim config-files/config-developer.yml 
vim config-files/config-user.yml 
vim config-files/config-developer.yml 
vim config-files/config-user.yml 
vim recipes/recipe_JASMIN.yml 
vim config-files/config-user.yml 
vim recipes/recipe_JASMIN.yml 
sci2
ls
mkdir out
ls
mv out figures
ls
cd /badc/
ls
cd cmip6
ls
cd data
ls
cd CMIP6/CMIP
ls
cd /
find type -d name AWI-CM1-1-MR
find AWI-CM1-1-MR
find ./ AWI-CM1-1-MR
cd gws/nopw/j04/cmip6_prep_vol1/cmip6_data_vols/
ls
cd cmip6-data-vols/
ls
cd /
ls
cd
esmval_env 
sh bash_scripts/exe
ls bash_scripts/h
ls bash_scripts/
cat bash_scripts/esmvaltool_job.sh 
mv bash_scripts/esmvaltool_job.sh bash_scripts/execute_recipe.sh
sh bash_scripts/execute_recipe.sh 
cd /badc/cmip6/data/CMIP6/
ls
CMIP/
ls
cd CMIP/
ls
cd CNRM-CERFACS/
ls
cd ..
cd NCAR/
ls
ls CESM2
cd CESM2
ls 1pctCO2/r1i1p1f1/Amon/tas
ls 1pctCO2/r1i1p1f1/Amon/tas/gn/
ls esm-hist/r1i1p1f1/Amon/tas/gn/
ls esm-hist/r1i1p1f1/Amon/tas/gn/latest/
cd ..
ls
cd ..
ls
cd MRI/MRI-ESM2-0/
ls
cd 
cd /badc/cmip6/data/CMIP6/HighResMIP/
cd NOAA-GFDL/GFDL-CM4C192/
ls
cd highresSST-future/r1i1p1f1/Amon/
ls
cd /gws/nopw/j04/
ls
cd /badc/cmip6/
ls
cd data/CMIP6/
ls
cd CMIP
ls
cd BCC
ls
cd BCC-CSM2-MR/
ls
cd ..
ls
ls AS-RCEC/TaiESM1/
ls CCCR-IITM/
ls CCCR-IITM/IITM-ESM/
ls CAS/FGOALS-g3/
ls
ls CCCma/CanESM5
ls E3SM-Project/E3SM-1-1-ECA/
ls EC-Earth-Consortium/EC-Earth3/
ls EC-Earth-Consortium/EC-Earth3-Veg/
ls EC-Earth-Consortium/EC-Earth3-Veg/historical/
ls EC-Earth-Consortium/EC-Earth3-Veg/historical/r1º
ls EC-Earth-Consortium/EC-Earth3-Veg/historical/r1i1p1f1/Amon
ls
ls FIO-QLNM/FIO-ESM-2-0/
ls IPSL/IPSL-CM6A-LR/
ls
ls MRI/MRI-ESM2-0/
ls NASA-GISS/GISS-E2-1-H/
ls NASA-GISS/GISS-E2-2-G/
ls 
ls NUIST/NESM3/
ls SNU/SAM0-UNICON/historical/ls
ls THU/CIESM/historical/
ls UA/MCM-UA-1-0/historical/
cd
sh bash_scripts/execute_recipe.sh 
ls SNU/SAM0-UNICON/historical/ls
sh bash_scripts/execute_recipe.sh 
cd ../jvegas/
ls
vim config-user.yml
vim config-developer.yml 
vim config-esmeval.yml 
cd /badc/cmip6/data/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Amon/tas/gn/latest/tas_Amon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc cd
cd
sh bash_scripts/execute_recipe.sh 
cd /gws/smf/j04/primavera/envs/esmvaltool/
ls
sh bash_scripts/execute_recipe.sh 
cd
sh bash_scripts/execute_recipe.sh 
cd /gws/smf/j04/
ls
conda deactivate
esmval_env 
cd /gws/smf/j04/
ls
cd primavera/
ls
cd envs/esmval_primavera/
ls
ls bin/
find ./config-developer.yml
cd
sh bash_scripts/execute_recipe.sh 
cd /gws/smf/j04/primavera/conf-files/config-user.yml 
vim /gws/smf/j04/primavera/conf-files/config-user.yml 
sh bash_scripts/execute_recipe.sh 
vim bash_scripts/execute_recipe.sh 
conda deactivate
conda activate /gws/smf/j04/primavera/envs/esmvaltool/
sh bash_scripts/execute_recipe.sh 
conda deactivate
esmval_env 
sh bash_scripts/execute_recipe.sh 
cd /gws/smf/j04/primavera/envs/esmval_primavera/
ls bin/
:wq
find ./ -name config-developper.yml
find . -name config-developper.yml
find . -name config-developer.yml
vim lib/python3.7/site-packages/esmvalcore/config-developer.yml 
cd
sh bash_scripts/execute_recipe.sh 
sci2
cd /badc/
ls
cd ecmwf-era-interim/
ls
cd data/
ls
cd monthly-means/
ls
ls gg
ls ha
ls hg
cd /
ls
ls bodc/
ls edc
ls neodc/
ls ngdc/
ls sparc/
cd badc/cru/
ls
cd data/
ls cru_cy/
ls cru_jra/cru_jra_2.0/
ls crutem/
ls PDSI/cd
dc
ls
vim datasets.txt 
exit
cd /badc/
cd cmip6
ls
cd data/
ls
CMIP6
cd CMIP6/HighResMIP/
ls
cd CNRM-CERFACS/
ls
cd CNRM-CM6-1/
ls
ls hist-1950/
exit
vim recipes/recipe_JASMIN.yml 
ls bash_scripts/
sh bash_scripts/execute_recipe.sh 
esmval_env 
sh bash_scripts/execute_recipe.sh 
vim ../jvegas/config-developer.yml 
vim config-files/config-user.yml 
cd /group_workspaces/jasmin4/esmeval/obsdata-v2
ls
cd..
cd ..
ls
cd ..
ls
cd ..
ls
ls jasmin
ls jasmin2
ls jasmin2/leicester/
ls jasmin3 
ls jasmin3/esgf_repl/
ls jasmin4
ls jasmin4/esmeval/
ls jasmin4/esmeval/esmvaltool-tutorial-docs/
ls jasmin4/esmeval/esmvaltool-tutorial-docs/ukesm-April-2020/
ls jasmin4/esmeval/esmvaltool-tutorial-docs/ukesm-April-2020/examples/
cd jasmin4/esmeval/
ls
vim recipe_hands_on_example.yml 
vim environment.sh 
vim config-user.yml 
ls obsdata-v2/
ls obsdata-v2/Tier1
ls obsdata-v2/Tier2
ls obsdata-v2/Tier3
ls obsdata-v2/Tier3/ERA-Interim/
ls obsdata-v2/Tier3/ERA-Interim/*tas*
ls obsdata-v2/Tier3/ERA-Interim/*pr*
wq
ls obsdata-v2/Tier2/CRU/OBS_CRU_reanaly_TS4.02_Amon_tas*
ls obsdata-v2/Tier2/CRU/OBS_CRU_reanaly_TS4.02_Amon_pr*
ls obsdata-v2/Tier2
ls obsdata-v2/Tier3
ls obsdata-v2/Tier1
cd /badc/
ls
cd cmip6
ls
cd data/
ls
cd CMIP6/CMIP/
ls
cd AWI/AWI-ESM-1-1-LR/
ls
cd ..
find ./ ssp585
ls
ls CCCma/
ls CCCma/CanESM5
ls CCCma/CanESM5/esm-piControl/
ls CCCma/CanESM5/amip/
cd ..
ls
cd ScenarioMIP/
ls
ls AWI/
ls AWI/AWI-CM-1-1-MR/
cd ..
ls
cd HighResMIP/
ls
find ./ ssp585
find ./ -name ssp585
ls
ls ECMWF/ECMWF-IFS-HR/highresSST-present/r1i1p1f1/Amon/*tas*
ls 
ls CAS/FGOALS-f3-H/highresSST-present/r1i1p1f1/
ls CMCC/CMCC-CM2-HR4/highresSST-future/r1i1p1f1/Amon/tas/gn/latest/tas_Amon_CMCC-CM2-HR4_highresSST-future_r1i1p1f1_gn_20
ls CMCC/CMCC-CM2-HR4/highresSST-future/r1i1p1f1/Amon/tas/gn/latest/tas_Amon_CMCC-CM2-HR4_highresSST-future_r1i1p1f1_gn_20*
ls
ls CMCC/CMCC-CM2-HR4/highres-future/r1i1p1f1/Amon/tas/gn/latest/tas_Amon_CMCC-CM2-HR4_highres-future_r1i1p1f1_gn_20
ls 
ls CMCC/CMCC-CM2-
ls EC-Earth-Consortium/EC-Earth3P/
ls EC-Earth-Consortium/EC-Earth3P-HR/
ls -lt EC-Earth-Consortium/EC-Earth3P/
ls 
ls MOHC/HadGEM3-GC31-LL
ls MOHC/HadGEM3-GC31-LM
ls MOHC/HadGEM3-GC31-MH
ls MOHC/HadGEM3-GC31-MM
ls
ls MRI/MRI-AGCM3-2-S
find ./ -name gn
find . -name gn
find . -name gr
find ./ -type f -iname gn
find . -type f -iname gn
find . -type f -name gn
find ./ -type f -name gn
find ./ -type d -name gn
find ./ -type d -name gr
find . -type d -name gr
find . -type d -name gn
ls 
ls CNRM-CERFACS/CNRM-CM6-1-HR/hist-1950/r1i1p1f2/Amon/tas/gr/
ls
ls MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/Amon/tas/gn/
ls EC-Earth-Consortium/EC-Earth3P/hist-1950/r1i1p2f1/Amon/tas/gr/
ls
ls MRI/MRI-AGCM3-2-H/highresSST-present/r1i1p1f1/Amon/tas/gn/
ls NOAA-GFDL/GFDL-CM4C192/hist-1950/r1i1p1f1/Amon/ts/gr3/cd
cd
sh bash_scripts/execute_recipe.sh 
vim /work/scratch-nompiio/pcos/recipe_JASMIN_20200525_094414/run/mean/mean/log.txt
sh bash_scripts/execute_recipe.sh 
cd /group_workspaces/jasmin4/
ls
exit
cd /group_workspaces/jasmin4/esmeval/obsdata-v2
ls
ls Tier1
ls Tier2
cd Tier2
ls
cd CRU/OBS_CRU_reanaly_TS4.02_Amon_:[B
cd
sh bash_scripts/execute_recipe.sh 
esmval_env 
sh bash_scripts/execute_recipe.sh 
cd /group_workspaces/jasmin4/esmeval/obsdata-v2
ls
ls Tier1
ls Tier2
ls Tier3
ls Tier3/ERA-Interim/OBS6_ERA-Interim_reanaly_1_tas*
ls Tier3/ERA-Interim/OBS6_ERA-Interim_reanaly_1_Amon_*
cd
sh bash_scripts/execute_recipe.sh 
exit
vim .bashrc 
vim recipes/recipe_JASMIN.yml 
sci2
exit
ls
git config --global use.name "jcos"
git config --global user.email "josep.cos@bsc.es"
git config --global --list
git remote -v
git clone git@earth.bsc.es:jcos/esmval_jasmin_pcos.git
git clone https://earth.bsc.es/gitlab/jcos/esmval_jasmin_pcos.git
git status
ls
rm -r esmval_jasmin_pcos/
ls
rm -r esmval_jasmin_pcos/
ls
ls exs/
rm -r exs
git init
git add .
git commit -m "initial commit"
git remote add origin git@earth.bsc.es:jcos/esmval_jasmin_pcos.git
git push origin master
git remote add origin https://earth.bsc.es/gitlab/jcos/esmval_jasmin_pcos.git
git push -u origin master
vim bash_scripts/extract_HRMIP_proj.sh 
vim bash_scripts/submit-all.sh 
vim bash_scripts/demo.sh 
rm bash_scripts/demo.sh 
ls
vim recipes/recipe_JASMIN.yml 
ls
git init
ls -la
rm -r .git
rm -rf .git
git init
git add .
git commit -m "initial commit"
git remote add origin https://github.com/sopa-de-venuts/esmval_JASMIN.git
git push -u origin master
git pull https://github.com/sopa-de-venuts/esmval_JASMIN.git
ls
ls -lt
git status
git add .
git commit -m "initial commit"
git push origin master
git add .gitignore
vim .gitignore
git add .
git commit -m ".gitignore"
git push origin master
vim recipes/recipe_JASMIN.yml 
vim bash_scripts/execute_recipe.sh 
vim recipes/recipe_JASMIN.yml 
vim scripts/blank
vim scripts/blank.py 
vim recipes/recipe_JASMIN.yml 
vim config-files/config-user.yml 
vim recipes/recipe_JASMIN.yml 
vim scripts/martin_plots.py 
vim recipes/recipe_JASMIN.yml 
vim bash_scripts/execute_recipe.sh 
vim recipes/recipe_JASMIN.yml 
vim scripts/martin_plots.py 
vim config-files/config-user.yml 
ls -lt scripts/
vim recipes/recipe_JASMIN.yml 
vim bash_scripts/execute_recipe.sh 
ls
ls latest_recipie 
cd /work/scratch-nom
vim config-files/config-user.yml 
ls
ls figures/
git add .
git commit -m "figures test"
git push origin master
exit
cd work
cd /work/scratch-nompiio/
ls
cd pcos/
ls
vim recipe_JASMIN_20200525_100733/run/mean/mean/log.txt recipe_JASMIN_2020052
top
cd /work/scratch-nompiio/
ls
cd pcos/
ls -lt
cd recipe_JASMIN_20200525_110456
ls
cd plots
cd mean/mean/
ls
cd ..
mv mean/ ~/figures/
exit
vim config-files/config-user.yml 
sci2
cd /group_workspaces/
ls
cd
cd /
ls
cd work/
ls
sci2
exit
cd /group_workspaces/
ls
cd jasmin4
cd esmeval/obsdata-v2/
ls
ls Tier1
ls Tier2
ls Tier3
cd Tier1
ls CALIOP/
ls
ls AIRS
ls ATSR
ls CERES-EBAF/
ls CloudSat/
ls 
ls ESACCI-GHG/
ls GPCP-SG/
ls ISCCP/
ls
ls JRA-55/
ls MODIS/
ls SSMI
ls SSMI-MERIS/
ls TRMM-L3/pr
ls TRMM-L3/
ls ghgcci/
cd ..
cd Tier2
ls
ls CALIPSO-GOCCP/OBS_CALIPSO-GOCCP_sat_3.1.2_cfMon_clcalipso_200701-201512.nc 
ls CERES-EBAF/
ls CRU/
ls Duveiller2018/
ls ESACCI-SOILMOISTURE/
ls Eppley-VGPM-MODIS/
ls
ls GCP/OBS_GCP_reanaly_2018_
ls GHCN/
las HadCRUT3
la HadCRUT3
ls HadCRUT3
ls HadCRUT4
ls ../Tier1/GPCP-SG/
ls HadCRUT4-
ls HadCRUT4-clim/
ls HadISST/
ls ISCCP-FH/
ls Landschuetzer2016/
ls
ls NCEP/
ls
ls OSI-450-
ls OSI-450-nh/
ls
ls OSI-450-sh/
ls
ls PATMOS-x/
ls PHC/
ls PIOMAS/
ls
ls WOA/
cd ../Tier3
ls
ls ERA-Interim/
ls ../Tier2/CRU/
ls ../Tier2/
cd /badc/cmip6/data/
ls
cd CMIP6/CMIP/
find . GISS*
ls
cd /group_workspaces/
cd /group_workspaces/jasmin4/esmeval/
ls
cd obsdata-v2/
ls
ls Tier3/ERA-Interim/OBS6_ERA-Interim_reanaly_1_Amon_pr_
cd
pwd
ls
mkdir job_output
cd job_output/
pwd
exit
vim recipes/recipe_martin_fig.yml 
ls recipes/
vim recipes/recipe_JASMIN.yml 
vim scripts/martin_plots.py 
sci2
exit
cd /badc/cmip6/data/CMIP6/HighResMIP/
ls
ls CMCC/
ls CNRM-CERFACS/
ls EC-Earth-Consortium/
ls ECMWF/
ls INM/
ls MOHC/
ls MPI-M/
ls NERC/
ls CMCC/CMCC-CM2-HR4/
ls -lt  CMCC/CMCC-CM2-HR4/
ls -lt  CMCC/CMCC-CM2-VHR4/
ls IPSL/
ls ECMWF/ECMWF-IFS-HR/
exi
exit
sci2
exit
cd /work/scratch-nompiio/pcos/
ls -lt
vim recipe_JASMIN_20200525_100733/
vim recipe_JASMIN_20200525_100733/run/mean/mean/log.txt 
ls
cd
ls
ls scripts/
rm scripts/martin_plots.py 
mv martin_plots.py scripts/
ls
ls scripts/
ls -lt
ls -lt scripts/
esmval_env 
vim bash_scripts/execute_recipe.sh 
sh bash_scripts/execute_recipe.sh 
vim /work/scratch-nompiio/pcos/recipe_JASMIN_20200526_074431/run/mean/mean/log.txt
ls -lt
ls -lt scripts/
vim /work/scratch-nompiio/pcos/recipe_JASMIN_20200526_074431/run/mean/mean/log.txt
cd /work/scratch-nompiio/pcos/recipe_JASMIN_20200526_074431/run/mean/mean; MPLBACKEND="Agg" /gws/smf/j04/primavera/envs/esmval_primavera/bin/python /home/users/pcos/scripts/martin_plots.py /work/scratch-nompiio/pcos/recipe_JASMIN_20200526_074431/run/mean/mean/settings.yml -f
ls
cd
ls /work/scratch-nompiio/pcos/recipe_JASMIN_20200526_074431/
ls /work/scratch-nompiio/pcos/recipe_JASMIN_20200526_074431/plots/mean/mean/ 
cp /work/scratch-nompiio/pcos/recipe_JASMIN_20200526_074431/plots/mean/mean/ figures/
cp -r /work/scratch-nompiio/pcos/recipe_JASMIN_20200526_074431/plots/mean/mean/ figures/
ls figures/
ls -lt figures/
ls -lt figures/mean/
git .add
git add .
git commit -m "martin_figures.py rework"
git push origin master
vim recipes/recipe_JASMIN.yml 
vim bash_scripts/execute_recipe.sh 
ls
ls -lt
vim bash_scripts/execute_recipe.sh 
sh bash_scripts/execute_recipe.sh 
wq
vim bash_scripts/execute_recipe.sh 
vim recipes/recipe_JASMIN.yml 
sh bash_scripts/execute_recipe.sh 
vim recipes/recipe_JASMIN.yml 
vim scripts/martin_plots.py 
vim recipes/recipe_
vim recipes/recipe_JASMIN.yml 
cd recipes/
ls
cp recipe_JASMIN.yml CRU.yml GHCN.yml 
cp recipe_JASMIN.yml CRU.yml
cp recipe_JASMIN.yml GHCN.yml
cp recipe_JASMIN.yml HadCRUT4.yml
cp recipe_JASMIN.yml NCEP.yml
ls 
vim CRU.yml 
vim GHCN.yml 
vim HadCRUT4.yml 
vim NCEP.yml 
cp recipe_JASMIN.yml ERAINTERIM.yml
vim ERAINTERIM.yml 
cd ..
cd bash_scripts/
cat execute_recipe.sh 
cat submit-all.sh 
vim job.sh
bjobs
vim execute_recipe.sh 
vim job.sh
ls ../recipes/
vim job.sh
bsub < job.sh 
vim job.sh
cd
vim .bashrc
vim job.sh
cd bash_scripts/
vim job.sh
bsub < job.sh 
bjobs
vim job.sh
bjobs
cd
ls
vim job_output/1587812.err 
vim bash_scripts/job.sh 
bsub < bash_scripts/job.sh 
bjobs
vim recipes/CRU.yml 
bjobs
vim recipes/recipe_JASMIN.yml 
ls /badc/cmip6/data/CMIP6/CMIP/AWI/AWI-CM-1-1-MR/piControl/r1i1p1f1/Amon/tas/gn/
bjobs
cd /gws/nopw/j04/primavera5/stream1/CMIP6/HighResMIP
ls
ls AWI/AWI-CM-1-0-HR/
ls AWI
ls CMCC/
ls CNRM-CERFACS/
ls EC-Earth-Consortium/
ls ECMWF/
ls MOHC/
ls MPI-M/
ls NERC/
ls CMCC/CMCC-CM2-HR4/
ls CMCC/CMCC-CM2-VHR4/
ls AWI/AWI-CM-1-0-LR/
ls AWI/AWI-CM-1-0-HR/hist-1950/r1i1p1f002/cd
cd
bjobs
ls
cd job_output/
ls
vim 1588487.out 
cd ..
vim recipes/CRU.yml 
vim recipes/GHCN.yml 
vim recipes/HadCRUT4.yml 
vim recipes/ERAINTERIM.yml 
vim recipes/NCEP.yml 
vim bash_scripts/job.sh 
vim job_output/1588487.err 
vim bash_scripts/job.sh 
bsub < bash_scripts/job.sh 
bjobs
exit
sci2
exit
