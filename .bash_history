exit
pwd
pdu
groups
cd ..
ls
cd ..
ls
cd ..
ls
cd open
cd var
ls
cd ..
ls media/
ls -lt
ls -lt datacentre/
ls -lt datacentre/archvol5
module available
ls
cd badc/
ls
ls cordex
cd cordex/
ls
cd data/
ls
ls -lt
cd CORDEX/
ls
cd output/
ls
cd EUR-11
ls
cd MOHC/
ls
cd MOHC-HadGEM2-ES/
ls
ls historical
cd rcp85/r1i1p1/MOHC-HadREM3-GA7-05/v1/mon/pr/latest/pr_EUR-11_MOHC-HadGEM2-ES_rcp85_r1i1p1_MOHC-HadREM3-GA7-05_v1_mon_20
cd 
vim datasets
cd /badc
ls
cd cru/
ls
cd data/
lsls
ls
vim 00README 
ls cru_cy
exit
pwd
cd ..
ls
cd ..
ls
cd ..
ls
cd group_workspaces/
ls
cd ..
cd gws
ls
cd ..
ssh -A pcos@jasmin-sci3.ceda.ac.uk
ls
vim datasets 
cd /badc
cd /
ls
ls work/
cd work/
ls
cd ..
ls media/
ls tmp/
cd ..
ls
ls gws
ls home/
ls var
exit
ssh -A pcos@jasmin-sci2.ceda.ac.uk
ls
cd /
ls
ls -lt
ls badc/
cd badc/cmip6/
ls
ls metadata/
cd data/
ls
cd CMIP6/
ls
cd HighResMIP/
ls
cd ECMWF/
s
ls
ECMWF-IFS-HR
cd ECMWF-IFS-HR
ls
ls hist-1950/
ls spinup-1950/
ls spinup-1950/r1i1p1f1/
ls spinup-1950/r1i1p1f1/Omon/
ls hist-1950/r1i1p1f1/+
ls hist-1950/r1i1p1f1/
ls hist-1950/r1i1p1f1/Amon/
cd /home/users/valeriu/roses
ls
cd u-bd684/
ls
cd
esmvaltool --help
xit
exit
exi
exit
cd /badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/
ls
cd /badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/latest/tas_Amon_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_20
cd /badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/latest/
ls
01cd
cd
pwd
ls exs/
ls bash_scripts/
$USER
ls exs
rm exs/i*
ls exs
vim exs/2015.1170414.err
rm exs/*
ls exs
ls exs/
vim exs/2015.1170459.err
cd /
ls
exit
cd /group_workspaces/
ls
ls jasmin
ls jasmin2
ls jasmin2/workshop/
ls jasmin2/workshop/users/
ls jasmin3
ls jasmin4
ls cems
ls cems2
cd
ls
exit
top
vim .bashrc 
source .bashrc 
xfer
exit
cp ../jvegas/config-developer.yml .
ls
cp ../jvegas/config-user.yml .
ls
ls -lt
vim config-developer.yml 
cp ../jvegas/launch_esmvaltool.sh .
ls
cd ..
ls
cd ..
ls
cd ..
ls
ls group_workspaces/
cd
ls
cp ../jvegas/miniconda3/ .
cp -r ../jvegas/miniconda3 .
ls
ls miniconda3/
rm -r miniconda3/
ls
mv config-developer-jv.yml 
mv config-developer.yml config-developer-jv.yml 
mv config-user.yml config-user-jv.yml 
ls
mv config-* config-files/
ls
ls config-files/
ls
mkdir bash_scripts
mv launch_esmvaltool.sh bash_scripts/
ls
ssh jasmin-sci1
vim .bashrc 
source .bashrc 
sci3
vim .bashrc 
source .bashrc 
xfer
sci1
sci2
ls 
ls ex04/
ncdump -h ou
module add NCO
sci1
cd /
ls
cd gws
ls
cd
sha256sum 957d2f0f0701c3d1335e3b39f235d197837ad69a944fa6f5d8ad2c686b69df3b
mkdir miniconda3
cd miniconda3/
sha256sum Miniconda3-latest-Linux-x86_64.sh
sha256sum 957d2f0f0701c3d1335e3b39f235d197837ad69a944fa6f5d8ad2c686b69df3b
sha256sum Miniconda3-latest-Linux-x86.sh
idl
ls
cd
vim .bashrc 
cd miniconda3/
ls
sh Miniconda3-latest-Linux-x86_64.sh 
-u
ls
cd ..
mv miniconda3/ miniconda/
cd miniconda/
sh Miniconda3-latest-Linux-x86_64.sh 
python --version
conda init
conda activate
exit
conda activate /gws/smf/j04/primavera/envs/esmvaltool
conda info --envs
cd
.vim bashrc
vim .bashrc
source .bashrc
esmval_env 
esmvaltool -v
conda activate /gws/smf/j04/primavera/envs/esmval_dev
esmvaltool -v
conda activate /gws/smf/j04/primavera/envs/esmval_primavera
esmvaltool -v
pwd
ls
mv config-files/config-user.yml config-files/config-user-bsc.yml
ls config-files/
ls
ls config-files/
vim config-files/config-user-bsc.yml 
ls
ls miniconda
rm -r miniconda/
ls
ls recipes/
ls bash_scripts/
ls -lt
lt scripts/
l scripts/
ls -lt scripts/
ls -ltbash
ls -lt bash_scripts/
chmod +x *.sh
chmod +x bash_scripts/*.sh
ls bash_scripts/
vim laun
vim bash_scripts/launch_esmvaltool.sh 
vim recipes/recipe_martin_fig.yml 
ls recipes/
rm -r recipes/
pwd
ls
rm *.yml
ls
mv datasets datasets.txt
ls recipes/
cd recipes/
cp recipe_martin_fig.yml recipe_JASMIN.yml
vim recipe_JASMIN.yml 
cd ../config-files/
vim config-user.yml 
cd ..
vim bash_scripts/exe
ls bash_scripts/
vim bash_scripts/esmvaltool_job.sh 
vim recipes/recipe_JASMIN.yml 
vim config-files/config-developer.yml 
vim recipes/recipe_JASMIN.yml 
cd bash_scripts/
ls
cd 
cd bash_scripts/esmvaltool_job.sh 
sh bash_scripts/esmvaltool_job.sh 
vim recipes/recipe_JASMIN.yml 
sh bash_scripts/esmvaltool_job.sh 
exit
python --version
conda activate /gws/smf/j04/primavera/envs/esmvaltool
sc1
sci1
exit
cd /
ls
cd badc/cmip6/
ls
cd data/
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
