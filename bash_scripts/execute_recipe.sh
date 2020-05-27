#!/bin/bash
recipe=recipe_JASMIN
configuser=config-user
echo "Starting $recipe.yml recipe"
echo "With config-user --> $configuser.yml"
esmvaltool -c /home/users/pcos/config-files/${configuser}.yml /home/users/pcos/recipes/${recipe}.yml --max-years 2 #--max-datasets 5 #--skip-nonexistent
ln -sfT $(ls -1d /work/scratch-nompiio/${USER}/${recipe}* | tail -1) /home/users/pcos/latest_recipie
