import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import iris.coord_categorisation
from iris.time import PartialDateTime
import cartopy.crs as ccrs

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata
# from esmvaltool.diag_scripts.shared import plot
from esmvalcore.preprocessor._area import area_statistics
from esmvalcore.preprocessor._time import climate_statistics, anomalies, seasonal_statistics

class Example(object):
    def __init__(self, config):
        # ---------------------------------------------------------------------
        # config is a dictionary containing metadata regarding input files and
        # overall, as the name suggests, configuration options.
        # ---------------------------------------------------------------------
        self.cfg = config
        # self.extra_parametre = self.cfg.get('extra_parametre')

    def compute(self):
        # ---------------------------------------------------------------------
        # Every dataset in the recipe is associated with an alias. We are going
        # to use the alias and the group_metadata shared function to loop over
        # the datasets.
        #----------------------------------------------------------------------
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        total = {}
        # Loop over the datasets.
        for alias in data:
            # -----------------------------------------------------------------
            # Use the group_metadata function again so that, for each dataset,
            # the metadata dictionary is organised by the variables'
            # short name.
            # -----------------------------------------------------------------
            variables = group_metadata(data[alias], 'short_name')
            # Returns the path to the preprocessed files.
            tas_file = variables['tas'][0]['filename']
#           tos_file = variables['tos'][0]['filename']

            # -----------------------------------------------------------------
            # Now it is up to you to load the data and work with it with your
            # preferred libraries and methods. We are going to continue the
            # example using Iris and ESMValTool functions.
            # Note that all preprocessor routines can be called inside the
            # diag_script by adding the corresponding imports.
            # -----------------------------------------------------------------
            tas = iris.load(tas_file)[0]
            tas.convert_units('degC')
            y_mn = seasonal_statistics(tas, 'mean')#[0]
            # Returns the path to the preprocessed files.
            mn = climate_statistics(tas,'mean','season')[0]
            std = climate_statistics(tas,'std_dev','season')[0]
            #anom = anomalies(y_mn, 'season')
            with open("/esarchive/scratch/jcos/esmvaltool/scripts/sample.txt", "w") as f:
                f.write(str(tas))
                f.write('******************************** \\ y_mean******')
                f.write(str(y_mn[0]))
                f.write('********************************\\ mn******')
                f.write(str(y_mn[1]))
                f.write('******************************\\ std*******')
                f.write(str(std))
                f.write('********************************')
                f.write(str(0))
                f.write('********************************')
            # Returns the path to the preprocessed files.
            #self.plot(mn, 'eur_tas_mn_19502012.png')
            #self.plot(std, 'eur_tas_std_19502012.png')
            #------------------------------------------------------------------
            # Substract the data and create a new cube with the result but
            # keeping tas' metadata
            #------------------------------------------------------------------
            #cube = tas

            # Calling a preproc function inside a diag_script
            #timeseries = area_statistics(cube, 'mean')

            #timeseries.standard_name = None
            #timeseries.long_name = 'global_mean_difference_between_tas_and_tos'
            #timeseries.short_name = 'tastos'
            #total[alias] = timeseries
            # Save the outputs for each dataset.
            #self.save(timeseries, alias, data)

        # Plot the results.
        #self.plot_1D(total, data)
 
    def plot(self, cube, name):
         plt.figure(figsize=(8, 6))
         proj = ccrs.EuroPP()
         plt.axes(projection=proj)
         contour = qplt.contourf(cube)
         plt.gca().coastlines()
         plt.clabel(contour, inline=False)
         plt.savefig(os.path.join('/esarchive/scratch/jcos/esmvaltool/output/figures/', name))
         plt.clf()

    def save(self, result, alias, data):
        # More examples to access the metadata.
        script = self.cfg[n.SCRIPT]
        project = data[alias][0]['project']
        dataset = data[alias][0]['dataset']
        exp = data[alias][0]['exp']
        start = data[alias][0]['start_year']
        end = data[alias][0]['end_year']

        # Define a filename for the  NetCDF output.
        output = '{project}_{dataset}_{exp}_{script}_{start}_{end}.nc'.format(
            project=project,
            dataset=dataset,
            exp=exp,
            script=script,
            start=start,
            end=end
            )

        # Get the paths to the work directory.
        work_dir = self.cfg[n.WORK_DIR]

        #----------------------------------------------------------------------
        # Again, up to you to decide how do you want to save the data. Note
        # that this will save the NetCDF file in the work directory.
        #----------------------------------------------------------------------
        #iris.save(result, os.path.join(work_dir, output))
        iris.save(result, os.path.join(iris.save('/home/Earth/jcos/es-esmvaltool/output/', output)))

    def plot_1D(self, results, data):
        # Plot the results for each dataset in the same plot
        for result in results:
            dataset = data[result][0]['dataset']
            qplt.plot(results[result], label=dataset)
        plot_name = 'Timeseries_difference_between_tas_and_tos.{out}'.format(
            out=self.cfg[n.OUTPUT_FILE_TYPE]
        )

        # Get the path to the plot directory
        plot_dir = self.cfg[n.PLOT_DIR]

        # Some tweaks using matplotlib
        plt.legend()
        plt.tick_params(axis='x', labelsize=8)

        # Save the plot
        # plt.savefig(os.path.join(plot_dir, plot_name))
        plt.savefig(os.path.join('/home/Earth/jcos/es-esmvaltool/output/', plot_name))
        plt.close()

    def plot_2D(self, results, data):
        # Plot the results for each dataset in the same plot
        for result in results:
            dataset = data[result][0]['dataset']
            qplt.contourf(results[result], label=dataset)
        plot_name = 'tas_mn_20002010_djf.{out}'.format(
            out=self.cfg[n.OUTPUT_FILE_TYPE]
        )

        # Get the path to the plot directory
        plot_dir = self.cfg[n.PLOT_DIR]

        # Some tweaks using matplotlib
        plt.legend()
        plt.tick_params(axis='x', labelsize=8)

        # Save the plot
        # plt.savefig(os.path.join(plot_dir, plot_name))
        plt.savefig(os.path.join('/home/Earth/jcos/es-esmvaltool/output/', plot_name))
        plt.close()

def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Example(config).compute()

if __name__ == "__main__":
    main()
