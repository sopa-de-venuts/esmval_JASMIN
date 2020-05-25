import os

import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import iris.coord_categorisation
from iris.time import PartialDateTime
import cartopy.crs as ccrs
import numpy as np
import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata
# from esmvaltool.diag_scripts.shared import plot
from esmvalcore.preprocessor._area import area_statistics
from esmvalcore.preprocessor._time import annual_statistics, climate_statistics, anomalies, seasonal_statistics

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
        with open("/esarchive/scratch/jcos/esmvaltool/scripts/sample.txt", "w") as f:
            f.write('******')
        # Loop over the datasets.
        for i, alias in enumerate(data):
            # -----------------------------------------------------------------
            # Use the group_metadata function again so that, for each dataset,
            # the metadata dictionary is organised by the variables'
            # short name.
            # -----------------------------------------------------------------
            variables = group_metadata(data[alias], 'short_name')
            # Returns the path to the preprocessed files.
            tas_file = variables['tas'][0]['filename']

            cube = iris.load(tas_file)[0]
            cube.convert_units('degC')
            window = 96
            wgts = self.low_pass_weights(window, 1. / 84.)
            filt = cube.rolling_window('time', iris.analysis.MEAN, len(wgts), weights=wgts)
            #tas = annual_statistics(cube, 'mean')
            #regional_month_mean = tas.collapsed(['longitude', 'latitude'], iris.analysis.MEAN) 
            #------------------------------------------------------------------
            # Substract the data and create a new cube with the result but
            # keeping tas' metadata
            #------------------------------------------------------------------
            tas = cube # = tas

            # Calling a preproc function inside a diag_script
            #timeseries = area_statistics(cube, 'mean')

            tas.standard_name = None
            tas.long_name = 'tas'
            tas.short_name = alias
            total[alias] = tas
            filt.standard_name = None
            filt.long_name = 'filtered tas'
            filt.short_name = alias+'_filt'
            total[alias+'filt'] = filt
            with open("/esarchive/scratch/jcos/esmvaltool/scripts/sample.txt", "a") as f:
                f.write(str(tas))
                f.write('******')
                f.write(str(total))
                f.write('******')

            # Save the outputs for each dataset.
            #self.save(timeseries, alias, data)

        # Plot the results
        self.plot1D(total, data, 'filt_proj_mean_tas.png')


    def low_pass_weights(self, window, cutoff):
         """Calculate weights for a low pass Lanczos filter.

         Args:

         window: int
             The length of the filter window.

         cutoff: float
            The cutoff frequency in inverse time steps.

         """
         order = ((window - 1) // 2) + 1
         nwts = 2 * order + 1
         w = np.zeros([nwts])
         n = nwts // 2
         w[n] = 2 * cutoff
         k = np.arange(1., n)
         sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
         firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
         w[n-1:0:-1] = firstfactor * sigma
         w[n+1:-1] = firstfactor * sigma
         return w[1:-1]

    def plot1D(self, results, data, name):
         plt.figure(figsize=(9, 4))
         colors = ['#ff0000','#ff0000','#0000ff','#0000ff','#00ff00','#00ff00']
         for result,c in zip(results,colors):
            if 'filt' in result:
                iplt.plot(results[result], color=c, linewidth=2., linestyle='-', alpha=1., label=dataset+' '+exp)
            else:
                dataset = data[result][0]['dataset']
                exp = data[result][0]['exp']
                iplt.plot(results[result], color=c, linewidth=1., linestyle='-', alpha=.5)
         plt.title('TAS 7-year filter applied')
         plt.xlabel('Time')
         plt.ylabel('T [K]')
         plt.legend(fontsize=10)
         plt.savefig(os.path.join('/esarchive/scratch/jcos/esmvaltool/output/figures/', name))

    def plot2D(self, cube, name):
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
