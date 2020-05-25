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
            with open("/esarchive/scratch/jcos/esmvaltool/scripts/sample.txt", "a") as f:
                f.write(str(cube))
                f.write('******')
            window = 96
            wgts = self.low_pass_weights(window, 1. / 84.)
            filt = cube.rolling_window('time', iris.analysis.MEAN, len(wgts), weights=wgts)
            #tas = annual_statistics(cube, 'mean')
            #regional_month_mean = tas.collapsed(['longitude', 'latitude'], iris.analysis.MEAN) 

            #------------------------------------------------------------------
            # Substract the data and create a new cube with the result but
            # keeping tas' metadata
            #------------------------------------------------------------------
            filt.standard_name = None
            filt.long_name = 'filtered tas'
            filt.short_name = alias
            total[alias] = filt

        # Plot the results
        self.plot1D(total, data, 'flt_cmip6_proj_ssp126.png')


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
         for result in results:
             dataset = data[result][0]['dataset']
             exp = data[result][0]['exp']
             iplt.plot(results[result], linewidth=2., linestyle='-', alpha=1., label=dataset+' '+exp)
         plt.title('TAS 7-year filter applied')
         plt.xlabel('Time')
         plt.ylabel('T [degC]')
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


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Example(config).compute()

if __name__ == "__main__":
    main()