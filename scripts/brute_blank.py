import os

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.coord_categorisation
from iris.time import PartialDateTime
import cartopy.crs as ccrs
import numpy as np
from iris.coords import DimCoord
from iris.cube import Cube

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata
#from esmvaltool.diag_scripts.shared import plot
from esmvalcore.preprocessor._mask import mask_landsea
from esmvalcore.preprocessor._time import climate_statistics, anomalies, seasonal_statistics
from esmvalcore.preprocessor._regrid import regrid
from esmvalcore.preprocessor._multimodel import multi_model_statistics

class Example(object):
    def __init__(self, config):
        # ---------------------------------------------------------------------
        # config is a dictionary containing metadata regarding input files and
        # overall, as the name suggests, configuration options.
        # ---------------------------------------------------------------------
        self.cfg = config
        # self.extra_parametre = self.cfg.get('extra_parametre')

    def compute(self):
        print('----------- COMPUTE ----------')
        # ---------------------------------------------------------------------
        # Every dataset in the recipe is associated with an alias. We are going
        # to use th:We alias and the group_metadata shared function to loop over
        # the datasets.
        #----------------------------------------------------------------------
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        ssp_ts = {}
        hist_ts = {} 
        rean_ts = {}
        hist = 0
        ssp = 0
        rean = 0
        # Loop over the datasets.
        for i, alias in enumerate(data):
            exp = data[alias][0]['exp']
            variables = group_metadata(data[alias], 'short_name')
            # Returns the path to the preprocessed files.
            tas_file = variables['tas'][0]['filename']
            tas = iris.load(tas_file)[0]
            tas.convert_units('degC')
            if i == 0:
                climatology = self.ref_clim(tas, 1960, 1962)
            #anomaly = tas - climatology
            #timeseries = anomaly.collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
            timeseries = tas.collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
            #timeseries.long_name = 'med_r_timeseries_tas'
  
            # Calculate Trends
            nlat = tas.coord('latitude').shape[0]
            nlon = tas.coord('longitude').shape[0]
            lat = tas.coord('latitude').points
            lon = tas.coord('longitude').points
            time_array = np.arange(1,tas.coord('time').shape[0]+1,1)
            regr = np.zeros([nlat, nlon])
            for j in range(nlat):
                for k in range(nlon):
                    p = np.polyfit(time_array, tas[:,j,k].data, 1)
                    regr[j, k] = p[0]*10 # the 10 is to convert to decadal
            latitude = DimCoord(lat, standard_name='latitude', units='degrees')
            longitude = DimCoord(lon, standard_name='longitude', units='degrees')
            regr_cube = Cube(regr, dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
            ### ---------- remask -------------- ###
            # finding the trends turns the remask usseless as pyplot doesn't care about masked arrays
            # Ergo another remask is needed.
            output_trend = mask_landsea(regr_cube, ['/blablabla/where/the/fx/at/'] ,'sea', True)   
            # Save the output trends in the cube dict
            output_trend.standard_name = None 
            output_trend.long_name = 'tas_trend_med'
            output_trend.short_name = 'tastrend'
        
            # Calculate Climatology
            output_clim = climate_statistics(tas)           
            output_clim.standard_name = None
            output_clim.long_name = 'tas_clim_med'
            output_clim.short_name = 'tasclim' 

            # Save diagnosed dataset to dict. TODO: what about averaging first? 
            if exp == 'historical':
                hist_ts[alias] = timeseries 
                if hist == 0:
                   mean_hist_trend = output_trend
                   mean_hist_clim = output_clim
                   hist += 1
                else:
                   mean_hist_trend = (mean_hist_trend + output_trend)
                   mean_hist_clim = (mean_hist_clim + output_clim)
                   hist += 1
            if exp == 'ssp585':
                ssp_ts[alias] = timeseries 
                if ssp == 0:
                   mean_ssp_trend = output_trend
                   mean_ssp_clim = output_clim
                   ssp += 1
                else:
                   mean_ssp_trend = (mean_ssp_trend + output_trend)
                   mean_ssp_clim = (mean_ssp_clim + output_clim)
                   ssp += 1
            if exp == 'reanaly':
                rean_ts[alias] = timeseries 
                if rean == 0:
                   mean_rean_trend = output_trend
                   mean_rean_clim = output_clim
                   rean += 1
                else:
                   mean_rean_trend = (mean_rean_trend + output_trend)
                   mean_rean_clim = (mean_rean_clim + output_clim)
                   rean += 1
            
        mean_hist_trend = mean_hist_trend/hist
        mean_hist_clim = mean_hist_clim/hist
        mean_ssp_trend = mean_ssp_trend/ssp
        mean_ssp_clim = mean_ssp_clim/ssp
        #mean_rean_trend = mean_rean_trend/ssp
        #mean_rean_clim = mean_rean_clim/ssp
        mean_ssp_trend.long_name = 'ssp_trend_Med'
        mean_ssp_clim.long_name = 'ssp_clim_Med'
        mean_hist_trend.long_name = 'hist_trend_Med'
        mean_hist_clim.long_name = 'hist_clim_Med'
       # mean_rean_trend.long_name = 'rean_trend_Med'
       # mean_rean_clim.long_name = 'rean_clim_Med'
                
        ##### Biases #####
        #trend_bias = mean_hist_trend - mean_rean_trend
        #clim_bias = mean_hist_clim - mean_rean_clim


            # Save the outputs for each dataset.
            #self.save(output, alias, data)
        # Plot the results.
        self.plot_2D(mean_ssp_trend)
        self.plot_1D(timeseries)   
        print(mean_ssp_trend)  
 
    def plot_2D(self, result):
        extent = [-10, 40, 25, 50]
        #ax = plt.axes(projection=ccrs.AlbersEqualArea(central_longitude=c_lon, central_latitude=c_lat))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extent)
        fill = iplt.pcolormesh(result, vmin=0, vmax=1, coords=('longitude','latitude'), cmap=plt.cm.Reds)
        plt.gca().coastlines()
        title = 'Trend ssp'
        plt.title(title)
        cb = plt.colorbar(fill, orientation='horizontal', format='%.2f')
        plot_name = 'trend_med_tas_ssp.{out}'.format(
             out=self.cfg[n.OUTPUT_FILE_TYPE],
            )

        # Get the path to the plot directory
        plot_dir = self.cfg[n.PLOT_DIR]

        # Some tweaks using matplotlib
        #plt.legend()
        plt.tick_params(axis='x', labelsize=8)

        # Save the plot
        plt.savefig(os.path.join(plot_dir, plot_name))
        # plt.savefig(os.path.join('/home/Earth/jcos/es-esmvaltool/output/', plot_name))
        plt.close()

    def plot_1D(self, result):
        #for result in results:
        iplt.plot(result)
        title = 'Timeseries ssp'
        plt.title(title)
        plot_name = 'ts_med_tas_ssp.{out}'.format(
             out=self.cfg[n.OUTPUT_FILE_TYPE],
            )

        # Get the path to the plot directory
        plot_dir = self.cfg[n.PLOT_DIR]

        # Some tweaks using matplotlib
        #plt.legend()
        plt.tick_params(axis='x', labelsize=8)

        # Save the plot
        plt.savefig(os.path.join(plot_dir, plot_name))
        # plt.savefig(os.path.join('/home/Earth/jcos/es-esmvaltool/output/', plot_name))
        plt.close()


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
        iris.save(result, os.path.join(work_dir, output))
        #iris.save(result, os.path.join(iris.save('/home/Earth/jcos/es-esmvaltool/output/', output)))

    def ref_clim(self, cube, start_y,end_y):
        t1 = PartialDateTime(year=start_y)
        t2 = PartialDateTime(year=end_y)
        constr_1 = iris.Constraint(time=lambda t: t1 <= t.point <= t2)
        cube_out = cube.extract(constr_1)
        cube_out = cube_out.collapsed('time', iris.analysis.MEAN)
        return cube_out

def main():
    print('----------- MAIN ----------')
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Example(config).compute()

if __name__ == "__main__":
    main()
