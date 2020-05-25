import os

import iris
# import iris.quickplot as qplt
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.coord_categorisation
import cartopy.crs as ccrs
import numpy as np

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata
# from esmvaltool.diag_scripts.shared import plot
from esmvalcore.preprocessor._mask import mask_landsea
# from esmvalcore.preprocessor._time import anomalies  # ,  climate_statistics, seasonal_statistics, timeseries_filter
# from esmvalcore.preprocessor._regrid import regrid
from esmvalcore.preprocessor._detrend import detrend
from esmvalcore.preprocessor._multimodel import _get_overlap, _assemble_full_data, _assemble_overlap_data

'''
SCRIPT STRUCTURE:

class:
|
--> config (main function inside class)
|
--> statistics and calc funtions
|
--> plotting functions (maps, timeseries, histogram)
'''


class Example(object):
    def __init__(self, config):
        # ---------------------------------------------------------------------
        # config is a dictionary containing metadata regarding input files and
        # overall, as the name suggests, configuration options.
        # ---------------------------------------------------------------------
        self.cfg = config

    def compute(self):
        '''
        Here is where all the interesting stuff begins.

        (1) The function calls the variables and stores all the cubes according
        to their experiment (ssp585, historical, reanalysis...)

        (2) Calls the functions that calculates the output data to be shown

        (3) Calls the plotting functions
        '''
        print('----------- COMPUTE ----------')
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        hist_ls = []
        ssp_ls = []
        rean_ls = []
        # iteration that sorts the dataset cubes into it's own experiment list.
        # (hist_ls, ssp_ls, rean_ls)
        for i, alias in enumerate(data):
            print(alias)
            variables = group_metadata(data[alias], 'short_name')
            tas_file = variables['tas'][0]['filename']

            # TODO:
            # pr_file = variables['pr'][0]['filename']
            # pr = iris.load(tas_file)[0]
            # pr.convert_units('mm')

            # LOAD DATA IN CUBES
            # OBS --> rean data, they don't have atribute 'exp'
            if alias != 'OBS':
                exp = data[alias][0]['exp']
                tas = iris.load(tas_file)[0]
                tas.convert_units('degC')
            else:
                tas = iris.load(tas_file)[0]
                tas.convert_units('degC')
                exp = 'a'
            # LIST THE CUBES ACCORDING TO THE 'EXP' THEY BELONG
            if exp == 'a':
                rean_ls.append(tas)
            if exp == 'ssp585':
                ssp_ls.append(tas)
            if exp == 'historical':
                hist_ls.append(tas)

        # Find the experiments mean and/or median (rean, ssp, hist)
        ssp_cube_stats = self.multimodel_stats(ssp_ls)
        hist_cube_stats = self.multimodel_stats(hist_ls)
        if len(rean_ls) < 2:  # safety blt in case only a rean dataset is avail.
            rean_cube_stats = {'mean': rean_ls[0]}
        else:
            rean_cube_stats = self.multimodel_stats(rean_ls)
        # extract the multimodel cube mean
        ssp_cube_mn = ssp_cube_stats['mean']
        hist_cube_mn = hist_cube_stats['mean']
        rean_cube_mn = rean_cube_stats['mean']
        # interval to calculate the trends [degC/deade
        interval = 10  # decadal
        # conversion factor from yearly to decadal
        hist_decades = len(hist_cube_mn.coord('time').points)/interval
        rean_decades = len(rean_cube_mn.coord('time').points)/interval
        ssp_decades = len(ssp_cube_mn.coord('time').points)/interval
        # Calculates the statistics needed for each multimodel experiment
        ssp_trends, ssp_anom, ssp_mean, ssp_ts = self.trend_anom_mean(
            ssp_cube_mn, ssp_decades, past=False)
        hist_trends, hist_anom, hist_mean, hist_ts = self.trend_anom_mean(
            hist_cube_mn, hist_decades)
        rean_trends, rean_anom, rean_mean, rean_ts = self.trend_anom_mean(
            rean_cube_mn, rean_decades)
        # finds the ssp585 anomaly timeseries, outside the function,
        # as it needs to be calculated with historical data
        ssp_ts_detrender = hist_mean.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
        ssp_ts = ssp_ts - ssp_ts_detrender.data
        # trend biases #
        trend_bias = hist_trends - rean_trends
        # climatological biases #
        clim_bias = hist_mean - rean_mean
        # Trend distributions
        rean_trend_distr = self.trend_distr(rean_ls, rean_trends, rean_decades)
        hist_trend_distr = self.trend_distr(hist_ls, hist_trends, hist_decades)
        ssp_trend_distr = self.trend_distr(ssp_ls, ssp_trends, ssp_decades)
        # plot the histogram
        self.plot_histogram([rean_trend_distr, hist_trend_distr, ssp_trend_distr], [
                            'ERA40 reanalysis', 'CMIP6 hist.', 'CMIP6 ssp585'])

        # TODO: filter the timeseries
        # ssp_ts = timeseries_filter(ssp_ts, 10, 5)
        self.timeseries_plot(ssp_ts, hist_ts, rean_ts)

        # prepare 2D data to plot (regrid to avoid matplotlib funny behaviour)
        cube_ls_out = self.regrid_longitude_coord([ssp_trends, trend_bias, clim_bias, rean_trends])
        ssp_trends, trend_bias, clim_bias, rean_trends = cube_ls_out
        # Plot the maps.
        deg = '$^o$C'
        deg_decade = '$^o$C decade$^{-1}$'
        title = 'CMIP6 SSP585 trends future period (2015-2050)'
        self.med_r_plot(
            mask_landsea(ssp_trends, ['/blablabla/where/the/fx/at/'], 'sea', True),
            [0, 0.8], title, 'ssp_trends',
            plt.cm.Reds, deg_decade)
        title = 'CMIP6 historical trend bias (1960-2014)'
        self.med_r_plot(
            mask_landsea(trend_bias, ['/blablabla/where/the/fx/at/'], 'sea', True),
            [-0.25, 0.25], title, 'hist_trend_bias',
            plt.cm.RdBu_r, deg_decade)
        title = 'CMIP6 historical climatological bias (1960-2014)'
        self.med_r_plot(
            mask_landsea(clim_bias, ['/blablabla/where/the/fx/at/'],
                         'sea', True),
            [-4, 4], title, 'hist_clim_bias',
            plt.cm.RdBu_r, deg)
        title = 'AVG. reanalysis trend (1960-2002)'
        self.med_r_plot(
            mask_landsea(rean_trends, ['/blablabla/where/the/fx/at/'], 'sea', True),
            [0, 0.5], title, 'era40_trend',
            plt.cm.Reds, deg_decade)

# ----------------------- stats, calc & util functions ------------------------

    def multimodel_stats(self, cubes):
        '''
        Finds mean or median from a list of cubes
        use: get multimodel and/or ensemble stats.
        '''
        # Get time overlaps
        interval = _get_overlap(cubes)
        # Define list of statistics to compute
        statistics = ['mean']
        # Define time-span, can be full or overlap
        span = 'full'  # 'overlap'
        # Compute statistics and store resulting cube in a dict
        # where each entry is an statistic
        results = {}
        for statistic in statistics:
            if span == 'overlap':
                statistic_cube = _assemble_overlap_data(cubes, interval, statistic)
            elif span == 'full':
                statistic_cube = _assemble_full_data(cubes, statistic)
            results[statistic] = statistic_cube
        return results

    def regrid_longitude_coord(self, cube_ls):
        '''
        Regrids a lis of cubes.
        Regrids each cube longitude from a 0/360deg shape to -180/180 deg
         - Useful for plotting regions that sit on top of the 0deg  meridian
        '''
        cube_ls_out = []
        for cube in cube_ls:
            # make a list with the 'longitude' coord in the form: 0/180/-180/0
            neg_lons = ((cube.coord('longitude').points + 180) % 360)-180
            # interpolates the cube data to the new 'longitude' dimensions
            cube = cube.interpolate([('longitude', neg_lons)], iris.analysis.Linear())
            # sorts the dimensions to the form: -180/0/180
            cube_ls_out.append(self.sorted_dim(cube))
        return cube_ls_out

    def sorted_dim(self, cube, coord='longitude'):
        '''
        Sorts the cube data according to the longitude coordinate values
        for example: 180/-180 --> -180/180
        (useful to avoid messy map plots)
        '''
        coord_to_sort = cube.coord(coord)
        assert coord_to_sort.ndim == 1, 'Coord should be 1-dimensional.'
        dim, = cube.coord_dims(coord_to_sort)
        index = [slice(None)] * cube.ndim
        index[dim] = np.argsort(coord_to_sort.points)
        return cube[tuple(index)]

    def cube_trend(self, cube, decades):
        '''
        Finds the trend of a cube by subtracting the detrended cube
        (ESMVal preproc) to the original and finding the slope of the
        resultant cube by subtracting the first position (0) to the last (-1)
     (data)
        |            (-1)
        |          /
        |        /
        |      /
        |    /
        | (0)
        ------------------(time)
        '''
        cube_detr = detrend(cube)
        cube_trend_val = cube - cube_detr
        cube_trends = (cube_trend_val[-1]-cube_trend_val[0])/decades
        return cube_trends

    def trend_distr(self, exp_ls, exp_trends, exp_decades):
        '''
        It calculates the domain mean trends of all the models of an experiment
        (degC/decade)
        and lists the values.
        '''
        exp_trend_distr = []
        for exp in exp_ls:
            # find the model trend
            e_trd = self.cube_trend(exp, exp_decades)
            # find the trend spatial average of the model
            exp_trend = e_trd.collapsed(['latitude', 'longitude'], iris.analysis.MEAN).data
            exp_trend_distr.append(float(exp_trend))
        return exp_trend_distr

    def trend_anom_mean(self, cube, decades, past=True):
        '''
        Finds the trend, anomaly, mean and timeseries of a given cube in the working domain/region.
        TODO: implement the anomalies ESMVal preproc and use CMIP6hist exp as the baseline period for CMIP6ssp585
        '''
        cube_mean = cube.collapsed('time', iris.analysis.MEAN)
        # Baseline currently not implementd
        # baseline_period = {'start_year': 1995,
        #                     'start_month': 1,
        #                     'start_day': 1,
        #                     'end_year': 2002,
        #                     'end_month': 1,
        #                     'end_day': 1}
        # if the experiment is not from the past period a different baseline period needs to be defined
        # TODO: imp. anomalies & past baseline for future projections
        if past:
            # cube_anom = anomalies(cube, 'full', baseline_period)
            cube_anom = cube - cube_mean  # crappy way of finding the anomaly, no baseline perid, the whole is used
        else:
            # if the exp is a projection the anomaly is computed inside self.config(...)
            cube_anom = cube
        cube_trends = self.cube_trend(cube, decades)
        # timeseries
        cube_timeseries = cube_anom.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
        return cube_trends, cube_anom, cube_mean, cube_timeseries

# ------------------------- ploting functions ----------------------------------

    def med_r_plot(self, cube, v_min_max, title, plot_name, cmap, units):
        '''
        Plots 2D data over the mediterranean region
        and saves it at {recipe_folder}/{plots}...
        '''
        plt.figure()
        extent = [-10, 40, 25, 50]
        vm, vx = v_min_max
        # TODO: use a sexier projection.
        # proj=ccrs.AlbersEqualArea(central_longitude=25, central_latitude=37.5) #ccrs.PlateCarree()
        proj = ccrs.PlateCarree()
        ax = plt.axes(projection=proj)
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        fill = iplt.pcolormesh(cube, vmin=vm, vmax=vx, coords=('longitude', 'latitude'), cmap=cmap)
        ax.coastlines('50m', linewidth=0.8)
        plt.title(title)
        cb = plt.colorbar(fill, orientation='horizontal', format='%.2f')
        cb.set_label(units)
        plt.tight_layout()
        self.save_fig(plot_name)

    def plot_histogram(self, exp_trend_distr, label):
        '''
        Plots a histogram with the trends of all the members of each experiments
        '''
        plt.figure()
        bins = 5
        for distribution, l in zip(exp_trend_distr, label):
            plt.hist(distribution, bins, alpha=0.5, label=l)
        plt.legend()
        plt.xlabel('ternd ($^o$C decade$^{-1}$)')
        plt.ylabel('number of models')
        self.save_fig('histogram')

    def timeseries_plot(self, ssp_ts, hist_ts, rean_ts):
        '''
        Timeseries plot.
        TODO: solve the x-axis appearence issue (coord.points)
        TODO: plot the variance shading
        '''
        plt.figure()
        plt.plot(ssp_ts.coord('time').points, ssp_ts.data, label='CMIP6 ssp585')
        plt.plot(hist_ts.coord('time').points, hist_ts.data, label='CMIP6 hist')
        plt.plot(rean_ts.coord('time').points, rean_ts.data, label='ERA40 reanalysis')
        plt.title('Mediterranean Temperature anomalies (baseline period: 1960-2002)')
        plt.ylabel('($^o$C)')
        plt.xlabel('iris cube time coord points. NEEDS FIX')
        plt.grid()
        plt.legend()
        self.save_fig('timeseries')

    def save_fig(self, pname):
        '''
        save a plot
        '''
        plot_name = 'martin_{pname}.{out}'.format(
            out=self.cfg[n.OUTPUT_FILE_TYPE],
            pname=pname
        )
        # Get the path to the plot directory
        plot_dir = self.cfg[n.PLOT_DIR]
        # Some tweaks using matplotlib
        # plt.legend()
        # plt.tick_params(axis='x', labelsize=8)
        plt.savefig(os.path.join(plot_dir, plot_name))
#        plt.savefig('/esarchive/scratch/jcos/esmvaltool/latest_recipie/plots/mean/mean/'+plot_name)
        plt.close()


def main():
    print('----------- MAIN ----------')
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Example(config).compute()


if __name__ == "__main__":
    main()
