import os

import iris
import datetime
import cf_units
from dateutil import relativedelta
# import iris.quickplot as qplt
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.colors as colors
import iris.plot as iplt
import iris.coord_categorisation
import cartopy.crs as ccrs
import numpy as np
# from iris.coords import AuxCoord
from iris.cube import Cube
# from iris.coords import DimCoord

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata
# from esmvaltool.diag_scripts.shared import plot
from esmvalcore.preprocessor._mask import mask_landsea
# anomalies,  climate_statistics, seasonal_statistics, timeseries_filter
from esmvalcore.preprocessor._time import timeseries_filter
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

    def retrieve_data(self):
        '''
        Here is where all the interesting stuff begins.

        (1) The function calls the variables and stores all the cubes according
        to their experiment (ssp585, historical, reanalysis...)

        (2) Calls the functions that calculates the output data to be shown

        (3) Calls the plotting functions
        '''
        print('----------- DATA RETREIVAL ----------')
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        hist_ls_tas, hist_ls_pr = [], []
        ssp_ls_tas, ssp_ls_pr = [], []
        rean_ls_tas, rean_ls_pr = [], []
        start_year_ssp, end_year_ssp, start_year_hist, end_year_hist, start_year_rean, end_year_rean = 0, 0, 0, 0, 0, 0
        w_density = iris.coords.AuxCoord(1000, long_name='water_density', units='kg m-3')
        # iteration that sorts the dataset cubes into it's own experiment list.
        # (hist_ls, ssp_ls, rean_ls)
        for i, alias in enumerate(data):
            print(alias)
            variables = group_metadata(data[alias], 'short_name')
            pr_exist, tas_exist = False, False
            if 'pr' in variables:
                pr_file = variables['pr'][0]['filename']
                pr = iris.load(pr_file)[0]
                pr = pr / w_density
                pr.convert_units('mm month-1')
                pr_exist = True
            if 'tas' in variables:
                tas_file = variables['tas'][0]['filename']
                tas = iris.load(tas_file)[0]
                tas.convert_units('degC')
                # coord_names = [coord.name() for coord in tas.coords()]
                tas_exist = True
            # LOAD DATA IN CUBES
            # OBS --> rean data, they don't have atribute 'exp'
            if 'OBS' not in alias:
                exp = data[alias][0]['exp']
                activity = data[alias][0]['activity']
            else:
                exp = 'a'
            # LIST THE CUBES ACCORDING TO THE 'EXP' THEY BELONG
            if exp == 'a':
                if (start_year_rean == 0) & (end_year_rean == 0):
                    start_year_rean = variables['tas'][0]['start_year']
                    end_year_rean = variables['tas'][0]['end_year']
                if tas_exist:
                    tas = self.regrid_time(tas, start_year_rean, end_year_rean)
                    rean_ls_tas.append(tas)
                if pr_exist:
                    pr = self.regrid_time(pr, start_year_rean, end_year_rean)
                    rean_ls_pr.append(pr)
            elif ((activity == 'HighResMIP')and(exp == 'highres-future')) or ((activity == 'cmip')and(exp == 'ssp585')):
                if (start_year_ssp == 0) & (end_year_ssp == 0):
                    start_year_ssp = variables['tas'][0]['start_year']
                    end_year_ssp = variables['tas'][0]['end_year']
                if tas_exist:
                    tas = self.regrid_time(tas, start_year_ssp, end_year_ssp)
                    ssp_ls_tas.append(tas)
                if pr_exist:
                    pr = self.regrid_time(pr, start_year_ssp, end_year_ssp)
                    ssp_ls_pr.append(pr)
            elif ((activity == 'HighResMIP')and(exp == 'hist-1950')) or ((activity == 'cmip')and(exp == 'historical')):
                if (start_year_hist == 0) & (end_year_hist == 0):
                    start_year_hist = variables['tas'][0]['start_year']
                    end_year_hist = variables['tas'][0]['end_year']
                if tas_exist:
                    tas = self.regrid_time(tas, start_year_hist, end_year_hist)
                    hist_ls_tas.append(tas)
                if pr_exist:
                    pr = self.regrid_time(pr, start_year_hist, end_year_hist)
                    hist_ls_pr.append(pr)

        # COMPUTE EVERYTHING & PLOT
        cube_ls_tas, ts_ls_tas = self.compute(rean_ls_tas, ssp_ls_tas, hist_ls_tas, 'tas')
#        cube_ls_pr, ts_ls_pr = self.compute(rean_ls_pr, ssp_ls_pr, hist_ls_pr, 'pr')
        self.tas_plot_caller(cube_ls_tas)
#        self.pr_plot_caller(cube_ls_pr)
        # # Compute standard deviation of the ensembles & plot
        # # TODO: Tidy the timeseries bit ---------
        STDs_tas = self.timeseries_std([ssp_ls_tas, hist_ls_tas, rean_ls_tas])
#        STDs_pr = self.timeseries_std([ssp_ls_pr, hist_ls_pr, rean_ls_pr])
        start_year_ls = [start_year_ssp, start_year_hist, start_year_rean]
        end_year_ls = [end_year_ssp, end_year_hist, end_year_rean]
        self.timeseries_plot(ts_ls_tas, start_year_ls, end_year_ls,
                             '_tas', STDs_tas, '($^o$C)', 'Temperature')
#        self.timeseries_plot(ts_ls_pr, start_year_ls, end_year_ls,
#                             '_pr', STDs_pr, '(mm month$^{-1}$)', 'Precipitation')

    def compute(self, rean_ls, ssp_ls, hist_ls, sht_nm):

        # Find the experiments mean and/or median (rean, ssp, hist)
        ssp_cube_stats = self.multimodel_stats(ssp_ls)
        hist_cube_stats = self.multimodel_stats(hist_ls)
        if len(rean_ls) < 2:  # safety blt in case only a rean dataset is avail.
            if len(rean_ls) != 0:
                rean_cube_stats = {'mean': rean_ls[0]}
            else:
                rean_cube_stats = {'mean': hist_ls[0]}
        else:
            rean_cube_stats = self.multimodel_stats(rean_ls)

        # extract the multimodel cube mean
        ssp_cube_mn = ssp_cube_stats['mean']
        hist_cube_mn = hist_cube_stats['mean']
        rean_cube_mn = rean_cube_stats['mean']
        # interval to calculate the trends [degC/deade
        interval = 10
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
        '''
        rean_trend_distr = self.trend_distr(rean_ls, rean_trends, rean_decades)
        hist_trend_distr = self.trend_distr(hist_ls, hist_trends, hist_decades)
        ssp_trend_distr = self.trend_distr(ssp_ls, ssp_trends, ssp_decades)
        # plot the histogram
        self.plot_histogram([rean_trend_distr, hist_trend_distr, ssp_trend_distr], [
                            'ERA40 reanalysis', 'CMIP6 hist.', 'CMIP6 ssp585'], sht_nm)
        '''

        # prepare 2D data to plot (regrid to avoid matplotlib funny behaviour)
        cube_ls_out = self.regrid_longitude_coord([ssp_trends, trend_bias, clim_bias, rean_trends])
        ssp_trends, trend_bias, clim_bias, rean_trends = cube_ls_out

        return cube_ls_out, [ssp_ts, hist_ts, rean_ts]

# ----------------------- stats, calc & util functions --------------------

    def regrid_time(self, cube, start_y, end_y):
        '''
        '''
        start = datetime.datetime(start_y, 1, 1)
        dt_array = np.array([start + relativedelta.relativedelta(years=i)
                             for i in range((end_y-start_y)+1)])

        new_t_unit_str = '{} since 1850-01-01 00:00:00'.format('days')
        new_t_unit = cf_units.Unit(new_t_unit_str, calendar=cf_units.CALENDAR_STANDARD)

        new_dt_points = [new_t_unit.date2num(new_dt) for new_dt in dt_array]
        new_t_coord = iris.coords.DimCoord(new_dt_points, standard_name='time', units=new_t_unit)

        t_coord_dim = cube.coord_dims('time')
        cube.remove_coord('time')
        cube.add_dim_coord(new_t_coord, t_coord_dim)
        return cube

    def multimodel_stats(self, cubes):
        '''
        Finds mean or median from a list of cubes
        use: get multimodel and/or ensemble stats.
        TODO: report issue with the 'ovelrap'
        '''
        # Get time overlaps
        interval = _get_overlap(cubes)
        # Define list of statistics to compute
        statistics = ['mean']
        # Define time-span, can be full or overlap
        span = 'full'  # 'full'  # 'overlap'
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
        Regrids a list of cubes.
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
            sorted_cube = self.sorted_dim(cube)
            cube_ls_out.append(sorted_cube)
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

    def timeseries_std(self, exp_ls_array):
        STDs = []
        for exp_ls in exp_ls_array:
            exp_ts_matrix = []
            for exp in exp_ls:
                exp_ts = exp.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
                exp_ts = timeseries_filter(
                    exp_ts, filter*2+1, filter, filter_type='lowpass', filter_stats='mean')
                exp_ts_matrix.append(exp_ts.data)
            # No iris function to compute STD between cubes?
            exp_ts_matrix = np.asarray(exp_ts_matrix)
            std_ts = np.std(exp_ts_matrix, axis=0)
            std_ts = std_ts.reshape(-1, )
            std_cube = Cube(std_ts, dim_coords_and_dims=[(exp_ts.coord('time'), 0)])
            STDs.append(std_cube)
        return STDs

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
        cube_timeseries = timeseries_filter(
            cube_timeseries, filter*2+1, filter, filter_type='lowpass', filter_stats='mean')
        coord_names_ts = [coord.name() for coord in cube_timeseries.coords()]
        if 'year' not in coord_names_ts:
            iris.coord_categorisation.add_year(cube_timeseries, 'time', name='year')
        return cube_trends, cube_anom, cube_mean, cube_timeseries

# ------------------------- ploting functions ----------------------------------

    def tas_plot_caller(self, tas_ls):
        ssp_trends_tas, trend_bias_tas, clim_bias_tas, rean_trends_tas = tas_ls
        # Plot the tas maps.
        deg = '$^o$C'
        deg_decade = '$^o$C decade$^{-1}$'
        title = 'CMIP6 SSP585 trends future period (2015-2050)'
        self.med_r_plot(
            mask_landsea(ssp_trends_tas, ['/blablabla/where/the/fx/at/'], 'sea', True),
            [0, 0.8, 0.05], title, 'ssp_trends_tas',
            plt.cm.Reds, deg_decade)
        title = 'CMIP6 historical trend bias (1960-2002)'
        self.med_r_plot(
            mask_landsea(trend_bias_tas, ['/blablabla/where/the/fx/at/'], 'sea', True),
            [-0.25, 0.25, 0.025], title, 'hist_trend_bias_tas',
            plt.cm.RdBu_r, deg_decade)
        title = 'CMIP6 historical climatological bias (1960-2014)'
        self.med_r_plot(
            mask_landsea(clim_bias_tas, ['/blablabla/where/the/fx/at/'],
                         'sea', True),
            [-4, 4, 0.4], title, 'hist_clim_bias_tas',
            plt.cm.RdBu_r, deg)
        title = 'AVG. reanalysis trend (1960-2002)'
        self.med_r_plot(
            mask_landsea(rean_trends_tas, ['/blablabla/where/the/fx/at/'], 'sea', True),
            [0, 0.5, 0.05], title, 'era40_trend_tas',
            plt.cm.Reds, deg_decade)

    def pr_plot_caller(self, pr_ls):
        ssp_trends_pr, trend_bias_pr, clim_bias_pr, rean_trends_pr = pr_ls
        # Plot the pr maps.
        pr = 'mm month$^{-1}$'
        pr_decade = 'mm month$^{-1}$ decade$^{-1}$'
        title = 'CMIP6 SSP585 trends future period (2015-2050)'
        self.med_r_plot(
            mask_landsea(ssp_trends_pr, ['/blablabla/where/the/fx/at/'], 'sea', True),
            [-10, 10, 1], title, 'ssp_trends_pr',
            plt.cm.BrBG, pr_decade)
        title = 'CMIP6 historical trend bias (1960-2014)'
        self.med_r_plot(
            mask_landsea(trend_bias_pr, ['/blablabla/where/the/fx/at/'], 'sea', True),
            [-10, 10, 1], title, 'hist_trend_bias_pr',
            plt.cm.BrBG, pr_decade)
        title = 'CMIP6 historical climatological bias (1960-2002)'
        self.med_r_plot(
            mask_landsea(clim_bias_pr, ['/blablabla/where/the/fx/at/'],
                         'sea', True),
            [-20, 20, 2], title, 'hist_clim_bias_pr',
            plt.cm.BrBG, pr)
        title = 'AVG. reanalysis trend (1960-2002)'
        self.med_r_plot(
            mask_landsea(rean_trends_pr, ['/blablabla/where/the/fx/at/'], 'sea', True),
            [-10, 10, 1], title, 'era40_trend_pr',
            plt.cm.BrBG, pr_decade)

    def region_to_square(self, region, dimension):
        if dimension == 'latitude':
            return([region['start_latitude'],
                    region['start_latitude'],
                    region['end_latitude'],
                    region['end_latitude'],
                    region['start_latitude']])
        elif dimension == 'longitude':
            return([region['start_longitude'],
                    region['end_longitude'],
                    region['end_longitude'],
                    region['start_longitude'],
                    region['start_longitude']])
        else:
            return('dimension unknown')

    def med_r_plot(self, cube, v_min_max, title, plot_name, cmap, units):
        '''
        Plots 2D data over the mediterranean region
        and saves it at {recipe_folder}/{plots}...
        '''
        vm, vx, stp = v_min_max
        region = {'start_longitude': -10,
                  'end_longitude': 40,
                  'start_latitude': 25,
                  'end_latitude': 50}
        projection = 'LambertConformal'
        plotextend = [region['start_longitude'],
                      region['end_longitude'],
                      region['start_latitude'],
                      region['end_latitude']]
        if projection == 'LambertConformal':
            # plotextend has to be a little larger so everything is on there
            plotextend = [plotextend[0]-1., plotextend[1]+1.,
                          plotextend[2]-1., plotextend[3]+1.]
            # path to cut out is exact though
            lons = self.region_to_square(region, 'longitude')
            lats = self.region_to_square(region, 'latitude')
            path_ext = [[lon, lat] for lon, lat in zip(lons, lats)]
            path_ext = mpath.Path(path_ext).interpolated(20)

        proj = ccrs.LambertConformal(central_longitude=np.sum(
            plotextend[: 2])/2., central_latitude=np.sum(plotextend[2:])/2.)
        plt.figure()
        ax = plt.axes(projection=proj)

        # now cut the plot out if its is LambertConformal
        if projection == 'LambertConformal':
            proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax) - ax.transData
            rect_in_target = proj_to_data.transform_path(path_ext)
            ax.set_boundary(rect_in_target, use_as_clip_path=True)

        bounds = np.arange(vm, vx+stp, stp)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        ax.set_extent(plotextend, crs=ccrs.PlateCarree())
        ax.coastlines('50m', linewidth=0.8)
        fill = iplt.pcolormesh(cube, norm=norm,
                               coords=('longitude', 'latitude'), cmap=cmap)
        cb = plt.colorbar(fill, orientation='horizontal', format='%.2f')
        cb.set_label(units)
        plt.title(title)
        plt.tight_layout()
        self.save_fig(plot_name)

    def plot_histogram(self, exp_trend_distr, label, sht_nm):
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
        self.save_fig('histogram'+sht_nm)

    # def rearrange_time(self, timeseries):
        # t = timeseries.coord('time')
        # cube_out = Cube(timeseries.data, dim_coords_and_dims=[(t, 0)])
        # return cube_out

    # def timeseries_plot(self, ssp_ts, hist_ts, rean_ts, sht_nm, STDs, units):
    def timeseries_plot(self, ts_ls, start_years, end_years, sht_nm, STDs, units, var):
        '''
        Timeseries plot.
        TODO: solve the x-axis appearence issue (coord.points)
        '''
        plt.figure()
        for ts, STD, start_y, end_y in zip(ts_ls, STDs, start_years, end_years):
            iplt.plot(ts.coord('year'), ts, label='.')
            # iplt.plot(ts, label='CMIP6 ssp585')
            start = datetime.datetime(start_y+filter, 1, 1)
            dt_array = np.array([start + relativedelta.relativedelta(years=i)
                                 for i in range((end_y-start_y)-filter*2+1)])
            y_dt_arr = [x.year for x in dt_array]
            minus, plus = ts.data-STD.data, ts.data+STD.data
            plt.fill_between(y_dt_arr, minus, plus, alpha=0.2)

        plt.title('Mediterranean '+var+' anomalies (baseline period: 1960-2002)')
        plt.ylabel(units)
        plt.xlabel('(iris cube time coord points. NEEDS FIX)')
        plt.grid()
        plt.legend()
        self.save_fig('timeseries'+sht_nm)

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
        Example(config).retrieve_data()


if __name__ == "__main__":
    filter = 5
    main()
