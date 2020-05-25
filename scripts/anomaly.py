import os
import numpy as np
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

            tas = iris.load(tas_file)[0]
            tas.convert_units('degC')
            iris.coord_categorisation.add_year(tas, 'time', name='year')
            start_year, end_year = self.extract_date(tas)
           
            if end_year-start_year >= 30:
                first_30 = tas.extract(iris.Constraint(year=lambda t: start_year < t.point < start_year+30)).collapsed('time', iris.analysis.MEAN)
            else: 
                raise ValueError("data series shorter than 30 years :(")
            last_10 = tas.extract(iris.Constraint(year=lambda t: end_year-11 < t.point < end_year+1)).collapsed('time', iris.analysis.MEAN)
            #year_mean = tas.aggregated_by('year', iris.analysis.MEAN)
            #year_mean = seasonal_statistics(tas, 'mean')
            #year_anom_sfc = year_mean - first_30
            #year_anom = year_anom_sfc.collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
            
            dec_anom = last_10 - first_30 

            #------------------------------------------------------------------
            # Substract the data and create a new cube with the result but
            # keeping tas' metadata
            #------------------------------------------------------------------
            #cube = tas

            # Calling a preproc function inside a diag_script
            #timeseries = area_statistics(cube, 'mean')

            dec_anom.standard_name = None
            dec_anom.long_name = 'medR_anom_tas'
            dec_anom.short_name = 'tas'+str(i)
            total[alias] = dec_anom
            with open("/esarchive/scratch/jcos/esmvaltool/scripts/sample.txt", "a") as f:
                f.write(str(dec_anom))
                f.write('******')
                f.write(str(total))
            # Save the outputs for each dataset.
            #self.save(timeseries, alias, data)

        # Plot the results
        self.plot2D(total, data, 'medRarea_anom_tas.png')
    
    def extract_date(self, cube):
        datetime = cube.coord('time')
        start_datetime = datetime.points[0]
        start_date = datetime.units.num2date(start_datetime)
        start_year = start_date.year
        end_datetime = datetime.points[-1]
        end_date = datetime.units.num2date(end_datetime)
        end_year = end_date.year
        return start_year, end_year

    def plot1D(self, results, data, name):
         plt.figure(figsize=(8, 6))
         for result in results:
            dataset = data[result][0]['dataset']
            qplt.plot(results[result], label=dataset)
         plt.legend()
         plt.grid()
         plt.savefig(os.path.join('/esarchive/scratch/jcos/esmvaltool/output/figures/', name))

    def plot2D(self, results, data, name):
      
         for result in results:
            dataset = data[result][0]['dataset']
            central_lon, central_lat = 20, 30
            extent = [-10, 40, 20, 60]
            plt.figure(figsize=(8,6))
            ax = plt.axes(projection=ccrs.Orthographic(central_lon, central_lat))
            ax.set_extent(extent)
            contour = qplt.contourf(results[result])
            plt.gca().coastlines()
            plt.clabel(contour, inline=False)
            plt.title('2000-05 period temperature anomaly, with respect to 1950-1980 climatology. '+dataset)
            plt.savefig(os.path.join('/esarchive/scratch/jcos/esmvaltool/output/figures/', dataset+name))
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
