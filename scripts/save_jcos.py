import os

import iris

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

class Example(object):
    def __init__(self, config):
        self.cfg = config
        self.latestversion = self.cfg.get('latestversion')

    def compute(self):
        rootpath = '/esarchive/exp/ecearth'
        data = group_metadata(self.cfg['input_data'].values(), 'project')

        for project in data:
            if 'ECEARTH' in project:
                self.drs_ecearth(rootpath, data, 'ECEARTH')

            if 'DCPP' in project:
                self.drs_dcpp(rootpath, data, 'DCPP')

    def extract_date(self, cube):
        datetime = cube.coord('time')

        start_datetime = datetime.points[0]
        start_date = datetime.units.num2date(start_datetime)
        start_year = start_date.year
        start_month = start_date.month
        start = str(start_year)+'{:02d}'.format(start_month)

        end_datetime = datetime.points[-1]
        end_date = datetime.units.num2date(end_datetime)
        end_year = end_date.year
        end_month = end_date.month
        end = str(end_year)+'{:02d}'.format(end_month)

        return start, end

    def save(self, cube, dirpath, filename):
        if not os.path.exists(dirpath):
            os.makedirs(dirpath, exist_ok = True)
        iris.save(cube, os.path.join(dirpath, filename))


    def drs_ecearth(self, rootpath, data, project):
        for var_index, var_entry in enumerate(data[project]):
            cube = iris.load_cube(data[project][var_index]['filename'])
            expid = data[project][var_index]['expid']
            activity = data[project][var_index]['activity']
            if isinstance(activity, list):
                activity = ''.join(activity)
            institute = data[project][var_index]['institute']
            if isinstance(institute, list):
                institute = ''.join(institute)
            dataset = data[project][var_index]['dataset']
            exp = data[project][var_index]['exp']
            ensemble = data[project][var_index]['ensemble']
            mip = data[project][var_index]['mip']
            shortname = (
                data[project][var_index]['short_name'] + self.cfg[n.SCRIPT]
                )
            grid = data[project][var_index]['grid']
            latestversion = self.latestversion
            start, end = self.extract_date(cube)

            jcospath = '/home/Earth/jcos'
            fullpath = os.path.join(
                jcospath,
                expid,
                shortname,
            )

            filename = ('{shortname}_'
                        '{mip}_'
                        '{dataset}_'
                        '{exp}_'
                        '{ensemble}_'
                        '{grid}_'
                        '{start}-'
                        '{end}.nc').format(
                            shortname=shortname,
                            mip=mip,
                            dataset=dataset,
                            exp=exp,
                            ensemble=ensemble,
                            grid=grid,
                            start=start,
                            end=end
                        )
            self.save(cube, fullpath, filename)


    def drs_dcpp(self, rootpath, data, project):
        for var_index, var_entry in enumerate(data[project]):
            cube = iris.load_cube(data[project][var_index]['filename'])
            expid = data[project][var_index]['expid']
            activity = data[project][var_index]['activity']
            if isinstance(activity, list):
                activity = ''.join(activity)
            institute = data[project][var_index]['institute']
            if isinstance(institute, list):
                institute = ''.join(institute)
            dataset = data[project][var_index]['dataset']
            exp = data[project][var_index]['exp']
            startdate = data[project][var_index]['startdate']
            ensemble = data[project][var_index]['ensemble']
            mip = data[project][var_index]['mip']
            shortname = (
                data[project][var_index]['short_name'] + self.cfg[n.SCRIPT]
                )
            grid = data[project][var_index]['grid']
            latestversion = self.latestversion
            start, end = self.extract_date(cube)

            jcospath = '/home/Earth/jcos'
            fullpath = os.path.join(
                jcospath,
                expid,
                shortname,
            )

            filename = ('{shortname}_'
                        '{mip}_'
                        '{dataset}_'
                        '{exp}_'
                        '{startdate}-'
                        '{ensemble}_'
                        '{grid}_'
                        '{start}-'
                        '{end}.nc').format(
                            shortname=shortname,
                            mip=mip,
                            dataset=dataset,
                            exp=exp,
                            startdate=startdate,
                            ensemble=ensemble,
                            grid=grid,
                            start=start,
                            end=end
                        )

            self.save(cube, fullpath, filename)

def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Example(config).compute()

if __name__ == "__main__":
    main()

