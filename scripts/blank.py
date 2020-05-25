import iris
import esmvaltool.diag_scripts.shared
from esmvaltool.diag_scripts.shared import group_metadata
from esmvalcore.preprocessor._multimodel import _get_overlap, _assemble_full_data, _assemble_overlap_data


class Example(object):
    def __init__(self, config):
        self.cfg = config

    def compute(self):
        print('----------- ACCESS ----------')
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        ls = []
        print('----------- DATA RETRIEVED -------')
        for i, alias in enumerate(data):
            variables = group_metadata(data[alias], 'short_name')
            if 'pr' in variables:
                pr_file = variables['pr'][0]['filename']
                pr = iris.load(pr_file)[0]
                print(pr)
            if 'tas' in variables:
                tas_file = variables['tas'][0]['filename']
                tas = iris.load(tas_file)[0]
                print(tas.coord('time'))
            ls.append(tas)
            print(variables)
        cube_stats = self.multimodel_stats(ls)
        print(cube_stats['mean'].coord('time'))
        print('------------ WELL DONE ----------')

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


def main():
    print('----------- MAIN ----------')
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Example(config).compute()


if __name__ == "__main__":
    main()
