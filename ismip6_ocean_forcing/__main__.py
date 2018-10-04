"""
Script for creating ISMIP6 ocean forcing
"""

import os
import argparse
import pkg_resources
from configparser import ConfigParser, ExtendedInterpolation

# make sure the 'Agg' backend is used before any local module get loaded
import matplotlib
matplotlib.use('Agg')

import ismip6_ocean_forcing
from ismip6_ocean_forcing.rignot2013.remap import rignot_to_ismip6_grid
from ismip6_ocean_forcing.bedmap2 import bedmap2_to_ismip6_grid
from ismip6_ocean_forcing.imbie import make_imbie_masks
from ismip6_ocean_forcing.woa.main import extrapolate_woa
from ismip6_ocean_forcing.model.extrap import extrapolate_model
from ismip6_ocean_forcing.model.anomaly import compute_anomaly_and_to_woa


def main():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('configFiles', metavar='CONFIG',
                        type=str, nargs='*', help='config file')
    parser.add_argument('-v', '--version',
                        action='version',
                        version='ismip6_ocean_forcing {}'.format(
                                ismip6_ocean_forcing.__version__),
                        help="Show version number and exit")
    args = parser.parse_args()

    for configFile in args.configFiles:
        if not os.path.exists(configFile):
            raise OSError('Config file {} not found.'.format(configFile))

    # add config.default to cover default not included in the config files
    # provided on the command line
    if pkg_resources.resource_exists('ismip6_ocean_forcing', 'config.default'):
        defaultConfig = pkg_resources.resource_filename('ismip6_ocean_forcing',
                                                        'config.default')
        configFiles = [defaultConfig] + args.configFiles
    else:
        print('WARNING: Did not find config.default.  Assuming other config '
              'file(s) contain a\n'
              'full set of configuration options.')
        defaultConfig = None
        configFiles = args.configFiles

    config = ConfigParser(interpolation=ExtendedInterpolation())
    config.read(configFiles)

    rignot_to_ismip6_grid(config)
    bedmap2_to_ismip6_grid(config)
    make_imbie_masks(config)
    extrapolate_woa(config)
    extrapolate_model(config)
    compute_anomaly_and_to_woa(config)


if __name__ == "__main__":
    main()
