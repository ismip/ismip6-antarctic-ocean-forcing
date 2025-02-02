"""
Script for creating ISMIP6 ocean forcing
"""

import os
import argparse
from importlib.resources import path
from configparser import ConfigParser, ExtendedInterpolation

from ismip6_ocean_forcing.version import __version__
from ismip6_ocean_forcing.rignot2013.remap import rignot_to_ismip6_grid
from ismip6_ocean_forcing.bedmap2 import bedmap2_to_ismip6_grid
from ismip6_ocean_forcing.imbie import make_imbie_masks
from ismip6_ocean_forcing.obs.main import process_obs
from ismip6_ocean_forcing.model.extrap import extrapolate_model
from ismip6_ocean_forcing.model.anomaly import compute_anomaly_and_to_obs


def main():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('config_files', metavar='CONFIG',
                        type=str, nargs='*', help='config file(s)')
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f'ismip6_ocean_forcing {__version__}',
                        help="Show version number and exit")
    parser.add_argument('--model_field', dest='model_field', type=str,
                        help="A field (temperature or salinity) to "
                             "extrapolate")
    parser.add_argument('--model_basin', dest='model_basin', type=str,
                        help="A basin (none, 1 to 16 or open_ocean) to "
                             "extrapolate")
    parser.add_argument('--model_combine_basins', dest='model_combine_basins',
                        action='store_true',
                        help="Whether to only combine basins that have "
                             "already been extrapolated")
    args = parser.parse_args()

    for config_file in args.config_files:
        if not os.path.exists(config_file):
            raise OSError(f'Config file {config_file} not found.')

    config_files = list()
    with path('ismip6_ocean_forcing', 'default.cfg') as default_config:
        config_files.append(str(default_config))
    config_files.extend(args.config_files)

    config = ConfigParser(interpolation=ExtendedInterpolation())
    config.read(config_files)
    if args.model_field is not None:
        config.set('model', 'field', args.model_field)
    if args.model_basin is not None:
        config.set('model', 'basin', args.model_basin)
    if args.model_combine_basins:
        config.set('model', 'combineBasins', 'True')
    elif args.model_basin is not None and args.model_basin != 'all':
        config.set('model', 'combineBasins', 'False')

    rignot_to_ismip6_grid(config)
    bedmap2_to_ismip6_grid(config)
    make_imbie_masks(config)
    process_obs(config)
    extrapolate_model(config)
    compute_anomaly_and_to_obs(config)


if __name__ == "__main__":
    main()
