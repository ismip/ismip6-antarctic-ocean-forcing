# make sure the 'Agg' backend is used before any local module get loaded
import matplotlib
matplotlib.use('Agg')

__version_info__ = (1, 1)
__version__ = '.'.join(str(vi) for vi in __version_info__)
