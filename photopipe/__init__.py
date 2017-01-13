"""
PhotoPipe - reduction and photometry pipeline developed at GSFC and UMD
"""
__all__ = ['__version__', 'version_info']

# Package version
import pkgutil
__version__ = pkgutil.get_data(__package__, 'VERSION').decode('ascii').strip()
version_info = tuple(int(v) if v.isdigit() else v
                     for v in __version__.split('.'))
del pkgutil

# Check minimum required Python version
import sys
if sys.version_info < (2, 7):
    print("PhotoPipe %s requires Python 2.7" % __version__)
    sys.exit(1)