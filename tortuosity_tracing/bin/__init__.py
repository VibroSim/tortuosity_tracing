import sys
import os 
import os.path

from ..svg_to_histogram import histogram_from_svgs,tortuosity_plots
from ..tiff2svg import convert_tiff_to_svg

try:
    # py2.x
    from urllib import pathname2url
    pass
except ImportError:
    # py3.x
    from urllib.request import pathname2url
    pass


class dummy(object):
    pass

def getstepurlpath():
    mypath = sys.modules[dummy.__module__].__file__
    mydir=os.path.split(mypath)[0]

    return [ pathname2url(os.path.join(mydir,"pt_steps")) ]
