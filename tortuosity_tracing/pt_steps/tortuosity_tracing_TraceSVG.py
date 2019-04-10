import sys

import limatix.timestamp

import tortuosity_tracing
import subprocess

import hashlib as hl

from limatix import dc_value

def eval_md5(filename):
    hasher=hl.md5()
    with open(filename,"rb") as svgdata:
        buf=svgdata.read()
        hasher.update(buf)
        pass
    return hasher.hexdigest()
    

def run(_xmldoc,_element,dc_tortuositysvg_href):
    
    svg_filename=dc_tortuositysvg_href.getpath()

    initialhash = eval_md5(svg_filename)
    status = subprocess.call([ "inkscape",svg_filename ] )

    finalhash = eval_md5(svg_filename)

    retval = {}

    if (status) :
        sys.stderr.write("Warning: Error return from inkscape on file \"%s\"\n" % (svg_filename))
        pass

    if finalhash == initialhash:
        sys.stderr.write("Warning: Did not modify tortuosity trace in file \"%s\"\n" % (svg_filename))
        pass

    if finalhash != initialhash:
       retval["dc:traced_svg"]=limatix.timestamp.now().isoformat()
       pass

    return retval
