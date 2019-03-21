import tortuosity_tracing
import subprocess

import hashlib

from limatix import dc_value

def eval_md5(filename):
    hasher=hashlib.md5()
    with open(filename,"rb") as svgdata:
        buf=svgdata.read()
        hasher.update(buf)
        pass
    return hasher.hexdigest()
    

def run(_xmldoc,_element,dc_tortuositysvg_hrefvalue):
    
    svg_filename=dc_tortuositysvg_hrefvalue.getpath()
    
    initialhash = eval_md5(svg_filename)
    status = subprocess.call([ "inkscape",svg_filename ] )

    finalhash = eval_md5(svg_filename)

    if (status) :
        sys.stderr.write("Warning: Error return from inkscape on file \"%s\"\n" % (svg_filename))
        pass

    if finalhash == initialhash:
        sys.stderr.write("Warning: Did not modify tortuosity trace in file \"%s\"\n" % (svg_filename))
        pass

    return 
