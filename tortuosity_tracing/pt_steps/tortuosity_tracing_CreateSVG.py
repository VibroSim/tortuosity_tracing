import posixpath

import tortuosity_tracing

from limatix import dc_value
from limatix.dc_value import hrefvalue as hrefv
from limatix.dc_value import numericunitsvalue as numericunitsv
from limatix.xmldoc import xmldoc

def run(_xmldoc,_element,dc_tortuositytiff_href,dc_XRes_numericunits,dc_YRes_numericunits):
    
    tiff_filename=dc_tortuositytiff_href.getpath()
    xResolution=dc_XRes_numericunits.value("m/px")
    yResolution=dc_YRes_numericunits.value("m/px")

    #(tiff_basename,tiff_ext)=os.path.splitext(tiff_filename)
    #assert(tiff_ext.lower()=='.tiff' or tiff_ext.lower()=='.tif')

    tiff_context = dc_tortuositytiff_href.leafless()

    tiff_fileurlname = dc_tortuositytiff_href.get_bare_quoted_filename()
    (tiff_baseurlname,tiff_urlext)=posixpath.splitext(tiff_fileurlname)
    assert(tiff_urlext.lower()=='.tiff' or tiff_urlext.lower()=='.tif')

    svg_fileurlname = tiff_baseurlname+".svg"
    
    svg_href = dc_value.hrefvalue(svg_fileurlname,contexthref=tiff_context)
    
    svg_filename=svg_href.getpath()

    
    
    tortuosity_tracing.convert_tiff_to_svg(tiff_filename,svg_filename,xResolution,yResolution)

    return { "dc:tortuositysvg": svg_href }
