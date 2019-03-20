import tortuosity_tracing

from limatix import dc_value

def run(_xmldoc,_element,dc_tortuositytiff_hrefvalue):
    
    tiff_filename=dc_tortuositytiff_hrefvalue.getpath()


    #(tiff_basename,tiff_ext)=os.path.splitext(tiff_filename)
    #assert(tiff_ext.lower()=='.tiff' or tiff_ext.lower()=='.tif')

    tiff_context = dc_tortuositytiff_hrefvalue.leafless()

    tiff_fileurlname = dc_tortuositytiff.get_bare_quoted_filename()
    (tiff_baseurlname,tiff_urlext)=posixpath.splitext(tiff_fileurlname)
    assert(tiff_urlext.lower()=='.tiff' or tiff_urlext.lower()=='.tif')

    svg_fileurlname = tiff_baseurlname+".svg"
    
    svg_href = dc_value.hrefvalue(svg_fileurlname,hrefcontext=tiff_context)
    
    svg_filename=svg_href.getpath()

    
    
    tortuosity_tracing.tiff2svg(tiff_filename,svg_filename)

    return { "dc_tortuositysvg": svg_href }
