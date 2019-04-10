import tortuosity_tracing
import subprocess

import hashlib

from limatix import dc_value


def run(_xmldoc,_element,fcutoff_numericunits):
    
    meastags = _xmldoc.xpath("dc:measurement[dc:traced_svg]")
    svg_filenames=[]

    fcutoff=fcutoff_numericunits.value(units="m^-1")

    for meastag in meastags:

        svgfile = _xmldoc.xpathsinglecontext(meastag,"dc:tortuositysvg")
        svg_href = dc_value.hrefvalue.fromxml(_xmldoc,svgfile)
        
        svg_filename = svg_href.getpath()
        
        svg_filenames.append(svg_filename)
        pass

    dest_el = _xmldoc.xpathsingle("dc:summary/dc:dest")
    dest_href = dc_value.hrefvalue.fromxml(_xmldoc,dest_el)

    (theta_final,
     thlength_final,
     filtered_final,
     eq_lengths_final,
     avg_mu,
     avg_filtered_mu,
     avg_sigma,
     avg_filtered_sigma) = tortuosity_tracing.histogram_from_svgs(svg_filenames,fcutoff)
    
    (unfiltered_filename,filtered_filename) = tortuosity_tracing.tortuosity_plots(
        theta_final,
        thlength_final,
        filtered_final,
        eq_lengths_final,
        avg_mu,
        avg_filtered_mu,
        avg_sigma,
        avg_filtered_sigma,savedir=dest_href.getpath())
    
    unfiltered_href = dc_value.hrefvalue(unfiltered_filename,contexthref=dest_href)
    filtered_href = dc_value.hrefvalue(filtered_filename,contexthref=dest_href)
    return { "dc_avg_filtered_mu": dc_value.numericunitsvalue(avg_filtered_mu,"radians"),
             "dc_avg_filtered_sigma": dc_value.numericunitsvalue(avg_filtered_sigma,"radians"),
             "dc_unfiltered_plot": unfiltered_href,
             "dc_filtered_plot": filtered_href
             }
