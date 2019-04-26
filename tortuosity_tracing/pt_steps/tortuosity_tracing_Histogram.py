import tortuosity_tracing
import subprocess

import hashlib

from limatix import dc_value


def run(_xmldoc,_element,fcutoff_numericunits,point_spacing_numericunits=dc_value.numericunitsvalue(1.0e-7,"m")):
    
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
     filtered_theta_final,
     eq_lengths_final,
     mu,
     mu_F,
     sigma,
     sigma_F) = tortuosity_tracing.histogram_from_svgs(svg_filenames,fcutoff,dest_href.getpath(),point_spacing_numericunits.value(units="m"))
    
    (unfiltered_filename,filtered_filename) = tortuosity_tracing.tortuosity_plots(
        theta_final,
        thlength_final,
        filtered_theta_final,
        eq_lengths_final,
        mu,
        mu_F,
        sigma,
        sigma_F,dest_href.getpath())
    
    unfiltered_href = dc_value.hrefvalue(unfiltered_filename,contexthref=dest_href)
    filtered_href = dc_value.hrefvalue(filtered_filename,contexthref=dest_href)
    return { "dc_filtered_mu": dc_value.numericunitsvalue(mu_F,"radians"),
             "dc_filtered_sigma": dc_value.numericunitsvalue(sigma_F,"radians"),
             "dc_unfiltered_plot": unfiltered_href,
             "dc_filtered_plot": filtered_href
             }
