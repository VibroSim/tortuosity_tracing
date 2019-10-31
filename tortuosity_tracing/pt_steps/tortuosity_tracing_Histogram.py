import tortuosity_tracing
import subprocess

import hashlib

from limatix import dc_value


def run(_xmldoc,_element,fcutoff_numericunits,point_spacing_numericunits=dc_value.numericunitsvalue(1.0e-7,"m")):
    
    meastags = _xmldoc.xpath("dc:measurement[dc:traced_svg]")
    svg_measnums = []
    svg_filenames=[]

    fcutoff=fcutoff_numericunits.value(units="m^-1")

    for meastag in meastags:

        svgfile = _xmldoc.xpathsinglecontext(meastag,"dc:tortuositysvg")
        svg_href = dc_value.hrefvalue.fromxml(_xmldoc,svgfile)
        
        svg_filename = svg_href.getpath()

        svg_measnum = int(_xmldoc.xpathsinglecontext(meastag,"string(dc:measnum)"))

        svg_filenames.append(svg_filename)
        svg_measnums.append(svg_measnum)

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
     sigma_F,
     tortuosity_path_filenames) = tortuosity_tracing.histogram_from_svgs(svg_filenames,svg_measnums,fcutoff,dest_href.getpath(),point_spacing_numericunits.value(units="m"))
    
    assert(len(tortuosity_path_filenames)==len(svg_filenames))
    summed_clicked_length=np.sum(thlength_final)
    summed_eq_length=np.sum(eq_lengths_final)

    (unfiltered_filename,filtered_filename) = tortuosity_tracing.tortuosity_plots(
	tortuosity_path_filenames,
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
    retval = [ ("dc:filtered_mu", dc_value.numericunitsvalue(mu_F,"radians")),
             ("dc:filtered_sigma", dc_value.numericunitsvalue(sigma_F,"radians")),
             ("dc:unfiltered_plot", unfiltered_href),
             ("dc:filtered_plot", filtered_href),
             ("dc:Full_clicked_length",dc_value.numericunitsvalue(summed_clicked_length,"meters"),
             ("dc:Full_eq_length",dc_value.numericunitsvalue(summed_eq_length,"meters")
        ]

    print("The length of the crack path based on the point clicks is {} m".format(summed_clicked_length))
    print("The length of the crack path based on the evenly spaced points is {} m".format(summed_eq_length))
    print("Check that these to numbers are consitent with each other and greater than the crack length in the specimen database.")
    for filecnt  in range(len(tortuosity_path_filenames)):
        tortuosity_path_filename = tortuosity_path_filenames[filecnt]
        retval.append((("dc:path_comparison", { "measnum" : str(svg_measnums[filecnt]) }), dc_value.hrefvalue(tortuosity_path_filename,contexthref=dest_href)))
        
        pass
    return retval
