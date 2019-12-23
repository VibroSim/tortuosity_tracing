import tortuosity_tracing
import subprocess
import numpy as np
from matplotlib import pyplot as pl

import hashlib

from limatix import dc_value


def run(_xmldoc,_element,fcutoff_numericunits,dc_specimen_str=None,point_spacing_numericunits=dc_value.numericunitsvalue(1.0e-7,"m")):

    if dc_specimen_str is None:
        dc_specimen = _xmldoc.xpathsingle("dc:summary/dc:specimen")
        dc_specimen_str = _xmldoc.xpathsinglecontextstr(dc_specimen,"string(.)")
        pass


    
    #meastags = _xmldoc.xpath("dc:measurement[dc:traced_svg]")
    meastags = _xmldoc.xpath("dc:measurement")
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
     tortuosity_path_filenames,
     tortuosity_plot_filenames,
     tortuosity_path_indexes,
     num_steps) = tortuosity_tracing.histogram_from_svgs(svg_filenames,svg_measnums,fcutoff,dc_specimen_str,dest_href.getpath(),point_spacing_numericunits.value(units="m"))
    


    #assert(len(tortuosity_path_filenames)==len(svg_filenames))
    summed_clicked_length=np.sum(thlength_final)
    summed_eq_length=np.sum(eq_lengths_final)

    (unfiltered_filename,filtered_filename) = tortuosity_tracing.tortuosity_plots(
	dc_specimen_str,
        theta_final,
        thlength_final,
        filtered_theta_final,
        eq_lengths_final,
        mu,
        mu_F,
        sigma,
        sigma_F,dest_href.getpath())
    
    mean=np.mean(thlength_final)*num_steps
    StdDv=np.std(thlength_final)*num_steps
    fig = pl.figure(str(svg_filenames[0]))
    pl.clf()
    #svg_to_histogram.py [49] divides spaces between clicks into 20 steps
    #So, to get distance between point clicks, we need to multiply back
    (N,B,P)=pl.hist(thlength_final*num_steps*10**6,bins=50)
    pl.title('Point Click Spacing',fontsize=30)
    pl.xlabel('Lengths [um]',fontsize=20)
    pl.figtext(0.55,0.75,('mu={}um\nsigma={}um'.format(round(mean*1e6,4),round(StdDv*1e6,4))),bbox={'facecolor':'white','alpha':0.8,'pad':10},fontsize=25)
    unfiltered_href = dc_value.hrefvalue(unfiltered_filename,contexthref=dest_href)
    filtered_href = dc_value.hrefvalue(filtered_filename,contexthref=dest_href)
    retval = [ ("dc:filtered_mu", dc_value.numericunitsvalue(mu_F,"radians")),
               ("dc:filtered_sigma", dc_value.numericunitsvalue(sigma_F,"radians")),
               ("dc:unfiltered_plot", unfiltered_href),
               ("dc:filtered_plot", filtered_href),
               ("dc:Full_clicked_length",dc_value.numericunitsvalue(summed_clicked_length,"meters")),
               ("dc:Full_eq_length",dc_value.numericunitsvalue(summed_eq_length,"meters")),
               ("dc:Mean_spacing",dc_value.numericunitsvalue(mean,"m")),
               ("dc:StdDv_spacing",dc_value.numericunitsvalue(StdDv,"m"))
           ]
    
    print("The length of the crack path based on the point clicks is {} m".format(summed_clicked_length))
    print("The length of the crack path based on the evenly spaced points is {} m".format(summed_eq_length))
    print("Check that these to numbers are consistent with each other and greater than the crack length in the specimen database.")
    for filecnt  in range(len(tortuosity_plot_filenames)):
        tortuosity_plot_filename = tortuosity_plot_filenames[filecnt]
        retval.append((("dc:path_comparison", { "measnum" : str(svg_measnums[filecnt]) }), dc_value.hrefvalue(tortuosity_plot_filename,contexthref=dest_href)))
        
        pass
    return retval
