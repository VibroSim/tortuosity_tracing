import tortuosity_tracing
import subprocess
import numpy as np
from matplotlib import pyplot as pl
import hashlib
import csv

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
     tortuosity_path_indexes) = tortuosity_tracing.histogram_from_svgs(svg_filenames,svg_measnums,fcutoff,dc_specimen_str,dest_href.getpath(),point_spacing_numericunits.value(units="m"))

    mean=np.mean(thlength_final)*20*10**6
    StDv=np.std(thlength_final)*20*10**6
    #print(mean,StDv)
    fig = pl.figure(str(dc_specimen))
    pl.clf()
    #svg_to_histogram.py [line 49] divides spaces between clicks into 20 steps
    #So, to get distance between point clicks, we need to multiply back
    (N,B,P)=pl.hist(thlength_final*20*10**6,bins=50)
    pl.title('Point Click Spacing',fontsize=30)
    pl.xlabel('Lengths [um]',fontsize=20)
    pl.figtext(0.55,0.75,('mu={}um\nsigma={}um'.format(round(mean*10**6,4),round(StDv*10**6,4))),bbox={'facecolor':'white','alpha':0.8,'pad':10},fontsize=25)
    retval = [ ("dc:Mean_spacing",dc_value.numericunitsvalue(mean,"microns")),
               ("dc:StDv_spacing",dc_value.numericunitsvalue(StDv,"microns"))
           ]
    for filecnt  in range(len(tortuosity_plot_filenames)):
        tortuosity_plot_filename = tortuosity_plot_filenames[filecnt]
        retval.append((("dc:path_comparison", { "measnum" : str(svg_measnums[filecnt]) }), dc_value.hrefvalue(tortuosity_plot_filename,contexthref=dest_href)))
        pass
    pl.show()
#    with open('./tortuosity_processed_output/tortuosity_statistics.csv', mode='w') as stats_file:
#        stats_writer = csv.writer(stats_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#
#        stats_writer.writerow([dc_specimen, mean, StDv])
#        pass
    return retval
    pass
