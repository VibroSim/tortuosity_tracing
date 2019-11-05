import sys
import os

import numpy as np
from matplotlib import pyplot as pl


from ..svg_to_histogram import histogram_from_svgs,tortuosity_plots


def main(args=None):
    if args is None:
        args=sys.argv
        pass
        
    if len(args) < 4:
        print("Usage: point_spacing_calc <fcutoff_in_1/m> <specimen> <svg_file.svg> ...")
        print("Creates tortuosity histograms from the given .svg files.")
        print("point_spacing is hard wired to 0.1 micron")
        sys.exit(0)
        pass
        
    fcutoff=float(args[1])
    specimen=str(args[2])
    svg_files=args[3:]
    point_spacing=0.1e-6
    

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
     tortuosity_path_indexes) = histogram_from_svgs(svg_files,None,fcutoff,specimen,None,point_spacing)

    mean=np.mean(thlength_final)
    StDv=np.std(thlength_final)
    #print(mean,StDv)
    fig = pl.figure(str(specimen))
    pl.clf()
    #svg_to_histogram.py [49] divides spaces between clicks into 20 steps
    #So, to get distance between point clicks, we need to multiply back
    (N,B,P)=pl.hist(thlength_final*20*10**6,bins=50)
    pl.title('Point Click Spacing',fontsize=30)
    pl.xlabel('Lengths [um]',fontsize=20)
    pl.figtext(0.55,0.75,('mu={}um\nsigma={}um'.format(round(mean*20*10**6,4),round(StDv*20*10**6,4))),bbox={'facecolor':'white','alpha':0.8,'pad':10},fontsize=25)


    pl.show()
    pass
