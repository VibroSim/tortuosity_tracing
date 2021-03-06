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
        print("Usage: svg_histograms <fcutoff_in_1/m> <frampwidth_in_1/m> <specimen> <svg_file.svg> ...")
        print("Creates tortuosity histograms from the given .svg files.")
        print("point_spacing is hard wired to 0.1 micron")
        sys.exit(0)
        pass
        
    fcutoff=float(args[1])
    frampwidth=float(args[2])
    specimen=args[3]
    svg_files=args[4:]
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
     tortuosity_path_indexes,num_steps) = histogram_from_svgs(svg_files,None,fcutoff,frampwidth,specimen,None,point_spacing)
    
    (unfiltered_filename,filtered_filename) = tortuosity_plots(
        specimen,
        theta_final,
        thlength_final,
        filtered_theta_final,
        eq_lengths_final,
        mu,
        mu_F,
        sigma,
        sigma_F,None)

    pl.show()
    pass


    
