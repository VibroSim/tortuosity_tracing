import xlrd
import sys
import os 
import os.path
import numpy as np
import itertools
import pandas as pd
import matplotlib.pyplot as pl
import pdb
from PIL import Image
from PIL import PngImagePlugin
from sklearn.linear_model import LinearRegression

from limatix import dc_value
from limatix.dc_value import hrefvalue as hrefv
from limatix.dc_value import numericunitsvalue as numericunitsv
from limatix.xmldoc import xmldoc

def run(_xmldoc,_element,
        dc_summarytable_href,
        dc_statsdir_href=hrefv("multiple_specimen_stats/",".")):
    #xls_file = ("./tortuosity_statistics.xls")
    #outputfile =_xmldoc.xpathcontext(_element,"/prx:inputfiles/dc:summaryspreadsheet")
    #xls_file = xmldoc.loadhref(hrefv.fromxml(_xmldoc,outputfile)) #grabs the .xlp href from the outputfile
    #xls_file = outputdoc.xpath("/dc:summaryspreadsheet") #look at the experiment log of each .xlp


    summarytable = pd.read_csv(dc_summarytable_href.getpath())

    specimens=np.array(summarytable["Specimen"])
    tortuosities = np.array(summarytable["Tortuosity Standard Deviation (deg)"],dtype='d')*np.pi/180.0 # tortuosity in radians
    means = np.array(summarytable["Mean Point-Click Spacing (um)"],dtype='d')/1e6 # spacing in meters
    spat_freqs = np.array(summarytable["Point Spacing Frequency (um^-1)"],dtype='d')*1e6 # Point spacing spatial frequency in meters^-1
    stddvs = np.array(summarytable["Standard Deviation of Point-Click Spacing (um)"],dtype='d')/1e6 # stddvs in meters 

    coefficients=np.polyfit(spat_freqs,tortuosities,1)
    freq_tort_func=np.poly1d(coefficients)

    n_bins=8
    tort_bins=np.linspace(np.min(tortuosities),np.max(tortuosities),n_bins)
    mean_bins=np.linspace(np.min(means),np.max(means),n_bins)
    stddv_bins=np.linspace(np.min(stddvs),np.max(stddvs),n_bins)
    
    if not os.path.exists(dc_statsdir_href.getpath()):
        os.mkdir(dc_statsdir_href.getpath())
        pass
    
    tort_hist_href=hrefv("tortuosity_histogram.png",dc_statsdir_href)
    mean_hist_href=hrefv("means_histogram.png",dc_statsdir_href)
    stddv_hist_href=hrefv("stddvs_histogram.png",dc_statsdir_href)
    tort_distribution_href=hrefv("tortuosity_distribution.png",dc_statsdir_href)
    mean_distribution_href=hrefv("means_distribution.png",dc_statsdir_href)
    stddv_distribution_href=hrefv("stddvs_distribution.png",dc_statsdir_href)
    inv_mean_vs_tort_href=hrefv("inv_means_vs_tortuosity.png",dc_statsdir_href)
    stddv_over_mean_href=hrefv("stddvs_over_means_distribution.png",dc_statsdir_href)
    
    ### Build the Histograms
    pl.figure()
    pl.clf()
    pl.title('Tortuosities Histogram')
    pl.hist(tortuosities*180.0/np.pi,bins=tort_bins*180.0/np.pi)
    pl.xlabel('Tortuosity [degrees]')
    pl.ylabel('Number of Specimens [unitless]')
    pl.savefig(tort_hist_href.getpath(),transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Average Point Spacing Histogram')
    pl.hist(means*1e6,bins=mean_bins*1e6)
    pl.xlabel('Average Point Spacing [um]')
    pl.ylabel('Number of Specimens [unitless]')
    pl.savefig(mean_hist_href.getpath(),transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Standard Deviation of Point Spacings Histogram')
    pl.hist(stddvs*1e6,bins=stddv_bins*1e6)
    pl.xlabel('Standard Deviations [um]')
    pl.ylabel('Number of Specimens [unitless]')
    pl.savefig(stddv_hist_href.getpath(),transparent=False)

    ### Build the distributions
    pl.figure()
    pl.clf()
    pl.title('Tortuosity Values',fontsize=20)
    pl.xticks(range(len(means)),specimens,rotation=90)
    pl.plot(range(len(tortuosities)) ,tortuosities*180.0/np.pi,'*')
    pl.ylabel('Tortuosity [degrees]',fontsize=15)
    pl.grid()
    pl.tight_layout()
    pl.savefig(tort_distribution_href.getpath(),transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Average Point Spacing Distribution (multiple specimens)')
    pl.xticks(range(len(means)),specimens,rotation=90)
    pl.errorbar(range(len(means)),means*1e6,fmt='o',yerr=stddvs*1e6)
    pl.ylabel('Average Point spacing [um]')
    pl.grid()
    pl.tight_layout()
    pl.savefig(mean_distribution_href.getpath(),transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Standard Deviation of Point Spacings Distribution')
    pl.xticks(range(len(means)),specimens,rotation=90)
    pl.plot(range(len(stddvs)),stddvs*1e6,'D')
    pl.ylabel('Point Spacing Standard Deviation [um]')
    pl.grid()
    pl.tight_layout()
    pl.savefig(stddv_distribution_href.getpath(),transparent=False)

    ### Other Figures
    pl.figure()
    pl.clf()
    pl.title('Point Spacing Frequency vs Tortuosity')
    pl.plot(spat_freqs/1e6,tortuosities*180.0/np.pi,'>')
    pl.plot(spat_freqs/1e6,freq_tort_func(spat_freqs)*180.0/np.pi)
    pl.legend(('True Values','Linear Regression: \ny=%.2fx+%.2f rad for x in m^-1' % (coefficients[0],coefficients[1])),prop={'size': 10},loc='best')
    pl.xlabel('Spatial Frequency [um^-1]')
    pl.ylabel('Tortuosity [degrees]')
    pl.tight_layout()
    pl.savefig(inv_mean_vs_tort_href.getpath(),transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Standard Deviation/Mean')
    pl.xticks(range(len(means)),specimens,rotation=90)
    pl.plot(range(len(stddvs)),np.array(stddvs)/np.array(means),'s')
    pl.ylabel('Standard Deviation / Mean [unitless]')
    pl.grid()
    pl.tight_layout()
    pl.savefig(stddv_over_mean_href.getpath(),transparent=False)
    
    #### Adapted from: https://stackoverflow.com/questions/10532614/can-matplotlib-add-metadata-to-saved-figures
    #METADATA = {"specimen,tortuosity[deg],mean[um],spatial frequency[um^-1],stddv[um]":str(all_data)}
    
    #tort_hist_image = Image.open(tort_hist_href)
    #mean_hist_image = Image.open(mean_hist_href)
    #stddv_hist_image = Image.open(stddv_hist_href)
    #tort_dist_image = Image.open(tort_distribution_href)
    #mean_dist_image = Image.open(mean_distribution_href)
    #stddv_dist_image = Image.open(stddv_distribution_href)
    #inv_mean_vs_tort_image = Image.open(inv_mean_vs_tort_href)
    #stddv_over_mean_image = Image.open(stddv_over_mean_href)
    
    #meta = PngImagePlugin.PngInfo()
    
    #for x in METADATA:
    #    meta.add_text(x, METADATA[x])
    #tort_hist_image.save(tort_hist_href, "png", pnginfo=meta)
    #mean_hist_image.save(mean_hist_href, "png", pnginfo=meta)
    #stddv_hist_image.save(stddv_hist_href, "png", pnginfo=meta)
    #tort_dist_image.save(tort_distribution_href, "png", pnginfo=meta)
    #mean_dist_image.save(mean_distribution_href, "png", pnginfo=meta)
    #stddv_dist_image.save(stddv_distribution_href, "png", pnginfo=meta)
    #inv_mean_vs_tort_image.save(inv_mean_vs_tort_href, "png", pnginfo=meta)
    #stddv_over_mean_image.save(stddv_over_mean_href, "png", pnginfo=meta)
    
    ###to read .png file:
    ###im2 = Image.open(PNG FILE)
    ###print im2.info

    #pl.show()
    return {
        "dc:tortuosity_vs_pointspacingfreq_slope":coefficients[0],
        "dc:tortuosity_vs_pointspacingfreq_offset":coefficients[1],
        "dc:tortuosity_histogram": tort_hist_href,
        "dc:mean_histogram": mean_hist_href,
        "dc:stddv_histogram": stddv_hist_href,
        "dc:tort_distribution": tort_distribution_href,
        "dc:mean_distribution": mean_distribution_href,
        "dc:inv_mean_vs_tort_plot": inv_mean_vs_tort_href,
        "dc:stddv_over_mean_plot": stddv_over_mean_href,
    }
    
