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

def run(_xmldoc,_element):
    xls_file = ("./tortuosity_statistics.xls")
    #outputfile =_xmldoc.xpathcontext(_element,"/prx:inputfiles/dc:summaryspreadsheet")
    #xls_file = xmldoc.loadhref(hrefv.fromxml(_xmldoc,outputfile)) #grabs the .xlp href from the outputfile
    #xls_file = outputdoc.xpath("/dc:summaryspreadsheet") #look at the experiment log of each .xlp

    worksheet = xlrd.open_workbook(xls_file)
    sheet = worksheet.sheet_by_index(0)
    specimens=[]
    tortuosities=[]
    means=[]
    stddvs=[]
    spat_freqs=[]

    for i in range(sheet.nrows):
        if sheet.cell_value(i, 2) == "":
            pass
        else:
            specimens.append(sheet.cell_value(i, 1))
            tortuosities.append(sheet.cell_value(i, 2))
            means.append(sheet.cell_value(i, 3))
            spat_freqs.append(sheet.cell_value(i, 5))
            stddvs.append(sheet.cell_value(i, 4))
            pass
        pass
    
    spat_freqs=np.array(spat_freqs[1:])
    specimens=specimens[1:]
    tortuosities=np.array(tortuosities[1:])
    means=np.array(means[1:])
    stddvs=np.array(stddvs[1:])
    all_data=zip(specimens,tortuosities,means,spat_freqs,stddvs)
    coefficients=np.polyfit(spat_freqs,tortuosities,1)
    freq_tort_func=np.poly1d(coefficients)

    n_bins=8
    tort_bins=np.linspace(np.min(tortuosities),np.max(tortuosities),n_bins)
    mean_bins=np.linspace(np.min(means),np.max(means),n_bins)
    stddv_bins=np.linspace(np.min(stddvs),np.max(stddvs),n_bins)
    
    tort_hist_href="./multiple_specimen_stats/tortuosity_histogram.png"
    mean_hist_href="./multiple_specimen_stats/means_histogram.png"
    stddv_hist_href="./multiple_specimen_stats/stddvs_histogram.png"
    tort_distribution_href="./multiple_specimen_stats/tortuosity_distribution.png"
    mean_distribution_href="./multiple_specimen_stats/means_distribution.png"
    stddv_distribution_href="./multiple_specimen_stats/stddvs_distribution.png"
    inv_mean_vs_tort_href="./multiple_specimen_stats/inv_means_vs_tortuosity.png"
    stddv_over_mean_href="./multiple_specimen_stats/stddvs_over_means_distribution.png"
    
    ### Build the Histograms
    pl.figure()
    pl.clf()
    pl.title('Tortuosities Histogram')
    pl.hist(tortuosities,bins=tort_bins)
    pl.xlabel('Tortuosity [degrees]')
    pl.ylabel('Number of Specimens [unitless]')
    pl.savefig(tort_hist_href,transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Average Point Spacing Histogram')
    pl.hist(means,bins=mean_bins)
    pl.xlabel('Average Point Spacing [um]')
    pl.ylabel('Number of Specimens [unitless]')
    pl.savefig(mean_hist_href,transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Standard Deviation of Point Spacings Histogram')
    pl.hist(stddvs,bins=stddv_bins)
    pl.xlabel('Standard Deviations [um]')
    pl.ylabel('Number of Specimens [unitless]')
    pl.savefig(stddv_hist_href,transparent=False)

    ### Build the distributions
    pl.figure()
    pl.clf()
    pl.title('Tortuosities Distribution')
    pl.xticks(range(len(means)),specimens,rotation=90)
    pl.plot(range(len(tortuosities)) ,tortuosities,'*')
    pl.ylabel('Tortuosity [degrees]')
    pl.grid()
    pl.savefig(tort_distribution_href,transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Average Point Spacing Distribution (multiple specimens)')
    pl.xticks(range(len(means)),specimens,rotation=90)
    pl.errorbar(range(len(means)),means,fmt='o',yerr=stddvs)
    pl.ylabel('Average Point spacing [um]')
    pl.grid()
    pl.savefig(mean_distribution_href,transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Standard Deviation of Point Spacings Distribution')
    pl.xticks(range(len(means)),specimens,rotation=90)
    pl.plot(range(len(stddvs)),stddvs,'D')
    pl.ylabel('Point Spacing Standard Deviation [um]')
    pl.grid()
    pl.savefig(stddv_distribution_href,transparent=False)

    ### Other Figures
    pl.figure()
    pl.clf()
    pl.title('Spacial Frequency vs Tortuosity')
    pl.plot(spat_freqs,tortuosities,'>')
    pl.plot(spat_freqs,freq_tort_func(spat_freqs))
    pl.legend(('True Values','Linear Regression: \ny=%.2fx+%.2f' % (coefficients[0],coefficients[1])),prop={'size': 10},loc='best')
    pl.xlabel('Spatial Frequency [um^-1]')
    pl.ylabel('Tortuosity [degrees]')
    pl.savefig(inv_mean_vs_tort_href,transparent=False)
    
    pl.figure()
    pl.clf()
    pl.title('Standard Deviation/Mean')
    pl.xticks(range(len(means)),specimens,rotation=90)
    pl.plot(range(len(stddvs)),np.array(stddvs)/np.array(means),'s')
    pl.ylabel('Standard Deviation / Mean [unitless]')
    pl.grid()
    pl.savefig(stddv_over_mean_href,transparent=False)
    
    ### Adapted from: https://stackoverflow.com/questions/10532614/can-matplotlib-add-metadata-to-saved-figures
    METADATA = {"specimen,tortuosity[deg],mean[um],spatial frequency[um^-1],stddv[um]":str(all_data)}
    
    tort_hist_image = Image.open(tort_hist_href)
    mean_hist_image = Image.open(mean_hist_href)
    stddv_hist_image = Image.open(stddv_hist_href)
    tort_dist_image = Image.open(tort_distribution_href)
    mean_dist_image = Image.open(mean_distribution_href)
    stddv_dist_image = Image.open(stddv_distribution_href)
    inv_mean_vs_tort_image = Image.open(inv_mean_vs_tort_href)
    stddv_over_mean_image = Image.open(stddv_over_mean_href)
    
    meta = PngImagePlugin.PngInfo()
    
    for x in METADATA:
        meta.add_text(x, METADATA[x])
    tort_hist_image.save(tort_hist_href, "png", pnginfo=meta)
    mean_hist_image.save(mean_hist_href, "png", pnginfo=meta)
    stddv_hist_image.save(stddv_hist_href, "png", pnginfo=meta)
    tort_dist_image.save(tort_distribution_href, "png", pnginfo=meta)
    mean_dist_image.save(mean_distribution_href, "png", pnginfo=meta)
    stddv_dist_image.save(stddv_distribution_href, "png", pnginfo=meta)
    inv_mean_vs_tort_image.save(inv_mean_vs_tort_href, "png", pnginfo=meta)
    stddv_over_mean_image.save(stddv_over_mean_href, "png", pnginfo=meta)
    
    ###to read .png file:
    ###im2 = Image.open(PNG FILE)
    ###print im2.info

    pl.show()
    return {"dc:Coefficient_m":coefficients[0],"dc:Coefficient_b":coefficients[1]}
    pass
