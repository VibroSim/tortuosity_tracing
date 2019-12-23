### Code adapted from
### "process_crack.prx"
import sys
import os 
import os.path
import numpy as np
#import cycler
import itertools
import pandas as pd
import matplotlib.pyplot as pl

import pdb

from limatix.dc_value import hrefvalue as hrefv
from limatix.dc_value import numericunitsvalue as numericunitsv
from limatix.xmldoc import xmldoc

def run(_xmldoc,_element):
  
  spreadsheetoutput_href = hrefv("tortuosity_statistics.xls",_xmldoc.getcontexthref())
  spreadsheetoutput_writer = pd.ExcelWriter(spreadsheetoutput_href.getpath())

  outputfiles = _xmldoc.xpathcontext(_element,"/prx:inputfiles/prx:inputfile/prx:outputfile")#all outputfiles in .pro (they hold the .xlp hrefs)

  tortuosity_table = pd.DataFrame(columns=["Specimen","Tortuosity Standard Deviation (deg)","Mean Point-Click Spacing (um)","Standard Deviation of Point-Click Spacing (um)"])
  specimens=[]
  tortuosities=[]
  pointspacing_means=[]
  pointspacing_stddvs=[]
  spat_freqs=[]
  for outputfile in outputfiles:
    outputdoc = xmldoc.loadhref(hrefv.fromxml(_xmldoc,outputfile)) #grabs the .xlp href from the outputfile
    crack = outputdoc.xpath("/dc:experiment") #look at the experiment log of each .xlp

    specimen=outputdoc.xpathsinglecontextstr(crack,"dc:summary/dc:specimen",default="UNKNOWN")
    tortuosity = outputdoc.xpathsinglecontextfloat(crack,"dc:filtered_sigma",default=np.nan)
    pointspacing_mean=outputdoc.xpathsinglecontextfloat(crack,"dc:Mean_spacing",default=np.nan)
    pointspacing_stddv=outputdoc.xpathsinglecontextfloat(crack,"dc:StdDv_spacing",default=np.nan)
    spatial_freq=1.0/pointspacing_mean

    spat_freqs.append(spatial_freq)
    specimens.append(specimen)
    tortuosities.append(tortuosity*180.0/np.pi)
    pointspacing_means.append(pointspacing_mean)
    pointspacing_stddvs.append(pointspacing_stddv)
    tortuosity_table = tortuosity_table.append({"Specimen": specimen,
                                              "Tortuosity Standard Deviation (deg)": tortuosity*180.0/np.pi,
                                              "Mean Point-Click Spacing (um)": pointspacing_mean,
                                              "Spatial Frequency (um^-1)": spatial_freq,
                                              "Standard Deviation of Point-Click Spacing (um)": pointspacing_stddv},ignore_index=True)


    pass
  '''
  markers=['o','v','^','<','>','s','p','+','x','D']
  colors=['b','g','r','c','m','y','k']
  all_markers=itertools.cycle(markers)
  all_colors=itertools.cycle(colors)
  all_torts = zip(specimens, tortuosities, all_markers, all_colors)
  all_stddvs = zip(specimens, pointspacing_stddvs, all_markers, all_colors)
  all_avgs = zip(specimens, pointspacing_means, all_markers, all_colors)


  fig_tort, ax0 =pl.subplots()
  for spec, tort, mark, color in all_torts:
    pl.plot(tort,linestyle='',color=color, marker=mark)
    pass
  pl.legend(specimens,prop={'size': 10})
  pl.ylabel('Standard Deviations of Tortuosities [degrees]')

  fig_avg, ax0 =pl.subplots()
  for spec, avgs, mark, color in all_avgs:
    pl.plot(avgs,linestyle='',color=color, marker=mark)
    pass
  pl.legend(specimens,prop={'size': 10})
  pl.ylabel('Average Lengths [um]')

  fig_stddv, ax0 =pl.subplots()
  for spec, stddev, mark, color in all_stddvs:
    pl.plot(stddev,linestyle='',color=color, marker=mark)
    pass
  pl.legend(specimens,prop={'size': 10})
  pl.ylabel('Standard Deviataions of Lengths [um]')

  pl.show()'''
  tortuosity_table.to_excel(spreadsheetoutput_writer,sheet_name="Tortuosities",float_format="%.2f")
  
  spreadsheetoutput_writer.close()
  return { 
    "dc:summaryspreadsheet": spreadsheetoutput_href,
  }  
  print ('Done!\nThe statistics are saved in "tortuosity_statistics.xls".')
  pass
