#! /usr/bin/env python

import sys
import os
import datetime

from limatix import xmldoc
from limatix import dc_value

try: 
    # py2.x
    from urllib import pathname2url
    from urllib import url2pathname
    import urlparse
    pass
except ImportError:
    import urllib.parse as urlparse # python3 
    from urllib.request import pathname2url
    from urllib.request import url2pathname
    pass

def main(args=None):
    if args is None:
        args=sys.argv
        pass

    if len(args) < 8 or os.path.splitext(args[4])[1] != '.xlg':
        print("Usage: tiffs_to_xlg <specimen> <date> <perfby> <xlg_file.xlg> <XResolution in m/px ("None" if you don't know)> <YResolution in m/px ("None" if you don't know)> <tiff_files.tif> ...")
        print(" Creates an experiment log from the given .tif files")
        sys.exit(0)
        pass
        
    specimen=args[1]
    date=args[2]
    perfby=args[3]
    xlg_file=args[4]
    xResolution=args[5]
    yResolution=args[6]
    tiff_files = args[7:]
    
    date_parsed=datetime.datetime.strptime(date,"%Y-%m-%d")  # Exception here if you don't give the date in a reasonable format
    
    if os.path.exists(xlg_file):
        raise ValueError("Experiment log file %s already exists; will not overwrite" % (xlg_file))
        
        
        
    xlg = xmldoc.xmldoc.newdoc("dc:experiment",nsmap={"dc": "http://limatix.org/datacollect",
                                                      "dcvalue": "http://limatix.org/dcvalue",
                                                      "lip":"http://limatix.org/provenance",
                                                      "xlink":"http://www.w3.org/1999/xlink"})

    # Give the experiment log a location so we can use relative URL's
    xlg.set_href(dc_value.hrefvalue(pathname2url(xlg_file),contexthref="."))

    summary = xlg.addelement(xlg.getroot(),"dc:summary")
    xlg.addsimpleelement(summary,"dc:specimen",(specimen,))
    xlg.addsimpleelement(summary,"dc:perfby",(perfby,))
    xlg.addsimpleelement(summary,"dc:date",(date,))
    xlg.addsimpleelement(summary,"dc:XRes",(xResolution,))
    xlg.addsimpleelement(summary,"dc:YRes",(yResolution,))
    dest_el=xlg.addelement(summary,"dc:dest")
    # dest_el will be the directory containing the first .tiff

    firsttiff_href = dc_value.hrefvalue(pathname2url(tiff_files[0]),contexthref=".")
    dest_href = firsttiff_href.leafless() # remove file part of url
    dest_href.xmlrepr(xlg,dest_el)
    

    
    measnum=1
    for tiff_file in tiff_files:
        meas_el=xlg.addelement(xlg.getroot(),"dc:measurement")
        measnum_el=xlg.addsimpleelement(meas_el,"dc:measnum",(str(measnum),))
        tiff_el = xlg.addelement(meas_el,"dc:tortuositytiff")
        tiff_href = dc_value.hrefvalue(pathname2url(tiff_file),contexthref=".")
        tiff_href.xmlrepr(xlg,tiff_el)
        measnum+=1
        pass
        
    xlg.close()

    print("Successfully generated experiment log %s" % (xlg_file))
    print(" ")
    print("Now you should probably create a limatix-git repository for the data")
    print(" e.g.  \"limatix-git init\"")
    print("       \"git commit\"")
    print("       \"limatix-git add -a\"")
    print("       \"git commit\"")


    pass
