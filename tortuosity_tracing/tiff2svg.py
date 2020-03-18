import os
import sys
from lxml import etree 
try:
    import ConfigParser
    pass
except ImportError:
    #python 3
    import configparser as ConfigParser
    pass
     
try: 
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from PIL import Image
from PIL.TiffTags import TAGS

try:
    from urllib import pathname2url
    pass
except ImportError:
    #Python 3 
    from urllib.request import pathname2url
    pass


def convert_tiff_to_svg(input_filename,output_filename,xResolution,yResolution):
    """The objective of this function is to read a given tif image and from it, create a small svg file.
    The svg file will have a locked image that conserves distances."""
    #    The conversion between physical distances and distances in the .svg
    #    is:  1 svg mm is equivalent to 1 physical micron. 
    #    This routine gets the scaling factor from the TIFF metadata stored 
    #in tag #34682 of TIFF files saved by an FEI SEM microscope, 

    #"""
    config=ConfigParser.ConfigParser()

    img=Image.open(input_filename)
    #(dirpath,filepart)=os.path.split(filename)
    #(filepartbase,filepartext)=os.path.splitext(filepart)
    if img.tag.keys()[-2]==34682:
        # For the FEI SEM microscope, TIFF Tag ID # 34682 contains 
        # the scale calibration as part of a .ini file structure, read with 
        # ConfigParser. 

        inidata=img.tag[34682][0] #print(inidata) shows the tags in readable way
        #inidata.get
        #use img.tag.keys() to view the dictionary
        #34682 corresponds to the entry from the microscope
        inidatafh=StringIO(inidata)
        config.readfp(inidatafh)
        sections=config.sections()
        ''' These are the elements in 'inidata' (Tag #34682)
[u'User',
 u'System',
 u'Beam',
 u'EBeam',
 u'GIS',
 u'Scan',
 u'EScan',
 u'Stage',
 u'Image',
 u'Vacuum',
 u'Specimen',
 u'Detectors',
 u'ETD',
 u'Accessories',
 u'EBeamDeceleration',
 u'PrivateFei',
 u'HiResIllumination']
    '''
        #Let's take out the useful information
        #DatabarHeight=float(config.get("PrivateFei","DatabarHeight"))
        ResolutionX=float(config.get("Image","ResolutionX"))
        ResolutionY=float(config.get("Image","ResolutionY"))
        PixelWidth=float(config.get("Scan","PixelWidth"))
        PixelHeight=float(config.get("Scan","PixelHeight"))

        assert PixelWidth==PixelHeight,'Pixels are not square. Check file.'  # Only support square pixels for now
        tif_width=PixelWidth*img.width*ResolutionX # physical width in mm of image as read from .tiff
        tif_height=PixelHeight*img.height*ResolutionY # physical height in mm of image as read from .tiff, including some non-physical height from bar at bottom of image

        # tif_width is in physical mm
        # svg_width needs to be in svg_mm
        # multply by 1 svg mm / 1 physical micron
        # multiply also by 1000 physical microns / 1 physical mmm
        svg_width= tif_width*1000.0
        svg_height=tif_height*1000.0
        pass
    elif xResolution== None or yResolution== None:
        print ("Resolution was not specified.")
        print ("Please Rerun 'tiffs_to_xlg' with x and y resolutions in m/px")
    else:
        inidata=img.tag
        #ResolutionX=float(inidata[282][0][0]/inidata[282][0][1])
        #ResolutionY=float(inidata[283][0][0]/inidata[283][0][1])
        ImageWidth=float(inidata[256][0])
        ImageHeight=float(inidata[257][0])
        ### ^ for reading the .TIFF files
        
        ResolutionX=xResolution#[m/px]
        ResolutionY=yResolution#[m/px]
        #width[m]=[px]*[m/px]
        tif_width=ImageWidth*ResolutionX # [m]
        tif_height=ImageHeight*ResolutionY # [m]

        ## 1 SVGum = 1.5436 um
        ## mm are the smallest metric inkscape unit
        ## so letting 1.5436 um = 1 SVGmm yields easily
        ## scaleable results
        svg_width= tif_width*1e6 #*1.5436
        svg_height=tif_height*1e6 #*1.5436
        #pdb.set_trace()
        pass
    ## Now, to build the svg file
    ## load in the svg file template:
    template=etree.XML("""<svg:svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   version="1.1"
   id="svg127"
   inkscape:version="0.92.2 (5c3e80d, 2017-08-06)">
  <metadata
     id="metadata133">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title></dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <defs
     id="defs131" />
  <sodipodi:namedview
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1"
     objecttolerance="10"
     gridtolerance="10"
     guidetolerance="10"
     inkscape:pageopacity="0"
     inkscape:pageshadow="2"
     inkscape:window-width="1327"
     inkscape:window-height="614"
     id="namedview129"
     showgrid="false"
     inkscape:zoom="0.41"
     inkscape:cx="121.7561"
     inkscape:cy="471.5"
     inkscape:window-x="0"
     inkscape:window-y="28"
     inkscape:window-maximized="1"
     inkscape:current-layer="svg127" />
  <image
     preserveAspectRatio="none"
     id="image135"
     x="0"
     y="0"
     sodipodi:insensitive="true"
/> 
</svg:svg>""")

    ## Next, Add in the calculated values and image:
    template.xpath("/svg:svg",namespaces={"svg":"http://www.w3.org/2000/svg"})[0].attrib["width"]=str(svg_width)+"mm"
    template.xpath("/svg:svg",namespaces={"svg":"http://www.w3.org/2000/svg"})[0].attrib["height"]=str(svg_height)+"mm"
    template.xpath("/svg:svg",namespaces={"svg":"http://www.w3.org/2000/svg"})[0].attrib["viewBox"]="0 0 "+str(svg_width)+" "+str(svg_height)
    template.xpath("/svg:svg",namespaces={"svg":"http://www.w3.org/2000/svg"})[0].attrib["{http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd}docname"]=os.path.split(output_filename)[1]
    template.xpath("/svg:svg/svg:image",namespaces={"svg":"http://www.w3.org/2000/svg"})[0].attrib["width"]=str(svg_width)
    template.xpath("/svg:svg/svg:image",namespaces={"svg":"http://www.w3.org/2000/svg"})[0].attrib["height"]=str(svg_height)
    template.xpath("/svg:svg/svg:image",namespaces={"svg":"http://www.w3.org/2000/svg"})[0].attrib["{http://www.w3.org/1999/xlink}href"]=pathname2url(os.path.split(input_filename)[1])
    viewdoc=etree.tostring(template)
    svgfile = etree.ElementTree(template)
    svgfile.write(output_filename)
    
    pass
