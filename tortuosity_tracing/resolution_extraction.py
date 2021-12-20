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


def convert_tiff_to_svg(input_filename):
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
    assert img.tag.keys()[-2]==34682, "Can not find parameters. File was probably not generated via SEM. Consider running 'svg_from_xlg()' with inputs manually entered in the xlg file."
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
    ResolutionX=float(config.get("Image","ResolutionX")) # also the float of img.width
    ResolutionY=float(config.get("Image","ResolutionY")) # also the float fo img.height
    PixelWidth=float(config.get("Scan","PixelWidth"))
    PixelHeight=float(config.get("Scan","PixelHeight"))
    assert PixelWidth==PixelHeight,'Pixels are not square. Check file.'  # Only support square pixels for now
    #width [um] = pixel size [m/px] * image size [px] * 1e6 [um] / 1 [m]
    svg_width=PixelWidth*ResolutionX*1e6 # physical width in mm of image as read from .tiff
    svg_height=PixelHeight*ResolutionY*1e6 # physical height in mm of image as read from .tiff, including SEM notations
    return (svg_width, svg_height)
    pass

# If images were collected without saved exif data, resolution was entered to the xlg file manually
def svg_from_xlg(input_filename,xResolution, yResolution):

    #images from the SEM should use the saved SEM data
    assert img.tag.keys()[-2]!=34682, 'SEM images should be evaluated using "convert_tiff_to_svg()".'
    #images missing the 305 tag might have a different set of exif data than the images from the Olympus Confocal Microscope
    assert img.tag.keys()[-2]==305, 'Check that exif data tags align with those generated from confocal microscope tiffs'
    img=Image.open(input_filename)
    inidata=img.tag
    ImageWidth=float(inidata[256][0]) #[px]
    ImageHeight=float(inidata[257][0]) #[px]
    ResolutionX=xResolution#[m/px]
    ResolutionY=yResolution#[m/px]
    #width[um]=[px]*[m/px]*1.0e6[um]/1.0[m]
    svg_width=ImageWidth*ResolutionX*1e6 # [um]
    svg_height=ImageHeight*ResolutionY*1e6 # [um]
    return (svg_width, svg_height)
    pass

def build_svg(output_filename,svg_width, svg_height):
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

    ## Next, Add in the calculated values and image
    ## let and 1 [svg mm] = 1 [physical um]
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
