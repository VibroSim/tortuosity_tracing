import numpy as np
from matplotlib import pyplot as pl
import scipy
from scipy import special
from scipy.signal import butter,lfilter
#from scipy import lfilter
import sys
from lxml import etree
import pdb
import glob
import os
import os.path


def get_num_paths(filename):
    """ Return number of paths in the specified file"""
    svgxml=etree.parse(filename)

    # Find <svg:path> elements
    path_els = svgxml.xpath("//svg:path",namespaces={"svg":"http://www.w3.org/2000/svg"})

    return len(path_els)

### inputs: [svgfiles],[i] <-- "i" will be assigned in the loop that this function is made to be used for
### outputs: [full_path_xout],[path_theta_out],[full_path_yout],[full_path_theta_out],
###          [full_path_thetax_out],[full_path_thetay_out],
###          [full_path_thlength_out],[path_thlength_out]
def get_thetas(filename,path_num):
    """ Read the specified .svg file, containing a single spline, 
as input. A scaling of 1000 between physical units and .svg units 
is assumed (i.e. 1 mm as recorded in the .svg is actually 1 micron). 
We break each segment into 20 steps, and return the angle
of each step and the length of each step, along with x and y coordinates.

    inputs: [filename] : name of .svg file to read
    outputs: 
    [full_path_xout]: concatenated x coordinates of the steps starting positions (meters)
    [path_theta_out]: list of array of angles, in radians, of the segments. 0 means horizontal 
    [full_path_yout]: concatenated y coordinates of the steps starting positions (meters)
    [full_path_theta_out], concatenated array of angles, in radians, of the segments
    [full_path_thetax_out], concatenated array of x coordinates of the steps center positions (meters)
    [full_path_thetay_out], concatenated array of y coordinates of the steps center positions (meters)
    [full_path_thlength_out], concatenated array of lengths of each step (meters)
    [path_thlength_out], list of arrays of step lengths for each segment (meters)
          """
    svgxml=etree.parse(filename)

    # Find <svg:path> elements
    path_els = svgxml.xpath("//svg:path",namespaces={"svg":"http://www.w3.org/2000/svg"})

    #assert(len(path_els)==1)  # Assume only one <path> in the .svg file

    # Extract "d=" attribute
    d_attr=path_els[path_num].attrib["d"]

    # Split d attribute by spaces and commas
    pathcmds=d_attr.replace(",",' ').split()

    # Initialize our position to the origin
    curpos = np.array((0,0),dtype='d')
    pathpos = 0  # index into pathcmds

    # Here is where we want to over sample and choose the value that is closest to the t that we want
    numsteps = 20 # number of steps through each curve segment
    seg_t=np.linspace(0.0,1.0,numsteps) # Times for coordinate evaluation within each segment
    seg_dt = 1.0/(numsteps-1.0)

    seg_deriv_t=seg_t[:-1]+seg_dt/2.0 # Times for derivative evaluation -- evaluate at centers of each segment

    seg_tsteps = [ 0.0 ]
    seg_tsteps.extend([ seg_dt ]*(numsteps-1))
    seg_tsteps=np.array(seg_tsteps,dtype='d')
    # seg_tsteps is now [ 0.0, seg_dt, seg_dt, seg_dt ... ]

    path_xout = []
    path_yout = []
    path_tstepout = []

    path_theta_out = []
    path_thlength_out = [] #the lengths of each segment
    d_seg_out =[] # the distance along the path
    path_thetax_out = []
    path_thetay_out = []

    # Step through the various path commands within the <svg:path> element's 'd=' attribute 
    while pathpos < len(pathcmds):
        pathcmd = pathcmds[pathpos]
        pathpos+=1
        if pathcmd == 'm':
            curpos += np.array([ float(pathparam) for pathparam in pathcmds[pathpos:pathpos+2] ],dtype='d')
            pathpos+=2

            more_params = pathpos < len(pathcmds) and not(pathcmds[pathpos].isalpha())
            if more_params: 

                raise ValueError("Path in file %s contains line segment(s) (implicit line after move command)")


            pass
        elif pathcmd == 'M':
            curpos = np.array([ float(pathparam) for pathparam in pathcmds[pathpos:pathpos+2] ],dtype='d')
            pathpos+=2

            more_params = pathpos < len(pathcmds) and not(pathcmds[pathpos].isalpha())
            if more_params: 

                raise ValueError("Path in file %s contains line segment(s) (implicit line after move command)")

            pass
        elif pathcmd == 'C' or pathcmd == "c":
            # Reference https://www.w3.org/TR/SVG/paths.html#PathElement

            more_params=True
            orig_curpos=curpos

            while more_params:
                cparams = pathcmds[pathpos:pathpos+6]
                cparams_float = [ float(cparam) for cparam in cparams ]
                (x1,y1,x2,y2,x,y) = cparams_float

                if pathcmd=='c':
                    # relative
                    if x1==0 and x2==0 and x==0 and y==0 and y1==0 and y2==0:
                        # Throws error if multiple nodes are in one spot. This would result in undefined slopes.
                        raise ValueError("Path in file "+str(filename)+" contains multiple nodes in the same position.\nPlease check the .svg file")
                        pass
                    x1=curpos[0]+x1
                    y1=curpos[1]+y1
                    x2=curpos[0]+x2
                    y2=curpos[1]+y2
                    x=curpos[0]+x
                    y=curpos[1]+y
                    pass


                #for tcnt in range(numsteps):
                xpos = ((1.0-seg_t)**3.0)*curpos[0] + (3.0*seg_t*(1-seg_t)**2.0)*x1 + (3.0*(seg_t**2.0)*(1-seg_t))*x2 + ((seg_t)**3.0)*x;
                ypos = ((1.0-seg_t)**3.0)*curpos[1] + (3.0*seg_t*(1-seg_t)**2.0)*y1 + (3.0*(seg_t**2.0)*(1-seg_t))*y2 + ((seg_t)**3.0)*y;

                path_xout.append(xpos)
                path_yout.append(ypos)
                path_tstepout.append(seg_tsteps)

                # dx/dt
                # (1-t)^3 -> -3*(1-t)^2
                # (3(1-t)^2)(t) -> -6(1-t)(t) + 3(1-t)^2
                # 3(1-t)(t^2) -> -3t^2 + 6(1-t)t
                # t^3  -> 3t^2

                dxpos_dt = (-3.0*(1.0-seg_deriv_t)**2.0)*curpos[0] + (-6.0*(1.0-seg_deriv_t)*seg_deriv_t + 3.0*(1.0-seg_deriv_t)**2.0)*x1 + (-3.0*seg_deriv_t**2.0 + 6.0*(1.0-seg_deriv_t)*seg_deriv_t)*x2 + (3.0*seg_deriv_t**2.0)*x;

                dypos_dt = (-3.0*(1.0-seg_deriv_t)**2.0)*curpos[1] + (-6.0*(1.0-seg_deriv_t)*seg_deriv_t + 3.0*(1.0-seg_deriv_t)**2.0)*y1 + (-3.0*seg_deriv_t**2.0 + 6.0*(1.0-seg_deriv_t)*seg_deriv_t)*y2 + (3.0*seg_deriv_t**2.0)*y;

                # theta = atan2(dy,dx)
                theta=np.arctan2(dypos_dt,dxpos_dt)
                th_length=np.sqrt(dxpos_dt**2 + dypos_dt**2)*seg_dt
                # Also store x and y positions at the segment midpoints (corresponding to theta)
                theta_x = ((1.0-seg_deriv_t)**3.0)*curpos[0] + (3.0*seg_deriv_t*(1-seg_deriv_t)**2.0)*x1 + (3.0*(seg_deriv_t**2.0)*(1-seg_deriv_t))*x2 + ((seg_deriv_t)**3.0)*x;
                theta_y = ((1.0-seg_deriv_t)**3.0)*curpos[1] + (3.0*seg_deriv_t*(1-seg_deriv_t)**2.0)*y1 + (3.0*(seg_deriv_t**2.0)*(1-seg_deriv_t))*y2 + ((seg_deriv_t)**3.0)*y;

                # analytical expressions for derivatives have been validated that they match
                # these numerical evaluations
                #theta_numerical = np.arctan2(ypos[1:]-ypos[:-1],xpos[1:]-xpos[:-1])
                #th_length_numerical = np.sqrt((ypos[1:]-ypos[:-1])**2.0 + (xpos[1:]-xpos[:-1])**2.0)

                path_theta_out.append(theta)
                path_thetax_out.append(theta_x)
                path_thetay_out.append(theta_y)
                path_thlength_out.append(th_length*10**-6) # 10**-6 converts from mm-based SVG with assumed 1000x scale to meters
                pathpos+=6
                # Does the spline continue? 
                more_params = pathpos < len(pathcmds) and not(pathcmds[pathpos].isalpha())
                curpos = np.array([ xpos[-1],ypos[-1] ],dtype='d')
                pass

            pass
        else:
            raise ValueError("Unknown path command %s" % (pathcmd))
            pass 
    full_path_xout=np.concatenate(path_xout)*10**-6 # 10**-6 converts from mm-based SVG with assumed 1000x scale to meters
    full_path_yout=np.concatenate(path_yout)*10**-6 # 10**-6 converts from mm-based SVG with assumed 1000x scale to meters
    full_path_theta_out=np.concatenate(path_theta_out)
    full_path_thetax_out=np.concatenate(path_thetax_out)*1e-6 # 10**-6 converts from mm-based SVG with assumed 1000x scale to meters
    full_path_thetay_out=np.concatenate(path_thetay_out)*1e-6 # 10**-6 converts from mm-based SVG with assumed 1000x scale to meters
    full_path_thlength_out=np.concatenate(path_thlength_out)
    path_thlength_out=np.array(path_thlength_out, dtype= 'd')
    return (full_path_xout,path_theta_out,full_path_yout,full_path_theta_out, full_path_thetax_out, full_path_thetay_out, full_path_thlength_out, path_thlength_out)
    pass

######################## funtion 2 ########################
### Inputs: [path_thlength_out][path_theta_out]
### Outputs: [sampling_thetas]
def evenly_spaced_thetas(path_thlength_out,path_theta_out,point_spacing):
    """
    Convert inconsistenly-spaced angles to consistently-spaced angles
    with the given point spacing. 

    point_spacing is the distance between the uniformly spaced
    points we will return, in meters"""

    ### We need to store the actual distances for each segment
    #d_seg=np.zeros(path_thlength_out.shape, dtype= 'd')
    #f=0.0
    #for i in range(path_thlength_out.shape[0]):
    #   for j in range(path_thlength_out.shape[1]):
    #           g=path_thlength_out[i,j]+f
    #           d_seg[i,j]=f
    #           f=g
    # This is done in 1 line where
    # d_seg[i,j] represents the jth element of the i'th point click
    # based on each segment being split into ~ 20 pieces (i)
    d_seg = np.cumsum(np.concatenate((np.array((0,),dtype='d'),path_thlength_out.reshape(np.prod(path_thlength_out.shape)))))[:-1].reshape(path_thlength_out.shape)#meters

    ### Pull in full crack length data ###
    crack_length= np.max(d_seg) # meters
    #d0=3.0e-7 # distance between evenly spaced points meters
    crack_dist = np.arange(np.max(d_seg)//point_spacing,dtype='d')*point_spacing
    path_theta_out=np.array(path_theta_out, dtype= 'd')#from touple to array
    sampling_thetas=np.zeros(crack_dist.shape)

    assert(crack_dist[-1] <= d_seg[-1,-1])  #crack_dist cannot go beyond data in d_seg

    # pick out elements closest to sample_dist
    rownumber=0
    for j in range(crack_dist.shape[0]):
        FoundMyIndex= False
        while not FoundMyIndex and rownumber < d_seg.shape[0]:  
            d_idx = np.argmin(np.abs(d_seg[rownumber,:]-crack_dist[j]))
            if d_idx< d_seg.shape[1]-1 or (rownumber<d_seg.shape[0]-1 and np.abs(d_seg[rownumber+1,0]-crack_dist[j])>np.abs(d_seg[rownumber,d_idx]-crack_dist[j])) :
                FoundMyIndex= True
                pass
            else:
                rownumber += 1
                pass
            pass
        if rownumber == d_seg.shape[0]:
            # overrun... use last element
            rownumber = d_seg.shape[0]-1
            d_idx=d_seg.shape[1]-1
            pass
        sampling_thetas[j]=path_theta_out[rownumber,d_idx]
        pass
    return sampling_thetas

######################## funtion 3 ########################
### Inputs: [sampling_thetas]
### Outputs: [filtered],[eq_lengths] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Need to change the filter!!!!!!!!!!!!!
def filtering(sampling_thetas,point_spacing,f_cutoff):
    """point_spacing is the uniform distance spacing of the sampling_thetas, in meters
    f_cutoff is the cutoff spatial frequency in m^-1"""
    N_samp =sampling_thetas.shape[0] #number of sample points
    #yf = np.fft.fft(sampling_thetas)
    #assert((np.linalg.norm(yf,2)/np.linalg.norm(yf.real,2))-1 <= 1.0)
    xf = np.arange(N_samp,dtype='d')/N_samp/point_spacing
    xf[N_samp//2:] -= 1.0/point_spacing # symmetric about 0 
    #adapted from:
    #https://stackoverflow.com/questions/12093594/how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter
    def butter_bandpass(lowcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        #high = highcut / nyq
        b, a = butter(order, low, btype='low')
        return b, a
        pass
    def butter_bandpass_filter(data, lowcut, fs, order=5):
        b, a = butter_bandpass(lowcut, fs, order=order)
        y = lfilter(b, a, data)
        return y
        pass
    
    low= f_cutoff
    fsamp=1.0/point_spacing
    order=5
    filtered_thetas=butter_bandpass_filter(sampling_thetas,low,fsamp,order=order)
    #pl.plot(xf,np.fft.ifft(filtered_thetas))
    eq_lengths=np.ones(filtered_thetas.shape[0])*point_spacing
    filtered_ypath=np.cumsum(np.sin(filtered_thetas))*point_spacing
    filtered_xpath=np.cumsum(np.cos(filtered_thetas))*point_spacing
    return (filtered_thetas,eq_lengths,filtered_xpath,filtered_ypath)


######################## function 4 ########################
### Inputs: [full_path_xout],[full_path_yout],[full_path_theta_out],

###         [full_path_thlength_out],[filtered],[eq_lengths],[filtered_xpath],[filtered_ypath],[draw_path=(T/F)]

### Outputs: [sigma],[sigma_F]
### calculate average and standard deviation

def calc_stdv(full_path_theta_out, full_path_thlength_out, filtered_thetas,eq_lengths):

        mu = np.average(full_path_theta_out,weights=full_path_thlength_out)
        sigma= np.sqrt(np.average((full_path_theta_out-mu)**2,weights=full_path_thlength_out))

        mu_F = np.average(filtered_thetas,weights=eq_lengths)
        sigma_F= np.sqrt(np.average((filtered_thetas-mu_F)**2,weights=eq_lengths))
        return (sigma, sigma_F, mu, mu_F)

def draw_path(full_path_xout,full_path_yout,filtered_xpath,filtered_ypath,measnum,path_index,specimen,savedir = None):
    pl.figure()
    pl.clf()
    pl.plot((full_path_xout-full_path_xout[0])*1.e6,(full_path_yout-full_path_yout[0])*1.e6,'-',
            filtered_xpath*1.e6,filtered_ypath*1.e6,'-')
    pl.xlabel('Horizontal position (um)')
    pl.ylabel('Vertical position (um)')
    if measnum is not None:
        pl.title('Original and filtered path; measnum=%d' % (measnum))
        pass
    else:
        pl.title('Original and filtered path')
        pass

    if savedir is not None:
        tortuosity_plot_filename="path_comparison_%s_meas%3.3d.%d.png" % (specimen,measnum,path_index)
        pl.savefig(os.path.join(savedir,tortuosity_plot_filename),dpi=300)
        # Interesting quantity: sum of theta*length, for left and right
        pass
    else:
        tortuosity_plot_filename=None
        pass
    return tortuosity_plot_filename
     

def tortuosity_plots(
        specimen,
        theta_final,
        thlength_final,
        filtered_thetas,
        eq_lengths_final,
        avg_mu,
        avg_filtered_mu,
        avg_sigma,avg_filtered_sigma,savedir):
    pl.figure()
    pl.clf()
    n_bins=50
    bins=np.linspace(-90.0,90.0,n_bins)
    dbin = (bins[-1]-bins[0])/(n_bins-1) # Theta stepsize per bin
    (n_01,b_01,p_01)=pl.hist(theta_final[:]*180.0/np.pi,bins=bins,weights=thlength_final*10**6)
    pl.plot(bins,(1.0/(np.sqrt(2.0*np.pi)*avg_sigma))*np.exp(-(bins*np.pi/180.0-avg_mu)**2.0/ (2.0*avg_sigma**2.0))*np.sum(thlength_final[:])*dbin*np.pi/180.0,'-')#fits to a gaussian curve
    pl.figtext(0.55,0.75,('mu=%.1fdeg\nsigma=%.1fdeg' %(avg_mu*180.0/np.pi,avg_sigma*180.0/np.pi)),bbox={'facecolor':'white','alpha':0.8,'pad':10},fontsize=25)
    pl.title('Angle Distribution',fontsize=30)
    pl.ylabel('Number of Instances',fontsize=20)
    pl.xlabel('Angle (degrees)',fontsize=20)
    if savedir is not None:
        unfiltered_filename="%s_histogram_unfiltered.png" % (specimen)
        pl.savefig(os.path.join(savedir,unfiltered_filename),dpi=300)
        pass
    else:
        unfiltered_filename=None
        pass

    pl.figure()
    pl.clf()
    (n_02,b_02,p_02)=pl.hist(filtered_thetas[:]*180.0/np.pi,bins=bins,weights=eq_lengths_final*10**6)
    pl.plot(bins,(1.0/(np.sqrt(2.0*np.pi)*avg_filtered_sigma))*np.exp(-(bins*np.pi/180.0-avg_filtered_mu)**2.0/(2.0*avg_filtered_sigma**2.0))*np.sum(eq_lengths_final[:])*dbin*np.pi/180.0,'-') #fits to a gaussian curve
    pl.figtext(0.55,0.75,('mu=%.1fdeg\nsigma=%.1fdeg' %(avg_filtered_mu*180.0/np.pi,avg_filtered_sigma*180.0/np.pi)),bbox={'facecolor':'white','alpha':0.8,'pad':10},fontsize=25)
    pl.title('Filtered Angle Distribution',fontsize=30)
    pl.ylabel('Number of Instances',fontsize=20)
    pl.xlabel('Angle (degrees)',fontsize=20)
    if savedir is not None:
        filtered_filename="%s_histogram_filtered.png" % (specimen)
        pl.savefig(os.path.join(savedir,filtered_filename),dpi=300)
        pass
    else:
        filtered_filename=None
        pass
        

    #if savedir is not None:
    #    pl.show()
    #    pass
    return (unfiltered_filename,filtered_filename)

### Now that the functions are all defined, time to use them.

def histogram_from_svgs(filenames,measnums,f_cutoff,specimen,savedir,point_spacing):
    all_path_xout=[]
    all_path_yout=[]
    all_theta=[]
    all_thlength=[]
    all_filtered_theta=[]
    all_eq_lengths=[]
    all_filtered_xpath=[]
    all_filtered_ypath=[]
    tortuosity_path_filenames=[]
    tortuosity_plot_filenames=[]
    tortuosity_path_indexes=[]

    for i in range(len(filenames)):
        
        for path_index in range(get_num_paths(filenames[i])):
            
            
            (full_path_xout, path_theta_out, full_path_yout, full_path_theta_out, full_path_thetax_out, full_path_thetay_out, full_path_thlength_out, path_thlength_out) = get_thetas(filenames[i],path_index)
            sampling_thetas=evenly_spaced_thetas(path_thlength_out,path_theta_out,point_spacing)
            (filtered_theta,eq_lengths,filtered_xpath,filtered_ypath)=filtering(sampling_thetas,point_spacing,f_cutoff)
            if measnums is None:
                measnum=None
                pass
            else :
                measnum=measnums[i]
                pass
                
            #(dirpath,filename) = os.path.split(filenames[i]])
            #(filebasename,ext) = os.path.splitext(filename)
            #specimen=filebasename
            tortuosity_plot_filename=draw_path(full_path_xout,full_path_yout,filtered_xpath,filtered_ypath,measnum,path_index,specimen,savedir)
            tortuosity_plot_filenames.append(tortuosity_plot_filename)
            tortuosity_path_filenames.append(filenames[i])
            tortuosity_path_indexes.append(path_index)
            all_path_xout.append(full_path_xout[:])
            all_path_yout.append(full_path_yout[:])
            all_theta.append(full_path_theta_out[:])
            all_thlength.append(full_path_thlength_out[:])
            all_filtered_theta.append(filtered_theta[:])
            all_eq_lengths.append(eq_lengths[:])
            all_filtered_xpath.append(filtered_xpath[:])
            all_filtered_ypath.append(filtered_ypath[:])
            pass
        pass
    theta_final=np.concatenate(all_theta)
    thlength_final=np.concatenate(all_thlength)
    filtered_theta_final=np.concatenate(all_filtered_theta)
    eq_lengths_final=np.concatenate(all_eq_lengths)
    (sigma, sigma_F, mu, mu_F)=calc_stdv(theta_final,thlength_final, filtered_theta_final, eq_lengths_final)

    return (
        theta_final,
        thlength_final,
        filtered_theta_final,
        eq_lengths_final,
        mu,
        mu_F,
        sigma,
        sigma_F,
        tortuosity_path_filenames,
        tortuosity_plot_filenames,
        tortuosity_path_indexes)
    pass

