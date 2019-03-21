import numpy as np
from matplotlib import pyplot as pl
import sys
from lxml import etree
import pdb
import glob
    
######################## funtion 1 ########################
### inputs: [svgfiles],[i] <-- "i" will be assigned in the loop that this function is made to be used for
### outputs: [full_path_xout],[path_theta_out],[full_path_yout],[full_path_theta_out],
###          [full_path_thetax_out],[full_path_thetay_out],
###          [full_path_thlength_out],[path_thlength_out]
def get_thetas(filenames,i):
    name=filenames[i]
    svgxml=etree.parse(name)
        
    # Find <svg:path> elements
    path_els = svgxml.xpath("//svg:path",namespaces={"svg":"http://www.w3.org/2000/svg"})
    
    assert(len(path_els)==1)  # Assume only one <path> in the .svg file
    
    # Extract "d=" attribute
    d_attr=path_els[0].attrib["d"]
    
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
            pass
        elif pathcmd == 'M':
            curpos = np.array([ float(pathparam) for pathparam in pathcmds[pathpos:pathpos+2] ],dtype='d')
            pathpos+=2
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
                path_thlength_out.append(th_length*10**-6) #meters
                pathpos+=6
                # Does the spline continue? 
                more_params = pathpos < len(pathcmds) and not(pathcmds[pathpos].isalpha())
                curpos = np.array([ xpos[-1],ypos[-1] ],dtype='d')
                pass
            
            pass
        else:
            raise ValueError("Unknown path command %s" % (pathcmd))
            pass 
    full_path_xout=np.concatenate(path_xout)*10**-6 # Converted into meters
    full_path_yout=np.concatenate(path_yout)*10**-6 
    full_path_theta_out=np.concatenate(path_theta_out)
    full_path_thetax_out=np.concatenate(path_thetax_out)
    full_path_thetay_out=np.concatenate(path_thetay_out)
    full_path_thlength_out=np.concatenate(path_thlength_out)
    path_thlength_out=np.array(path_thlength_out, dtype= 'd')
    return (full_path_xout,path_theta_out,full_path_yout,full_path_theta_out, full_path_thetax_out, full_path_thetay_out, full_path_thlength_out, path_thlength_out)
    pass

######################## funtion 2 ########################
### Inputs: [path_thlength_out]
### Outputs: [sampling_thetas]
def evenly_spaced_thetas(path_thlength_out,path_theta_out):
    ### We need to store the actual distances for each segment
    #d_seg=np.zeros(path_thlength_out.shape, dtype= 'd')
    #f=0.0
    #for i in range(path_thlength_out.shape[0]):
    #	for j in range(path_thlength_out.shape[1]):
    #		g=path_thlength_out[i,j]+f
    #		d_seg[i,j]=f
    #		f=g
    # This is done in 1 line where
    # d_seg[i,j] represents the jth element of the i'th point click
    # based on each segment being split into ~ 20 pieces (i)
    d_seg = np.cumsum(np.concatenate((np.array((0,),dtype='d'),path_thlength_out.reshape(np.prod(path_thlength_out.shape)))))[:-1].reshape(path_thlength_out.shape)#meters
    
    ### Pull in full crack length data ###
    crack_length= np.max(d_seg) # meters
    d0=3.0e-7 # distance between evenly spaced points meters
    crack_dist = np.arange(np.max(d_seg)//d0,dtype='d')*d0
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
    return (sampling_thetas,d0)

######################## funtion 3 ########################
### Inputs: [sampling_thetas]
### Outputs: [filtered],[eq_lengths]
def filtering(sampling_thetas,d0,f_cutoff):
        N_samp =sampling_thetas.shape[0] #number of sample points
        yf = np.fft.fft(sampling_thetas)
        assert((np.linalg.norm(yf,2)/np.linalg.norm(yf.real,2))-1 <= 1.0)
        #This assertion compares the 2-norm of complex yf with real yf
        #this ensures that the complex part is never larger than the real
        xf = np.arange(N_samp,dtype='d')/N_samp/d0
        xf[N_samp//2:] -= 1/d0 # symmetric about 0    
        ### Apply raised cosine/hanning window
        f_filter_region = (xf > -f_cutoff) & (xf < f_cutoff) #recall: f_cutoff specified by user
        raised_cos=np.cos((np.pi/4.0)*(2.0/f_cutoff)*xf[f_filter_region])**2.0
        #creates raised cosine with slope inflection point at f_cutoff/2.0
        yf[f_filter_region]=yf[f_filter_region]*raised_cos
        yf[~f_filter_region]=0.0
        filtered=np.real(np.fft.ifft(yf)) #can be done because of assertion
        eq_lengths=np.ones(filtered.shape[0])*d0
        fil_y=np.cumsum(np.sin(filtered))*d0
        fil_x=np.cumsum(np.cos(filtered))*d0
        return (filtered,eq_lengths,fil_x,fil_y)
        pass

######################## funtion 4 ########################
### Inputs: [full_path_xout],[full_path_yout],[full_path_theta_out],
###         [full_path_thlength_out],[filtered],[eq_lengths],[fil_x],[fil_y],[draw_path=(T/F)]
### Outputs: [sigma],[sigma_F]
### calculate average and standard deviation
def calc_stdv(full_path_xout,
                  full_path_yout,
                  full_path_theta_out,
                  full_path_thlength_out,
                  filtered,eq_lengths,
                  fil_x,
                  fil_y,
                  draw_path):
        mu = np.average(full_path_theta_out,weights=full_path_thlength_out)
        sigma= np.sqrt(np.average((full_path_theta_out-mu)**2,weights=full_path_thlength_out))
        mu_F = np.average(filtered,weights=eq_lengths)
        sigma_F= np.sqrt(np.average((filtered-mu_F)**2,weights=eq_lengths))

        if draw_path:
            pl.figure(20)
            pl.clf()
            pl.plot(full_path_xout-full_path_xout[0],full_path_yout-full_path_yout[0],'-',fil_x,fil_y,'-')
            pl.xlabel('Horizontal position (um)')
            pl.ylabel('Vertical position (um)')
            pl.show()
            # Interesting quantity: sum of theta*length, for left and right
            pass
        return (sigma, sigma_F, mu, mu_F)
        pass

def tortuosity_plots(
        theta_final,
        thlength_final,
        filtered_final,
        eq_lengths_final,
        avg_mu,
        avg_filtered_mu,
        avg_sigma,
        avg_filtered_sigma,
        savedir = None):

    pl.figure(1)
    pl.clf()
    n_bins=50
    bins=np.linspace(-90.0,90.0,n_bins)
    dbin = (bins[-1]-bins[0])/(n_bins-1) # Theta stepsize per bin
    (n_01,b_01,p_01)=pl.hist(theta_final[:]*180.0/np.pi,bins=bins,weights=thlength_final)
    #pl.plot(bins,(1.0/(np.sqrt(2.0*np.pi)*avg_sigma))*np.exp(-(bins*np.pi/180.0-avg_mu)**2.0/ (2.0*avg_sigma**2.0))*np.sum(thlength_final[:])*dbin*np.pi/180.0,'-')#fits to a gaussian curve
    pl.figtext(0.55,0.75,('mu=%.1fdeg\nsigma=%.1fdeg' %(avg_mu*180.0/np.pi,avg_sigma*180.0/np.pi)),bbox={'facecolor':'white','alpha':0.8,'pad':10},fontsize=25)
    pl.title('Angle Distribution',fontsize=30)
    pl.ylabel('Number of Instances',fontsize=20)
    pl.xlabel('Angle (degrees)',fontsize=20)
    if savedir:
        unfiltered_filename="histogram_unfiltered.png"
        pl.savefig(os.path.join(savedir,unfiltered_filename),dpi=300)
        pass

    pl.figure(2)
    pl.clf()
    (n_02,b_02,p_02)=pl.hist(filtered_final[:]*180.0/np.pi,bins=bins,weights=eq_lengths_final)
    #pl.plot(bins,(1.0/(np.sqrt(2.0*np.pi)*avg_filtered_sigma))*np.exp(-(bins*np.pi/180.0-avg_filtered_mu)**2.0/(2.0*avg_filtered_sigma**2.0))*np.sum(eq_lengths_final[:])*dbin*np.pi/180.0,'-') #fits to a gaussian curve
    pl.figtext(0.55,0.75,('mu=%.1fdeg\nsigma=%.1fdeg' %(avg_filtered_mu*180.0/np.pi,avg_filtered_sigma*180.0/np.pi)),bbox={'facecolor':'white','alpha':0.8,'pad':10},fontsize=25)
    pl.title('Filtered Angle Distribution',fontsize=30)
    pl.ylabel('Number of Instances',fontsize=20)
    pl.xlabel('Angle (degrees)',fontsize=20)
    if savedir:
        filtered_filename="histogram_filtered.png"
        pl.savefig(os.path.join(savedir,filtered_filename),dpi=300)
        pass

    #pl.savefig("./Histograms/"+base_name+"_filtered.png",dpi=300)

    if savedir is not None:
        pl.show()
        pass
    return (unfiltered_filename,filtered_filename)

### Now that the functions are all defined, time to test them.


def histogram_from_svgs(filenames,f_cutoff):
    unfiltered_mu=[]
    filtered_mu=[]
    unfiltered_sigma=[]
    filtered_sigma=[]
    all_theta=[]
    all_thlength=[]
    all_filtered=[]
    all_eq_lengths=[]
    for i in range(len(filenames)):
	(full_path_xout, path_theta_out, full_path_yout, full_path_theta_out, full_path_thetax_out, full_path_thetay_out, full_path_thlength_out, path_thlength_out) = get_thetas(filenames,i)
	(sampling_thetas,d0)=evenly_spaced_thetas(path_thlength_out,path_theta_out)
	(filtered,eq_lengths,fil_x,fil_y)=filtering(sampling_thetas,d0,f_cutoff)
	(sigma, sigma_F, mu, mu_F)=calc_stdv(full_path_xout,
                                             full_path_yout,
                                             full_path_theta_out,
                                             full_path_thlength_out,
                                             filtered,eq_lengths,
                                             fil_x,
                                             fil_y,
                                             draw_path=True)
        unfiltered_mu.append(mu) # mu and sigma in radians
	filtered_mu.append(mu_F)
	unfiltered_sigma.append(sigma)
	filtered_sigma.append(sigma_F)
	all_theta.append(full_path_theta_out[:])
	all_thlength.append(full_path_thlength_out[:])
	all_filtered.append(filtered[:])
	all_eq_lengths.append(eq_lengths[:])
	pass

    # !!!*** NOTE: NEED TO DO AVERAGE/STD OVER ALL TRACES TOGETHER
    theta_final=np.concatenate(all_theta)
    thlength_final=np.concatenate(all_thlength)
    filtered_final=np.concatenate(all_filtered)
    eq_lengths_final=np.concatenate(all_eq_lengths)
    avg_mu=np.sum(unfiltered_mu)/len(unfiltered_mu)
    avg_filtered_mu=np.sum(filtered_mu)/len(filtered_mu)
    avg_sigma=np.sum(unfiltered_sigma)/len(unfiltered_sigma)
    avg_filtered_sigma=np.sum(filtered_sigma)/len(filtered_sigma)

    return (
        theta_final,
        thlength_final,
        filtered_final,
        eq_lengths_final,
        avg_mu,
        avg_filtered_mu,
        avg_sigma,
        avg_filtered_sigma)

    #doplots(
    #    theta_final,
    #    thlength_final,
    #    filtered_final,
    #    eq_lengths_final,
    #    avg_mu,
    #    avg_filtered_mu,
    #    avg_sigma,
    #    avg_filtered_sigma)

