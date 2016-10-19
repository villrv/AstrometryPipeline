import alipy
import glob
import numpy as np
from photutils.morphology import centroid_com
from astropy.io import fits
from astropy import wcs
import sep
import matplotlib.pyplot as plt
from photutils import CircularAperture,EllipticalAperture,aperture_photometry
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import warnings
from astropy import cosmology

warnings.simplefilter('ignore', UserWarning)

MY_SN = np.loadtxt('input.dat',usecols=(0,),dtype='str')

Z, PIXEL_SCALE,THRES_REF,THRES_TRANS,GAIN_REF,GAIN_TRANS,DEBLEND_REF,DEBLEND_TRANS = \
        np.loadtxt('input.dat',usecols=(1,2,3,4,5,6,7,8),unpack=True)



my_table = Table.read('master_list.txt',format='ascii.csv',guess=False)

cut = my_table['SN name'] == MY_SN
sncoord_wcs = SkyCoord(my_table['SN RA'][cut][0]+' '+my_table['SN Dec'][cut][0], unit=(u.hourangle, u.deg))
hostcoord_wcs = SkyCoord(my_table['Host RA'][cut][0]+' '+my_table['Host Dec'][cut][0], unit=(u.hourangle, u.deg))



def loadimage(filename):
    hdulist = fits.open(filename)
    head = hdulist[0].header
    w = wcs.WCS(head)
    imagedata = hdulist[0].data
    hdulist.close()
    return imagedata, w, head

def match_object(image,thresh,gain,deblend,wcs_coord,plotting=True):

    #run sep to detect objects
    bkg = sep.Background(image.copy(order='C'))

    #bkg = sep.Background(image)
    thresh = thresh * bkg.globalrms
    objects,segmap = sep.extract(image, thresh, err=bkg.rms(),gain=gain,segmentation_map=True,deblend_cont = deblend)
    pos = np.column_stack((objects['x'],objects['y']))
    pos_wcs = w.wcs_pix2world(pos,1)
    pos_wcs_sc = SkyCoord(pos_wcs, unit='deg')
    idx, d2d, d3d = wcs_coord.match_to_catalog_sky(pos_wcs_sc)
    matched_pixel_coord = pos[idx]
    x_uncertainty, y_uncertainty = np.sqrt(objects[idx]['errx2']),np.sqrt(objects[idx]['erry2'])

    #Figure out the brightest point in the galaxy, and the offsets from the center and the transient
    x_peak = objects[idx]['xcpeak']
    y_peak = objects[idx]['ycpeak']

    if plotting == True:

        apertures = CircularAperture(pos,r=7.)
        mean, median, std = sigma_clipped_stats(image[200:800,200:800], sigma=3.0, iters=5)
        plt.figure()
        plt.imshow(image, vmin= mean-3.*std, vmax=mean+3.*std, cmap='Greys', origin='lower',interpolation="none")
        apertures.plot(color='red',lw=1.5)
        #use phot utils to find the nearest object to the host in background
        plt.plot(x_peak,y_peak,'r+')
        match_aper = CircularAperture(matched_pixel_coord, r=9.)
        match_aper.plot(color='green',lw=3.0)
        plt.show()


    object_to_peak = np.sqrt((x_peak-matched_pixel_coord[0])**2 + (y_peak - matched_pixel_coord[1])**2)



    r, flag = sep.flux_radius(image, matched_pixel_coord[0], matched_pixel_coord[1],
                          3.*objects['a'][idx], 0.5, subpix=5)


    radius = r * PIXEL_SCALE



    #Calculate total flux
    total_flux = objects[idx]['flux']
    #get the fluxes in the object from segmap
    all_fluxes = image[np.where(segmap == (idx+1))]



    return matched_pixel_coord,x_uncertainty,y_uncertainty,object_to_peak,segmap,radius,total_flux,all_fluxes

images_to_align = ["trans.fits"]
ref_image = "ref.fits"


idn = alipy.ident.run(ref_image, images_to_align, visu=False)
idn = idn[0]
outputshape = alipy.align.shape(ref_image)


alipy.align.irafalign(idn.ukn.filepath, idn.uknmatchstars, idn.refmatchstars, shape=outputshape, makepng=False)

#Uncomment if you want affine transform instead
alipy.align.affineremap(idn.ukn.filepath, idn.trans, shape=outputshape, makepng=True)

print "FINISHED ORIGINAL IDs"

aligned_image = "./alipy_out/trans_gregister.fits"

idn2 = alipy.ident.run(ref_image, [aligned_image], visu=False)
idn2 = idn2[0]
rms_x = 0
rms_y = 0
counter = 0

for i,ref_star in enumerate(idn2.refmatchstars):
    ukn_star = idn2.uknmatchstars[i]

    rms_x = rms_x + (ref_star.x - ukn_star.x)**2
    rms_y = rms_y + (ref_star.y - ukn_star.y)**2
    counter += 1


rms_in_x = np.sqrt(rms_x)/counter
rms_in_y = np.sqrt(rms_y)/counter

astrometric_error_x,astrometric_error_y = rms_in_x,rms_in_y


total_rms = np.sqrt(rms_in_x**2 + rms_in_y**2) * PIXEL_SCALE


#take the image data
data_ref,w,head = loadimage(ref_image)
data_sep = data_ref.byteswap().newbyteorder()

host_rms = sep.Background(data_sep).globalrms

print "HOST RMS",host_rms

#run match_object on host
match_pix_host,host_x_stat_unc,host_y_stat_unc,host_dist_peak,host_segmap,host_radius,host_total_flux,host_all_fluxes =  \
        match_object(data_sep,THRES_REF,GAIN_REF,DEBLEND_REF,hostcoord_wcs)


#use phot utils to find the nearest object to the transient in transient
data_trans,w,head = loadimage(aligned_image)
data_trans = np.asarray(data_trans,dtype=float)
trans_bkg = sep.Background(data_trans)
data_trans = data_trans - trans_bkg

#run match_object on transient
match_pix_sn,trans_x_stat_unc,trans_y_stat_unc,trans_dist_peak,trans_segmap,trans_radius,trans_total_flux,\
        trans_all_fluxes = match_object(data_trans,THRES_TRANS,GAIN_TRANS,DEBLEND_TRANS,sncoord_wcs)

trans_stat_unc = np.sqrt(trans_x_stat_unc**2 + trans_y_stat_unc**2) * PIXEL_SCALE

#find systematic error due to thres value
thresholds = np.linspace(2,10,20) * host_rms
pix_locations = np.zeros((len(thresholds),2))
for i,threshold in enumerate(thresholds):
    objects = sep.extract(data_sep, threshold)
    pos = np.column_stack((objects['x'],objects['y']))
    pos_wcs = w.wcs_pix2world(pos,1)
    pos_wcs_sc = SkyCoord(pos_wcs, unit='deg')
    idx, d2d, d3d = hostcoord_wcs.match_to_catalog_sky(pos_wcs_sc)
    match_pix_temp = pos[idx]
    pix_locations[i,:] = match_pix_temp


systematic_error = np.sqrt(np.var(pix_locations[:,0])+np.var(pix_locations[:,1]))

total_host_error = np.sqrt(systematic_error**2 + host_x_stat_unc**2 + host_y_stat_unc**2) * PIXEL_SCALE

#compute offset! woo
offset_arcsec = np.sqrt(np.sum((match_pix_host - match_pix_sn)**2)) * PIXEL_SCALE
normalized_offset = offset_arcsec / host_radius


#Get error position for the transient on the host image - with astrometric and transient uncertainty
aperture_size_x,aperture_size_y = np.sqrt(trans_x_stat_unc**2+astrometric_error_x**2),np.sqrt(trans_y_stat_unc**2+astrometric_error_y**2)

#Create a aperture on this position with this error size and take average
aperture = EllipticalAperture(match_pix_sn,a=aperture_size_x,b=aperture_size_y,theta=0)
aperture_sum = aperture_photometry(data_ref, aperture)['aperture_sum']
SNpos_flux = aperture_sum/aperture.area()

print "APERTURE SIZE",aperture_size_y,aperture_size_x

#plot this segmap to make sure there aren't any stars that are included
plt.figure()
plt.imshow(data_ref[match_pix_host[0]-100:match_pix_host[0]+100,match_pix_host[1]-100:match_pix_host[1]+100],cmap='Greys')
plt.imshow(host_segmap[match_pix_host[0]-100:match_pix_host[0]+100,match_pix_host[1]-100:match_pix_host[1]+100]/
                host_segmap[match_pix_host[0]-100:match_pix_host[0]+100,match_pix_host[1]-100:match_pix_host[1]+100],alpha=0.4)
plt.show()


#Generate many fake hosts and calculate the fractional flux
num_trials = 5000
my_fractional_fluxes = []
for trial in np.arange(num_trials):
    new_host = np.random.normal(0.0,host_rms,np.shape(host_all_fluxes)) + np.random.poisson(np.abs(host_all_fluxes),np.shape(host_all_fluxes))
    new_total_flux = np.sum(new_host,axis=None)
    supernova_flux = np.random.poisson(SNpos_flux) + np.random.normal(0.0,host_rms)
    my_fractional_fluxes = np.append(my_fractional_fluxes,np.sum(new_host[np.where(new_host<=supernova_flux)])/new_total_flux)

plt.hist(my_fractional_fluxes)
plt.show()

print "Fractinoal Flux error",np.mean(my_fractional_fluxes),np.std(my_fractional_fluxes)

#take array from segmap and do percentile calculation

print "FRACTIONAL FLUX",np.sum(host_all_fluxes[np.where(host_all_fluxes<=SNpos_flux)])/host_total_flux
print "FRACTIONAL FLUX",np.sum(host_all_fluxes[np.where(host_all_fluxes<=330)])/host_total_flux
print "FRACTIONAL FLUX",np.sum(host_all_fluxes[np.where(host_all_fluxes<=390)])/host_total_flux

print "SN_pos flux and app size",SNpos_flux,aperture_size_x,aperture_size_y
print match_pix_sn

#calculate a distance
cosmo = cosmology.FlatLambdaCDM(H0=70.,Om0=0.3)
R_kpc = cosmo.angular_diameter_distance(Z) * 1000. * offset_arcsec * 4.84814e-6

offset_error = np.sqrt(total_rms**2+trans_stat_unc**2+total_host_error**2)

f = open('output.dat','w')

print np.shape(R_kpc),R_kpc.value
f.write(str(MY_SN) + ' ' + str(total_rms) + ' ' + str(trans_stat_unc) + ' ' + str(total_host_error) + ' ' + \
        str(offset_arcsec) + ' ' + str(offset_error) + ' ' + str(host_radius) + ' ' + str(R_kpc.value)\
     + ' ' + str(np.mean(my_fractional_fluxes)) + ' ' + str(np.std(my_fractional_fluxes)))

f.close()
