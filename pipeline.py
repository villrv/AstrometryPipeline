import alipy
import glob
import numpy as np
from photutils.morphology import centroid_com
from astropy.io import fits
from astropy import wcs
import sep
import matplotlib.pyplot as plt
from photutils import CircularAperture
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u


MY_SN = 'sn2007i'

PIXEL_SCALE = 0.26

def loadimage(filename):
    hdulist = fits.open(filename)
    head = hdulist[0].header
    w = wcs.WCS(head)
    imagedata = hdulist[0].data
    hdulist.close()
    return imagedata, w, head

images_to_align = sorted(glob.glob("trans.fits"))
ref_image = "ref.fits"

identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
# That's it !
# Put visu=True to get visualizations in form of png files (nice but much slower)
# On multi-extension data, you will want to specify the hdu (see API doc).

'''# The output is a list of Identification objects, which contain the transforms :
for id in identifications: # list of the same length as images_to_align.
        if id.ok == True: # i.e., if it worked

                print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
                # id.trans is a alipy.star.SimpleTransform object. Instead of printing it out as a string,
                # you can directly access its parameters :
                #print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
                #print id.trans.matrixform()
                #print id.trans.inverse() # this returns a new SimpleTransform object

        else:
                print "%20s : no transformation found !" % (id.ukn.name)

# Minimal example of how to align images :
'''

outputshape = alipy.align.shape(ref_image)

for id in identifications:
        if id.ok == True:

                # Variant 1, using only scipy and the simple affine transorm :
                alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=True)

                # Variant 2, using geomap/gregister, correcting also for distortions :
                alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape, makepng=False)
                # id.uknmatchstars and id.refmatchstars are simply lists of corresponding Star objects.

                # By default, the aligned images are written into a directory "alipy_out".


print "FINISHED ORIGINAL IDs"

#Print the calculated RMS values from alignment...
aligned_image = "./alipy_out/trans_gregister.fits"

identifications = alipy.ident.run(ref_image, [aligned_image], visu=False)

rms_x = 0
rms_y = 0
counter = 0
for id in identifications:
    if id.ok == True:
        for i,ref_star in enumerate(id.refmatchstars):
            ukn_star = id.uknmatchstars[i]

            rms_x = rms_x + (ref_star.x - ukn_star.x)**2
            rms_y = rms_y + (ref_star.y - ukn_star.y)**2
            counter += 1


print "RMS in X:",np.sqrt(rms_x)/counter
print "RMS in Y:",np.sqrt(rms_y)/counter


#take the image data
data_ref,w,head = loadimage(ref_image)
data_sep = data_ref.byteswap().newbyteorder()

#run sep to detect objects on the background *and* transient image
bkg = sep.Background(data_sep)
thresh = 3. * bkg.globalrms
objects = sep.extract(data_sep, thresh)
pos = np.column_stack((objects['x'],objects['y']))

pos_wcs = w.wcs_pix2world(pos,1)
pos_wcs_sc = SkyCoord(pos_wcs, unit='deg')

apertures = CircularAperture(pos,r=7.)

mean, median, std = sigma_clipped_stats(data_ref[200:800,200:700], sigma=3.0, iters=5)
plt.figure()
plt.imshow(data_ref, vmin= mean-3.*std, vmax=mean+3.*std, cmap='Greys', origin='lower')
apertures.plot(color='red',lw=1.5)
#use phot utils to find the nearest object to the host in background

plt.show()

my_table = Table.read('master_list.txt',format='ascii.csv',guess=False)

cut = my_table['SN name'] == MY_SN
sncoord_wcs = SkyCoord(my_table['SN RA'][cut][0]+' '+my_table['SN Dec'][cut][0], unit=(u.hourangle, u.deg))
hostcoord_wcs = SkyCoord(my_table['Host RA'][cut][0]+' '+my_table['Host Dec'][cut][0], unit=(u.hourangle, u.deg))


idx, d2d, d3d = hostcoord_wcs.match_to_catalog_sky(pos_wcs_sc)
match_wcs_host = pos_wcs_sc[idx]
match_pix_host = pos[idx]

match_aper = CircularAperture(match_pix_host, r=9.)
match_aper.plot(color='green')


r, flag = sep.flux_radius(data_sep, match_pix_host[0], match_pix_host[1],
                          3.*objects['a'][idx], 0.5, subpix=5)


host_radius = r * PIXEL_SCALE

#use phot utils to find th enearest object to the transient in transient
data_trans,w,head = loadimage(aligned_image)
data_trans = np.asarray(data_trans,dtype=float)



#data_trans = data_trans.byteswap().newbyteorder()
bkg = sep.Background(data_trans.copy(order='C'))

print bkg.globalrms

thresh = 1.5 * bkg.globalrms
objects = sep.extract(data_trans - bkg, thresh,deblend_cont=1e-3)
pos = np.column_stack((objects['x'],objects['y']))
pos_wcs = w.wcs_pix2world(pos,1)
pos_wcs_sc = SkyCoord(pos_wcs, unit='deg')
apertures = CircularAperture(pos,r=7.)



idx, d2d, d3d = sncoord_wcs.match_to_catalog_sky(pos_wcs_sc)
match_wcs_sn = pos_wcs_sc[idx]
match_pix_sn = pos[idx]

match_aper = CircularAperture(match_pix_sn, r=3.)

#compute offset! woo
offset_arcsec = np.sqrt(np.sum((match_pix_host - match_pix_sn)**2)) * PIXEL_SCALE
print "THE OFFSET!!!",offset_arcsec
print "normalized offset",offset_arcsec / host_radius
#accept that our offsets are meaningless until we compute the half-light radii :]


mean, median, std = sigma_clipped_stats(data_trans[200:800,200:700], sigma=3.0, iters=5)
plt.figure()
plt.imshow(data_trans,interpolation='none', vmin= mean-3.*std, vmax=mean+3.*std, cmap='Greys', origin='lower')
apertures.plot(color='red',lw=1.5)
match_aper.plot(color='green')


plt.show()

