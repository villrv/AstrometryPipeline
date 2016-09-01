import alipy
import glob
import numpy as np
from photutils.morphology import centroid_com
from astropy.io import fits
from astropy import wcs

images_to_align = sorted(glob.glob("trans.fits"))
ref_image = "ref.fits"

identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
# That's it !
# Put visu=True to get visualizations in form of png files (nice but much slower)
# On multi-extension data, you will want to specify the hdu (see API doc).

# The output is a list of Identification objects, which contain the transforms :
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


#Next, get the host galaxy using WCS and compute its (and its transient's) centroid
hdulist = fits.open(ref_image)
w = wcs.WCS(hdulist[0].header)
print(w.wcs.name)

world_coord = np.array([[179.8042,-1.6053]], np.float_)

pixcrd = w.wcs_world2pix(world_coord,1)

print pixcrd

#The centroid command is just centroid_com(data)