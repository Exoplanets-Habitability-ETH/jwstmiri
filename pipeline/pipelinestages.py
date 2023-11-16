############
# Define the pipeline stages beginning from stage 2 to stage 3
############

import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.stats import sigma_clip
import jwst
# JWST pipelines (encompassing many steps)
from jwst.pipeline import Detector1Pipeline
from jwst.pipeline import Spec2Pipeline
from jwst.pipeline import Spec3Pipeline

# JWST pipeline utilities
from jwst import datamodels # JWST datamodels
from jwst.associations import asn_from_list as afl # Tools for creating association files
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase # Definition of a Lvl2 association file
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base # Definition of a Lvl3 association file
from astropy.io import fits


def runspec2(filename, output):
    if not os.path.exists(output):
        os.mkdir(output)
    spec2 = Spec2Pipeline()

    spec2.output_dir = output

    spec2.assign_wcs.skip = False  # assign world coordinate system
    spec2.bkg_subtract.skip = True
    spec2.flat_field.skip = False  # every pixel has same sensitivity
    spec2.srctype.skip = False  # source type
    spec2.straylight.skip = False  # scattered light inside detector corrected for in spatial dim
    spec2.fringe.skip = False  # coherent intereference inside detector in wvl dim, approximation
    spec2.photom.skip = False  # DN/s to MJy/sr
    spec2.residual_fringe.skip = True  # improve the fringe correction, fitting sin/cos
    spec2.cube_build.skip = True  # build cubes (but before dithering)
    spec2.extract_1d.skip = True
    spec2.save_results = True
    spec2(filename)


# Define a useful function to write out a Lvl3 association file from an input list
# Note that any background exposures will have to be of type x1d, so we'll convert any _cal.fits that were backgrounds to look for the corresponding x1d.
def writel3asn(files, asnfile, prodname, **kwargs):
    # Read the headers of each file to identify dedicated bg
    nfiles = len(files)
    scifiles = files
    # Define the basic association of science files
    asn = afl.asn_from_list(scifiles, rule=DMS_Level3_Base, product_name=prodname)
    # Write the association to a json file
    _, serialized = asn.dump()
    with open(asnfile, 'w') as outfile:
        outfile.write(serialized)


def runspec3(filename, out_dir):
    # This initial setup is just to make sure that we get the latest parameter reference files
    # pulled in for our files.  This is a temporary workaround to get around an issue with
    # how this pipeline calling method works.
    crds_config = Spec3Pipeline.get_config_from_reference(filename)
    spec3 = Spec3Pipeline.from_config_section(crds_config)

    spec3.output_dir = out_dir
    spec3.save_results = True

    # Cube building configuration
    spec3.cube_build.weighting = 'drizzle'  # 'emsm' or 'drizzle' #algorithm interpolating from point cloud to grid cube
    spec3.cube_build.coord_system = 'skyalign'  # 'ifualign', 'skyalign', or 'internal_cal' #which coordinate system for cube (often, ifualign (orthogonal to IFU), skyalign (in alpha dec, orthogonal to sky coord.))

    spec3.assign_mtwcs.skip = False  # world coord system to mosaic
    spec3.master_background.skip = True  # no master background available (otherwise in writel3asn)
    spec3.outlier_detection.skip = False  # bad pixels (stable), hot pixels (time dependent), cosmic ray showers
    spec3.mrs_imatch.skip = True  # background gets matched
    spec3.cube_build.skip = False  # build cube
    spec3.extract_1d.skip = False  # average of pixels as 1d spectrum (fast)
    # new parameters
    # spec3.extract_1d.ifu_autocen = True #, autocenter of circle where to take the mean
    spec3.extract_1d.center_xy = 24, 29
    spec3.extract_1d.ifu_rfcorr = True  # , residual fringe correction instead in spec2
    spec3.extract_1d.subtract_background = False  # , take ring around as background and subtract, only do this the first time
    spec3.extract_1d.ifu_rscale = 1  # set number of FWHMs fro radius
    spec3(filename)


def checkdata()
    data_str = [d for d in glob.glob(output_cubes + "*_x1d.fits")]
    data_str = sorted(data_str)
    # data = fits.getdata('/Users/helenakuehnle/Dateien/PhD/Analysis/WISE_J0855/pipelined/stage3/cubes/Level3_ch2-short_x1d.fits')
    data_all = []

    for i in range(len(data_str)):
        data = data_str[i]
        print(data)
        data = fits.getdata(data)
        dataz = np.array(list(zip(*data)))
        data_all.append(dataz)

    data_str3d = [d for d in glob.glob(output_cubes + "*_s3d.fits")]
    data_str3d = sorted(data_str3d)
    # data = fits.getdata('/Users/helenakuehnle/Dateien/PhD/Analysis/WISE_J0855/pipelined/stage3/cubes/Level3_ch2-short_x1d.fits')
    data_all3d = []

    for i in range(len(data_str3d)):
        data = data_str3d[i]
        print(data)
        data = fits.getdata(data)
        dataz = np.array(list(zip(*data)))
        data_all3d.append(dataz)

    # plot the model
    plt.figure(figsize=(12, 4))
    for i in range(len(data_str)):
        plt.plot(data_all[i][0], data_all[i][1])

    plt.plot(data_nirspec_prism[0], data_nirspec_prism[1], color='k')
    plt.plot(data_nirspec_395m[0], data_nirspec_395m[1], color='r')
    plt.ylim([-0.0001, 0.002])
    plt.xlim([3, 24])
    # plt.plot(x,y,linewidth=0.5,color='k')
    plt.xlabel('wavelength [mum]')
    plt.ylabel('Flux [Jy]')

if __name__ == '__main__':
    input_dir = '/Users/helenakuehnle/Dateien/PhD/Analysis/WISE_J0855/'

    # Point to where you want the output science results to go
    output_dir = os.path.join(input_dir, 'pipelined/')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    det1_dir = os.path.join(output_dir, 'stage1/')  # Detector1 pipeline outputs will go here
    spec2_dir = os.path.join(output_dir, 'stage2/')  # Spec2 pipeline outputs will go here
    spec3_dir = os.path.join(output_dir, 'stage3/')  # Spec3 pipeline outputs will go here

    # We need to check that the desired output directories exist, and if not create them
    if not os.path.exists(det1_dir):
        os.makedirs(det1_dir)
    if not os.path.exists(spec2_dir):
        os.makedirs(spec2_dir)
    if not os.path.exists(spec3_dir):
        os.makedirs(spec3_dir)

    # reordering of files
    files = [f for f in glob.glob(input_dir + "*/*rate.fits")]
    files = sorted(files)
    print(f"Found {len(files)} files to process")

    # run stage 2
    for f in files:
        runspec2(f, spec2_dir + f.split("/")[-2] + "/")

    # Find and sort all of the input files
    output_cubes = spec3_dir + "cubes_man1b/"
    if not os.path.exists(output_cubes):
        os.mkdir(output_cubes)

    # run stage 3
    for folder in ["obs5_MRS_MEDIUM"]:  # ["obs5_MRS_SHORT", "obs5_MRS_MEDIUM", "obs5_MRS_LONG"]:
        sstring = spec2_dir + folder + '/*calnodsub.fits'
        calfiles = np.array(sorted(glob.glob(sstring)))
        print('Found ' + str(len(calfiles)) + f' input files to process for folder {folder}')
        # Make an association file that includes all of the different exposures
        asnfile = os.path.join(spec2_dir + folder, 'l3asn.json')
        writel3asn(calfiles, asnfile, 'Level3')
        runspec3(asnfile, output_cubes)
