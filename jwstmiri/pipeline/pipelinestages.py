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


def pipeline(input_dir, pid, obs, res, det, stages, output_dir, check, cores): # define in and output directory (strings), which stages (list with strings included) should be preocessed, which pid, obs, resolution (MRS or LRS or IMA) and detector is going to be processed (strings), and whether the data should be loaded and plotted as a check (bool)
    print(f"input_dir: {input_dir}, output_dir: {output_dir}, pid: {pid}, obseravtion: {obs}, resolution: {res}, detector: {det}, stages: {stages}, load data and plot: {check}, number of cores: {cores}")
    # TODO: add to run the code in parallel of n cores, specify cores as input in function
    # Point to where you want the output science results to go
    output_dir = os.path.join(output_dir, 'pipelined/')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # We need to check that the desired output directories exist, and if not create them
    det1_dir = os.path.join(output_dir, 'stage1/')  # Detector1 pipeline outputs will go here
    spec2_dir = os.path.join(output_dir, 'stage2/')  # Spec2 pipeline outputs will go here
    spec3_dir = os.path.join(output_dir, 'stage3/')  # Spec3 pipeline outputs will go here

    if 1 in stages:
        if not os.path.exists(det1_dir):
            os.makedirs(det1_dir)
    elif 2 in stages:
        if not os.path.exists(spec2_dir):
            os.makedirs(spec2_dir)
    elif 3 in stages:
        if not os.path.exists(spec3_dir):
            os.makedirs(spec3_dir)


    folders = []
    for o in obs:
        for r in res:
            for d in det:
                folders.append(f"obs{o}_{r}_{d}") # ["obs5_MRS_SHORT", "obs5_MRS_MEDIUM", "obs5_MRS_LONG"]
    print(folders)
    print(f"Use {len(folders)} to process")

    # reordering of files
    files = []
    for fol in folders:
        files.append([f for f in glob.glob(input_dir + fol + "/*rate.fits")])
    files = [item for sublist in files for item in sublist]
    files = sorted(files)
    print(files)
    print(f"Found {len(files)} files to process")

    stages_post = []

    # setup and run stage 2
    if 2 in stages:
        print(f"Start stage 2 processing")
        setupstage2(files, spec2_dir, folders)
        stages_post.append(2)
    else:
        print(f"Skip stage 2")

    # setup and run stage 3
    if (3 in stages) and (os.path.exists(spec2_dir)):
        print(f"Start stage 3 processing")
        out3 = setupstage3(spec2_dir, spec3_dir, folders)
        stages_post.append(3)
    else:
        print(f"Skip stage 3")

    # finish
    print(f"Data reduction of stages {stages_post} successful.")

    if check:
        checkdata(out3)


def setupstage2(files, spec2_dir, folders):
    # run stage 2 for all detector files
    for f in files:
        output2 = spec2_dir + f.split("/")[-2] + "/"
        # run stage 2
        runspec2(f, output2)

    # do nod subtraction
    nodsubtraction(spec2_dir, folders)


def setupstage3(spec2_dir, spec3_dir, folders):
    # Find and sort all of the input files
    output3 = spec3_dir + "cubes/"
    if not os.path.exists(output3):
        os.mkdir(output3)

    # run stage 3
    for folder in folders:
        sstring = spec2_dir + folder + '/*calnodsub.fits'
        calfiles = np.array(sorted(glob.glob(sstring)))
        print('Found ' + str(len(calfiles)) + f' input files to process for folder {folder}')

        # Make an association file that includes all of the different exposures
        asnfile = os.path.join(spec2_dir + folder, 'l3asn.json')
        writel3asn(calfiles, asnfile, 'Level3')
        runspec3(asnfile, output3)
    return output3


def runspec2(filename, output): # TODO: add dictionary with all parameters as inputs into the pipeline
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



def runspec3(filename, out_dir): # TODO: add dictionary with all parameters as inputs into the pipeline
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
    # TODO: implement switch between channels or rerun entire stage and replace ceratin files: use autocen for default, need to correct e.g. channel 1b manually as flux is to small
    spec3.extract_1d.ifu_autocen = True #, autocenter of circle where to take the mean
    #spec3.extract_1d.center_xy = 24, 29
    spec3.extract_1d.ifu_rfcorr = True  # , residual fringe correction instead in spec2
    spec3.extract_1d.subtract_background = False  # , take ring around as background and subtract, only do this the first time
    spec3.extract_1d.ifu_rscale = 1  # set number of FWHMs fro radius
    spec3(filename)

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

def nodsubtraction(spec2_dir, folders):
    # nod subtraction on all channels, improves channel 4B and 4C the most
    for folder in folders:
        files = [f for f in glob.glob(spec2_dir + folder + "/*cal.fits")]
        files = sorted(files)
        print(f"Found {len(files)} files to process")
        # find pairs
        filesshort = []
        fileslong = []
        for f in files:
            if f.split("_")[-2] == "mirifushort":
                filesshort.append(f)
            else:
                fileslong.append(f)

        fnew1s = filesshort[0].replace("cal.fits", "calnodsub.fits")
        fnew2s = filesshort[1].replace("cal.fits", "calnodsub.fits")

        im1 = jwst.datamodels.open(filesshort[0])
        im2 = jwst.datamodels.open(filesshort[1])

        temp1 = im1.data - im2.data
        temp2 = im2.data - im1.data

        im1.data = temp1
        im2.data = temp2

        im1.to_fits(fnew1s, overwrite=True)
        im2.to_fits(fnew2s, overwrite=True)

        fnew1l = fileslong[0].replace("cal.fits", "calnodsub.fits")
        fnew2l = fileslong[1].replace("cal.fits", "calnodsub.fits")

        im1 = jwst.datamodels.open(fileslong[0])
        im2 = jwst.datamodels.open(fileslong[1])

        temp1 = im1.data - im2.data
        temp2 = im2.data - im1.data

        im1.data = temp1
        im2.data = temp2

        im1.to_fits(fnew1l, overwrite=True)
        im2.to_fits(fnew2l, overwrite=True)

def checkdata(output_cubes):
    data_str = [d for d in glob.glob(output_cubes + "*_x1d.fits")]
    data_str = sorted(data_str)
    data_all = []

    for i in range(len(data_str)):
        data = data_str[i]
        print(data)
        data = fits.getdata(data)
        dataz = np.array(list(zip(*data)))
        data_all.append(dataz)

    # data_str3d = [d for d in glob.glob(output_cubes + "*_s3d.fits")]
    # data_str3d = sorted(data_str3d)
    # data_all3d = []
    #
    # for i in range(len(data_str3d)):
    #     data = data_str3d[i]
    #     print(data)
    #     data = fits.getdata(data)
    #     dataz = np.array(list(zip(*data)))
    #     data_all3d.append(dataz)

    # plot the model
    plt.figure(figsize=(12, 4))
    for i in range(len(data_str)):
        plt.plot(data_all[i][0], data_all[i][1])

    plt.ylim([-0.0001, 0.002])
    plt.xlim([3, 24])
    # plt.plot(x,y,linewidth=0.5,color='k')
    plt.xlabel('wavelength [mum]')
    plt.ylabel('Flux [Jy]')
    plt.show()

