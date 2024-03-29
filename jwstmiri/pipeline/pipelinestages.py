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


def pipeline(input_dir, pid, obs, res, det, stages, input_vars, output_dir, check, cores): # define in and output directory (strings), which stages (list with strings included) should be preocessed, which pid, obs, resolution (MRS or LRS or IMA) and detector is going to be processed (strings), and whether the data should be loaded and plotted as a check (bool)
    print(f"input_dir: {input_dir}, output_dir: {output_dir}, pid: {pid}, obseravtion: {obs}, resolution: {res}, detector: {det}, stages: {stages}, input variables: {input_vars}, load data and plot: {check}, number of cores: {cores}")
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
    if 2 in stages:
        if not os.path.exists(spec2_dir):
            os.makedirs(spec2_dir)
    if 3 in stages:
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
        #for 1275, obs4 test if showercorrected files show less offsets
        #files.append([f for f in glob.glob(input_dir + fol + "/*rate_showers_masked.fits")])
    files = [item for sublist in files for item in sublist]
    files = sorted(files)
    print(files)
    print(f"Found {len(files)} files to process")

    stages_post = []

    # setup and run stage 2
    if 2 in stages:
        print(f"Start stage 2 processing")
        setupstage2(files, spec2_dir, folders, input_vars)
        stages_post.append(2)
    else:
        print(f"Skip stage 2")

    # setup and run stage 3
    if (3 in stages) and (os.path.exists(spec2_dir)):
        print(f"Start stage 3 processing")
        out3 = setupstage3(spec2_dir, spec3_dir, folders, input_vars)
        stages_post.append(3)
    else:
        print(f"Skip stage 3")

    # finish
    print(f"Data reduction of stages {stages_post} successful.")

    if check:
        checkdata(out3)


def setupstage2(files, spec2_dir, folders, input_vars):
    # run stage 2 for all detector files
    for f in files:
        output2 = spec2_dir + f.split("/")[-2] + "/"
      #   run stage 2
        runspec2(f, output2, input_vars)

    # do nod subtraction
    if input_vars['stage2']['nodsub']:
        nodsubtraction(spec2_dir, folders)


def setupstage3(spec2_dir, spec3_dir, folders, input_vars):
    # Find and sort all of the input files
    output3 = spec3_dir + "cubes/"
    if not os.path.exists(output3):
        os.mkdir(output3)

    # run stage 3
    for folder in folders:
        if input_vars['stage2']['nodsub']:
            sstring = spec2_dir + folder + '/*calnodsub.fits'
        else:
            #test if this changes the spectrum significantly for 1275, obs4
            sstring = spec2_dir + folder + '/*cal.fits'
        calfiles = np.array(sorted(glob.glob(sstring)))
        print('Found ' + str(len(calfiles)) + f' input files to process for folder {folder}')

        # Make an association file that includes all of the different exposures
        asnfile = os.path.join(spec2_dir + folder, 'l3asn.json')
        writel3asn(calfiles, asnfile, 'Level3')
        runspec3(asnfile, output3, input_vars)

        man_bckg = False #try manual background subtraction for pid1278 obs40 -> not better, ignored
        if man_bckg:
            print('#############################')
            print('#### add manual background subtraction')
            print('#############################')
            man_bckg_rm(spec3_dir, folder, input_vars)
            print('manual background subtraction done')
    return output3


def runspec2(filename, output, input_vars):
    if not os.path.exists(output):
        os.mkdir(output)
    spec2 = Spec2Pipeline()

    spec2.output_dir = output

    vars = input_vars['stage2']
    spec2.assign_wcs.skip = vars['assign_wcs_skip'] #False  # assign world coordinate system
    spec2.bkg_subtract.skip = vars['bkg_subtract_skip'] #True
    spec2.flat_field.skip = vars['flat_field_skip'] #False  # every pixel has same sensitivity
    spec2.srctype.skip = vars['srctype_skip'] #False  # source type
    spec2.straylight.skip = vars['straylight_skip'] #False  # scattered light inside detector corrected for in spatial dim
    spec2.fringe.skip = vars['fringe_skip'] #False  # coherent intereference inside detector in wvl dim, approximation
    spec2.photom.skip = vars['photom_skip'] #False  # DN/s to MJy/sr
    spec2.residual_fringe.skip = vars['residual_fringe_skip'] #True  # improve the fringe correction, fitting sin/cos
    spec2.cube_build.skip = vars['cube_build_skip'] #True  # build cubes (but before dithering)
    spec2.extract_1d.skip = vars['extract_1d_skip'] #True
    spec2.save_results = vars['save_results'] #True
    spec2(filename)



def runspec3(filename, out_dir, input_vars):
    # This initial setup is just to make sure that we get the latest parameter reference files
    # pulled in for our files.  This is a temporary workaround to get around an issue with
    # how this pipeline calling method works.
    crds_config = Spec3Pipeline.get_config_from_reference(filename)
    spec3 = Spec3Pipeline.from_config_section(crds_config)

    vars = input_vars['stage3']

    if vars['ifu_autocen']:
        spec3.output_dir = out_dir
    else:
        out_alt = os.path.join(out_dir, 'manualcent_'+str(vars['center_x'])+'_'+str(vars['center_y'])+'/')
        if not os.path.exists(out_alt):
            os.mkdir(out_alt)
        spec3.output_dir = out_alt

    spec3.save_results = vars['save_results'] #True

    # Cube building configuration
    spec3.cube_build.weighting = vars['weighting'] #'drizzle'  # 'emsm' or 'drizzle' #algorithm interpolating from point cloud to grid cube
    spec3.cube_build.coord_system = vars['coord_system'] #'skyalign'  # 'ifualign', 'skyalign', or 'internal_cal' #which coordinate system for cube (often, ifualign (orthogonal to IFU), skyalign (in alpha dec, orthogonal to sky coord.))

    spec3.assign_mtwcs.skip = vars['assign_mtwcs_skip'] #False  # world coord system to mosaic
    spec3.master_background.skip = vars['master_background_skip'] #True  # no master background available (otherwise in writel3asn)
    spec3.outlier_detection.skip = vars['outlier_detection_skip'] #False  # bad pixels (stable), hot pixels (time dependent), cosmic ray showers
    spec3.mrs_imatch.skip = vars['mrs_imatch_skip'] #True  # background gets matched
    spec3.cube_build.skip = vars['cube_build_skip'] #False  # build cube
    spec3.extract_1d.skip = vars['extract_1d_skip'] #False  # average of pixels as 1d spectrum (fast)
    if vars['ifu_autocen']:
        spec3.extract_1d.ifu_autocen = True #, autocenter of circle where to take the mean
    else:
        spec3.extract_1d.center_xy = vars['center_x'], vars['center_y'] #24, 29
    spec3.extract_1d.ifu_rfcorr = vars['ifu_rfcorr'] #True  # , residual fringe correction instead in spec2
    spec3.extract_1d.subtract_background = vars['subtract_background'] #False  # , take ring around as background and subtract, only do this the first time
    spec3.extract_1d.ifu_rscale = vars['ifu_rscale'] #1  # set number of FWHMs fro radius
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
        print('#############################')
        print('start nod subtraction')
        print('#############################')
        files = [f for f in glob.glob(spec2_dir + folder + "/*cal.fits")]
        files = sorted(files)
        print(f"Found {len(files)} files to process")
        # find pairs
        filesshort = []
        fileslong = []
        for f in files:
            # 1275, obs5 change input files, from -2 to -5
            if f.split("_")[-2] == "mirifushort":
            #if f.split("_")[-5] == "mirifushort":
                filesshort.append(f)
            else:
                fileslong.append(f)

        dither = len(files)/2

        if dither == 2:
            fnew1s = filesshort[0].replace("cal.fits", "calnodsub.fits")
            fnew2s = filesshort[1].replace("cal.fits", "calnodsub.fits")
            im1 = jwst.datamodels.open(filesshort[0])
            im2 = jwst.datamodels.open(filesshort[1])
        
            temp1 = im1.data - im2.data
            temp2 = im2.data - im1.data
        
        elif dither == 4:
            fnew1s = filesshort[0].replace("cal.fits", "calnodsub.fits")
            fnew2s = filesshort[1].replace("cal.fits", "calnodsub.fits")
            fnew3s = filesshort[2].replace("cal.fits", "calnodsub.fits")
            fnew4s = filesshort[3].replace("cal.fits", "calnodsub.fits")
            im1 = jwst.datamodels.open(filesshort[0])
            im2 = jwst.datamodels.open(filesshort[1])
            im3 = jwst.datamodels.open(filesshort[2])
            im4 = jwst.datamodels.open(filesshort[3])

            temp1 = im1.data - im2.data
            temp2 = im2.data - im1.data
            temp3 = im3.data - im4.data
            temp4 = im4.data - im3.data

        # for 1189obs16, add manual background subtraction for channels 1a and 2a 
        man_bgs = False
        if man_bgs:
            print('#############################')
            print('manual background subtraction')
            print('#############################')
            #if folder.split("_")[-1] == "SHORT":
            bkg1 = np.nanmedian(temp1[800:1000, :516], axis=0)
            bkg2 = np.nanmedian(temp1[20:200, 516:], axis=0)
            temp1[:, :516] = np.subtract(temp1[:, :516], bkg1)
            temp1[:, 516:] = np.subtract(temp1[:, 516:], bkg2)

            bkg1 = np.nanmedian(temp2[800:1000, :516], axis=0)
            bkg2 = np.nanmedian(temp2[20:200, 516:], axis=0)
            temp2[:, :516] = np.subtract(temp2[:, :516], bkg1)
            temp2[:, 516:] = np.subtract(temp2[:, 516:], bkg2)

        
        if dither == 2:
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
        
        elif dither == 4:
            im1.data = temp1
            im2.data = temp2
            im3.data = temp3
            im4.data = temp4

            im1.to_fits(fnew1s, overwrite=True)
            im2.to_fits(fnew2s, overwrite=True)
            im3.to_fits(fnew3s, overwrite=True)
            im4.to_fits(fnew4s, overwrite=True)

            fnew1l = fileslong[0].replace("cal.fits", "calnodsub.fits")
            fnew2l = fileslong[1].replace("cal.fits", "calnodsub.fits")
            fnew3l = fileslong[2].replace("cal.fits", "calnodsub.fits")
            fnew4l = fileslong[3].replace("cal.fits", "calnodsub.fits")

            im1 = jwst.datamodels.open(fileslong[0])
            im2 = jwst.datamodels.open(fileslong[1])
            im2 = jwst.datamodels.open(fileslong[2])
            im3 = jwst.datamodels.open(fileslong[3])

            temp1 = im1.data - im2.data
            temp2 = im2.data - im1.data
            temp3 = im3.data - im4.data
            temp4 = im4.data - im3.data

            im1.data = temp1
            im2.data = temp2
            im3.data = temp3
            im4.data = temp4

            im1.to_fits(fnew1l, overwrite=True)
            im2.to_fits(fnew2l, overwrite=True)
            im3.to_fits(fnew3l, overwrite=True)
            im4.to_fits(fnew4l, overwrite=True)

def man_bckg_rm(spec3_dir, folder, input_vars): #this function only for pid 1278 obs40, only works with coord_system = ifualign as input
    
    print('#############################')
    print('start manual background subtraction')
    print('#############################')
    files = [f for f in glob.glob(spec3_dir + "cubes/*3d.fits")]
    print(files)
    print(spec3_dir + "cubes/*3d.fits")
    files = sorted(files)
    print(f"Found {len(files)} files to process")
        
    for f in files:
        da1 = jwst.datamodels.open(f)

        corr = np.median(da1[10:25,:,5:15],2)

        da1_corr = np.copy(da1)
        corr_3d = np.concatenate([[corr]] * np.shape(da1[0,0,:])[0], axis=0)
        corr_3d = np.transpose(corr_3d, (1, 2, 0))

        da1[10:25,:,:] = da1[10:25,:,:]-corr_3d
        da1.to_fits(da1, overwrite=True)
        
        print('#############################')
        print('start 1d extraction')
        print('#############################')

        out_dir = os.path.join(spec3_dir, folder)
        extract_1d(f, out_dir, input_vars)

    print('#############################')
    print('manual background subtraction done')
    print('#############################')
        
def extract_1d(file, out_dir, input_vars): #only extract 1d array from cube
    print('start 1d extraction')
    im = datamodels.IFUCubeModel(file)
    step = pipeline.calwebb_spec3.extract_1d_step.Extract1dStep()
    step.output_dir = out_dir
    
    vars = input_vars['stage3']
    step.save_results = vars['save_results']

    if vars['ifu_autocen']:
        step.ifu_autocen = True #, autocenter of circle where to take the mean
    else:
        step.center_xy = vars['center_x'], vars['center_y'] #24, 29
    
    step.ifu_rfcorr = vars['ifu_rfcorr'] #True  # , residual fringe correction instead in spec2
    step.subtract_background = vars['subtract_background'] #False  # , take ring around as background and subtract, only do this the first time
    step.ifu_rscale = vars['ifu_rscale'] #1  # set number of FWHMs fro radius
    
    print('run step')
    step(im)
    print('done')

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

