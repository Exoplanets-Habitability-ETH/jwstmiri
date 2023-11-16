############
# Define the pipeline stages beginning from stage 2
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