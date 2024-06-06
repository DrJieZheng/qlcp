# -*- coding: utf-8 -*-
"""
    v1 201901, Dr. Jie Zheng, Beijing & Xinglong, NAOC
    v2 202101, Dr. Jie Zheng & Dr./Prof. Linqiao Jiang
    v3 202201, Zheng & Jiang
    v4 202304, Upgrade, restructure, Zheng & Jiang
    Quick_Light_Curve_Pipeline
"""


# import os
import numpy as np
import astropy.io.fits as fits
from .u_conf import config, workmode
from .u_log import init_logger
from .u_utils import loadlist, rm_ix, uttimestr


def biascomb(
        conf:config,
        raw_dir:str,
        red_dir:str,
        mode:workmode=workmode(workmode.MISS_SKIP+workmode.EXIST_APPEND),
):
    """
    bias combine
    Combine bias frames to a master bias frame, by median.
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param mode: input files missing or output existence mode
    :returns: Nothing
    """
    logf = init_logger("biascomb", f"{red_dir}/log/biascomb.log", conf)
    mode.reset_append(workmode.EXIST_OVER)

    # list file, and load list
    listfile = f"{red_dir}/lst/bias.lst"
    if mode.missing(listfile, "bias list", logf):
        return
    bias_list = loadlist(listfile, base_path=raw_dir)

    # output file
    bias_fits = f"{red_dir}/bias.fits"

    # check file exists
    if mode.exists(bias_fits, "master bias", logf):
        return

    # check file missing
    ix = []
    for i, f in enumerate(bias_list):
        if mode.missing(f, "raw bias", logf):
            ix.append(i)
    # remove missing file
    rm_ix(ix, bias_list)
    nf = len(bias_list)

    ###############################################################################

    # get size of images
    hdr = fits.getheader(bias_list[0])
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    logf.debug(f"{nf:02d} bias files, image sizes {nx:4d}x{ny:4d}")

    # load images
    data_cube = np.empty((nf, ny, nx), dtype=np.float32)
    for i, f in enumerate(bias_list):
        logf.debug(f"Loading {i+1:02d}/{nf:02d}: {f:40}")
        data_cube[i, :, :] = fits.getdata(f)

    # get median
    data_med = np.median(data_cube, axis=0)

    # add process time to header
    hdr.append(('COMBTIME', uttimestr()))
    _ = hdr.tostring()  # force check the header

    # save new fits
    fits.writeto(bias_fits, data=data_med, header=hdr, overwrite=True)
    logf.info(f"Writing to: {bias_fits}")
