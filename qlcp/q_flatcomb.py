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


def flatcomb(
        conf:config,
        raw_dir:str,
        red_dir:str,
        band:str,
        use_bias:str=None,
        alt_bias:str=None,
        mode:workmode=workmode(),
):
    """
    flat combine
    Compose master flat from raw flats, by subtracting bias and norming by median
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param band: band of flats
    :param use_bias: alternative bias filename
    :param alt_bias: backup bias filename
    :param mode: input files missing or output existence mode
    :returns: Nothing
    """
    logf = init_logger("flatcomb", f"{red_dir}/log/flatcomb.log", conf)
    mode.reset_append(workmode.EXIST_OVER)

    # list file, and load list
    listfile = f"{red_dir}/lst/flat_{band}.lst"
    if mode.missing(listfile, f"{band}-flat list", logf):
        return
    flat_list = loadlist(listfile, base_path=raw_dir)

    # bias file, prefer use_bias, then bias today, then alt_bias
    bias_fits = use_bias if use_bias else f"{red_dir}/bias.fits"
    if workmode(workmode.MISS_SKIP).missing(bias_fits, "master bias", logf):
        bias_fits = alt_bias
    # output file
    flat_fits = f"{red_dir}/flat_{band}.fits"

    # check file exists
    if mode.exists(flat_fits, f"master {band}-flat", logf):
        return

    # check file missing
    if mode.missing(bias_fits, "master bias", logf):
        return
    ix = []
    for i, f in enumerate(flat_list):
        if mode.missing(f, "raw flat", logf):
            ix.append(i)
    # remove missing file
    rm_ix(ix, flat_list)
    nf = len(flat_list)

    ###############################################################################

    # get size of images
    hdr = fits.getheader(bias_fits)
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    logf.debug(f"{nf:02d} flat files, image sizes {nx:4d}x{ny:4d}")

    # load bias
    logf.debug(f"Loading Bias: {bias_fits}")
    data_bias = fits.getdata(bias_fits)

    # load images
    data_cube = np.empty((nf, ny, nx), dtype=np.float32)
    for i, f in enumerate(flat_list):
        data_tmp = fits.getdata(f) - data_bias
        data_tmp_med = np.median(data_tmp)
        # check flat level, if too low or too high, discard
        if conf.flat_limit_low < data_tmp_med < conf.flat_limit_high:
            data_cube[i, :, :] = data_tmp / data_tmp_med
            logf.debug(f"Loading {i+1:02d}/{nf:02d}: {f:40s} / Scaled by {data_tmp_med:7.1f}")
        else:
            data_cube[f, :, :] = np.nan
            logf.warning(f"Ignore  {i+1:02d}/{nf:02d}: {f:40s} / XXX MED = {data_tmp_med:7.1f}")

    # get median
    data_med = np.median(data_cube, axis=0)

    # add process time to header
    hdr.append(('COMBTIME', uttimestr()))
    _ = hdr.tostring()  # force check the header

    # save new fits
    fits.writeto(flat_fits, data=data_med, header=hdr, overwrite=True)
    logf.info(f"Writing to: {flat_fits}")
