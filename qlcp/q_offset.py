# -*- coding: utf-8 -*-
"""
    v1 201901, Dr. Jie Zheng, Beijing & Xinglong, NAOC
    v2 202101, Dr. Jie Zheng & Dr./Prof. Linqiao Jiang
    v3 202201, Zheng & Jiang
    v4 202304, Upgrade, restructure, Zheng & Jiang
    Quick_Light_Curve_Pipeline
"""


import numpy as np
import astropy.io.fits as fits
from qmatch import mean_xy, mean_offset1d
import matplotlib.pyplot as plt
from .u_conf import config, workmode
from .u_log import init_logger
from .u_utils import loadlist, rm_ix, hdr_dt, zenum, str2mjd, pkl_dump


def offset(
        conf:config,
        raw_dir:str,
        red_dir:str,
        obj:str,
        band:str,
        base_img:int|str=0,
        mode:workmode=workmode(),
):
    """
    Calculate offset of images
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param obj: object
    :param band: band
    :param base_img: the offset base index or filename
    :param mode: input files missing or output existence mode
    :returns: Nothing
    """
    logf = init_logger("offset", f"{red_dir}/log/offset.log", conf)
    mode.reset_append(workmode.EXIST_OVER)

    # list file, and load list
    listfile = f"{red_dir}/lst/{obj}_{band}.lst"
    if mode.missing(listfile, f"{obj} {band} list", logf):
        return
    raw_list = loadlist(listfile, base_path=raw_dir)
    offset_pkl = f"{red_dir}/offset_{obj}_{band}.pkl"
    offset_txt = f"{red_dir}/offset_{obj}_{band}.txt"
    offset_png = f"{red_dir}/offset_{obj}_{band}.png"

    # check file exists
    if mode.exists(offset_pkl, f"offset {obj} {band}", logf):
        return

    # check file missing
    ix = []
    for i, (f,) in zenum(raw_list):
        if mode.missing(f, "raw image", logf):
            ix.append(i)
    # remove missing file
    rm_ix(ix, raw_list)
    nf = len(raw_list)

    if nf == 0:
        logf.info(f"SKIP {obj} {band} Nothing")
        return

    # base image, type check, range check, existance check
    if isinstance(base_img, int):
        if 0 > base_img or base_img >= len(raw_list):
            base_img = 0
        base_img = [base_img]
    elif not isinstance(base_img, str):
        base_img = raw_list[0]
    # if external file not found, use 0th
    # special, fixed mode
    if workmode(workmode.MISS_SKIP).missing(base_img, "offset base image", logf):
        base_img = raw_list[0]

    ###############################################################################

    logf.debug(f"{nf:3d} files")

    # load base image
    logf.debug(f"Loading base image: {base_img}")
    base_x, base_y = mean_xy(fits.getdata(base_img))

    # xy offset array
    offset_x = np.empty(nf, int)
    offset_y = np.empty(nf, int)
    obs_mjd = np.empty(nf)

    # load images and process
    for i, (rawf,) in zenum(raw_list):
        logf.debug(f"Loading {i+1:03d}/{nf:03d}: {rawf:40s}")

        # process data
        raw_x, raw_y = mean_xy(fits.getdata(rawf))
        offset_x[i] = int(mean_offset1d(base_x, raw_x, max_d=conf.offset_max_dis))
        offset_y[i] = int(mean_offset1d(base_y, raw_y, max_d=conf.offset_max_dis))

        # mjd of obs
        hdr = fits.getheader(rawf)
        obs_dt = hdr_dt(hdr)[:19]
        obs_mjd[i] = str2mjd(obs_dt) + hdr.get("EXPTIME", 0.0) / 2 / 86400

        logf.debug(f"{'':10s}{obs_mjd[i]:12.7f} | {offset_x[i]:+5d} {offset_y[i]:+5d}")

    # save new fits
    with open(offset_txt, "w") as ff:
        for d, x, y, rawf in zip(obs_mjd, offset_x, offset_y, raw_list):
            ff.write(f"{d:12.7f}  {x:+5d} {y:+5d}  {rawf}\n")
    pkl_dump(offset_pkl, obs_mjd, offset_x, offset_y, raw_list)
    logf.debug(f"Writing {offset_pkl}")

    # draw offset figure
    fig = plt.figure(figsize=(8, 8))
    ax_xy = fig.add_axes([0.05, 0.05, 0.60, 0.60])
    ax_xt = fig.add_axes([0.05, 0.70, 0.60, 0.25])
    ax_ty = fig.add_axes([0.70, 0.05, 0.25, 0.60])
    ax_xy.plot(offset_x, offset_y, "k.:")
    ax_xt.plot(offset_x, obs_mjd, "k.:")
    ax_ty.plot(obs_mjd, offset_y, "k.:")
    ax_xy.set_xlabel("X offset")
    ax_xy.set_ylabel("Y offset")
    ax_xt.set_ylabel("MJD")
    ax_ty.set_xlabel("MJD")
    ax_xt.set_title(f"Offset {obj} {band}")
    fig.savefig(offset_png)
