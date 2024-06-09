# -*- coding: utf-8 -*-
"""
    v1 201901, Dr. Jie Zheng, Beijing & Xinglong, NAOC
    v2 202101, Dr. Jie Zheng & Dr./Prof. Linqiao Jiang
    v3 202201, Zheng & Jiang
    v4 202304, Upgrade, restructure, Zheng & Jiang
    Quick_Light_Curve_Pipeline
"""


import os
import numpy as np
import astropy.io.fits as fits
from qmatch import match2d
from .u_conf import config, workmode
from .u_log import init_logger
from .u_utils import loadlist, rm_ix, zenum, pkl_load, pkl_dump, meanclip, fnbase


def pick(
        conf:config,
        raw_dir:str,
        red_dir:str,
        obj:str,
        band:str,
        mode:workmode=workmode(),
) -> tuple[list, list, list, list]:
    """
    try to guess the variables and ref stars from general catalog
    match stars, create a star cube, pick the best ones
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param obj: object
    :param band: band
    :param mode: input files missing or output existence mode
    :returns: a tuple with 4 lists
    """
    logf = init_logger("pick", f"{red_dir}/log/pick.log", conf)
    mode.reset_append(workmode.EXIST_OVER)

    # list file, and load list
    listfile = f"{red_dir}/lst/{obj}_{band}.lst"
    if mode.missing(listfile, f"{obj} {band} list", logf):
        return [None]*4
    raw_list = loadlist(listfile, base_path=raw_dir)
    cat_fits_list = loadlist(listfile, base_path=red_dir,
                        suffix="cat.fits", separate_folder=True)
    offset_pkl = f"{red_dir}/offset_{obj}_{band}.pkl"
    pick_pkl = f"{red_dir}/pick_{obj}_{band}.pkl"
    pick_txt = f"{red_dir}/pick_{obj}_{band}.txt"

    # check file exists
    if mode.missing(offset_pkl, f"offset {obj} {band}", logf):
        return [None]*4
    # if mode.exists(pick_pkl, f"picking result catalog {obj} {band}", logf):
        # return [None]*4

    # check file missing
    ix = []
    for f, (f,) in zenum(cat_fits_list):
        if mode.missing(f, "image catalog", logf):
            ix.append(f)
    # remove missing file
    rm_ix(ix, raw_list, cat_fits_list)
    nf = len(cat_fits_list)

    if nf == 0:
        logf.info(f"SKIP {obj} {band} No File")

    ###############################################################################

    # load offset result, and transfer to dict
    _, offset_x, offset_y, raw_list = pkl_load(offset_pkl)
    offset_x = dict(zip(raw_list, offset_x))
    offset_y = dict(zip(raw_list, offset_y))
    fn_len_max = max([len(f) for f in raw_list])

    nx = fits.getval(cat_fits_list[0], "IMNAXIS1")
    ny = fits.getval(cat_fits_list[0], "IMNAXIS2")

    # load all catalogs, filter good stars
    cat_all = []
    # all star count and good star count
    cat_len = np.zeros(nf, int)
    cat_good = np.zeros(nf, int)
    for f, (catf,) in zenum(cat_fits_list):
        cat = fits.getdata(catf)
        cat_len[f] = len(cat)
        ix_good = np.where(cat["ErrAUTO"] < conf.pick_err_max)[0]
        cat_good[f] = len(ix_good)
        cat_all.append(cat[ix_good])

    # use image with most good stars as the base
    base_img = np.argmax(cat_good)
    n_good = cat_good[base_img]
    logf.info(f"Base image: No {base_img} ({n_good:4d})")

    # create the x, y, mag, cali-mag matrix
    magi = np.empty((nf, n_good), float)
    magi[:, :] = np.nan
    magc = magi.copy()
    ximg = magi.copy()
    yimg = magi.copy()
    # mag calibration const & std
    cali_cst = np.zeros(nf, float)
    cali_std = np.zeros(nf, float)
    # load base image data
    magi[base_img] = cat_all[base_img]["MagAUTO"]
    magc[base_img] = cat_all[base_img]["MagAUTO"]
    ximg[base_img] = cat_all[base_img]["X"] + offset_x[raw_list[base_img]]
    yimg[base_img] = cat_all[base_img]["Y"] + offset_y[raw_list[base_img]]

    # match with base image, and fill the matrix
    for f, (catf, rawf) in zenum(cat_fits_list, raw_list):
        if f == base_img:
            continue
        # align i-th image with base
        xf = cat_all[f]["X"] + offset_x[rawf]
        yf = cat_all[f]["Y"] + offset_y[rawf]
        kbase, ki = match2d(ximg[base_img], yimg[base_img], xf, yf, conf.match_max_dis)
        magi[f, kbase] = cat_all[f]["MagAUTO"][ki]
        ximg[f, kbase] = xf[ki]
        yimg[f, kbase] = yf[ki]
        cali_cst[f], cali_std[f] = meanclip(magi[f] - magi[base_img])
        magc[f] = magi[f] - cali_cst[f]
        logf.debug(f"Picking {f:3d}/{nf:3d}: "
                   f"N={cat_len[f]:4d}->{cat_good[f]:4d}  "
                   f"Cali-Const={cali_cst[f]:+6.3f}+-{cali_std[f]:5.3f} | "
                   f"{fnbase(catf)}")

    # the mean position of star throught all images
    x_mean = np.nanmean(ximg, axis=0)
    y_mean = np.nanmean(yimg, axis=0)

    # calc std of all stars, and then find the good enough stars
    # median and std of each star between all images
    magmed = np.nanmedian(magc, axis=0)
    magstd = np.nanstd(magc, axis=0)
    # diff between min and max of each star
    magdif = np.nanmax(magc, axis=0) - np.nanmin(magc, axis=0)
    # percent of bad (nan) of each star
    magbad = np.sum(np.isnan(magc), axis=0) / nf

    # pick variable stars, by mag std, and distance to center
    ix_var = np.where((magstd > conf.pick_var_std) &
                      (np.abs(x_mean / nx - 0.5) < conf.pick_var_rad) &
                      (np.abs(y_mean / ny - 0.5) < conf.pick_var_rad) &
                      (magbad < conf.pick_bad_max) )[0]

    # pick reference stars, by a error limit or a number limit or both
    ix_ref = np.where((magstd < conf.pick_ref_std) &
                      (magdif < conf.pick_ref_dif) &
                      (magbad < conf.pick_bad_max) )[0]
    ix_ref = ix_ref[np.argsort(magstd[ix_ref])][:conf.pick_ref_n]

    logf.info(f"Pick {len(ix_ref)} ref stars and {len(ix_var)} var stars")
    # print out the result
    for i, k in enumerate(ix_var):
        logf.info(f"  VAR {i:2d}: [{k:3d}] ({x_mean[k]:6.1f} {y_mean[k]:6.1f})"
                  f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f} !{magbad[k]:4.2f}%")
    for i, k in enumerate(ix_ref):
        logf.info(f"  REF {i:2d}: [{k:3d}] ({x_mean[k]:6.1f} {y_mean[k]:6.1f})"
                  f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f} !{magbad[k]:4.2f}%")

    # save the result to array
    xy_ref = [(x_mean[k], y_mean[k]) for k in ix_ref]
    xy_var = [(x_mean[k], y_mean[k]) for k in ix_var]

    if len(xy_ref) == 0 or len(xy_var) == 0:
        logf.error("No good star found!")
        return None, None, None, None
    elif len(xy_var) >= 0 and len(xy_ref) == 0:
        xy = xy_var
    elif len(xy_var) == 0 and len(xy_ref) >= 0:
        xy = xy_ref
    else:
        xy = np.vstack([xy_ref, xy_var])
    ind_var = np.arange(len(xy_var), dtype=int)
    ind_ref = np.arange(len(xy_var), len(xy), dtype=int)
    ind_chk = np.arange(len(xy_var), len(xy), dtype=int)

    # save all stars to txt file
    with open(pick_txt[:-3]+"test.txt", "w") as ff:
        for k in range(n_good):
            ff.write(f"[{k:3d}] ({x_mean[k]:6.1f} {y_mean[k]:6.1f})"
                     f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f} !{magbad[k]:4.2f}%\n")

    # save the result to text file
    with open(pick_txt, "w") as ff:
        for f, k in enumerate(ix_var):
            ff.write(f"VAR {f:2d}: [{k:3d}] ({x_mean[k]:6.1f} {y_mean[k]:6.1f})"
                     f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f} !{magbad[k]:4.2f}%\n")
        for f, k in enumerate(ix_ref):
            ff.write(f"REF {f:2d}: [{k:3d}] ({x_mean[k]:6.1f} {y_mean[k]:6.1f})"
                     f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f} !{magbad[k]:4.2f}%\n")
    # dump result to pickle file
    pkl_dump(pick_pkl, magi, magc, ximg, yimg, xy, ind_ref, ind_var, ind_chk)

    return xy, ind_var, ind_ref, ind_chk

