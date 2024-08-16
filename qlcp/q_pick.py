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
from astropy.stats import sigma_clipped_stats
from qmatch import match2d
from .u_conf import config
from .u_workmode import workmode
from .u_log import init_logger
from .u_utils import loadlist, rm_ix, zenum, pkl_load, pkl_dump,\
    meanclip, fnbase, tqdm_bar


def pick(
        conf:config,
        raw_dir:str,
        red_dir:str,
        obj:str,
        band:str,
        base_img:int|str=None,
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
    :param base_img: base image, if str (external file), treet as none
    :param mode: input files missing or output existence mode
    :returns: a tuple with 4 lists
    """
    logf = init_logger("pick", f"{red_dir}/log/pick.log", conf)
    mode.reset_append(workmode.EXIST_OVER)

    # list file, and load list
    listfile = f"{red_dir}/lst/{obj}_{band}.lst"
    if mode.missing(listfile, f"{obj} {band} list", logf):
        return [None]*4
    bf_fits_list = loadlist(listfile, base_path=red_dir,
                        suffix="bf.fits", separate_folder=True)
    cat_fits_list = loadlist(listfile, base_path=red_dir,
                        suffix="cat.fits", separate_folder=True)
    offset_pkl = f"{red_dir}/offset_{obj}_{band}.pkl"
    pick_pkl = f"{red_dir}/pick_{obj}_{band}.pkl"
    pick_txt = f"{red_dir}/pick_{obj}_{band}.txt"

    # check file exists
    if mode.missing(offset_pkl, f"offset {obj} {band}", logf):
        return [None]*4

    # check file missing
    mode.start_lazy()
    ix = []
    for f, (catf,) in zenum(cat_fits_list):
        if mode.missing(catf, "image catalog", None):
            ix.append(f)
    # remove missing file
    mode.end_lazy(logf)
    rm_ix(ix, bf_fits_list, cat_fits_list)
    nf = len(cat_fits_list)

    # base image, if str (external file), treet as none
    if isinstance(base_img, int):
        if 0 > base_img or base_img >= nf:
            base_img = 0
    else:
        base_img = 0
    
    # if load_last:
    #     # if load last result, check the file exists
    #     if os.path.isfile(pick_pkl):
    #         logf.info(f"LOAD saved picking result: {obj} {band}")
    #         _, _, _, _, xy, ind_ref, ind_var, ind_chk = pkl_load(pick_pkl)
    #         return xy, ind_ref, ind_var, ind_chk

    if nf == 0:
        logf.info(f"SKIP {obj} {band} No File")

    ###############################################################################

    # load offset result, and transfer to dict
    _, offset_x, offset_y, offset_bf_fits_list = pkl_load(offset_pkl)
    offset_bf_fits_list = [fnbase(f) for f in offset_bf_fits_list]
    offset_x = dict(zip(offset_bf_fits_list, offset_x))
    offset_y = dict(zip(offset_bf_fits_list, offset_y))
    fn_len_max = max(len(f) for f in offset_bf_fits_list)

    nx = fits.getval(bf_fits_list[0], "NAXIS1")
    ny = fits.getval(bf_fits_list[0], "NAXIS2")
    # nx  = ny = 2048

    # load base image catalog and choose good stars
    cat_base = fits.getdata(cat_fits_list[base_img])
    # ix_good = np.where(cat_base["ErrAUTO"] < conf.pick_err_max)[0]
    # ix_good = ix_good[np.argsort(cat_base["MagAUTO"][ix_good])][:conf.pick_star_n]
    ix_good = np.argsort(cat_base["ErrAUTO"])[:conf.pick_star_n]
    ix_good = ix_good[cat_base[ix_good]["ErrAUTO"] < conf.pick_err_max]
    cat_base = cat_base[ix_good]
    n_base_good = len(ix_good)
    # exctract x, y, mag
    bbff = fnbase(bf_fits_list[base_img])
    x_base = cat_base["X"] + offset_x[bbff]
    y_base = cat_base["Y"] + offset_y[bbff]
    m_base = cat_base["MagAUTO"]

    # create the x, y, mag, cali-mag matrix
    magi = np.empty((nf, n_base_good), float) + np.nan
    magc = magi.copy()
    ximg = magi.copy()
    yimg = magi.copy()
    # mag calibration const & std, and stars used
    cali_cst = np.zeros(nf, float)
    cali_std = np.zeros(nf, float)
    cali_n = np.zeros(nf, int)

    # match with base image, and fill the matrix
    pbar = tqdm_bar(nf, f"PICK {obj} {band}")
    for f, (catf, bff) in zenum(cat_fits_list, bf_fits_list):
        if f == base_img:
            magi[f] = m_base
            magc[f] = m_base
            ximg[f] = x_base
            yimg[f] = y_base
            cali_cst[f] = 0.0
            cali_std[f] = 0.0
            cali_n[f] = n_base_good
            n_good = n_match = n_base_good
        else:
            # align i-th image with base
            bbff = fnbase(bff)
            # load catalog and choose good stars
            cat_f = fits.getdata(catf)
            # ix_good = np.where(cat_f["ErrAUTO"] < conf.pick_err_max)[0]
            # ix_good = ix_good[np.argsort(cat_f["MagAUTO"][ix_good])][:conf.pick_star_n]
            ix_good = np.argsort(cat_f["ErrAUTO"])[:conf.pick_star_n]
            ix_good = ix_good[cat_f[ix_good]["ErrAUTO"] < conf.pick_err_max]
            # cat_f = cat_f[ix_good]
            n_good = len(ix_good)
            # extract x, y, mag
            xf = cat_f["X"] + offset_x[bbff]
            yf = cat_f["Y"] + offset_y[bbff]
            mf = cat_f["MagAUTO"]
            # match with base image
            ix_base, ix_f = match2d(x_base, y_base, xf, yf, conf.match_max_dis)
            n_match = len(ix_f)
            # fill the matrix
            magi[f, ix_base] = mf[ix_f]
            ximg[f, ix_base] = xf[ix_f]
            yimg[f, ix_base] = yf[ix_f]
            # calibration
            cali_cst[f], cali_std[f] = meanclip(mf[ix_f] - m_base[ix_base])
            cali_n[f] = sum(np.abs(mf[ix_f] - m_base[ix_base] - cali_cst[f]) <= 3 * cali_std[f])
            # if the count of calibrate stars is too small, then discard this image
            if cali_n[f] / n_base_good < conf.pick_bad_img:
                magc[f] = np.nan
            else:
                magc[f] = magi[f] - cali_cst[f]
        logf.debug(f"Picking {f:3d}/{nf:3d}: "
                f"N={n_good:4d}->{n_match:3d}  "
                f"Cali-Const={cali_cst[f]:+6.3f}+-{cali_std[f]:5.3f}/{cali_n[f]:3d} | "
                f"{bbff}")
        pbar.update(1)
    pbar.close()

    # the mean position of star throught all images
    x_mean = np.nanmean(ximg, axis=0)
    y_mean = np.nanmean(yimg, axis=0)
    x_std = np.nanstd(ximg, axis=0)
    y_std = np.nanstd(yimg, axis=0)

    # calc std of all stars, and then find the good enough stars
    # median and std of each star between all images
    magmed = np.nanmedian(magc, axis=0)
    magstd = np.nanstd(magc, axis=0)
    # diff between min and max of each star
    magdif = np.nanmax(magc, axis=0) - np.nanmin(magc, axis=0)
    # bad image count
    n_bad_image = sum(cali_n / n_base_good < conf.pick_bad_img)
    bad_image_id = ",".join(f"{i:d}" for i in np.where(cali_n / n_base_good < conf.pick_bad_img)[0])
    logf.debug(f"Bad image count: {n_bad_image}/{nf}: {bad_image_id}")
    # percent of bad (nan) of each star
    magbad = (np.sum(np.isnan(magc), axis=0) - n_bad_image) / nf

    # pick variable stars, by mag std, and distance to center
    ix_var = np.where((magstd > conf.pick_var_std) &
                      (np.abs(x_mean / nx - 0.5) < conf.pick_var_rad) &
                      (np.abs(y_mean / ny - 0.5) < conf.pick_var_rad) &
                      (magbad < conf.pick_bad_max) )[0]
    ix_var = ix_var[np.argsort(magstd[ix_var])][::-1][:conf.pick_var_n]

    # pick reference stars, by a error limit or a number limit or both
    ix_ref = np.where((magstd < conf.pick_ref_std) &
                      (magdif < conf.pick_ref_dif) &
                      (magbad < conf.pick_bad_max) )[0]
    ix_ref = ix_ref[np.argsort(magstd[ix_ref])][:conf.pick_ref_n]

    logf.info(f"Pick {len(ix_ref)} ref stars and {len(ix_var)} var stars {obj} {band}")
    # print out the result
    for i, k in enumerate(ix_var):
        logf.info(f"  VAR {i:2d}: [{k:3d}] ({x_mean[k]:6.1f} {y_mean[k]:6.1f})"
                  f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f} !{magbad[k]*100:4.1f}%")
    for i, k in enumerate(ix_ref):
        logf.info(f"  REF {i:2d}: [{k:3d}] ({x_mean[k]:6.1f} {y_mean[k]:6.1f})"
                  f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f} !{magbad[k]*100:4.1f}%")

    # save the result to array
    xy_ref = [(x_mean[k], y_mean[k]) for k in ix_ref]
    xy_var = [(x_mean[k], y_mean[k]) for k in ix_var]

    # combind the result, and make indice
    xy = xy_ref + xy_var
    ind_var = np.arange(len(xy_var), dtype=int)
    ind_ref = np.arange(len(xy_var), len(xy), dtype=int)
    ind_chk = np.arange(len(xy_var), len(xy), dtype=int)

    # save all stars to txt file
    with open(pick_txt[:-3]+"test.txt", "w") as ff:
        for k in range(n_base_good):
            ff.write("V " if k in ix_var else "R " if k in ix_ref else "  ")
            ff.write(f"[{k:3d}] "
                     f"({x_mean[k]:6.1f}~{x_std[k]:4.2f} {y_mean[k]:6.1f}~{y_std[k]:4.2f})"
                     f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f}/{magdif[k]:4.2f}"
                     f" !{magbad[k]*100:4.1f}%\n")

    # save the result to text file
    with open(pick_txt, "w") as ff:
        for f, k in enumerate(ix_var):
            ff.write(f"VAR {f:2d}: [{k:3d}] ({x_mean[k]:6.1f} {y_mean[k]:6.1f})"
                     f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f} !{magbad[k]*100:4.1f}%\n")
        for f, k in enumerate(ix_ref):
            ff.write(f"REF {f:2d}: [{k:3d}] ({x_mean[k]:6.1f} {y_mean[k]:6.1f})"
                     f"  {magmed[k]:5.2f}+-{magstd[k]:5.3f} !{magbad[k]*100:4.1f}%\n")
    # dump result to pickle file
    pkl_dump(pick_pkl, magi, magc, ximg, yimg, 
                cali_cst, cali_std, x_mean, x_std, y_mean, y_std,
                xy, ind_ref, ind_var, ind_chk,
                magmed, magstd, magdif, magbad)

    return xy, ind_var, ind_ref, ind_chk


def pick_last(
        conf:config,
        raw_dir:str,
        red_dir:str,
        obj:str,
        band:str,
        base_img:int|str=None,
        mode:workmode=workmode(),
) -> tuple[list, list, list, list]:
    """
    try to load last picking result, if not exists, do a new picking
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param obj: object
    :param band: band
    :param base_img: base image, if str (external file), treet as none
    :param mode: input files missing or output existence mode
    :returns: a tuple with 4 lists
    """
    logf = init_logger("pick", f"{red_dir}/log/pick.log", conf)
    pick_pkl = f"{red_dir}/pick_{obj}_{band}.pkl"
    if os.path.isfile(pick_pkl):
        logf.info(f"LOAD saved picking result: {obj} {band}")
        (
            magi, magc, ximg, yimg, 
            cali_cst, cali_std, x_mean, x_std, y_mean, y_std,
            xy, ind_ref, ind_var, ind_chk,
            magmed, magstd, magdif, magbad
        ) = pkl_load(pick_pkl)
    else:
        logf.info(f"NO saved picking result: {obj} {band}, picking")
        xy, ind_ref, ind_var, ind_chk = pick(
            conf, raw_dir, red_dir, obj, band, base_img, mode)
    return xy, ind_ref, ind_var, ind_chk