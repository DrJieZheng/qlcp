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
# import astropy.io.fits as fits
from .u_conf import config, workmode
from .u_log import init_logger
from .u_utils import loadlist, pkl_load, pkl_dump


def cali(
        conf:config,
        raw_dir:str,
        red_dir:str,
        obj:str,
        band:str,
        ind_tgt:int|list[int]=None,
        ind_ref:int|list[int]=None,
        ind_chk:int|list[int]=None,
        mode:workmode=workmode(),
):
    """
    chosen star info on given xy, draw finding chart
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param obj: object
    :param band: band
    :param ind_tgt: index of target star, if not, use 0
    :param ind_ref: index of reference star, if not, use 1:n-1
    :param ind_chk: index of check star, if not, use 1:n-1
    :param mode: input files missing or output existence mode
    :returns: Nothing
    """
    logf = init_logger("cali", f"{red_dir}/log/cali.log", conf)
    mode.reset_append(workmode.EXIST_OVER)

    # list file, and load list
    listfile = f"{red_dir}/lst/{obj}_{band}.lst"
    if mode.missing(listfile, f"{obj} {band} list", logf):
        return
    # raw_list = loadlist(listfile, base_path=raw_dir)
    # bf_fits_list = loadlist(listfile, base_path=red_dir,
    #                     suffix="bf.fits", separate_folder=True)
    # cat_fits_list = loadlist(listfile, base_path=red_dir,
    #                     suffix="cat.fits", separate_folder=True)
    # offset_pkl = f"{red_dir}/offset_{obj}_{band}.pkl"
    cata_pkl = f"{red_dir}/cata_{obj}_{band}.pkl"
    cali_pkl = f"{red_dir}/cali_{obj}_{band}.pkl"

    # check file exists
    if mode.missing(cata_pkl, f"general catalog {obj} {band}", logf):
        return
    if mode.exists(cali_pkl, f"calibrated catalog {obj} {band}", logf):
        return

    ###############################################################################

    # load final catalog
    cat_final, starxy, apstr = pkl_load(cata_pkl)
    nf = len(cat_final)
    ns = len(starxy)

    # process default value of indices
    if ind_tgt is None:
        ind_tgt = [0]
    elif isinstance(ind_tgt, int):
        ind_tgt = [ind_tgt]

    if ind_ref is None:
        # use all stars except the target star
        ind_ref = [i for i in range(ns) if i not in ind_tgt]
    elif isinstance(ind_ref, int):
        ind_ref = [ind_ref]

    if ind_chk is None:
        ind_chk = [i for i in range(ns) if i not in ind_tgt]
    elif isinstance(ind_chk, int):
        ind_chk = [ind_chk]

    n_tgt = len(ind_tgt)
    n_ref = len(ind_ref)
    n_chk = len(ind_chk)
    # ref string, used in filenames
    refs = "_".join([f"{i:02d}" for i in ind_ref])

    # stru of calibrated results, only results
    cat_cali_dt = [[
        (f"CaliTarget{a}", np.float32, (n_tgt,)),
        (f"ErrTarget{a}",  np.float32, (n_tgt,)),
        (f"CaliCheck{a}",  np.float32, (n_chk,)),
        (f"ErrCheck{a}",   np.float32, (n_chk,)),
        (f"CaliConst{a}",  np.float32,         ),
        (f"CaliStd{a}",    np.float32,         ),
    ] for a in apstr]
    cat_cali_dt = [b for a in cat_cali_dt for b in a]
    cat_cali = np.empty(nf, cat_cali_dt)

    # calibration each aperture
    for a in apstr:
        # create the dir for each aper
        aper_dir = f"{red_dir}/cali_{obj}_AP{a}"
        os.makedirs(aper_dir, exist_ok=True)

        # load mag and error of specified aperature
        mags = cat_final[f"Mag{a}"]
        errs = cat_final[f"Err{a}"]
        mags[mags == -1] = np.nan

        # compute the mean of ref stars as calibration const
        ref_mean = np.mean(mags[:, ind_ref], axis=1)
        ref_std = np.std(mags[:, ind_ref], axis=1)
        cat_cali[f"CaliConst{a}" ] = ref_mean
        cat_cali[f"CaliStd{a}"   ] = ref_std

        # calibration, each target and check star, substract ref_mean
        for i, k in enumerate(ind_tgt):
            cat_cali[f"CaliTarget{a}"][:, i] = mags[:, k] - ref_mean
            cat_cali[f"ErrTarget{a}" ][:, i] = errs[:, k]
        for i, k in enumerate(ind_chk):
            cat_cali[f"CaliCheck{a}"][:, i] = mags[:, k] - ref_mean
            cat_cali[f"ErrCheck{a}" ][:, i] = errs[:, k]

        # dump text files
        for i, k in enumerate(ind_tgt):
            with open(f"{aper_dir}/{obj}_{band}_vc{k:02d}_{refs}.txt", "w") as ff:
                for f in range(nf):
                    ff.write(f"{cat_final['BJD'][f]:15.7f} "
                             f"{cat_cali[f'CaliTarget{a}'][f, i]:6.3f} "
                             f"{cat_cali[f'ErrTarget{a}' ][f, i]:6.3f}\n")
        for i, k in enumerate(ind_chk):
            with open(f"{aper_dir}/{obj}_{band}_chk{k:02d}_{refs}.txt", "w") as ff:
                for f in range(nf):
                    ff.write(f"{cat_final['BJD'][f]:15.7f} "
                             f"{cat_cali[f'CaliCheck{a}'][f, i]:6.3f} "
                             f"{cat_cali[f'ErrCheck{a}' ][f, i]:6.3f}\n")

    # save to file
    pkl_dump(cali_pkl, cat_final, cat_cali, apstr, starxy, ind_tgt, ind_ref, ind_chk)

    logf.info(f"{n_tgt} target and {n_chk} check calibrated by {n_ref} ref.")
