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
from .u_conf import config
from .u_workmode import workmode
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

    logf.debug(f"Target: {ind_tgt}")
    logf.debug(f"Reference: {ind_ref}")
    logf.debug(f"Check: {ind_chk}")

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
    refs = "_".join([f"c{i:02d}" for i in ind_ref])

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
            with open(f"{aper_dir}/{obj}_{band}_v{k:02d}_{refs}.txt", "w") as ff:
                for f in range(nf):
                    ff.write(f"{cat_final['BJD'][f]:15.7f} "
                             f"{cat_cali[f'CaliTarget{a}'][f, i]:6.3f} "
                             f"{cat_cali[f'ErrTarget{a}' ][f, i]:6.3f}\n")
        for i, k in enumerate(ind_chk):
            with open(f"{aper_dir}/{obj}_{band}_ch{k:02d}_{refs}.txt", "w") as ff:
                for f in range(nf):
                    ff.write(f"{cat_final['BJD'][f]:15.7f} "
                             f"{cat_cali[f'CaliCheck{a}'][f, i]:6.3f} "
                             f"{cat_cali[f'ErrCheck{a}' ][f, i]:6.3f}\n")

    # save to file
    pkl_dump(cali_pkl, cat_final, cat_cali, apstr, starxy, ind_tgt, ind_ref, ind_chk)
    logf.info(f"{n_tgt} target and {n_chk} check calibrated by {n_ref} ref.")

    ###############################################################################

    # draw graph
    bjd = cat_final["BJD"]
    # space between each curve
    curve_space = 0.01
    # tgt,chk,ref in filename
    fn_part = "_".join(
        [f"v{i:02d}" for i in ind_tgt] +
        [f"ch{i:02d}" for i in ind_chk] +
        [f"c{i:02d}" for i in ind_ref]
    )
    # color and marker
    cm_tgt = lambda i: ("rs",)[i % 1]
    cm_chk = lambda i: ("b*", "g*", "m*", "c*")[i % 4]

    # draw graph for each aperture
    for a in apstr:
        # load mag and error of specified aperature
        mtgt = cat_cali[f"CaliTarget{a}"]
        mchk = cat_cali[f"CaliCheck{a}"]

        # compute std of each check star
        std_chk = np.nanstd(mchk, axis=0)
        # the mean of each target and check star
        mean_tgt = np.empty(n_tgt)
        mean_chk = np.empty(n_chk)

        # the half height of stars
        half_range_tgt = np.empty(n_tgt)
        half_range_chk = np.empty(n_chk)
        if nf >= 4:
            for i in range(n_tgt):
                mean_tgt[i], _, s = sc.sigma_clipped_stats(mtgt[:, i], sigma=3)
                half_range_tgt[i] = s * 4
            for i in range(n_chk):
                mean_chk[i], _, s = sc.sigma_clipped_stats(mchk[:, i], sigma=3)
                half_range_chk[i] = s * 4
        elif 2 <= nf <= 3:
            for i in range(n_tgt):
                mean_tgt[i] = np.nanmedian(mtgt[:, i])
                half_range_tgt[i] = (max(mtgt[:, i]) - min(mtgt[:, i])) / 2
            for i in range(n_chk):
                mean_chk[i] = np.nanmedian(mchk[:, i])
                half_range_chk[i] = (max(mchk[:, i]) - min(mchk[:, i])) / 2
        else:
            mean_tgt[:] = mtgt[0]
            mean_chk[:] = mchk[0]
            half_range_tgt[:] = curve_space
            half_range_chk[:] = curve_space

        # the base position of each target and check star
        base_tgt = np.zeros(n_tgt)
        for i in range(1, n_tgt):
            base_tgt[i] = base_tgt[i-1] - half_range_tgt[i-1] + half_range_tgt[i] - curve_space
        base_chk = np.zeros(n_chk)
        base_chk[0] = half_range_tgt[0] + half_range_chk[0] + curve_space
        for i in range(1, n_chk):
            base_chk[i] = base_chk[i-1] + half_range_chk[i-1] + half_range_chk[i] + curve_space

        # the y size of graph
        ysize = (n_tgt + n_chk) * curve_space + sum(half_range_tgt) + sum(half_range_chk)

        # draw graph
        fig, ax = plt.subplots(figsize=(10, ysize*20))
        for i, k in enumerate(ind_tgt):
            ax.plot(bjd, mtgt[:, i] - mean_tgt[i] + base_tgt[i],
                    cm_tgt(i), label=f"Target{k:2d}")
            # ax.axhline(y=base_tgt[i], color="k", linestyle="--")
            # ax.axhline(y=base_tgt[i] + range_tgt[i]/2, color="k", linestyle=":")
            # ax.axhline(y=base_tgt[i] - range_tgt[i]/2, color="k", linestyle=":")
        for i, k in enumerate(ind_chk):
            ax.plot(bjd, mchk[:, i] - mean_chk[i] + base_chk[i],
                    cm_chk(i), label=f"Check{k:2d} $\\sigma$={std_chk[i]:6.4f}")
            # ax.axhline(y=base_chk[i], color="k", linestyle="--")
            # ax.axhline(y=base_chk[i] + range_chk[i]/2, color="k", linestyle=":")
            # ax.axhline(y=base_chk[i] - range_chk[i]/2, color="k", linestyle=":")
        ax.legend()
        # ax.invert_yaxis()
        ax.set_ylim(base_chk[-1] + half_range_chk[-1] + curve_space,
                    base_tgt[-1] - half_range_tgt[-1] - curve_space)
        ax.set_xlabel("BJD")
        ax.set_ylabel("Relative Magnitude")
        ax.set_title(f"{obj}-{band} (Aperture={a})")
        fig.savefig(f"{red_dir}/lc_{obj}_{band}_AP{a}_{fn_part}.png", bbox_inches="tight")
        # plt.close()
        logf.info(f"Light-curve saved for {obj} {band}")
