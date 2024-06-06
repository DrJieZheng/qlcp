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
# import astropy.io.fits as fits
import matplotlib.pyplot as plt
import astropy.stats.sigma_clipping as sc
from .u_conf import config, workmode
from .u_log import init_logger
from .u_utils import loadlist, zenum, pkl_load#, pkl_dump


def graph(
        conf:config,
        raw_dir:str,
        red_dir:str,
        obj:str,
        band:str,
        mode:workmode=workmode(),
):
    """
    chosen star info on given xy, draw finding chart
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param obj: object
    :param band: band
    :param mode: input files missing or output existence mode
    :returns: Nothing
    """
    logf = init_logger("graph", f"{red_dir}/log/graph.log", conf)
    mode.reset_append(workmode.EXIST_OVER)

    # list file, and load list
    listfile = f"{red_dir}/lst/{obj}_{band}.lst"
    if mode.missing(listfile, f"{obj} {band} list", logf):
        return
    # raw_list = loadlist(listfile, base_path=raw_dir)
    cali_pkl = f"{red_dir}/cali_{obj}_{band}.pkl"

    # check file exists
    if mode.missing(cali_pkl, f"calibrated catalog {obj} {band}", logf):
        return

    ###############################################################################

    # load final catalog
    cat_final, cat_cali, apstr, starxy, ind_tgt, ind_ref, ind_chk = pkl_load(cali_pkl)
    nf = len(cat_final)
    n_tgt = len(ind_tgt)
    n_ref = len(ind_ref)
    n_chk = len(ind_chk)
    bjd = cat_final["BJD"]
    # space between each curve
    curve_space = 0.01

    # draw graph for each aperture
    for a in apstr:

        # load mag and error of specified aperature
        mtgt = cat_cali[f"CaliTarget{a}"]
        mchk = cat_cali[f"CaliCheck{a}"]

        # compute std of each check star
        std_chk = np.std(mchk, axis=1)

        # the height of mag, use th 2nd larggest and smallest
        range_tgt = np.empty(n_tgt)
        range_chk = np.empty(n_chk)
        if nf >= 4:
            for i in range(n_tgt):
                # ix = np.argsort(mtgt[:, i])
                # range_tgt[i] = mtgt[ix[-2], i] - mtgt[ix[1], i]
                _, _, s = sc.sigma_clipped_stats(mtgt[:, i], sigma=3)
                range_tgt[i] = s * 4
            for i in range(n_chk):
                # ix = np.argsort(mchk[:, i])
                # range_chk[i] = mchk[ix[-2], i] - mchk[ix[1], i]
                _, _, s = sc.sigma_clipped_stats(mchk[:, i], sigma=3)
                range_chk[i] = s * 4
        elif 2 <= nf <= 3:
            for i in range(n_tgt):
                ix = np.argsort(mtgt[:, i])
                range_tgt[i] = mtgt[ix[-1], i] - mtgt[ix[0], i]
            for i in range(n_chk):
                ix = np.argsort(mchk[:, i])
                range_chk[i] = mchk[ix[-1], i] - mchk[ix[0], i]
        else:
            range_tgt[:] = curve_space
            range_chk[:] = curve_space

        # the base position of each target and check star
        base_tgt = np.zeros(n_tgt)
        for i in range(1, n_tgt):
            base_tgt[i] = base_tgt[i-1] - (range_tgt[i-1] + range_tgt[i]) / 2 - curve_space
        base_chk = np.zeros(n_chk)
        base_chk[0] = (range_tgt[0] + range_chk[0]) / 2 + curve_space
        for i in range(1, n_chk):
            base_chk[i] = base_chk[i-1] + (range_chk[i-1] + range_chk[i]) / 2 + curve_space

        # the mean of each target and check star
        mean_tgt = np.mean(mtgt, axis=0)
        mean_chk = np.mean(mchk, axis=0)

        # the y size of graph
        ysize = (n_tgt + n_chk) * curve_space + sum(range_tgt) + sum(range_chk)

        # draw graph
        fig, ax = plt.subplots(figsize=(10, ysize*20))
        for i, k in enumerate(ind_tgt):
            ax.plot(bjd, mtgt[:, i] - mean_tgt[i] + base_tgt[i],
                    "*", label=f"Target{k:2d}")
            # ax.axhline(y=base_tgt[i], color="k", linestyle="--")
            # ax.axhline(y=base_tgt[i] + range_tgt[i]/2, color="k", linestyle=":")
            # ax.axhline(y=base_tgt[i] - range_tgt[i]/2, color="k", linestyle=":")
        for i, k in enumerate(ind_chk):
            ax.plot(bjd, mchk[:, i] - mean_chk[i] + base_chk[i],
                    "+", label=f"Check{k:2d} $\\sigma$={std_chk[i]:6.4f}")
            # ax.axhline(y=base_chk[i], color="k", linestyle="--")
            # ax.axhline(y=base_chk[i] + range_chk[i]/2, color="k", linestyle=":")
            # ax.axhline(y=base_chk[i] - range_chk[i]/2, color="k", linestyle=":")
        ax.legend()
        # ax.invert_yaxis()
        ax.set_ylim(base_chk[-1] + range_chk[-1], base_tgt[-1] - range_tgt[-1] )
        ax.set_xlabel("BJD")
        ax.set_ylabel("Relative Magnitude")
        ax.set_title(f"{obj}-{band} (Aperture={a})")
        fig.savefig(f"{red_dir}/lc_{obj}_{band}_AP{a}.png", bbox_inches="tight")
        plt.close()
        logf.info(f"Light-curve saved as {red_dir}/lc_{obj}_{band}_{a}.png")
