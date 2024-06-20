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
import matplotlib.pyplot as plt
from tqdm import tqdm
from qmatch import match2d
from .u_conf import config
from .u_workmode import workmode
from .u_log import init_logger
from .u_utils import loadlist, rm_ix, zenum, pkl_load, pkl_dump, \
    cat2txt, fnbase, meanclip, tqdm_bar


def cata(
        conf:config,
        raw_dir:str,
        red_dir:str,
        obj:str,
        band:str,
        starxy:list[list[float]]|np.ndarray,
        base_img:int|str=0,
        mode:workmode=workmode(),
):
    """
    chosen star info on given xy, draw finding chart
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param obj: object
    :param band: band
    :param starxy: list of star xy, each star is a tuple/list of x, y
        usually 0th is the target and the others as ref/chk
    :param base_img: the offset base index or filename
    :param mode: input files missing or output existence mode
    :returns: Nothing
    """
    logf = init_logger("cata", f"{red_dir}/log/cata.log", conf)
    mode.reset_append(workmode.EXIST_OVER)

    # list file, and load list
    listfile = f"{red_dir}/lst/{obj}_{band}.lst"
    if mode.missing(listfile, f"{obj} {band} list", logf):
        return
    # raw_list = loadlist(listfile, base_path=raw_dir)
    bf_fits_list = loadlist(listfile, base_path=red_dir,
                        suffix="bf.fits", separate_folder=True)
    cat_fits_list = loadlist(listfile, base_path=red_dir,
                        suffix="cat.fits", separate_folder=True)
    offset_pkl = f"{red_dir}/offset_{obj}_{band}.pkl"
    cata_pkl = f"{red_dir}/cata_{obj}_{band}.pkl"
    cata_fits = f"{red_dir}/cata_{obj}_{band}.fits"
    cata_txt = f"{red_dir}/cata_{obj}_{band}.txt"
    cata_png = f"{red_dir}/cata_{obj}_{band}.png"

    # check file exists
    if mode.missing(offset_pkl, f"offset {obj} {band}", logf):
        return
    if mode.exists(cata_pkl, f"general catalog {obj} {band}", logf):
        return

    # check file missing
    mode.start_lazy()
    ix = []
    for i, (f,) in zenum(cat_fits_list):
        if mode.missing(f, "image catalog", None):
            ix.append(i)
    # remove missing file
    rm_ix(ix, bf_fits_list, cat_fits_list)
    mode.end_lazy(logf)
    nf = len(cat_fits_list)

    # base image, type check, range check, existance check
    if isinstance(base_img, int):
        # if int, check range, then use nth image
        if 0 > base_img or base_img >= nf:
            base_img = 0
        base_img = bf_fits_list[base_img]
    elif not isinstance(base_img, str):
        # if not int nor str, use 0th
        base_img = bf_fits_list[0]
    # special, fixed mode
    elif workmode(workmode.MISS_SKIP).missing(base_img, "offset base image", logf):
        # if str, is external file, if not found, use 0th
        base_img = bf_fits_list[0]

    if nf == 0:
        logf.info(f"SKIP {obj} {band} No Star")

    ###############################################################################

    # confirm starxy is 2d np array, and split x, y, and star count
    starxy = np.array(starxy)
    starx = starxy[:, 0]
    stary = starxy[:, 1]
    ns = len(starxy)

    # exit if no star
    if ns == 0:
        logf.error(f"SKIP {obj} {band} No Star")

    # load offset result, and transfer to dict
    _, offset_x, offset_y, offset_bf_fits_list = pkl_load(offset_pkl)
    offset_bf_fits_list = [fnbase(f) for f in offset_bf_fits_list]
    offset_x = dict(zip(offset_bf_fits_list, offset_x))
    offset_y = dict(zip(offset_bf_fits_list, offset_y))
    fn_len_max = max(len(f) for f in offset_bf_fits_list)

    # aper info, including 0.0
    apstr = ('AUTO,' + fits.getval(cat_fits_list[0], "APERS")).split(",")

    # init the catalog structure
    cat_inst_magflux_dt = [[
        (f"Mag{a}",  np.float32, (ns,)),
        (f"Err{a}",  np.float32, (ns,)),
        (f"Flux{a}", np.float32, (ns,)),
        (f"FErr{a}", np.float32, (ns,)),
    ] for a in apstr]
    cat_inst_magflux_dt = [b for a in cat_inst_magflux_dt for b in a]
    cat_inst_dt =  [
        ("File",     (str, fn_len_max)),
        ("Band",     (str, 10),),
        ("Expt",     np.float32),
        ("DT",       (str, 22),),
        ("JD",       np.float64),
        ("BJD",      np.float64),
        ("HJD",      np.float64),
        ("ID" ,      np.uint16 , (ns,)),
        ("X",        np.float64, (ns,)),
        ("Y",        np.float64, (ns,)),
        ("FWHM",     np.float32, (ns,)),
        ("Elong",    np.float32, (ns,)),
    ] + cat_inst_magflux_dt + [
        ("Flags",    np.uint16 , (ns,)),
    ]
    cat_inst = np.empty(nf, cat_inst_dt)

    pbar = tqdm_bar(nf, f"CATA {obj} {band}")
    # load stars from images into the array, by matching x,y
    for i, (catf, bff) in zenum(cat_fits_list, bf_fits_list):
        bbff = fnbase(bff)
        # load image info
        hdr = fits.getheader(catf)
        # carry image global info to catalog
        cat_inst[i]["File"] = bbff
        cat_inst[i]["DT"  ] = hdr["DATE-OBS"]
        cat_inst[i]["Band"] = hdr["FILTER"]
        cat_inst[i]["Expt"] = hdr["EXPTIME"]
        cat_inst[i]["JD"  ] = hdr["JD"]
        cat_inst[i]["BJD" ] = hdr["BJD"]
        cat_inst[i]["HJD" ] = hdr["HJD"]

        # load stars from image
        cat_i = fits.getdata(catf, 1)
        # move to the same coordinate system
        x_i = cat_i["X"] + offset_x[bbff]
        y_i = cat_i["Y"] + offset_y[bbff]
        # match to the catalog
        ix_s, ix_k = match2d(starx, stary, x_i, y_i, conf.match_max_dis)
        # dump the matched stars
        cat_inst[i]["ID"][ix_s] = ix_k
        cat_inst[i]["X" ][ix_s] = cat_i[ix_k]["X"]
        cat_inst[i]["Y" ][ix_s] = cat_i[ix_k]["Y"]
        for a in apstr:
            cat_inst[i][f"Mag{a}" ][ix_s] = cat_i[ix_k][f"Mag{a}" ]
            cat_inst[i][f"Err{a}" ][ix_s] = cat_i[ix_k][f"Err{a}" ]
            cat_inst[i][f"Flux{a}"][ix_s] = cat_i[ix_k][f"Flux{a}"]
            cat_inst[i][f"FErr{a}"][ix_s] = cat_i[ix_k][f"FErr{a}"]
        cat_inst[i]["FWHM"][ix_s] = cat_i[ix_k]["FWHM"]
        cat_inst[i]["Elong"][ix_s] = cat_i[ix_k]["Elong"]
        logf.debug(f"Add {i+1:3d}/{nf:3d}: {len(cat_i):4d}->{len(ix_k):4d} {bbff}")
        pbar.update(1)
    pbar.close()

    # save catalog to bintable fits, and pickle
    pri_hdu = fits.PrimaryHDU()
    tb_hdu = fits.BinTableHDU(data=cat_inst)
    new_fits = fits.HDUList([pri_hdu, tb_hdu])
    new_fits.writeto(cata_fits, overwrite=True)
    pkl_dump(cata_pkl, cat_inst, starxy, apstr)
    logf.info(f"Result save to {cata_pkl}")

    # dump catalog to text file
    # cat2txt(cata_txt, cat_inst)

    # dump catalog to text file, each star in a line
    with open(cata_txt, "w") as ff:
        ff.write(f"# {'DateTime':^19s} {'Mag':^6s}{'ExpTime':7s}"
                 f"{'Band':4s}{'File':^{fn_len_max}s}{'BJD':^15s} "
                 f"{'Err':^6s} {'FWHM':5s} "
                 f"{'X':^6s} {'Y':^6s}\n")
        for c in cat_inst:
            for s in range(ns):
                ff.write(f"{c['DT']:21s} {c['MagAUTO'][s]:6.3f} {c['Expt']:5.1f} "
                         f"{c['Band']:2s} {c['File']:{fn_len_max}s} {c['BJD']:15.7f}  "
                         f"{c['ErrAUTO'][s]:6.4f} {s:2d} {c['FWHM'][s]:5.2f} "
                         f"{c['X'][s]:6.1f} {c['Y'][s]:6.1f}\n")

    # plot finding chart
    if os.path.isfile(base_img):
        img = fits.getdata(base_img)
        imtitle = os.path.basename(base_img)
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        imm, ims = meanclip(img)
        ax.imshow(img-imm, vmin=-1*ims, vmax=5*ims, origin="lower", cmap="gray")
        ax.scatter(starx-1, stary-1, marker="o", s=50, c="none", edgecolors="r")
        for i, (x, y) in zenum(starx, stary):
            ax.text(x+10, y+10, f"{i:d}", color="r")
        ax.set_title(f"{imtitle} ({len(starx)} stars)")
        ax.set_xlabel("X (pixel)")
        ax.set_ylabel("Y (pixel)")
        fig.savefig(cata_png, bbox_inches="tight")
        plt.close(fig)
