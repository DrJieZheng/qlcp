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
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from matplotlib import pyplot as plt
from .u_conf import config, workmode
from .u_log import init_logger
from .u_utils import loadlist, rm_ix, zenum


def phot(
        conf:config,
        red_dir:str,
        obj:str,
        band:str,
        se_cmd:str="source-extractor",
        aper:float|list[float]=None,
        mode:workmode=workmode(),
):
    """
    image bias & flat correction shell, check files and finally call the processing method
    :param conf: config object
    :param red_dir: red files dir
    :param obj: object
    :param band: band
    :param se_cmd: command of source-extractor, if empty or None, use sep package
    :param aper: one or multi aperture(s) for photometry
    :param mode: input files missing or output existence mode
    :returns: Nothing
    """
    logf = init_logger("phot", f"{red_dir}/log/phot.log", conf)
    mode.reset_append(workmode.EXIST_SKIP)

    # check se cmd
    if os.system(f"which {se_cmd} > /dev/null"):
        raise OSError(f"Wrong Source-Extractor Command: {se_cmd}")

    # list file, and load list
    listfile = f"{red_dir}/lst/{obj}_{band}.lst"
    if mode.missing(listfile, f"{obj} {band} list", logf):
        return
    bf_fits_list = loadlist(listfile, base_path=red_dir,
                        suffix="bf.fits", separate_folder=True)
    se_fits_list = loadlist(listfile, base_path=red_dir,
                        suffix="se.fits", separate_folder=True)
    cat_fits_list = loadlist(listfile, base_path=red_dir,
                        suffix="cat.fits", separate_folder=True)
    cat_txt_list = loadlist(listfile, base_path=red_dir,
                        suffix="cat.txt", separate_folder=True)
    cat_png_list = loadlist(listfile, base_path=red_dir,
                        suffix="cat.png", separate_folder=True)

    # check file missing
    ix = []
    for i, (scif, catf) in zenum(bf_fits_list, cat_fits_list):
        if mode.missing(scif, "corrected image", logf) or \
           mode.exists(catf, "catalog", logf):
            ix.append(i)
    # remove missing file
    rm_ix(ix, bf_fits_list, se_fits_list, cat_fits_list, cat_txt_list)
    nf = len(bf_fits_list)

    if nf == 0:
        logf.info(f"SKIP {obj} {band} Nothing")
        return

    ###############################################################################

    logf.debug(f"{nf} files for {obj} {band}")

    # se_command and parameters, use local if exists
    if os.path.isfile("default.sex"):
        se_sex = "default.sex"
    else:
        se_sex = conf.here + "default.sex"
    se_par = conf.here + "default.param"
    # apertures
    if not aper:
        aper = [5.0]
    aper = aper if isinstance(aper, (list, tuple)) else [aper]
    apstr  = [f"{a:04.1f}" for a in aper]
    apcomma = ",".join([f"{a:.1f}" for a in aper])
    # se command
    secall = f"{se_cmd} -c {se_sex} -parameters_name " \
             f"{se_par} {{}} -CATALOG_NAME {{}} "\
             f"-PHOT_APERTURES {apcomma} 2> /dev/null"

    # catalog datatypes
    mycatdt = [
        ("Num",      np.uint16 ),
        ("X",        np.float64),
        ("Y",        np.float64),
        ("Elong",    np.float32),
        ("FWHM",     np.float32),
        ("Mag00.0",  np.float32),
        ("Err00.0",  np.float32),
        ("Flux00.0", np.float32),
        ("FErr00.0", np.float32),
    ] + [
        (f"Mag{a}",  np.float32) for a in apstr] + [
        (f"Err{a}",  np.float32) for a in apstr] + [
        (f"Flux{a}", np.float32) for a in apstr] + [
        (f"FErr{a}", np.float32) for a in apstr] + [
        ("Flags",    np.uint16 ),
        ("Alpha",    np.float64),
        ("Delta",    np.float64),
    ]

    # image size
    hdr = fits.getheader(bf_fits_list[0])
    nx = hdr["NAXIS1"]
    ny = hdr["NAXIS2"]
    # border cut as 2 * max aperture size
    bc = max(aper) * 2

    # load images and process
    for i, (scif, sef, catf, txtf, pngf) in zenum(
            bf_fits_list, se_fits_list, cat_fits_list, cat_txt_list, cat_png_list):

        # call se
        logf.debug(f"SE {i+1:03d}/{nf:03d}: "+secall.format(scif, sef))
        os.system(secall.format(scif, sef))

        # load and processor catalog
        secat = fits.getdata(sef, 2)
        # remove stars at border
        ix = np.where(
            (bc < secat["X_IMAGE_DBL"] < bc) & (secat["X_IMAGE_DBL"] < nx - bc) &
            (bc < secat["Y_IMAGE_DBL"] < bc) & (secat["Y_IMAGE_DBL"] < ny - bc) )[0]

        ns = len(ix)
        mycat = np.empty(ns, mycatdt)
        mycat["Num"     ] = secat["NUMBER"]      [ix]
        mycat["X"       ] = secat["X_IMAGE_DBL"] [ix]
        mycat["Y"       ] = secat["Y_IMAGE_DBL"] [ix]
        mycat["Elong"   ] = secat["ELONGATION"]  [ix]
        mycat["FWHM"    ] = secat["FWHM_IMAGE"]  [ix]
        mycat["Mag00.0" ] = secat["MAG_AUTO"]    [ix]
        mycat["Err00.0" ] = secat["MAGERR_AUTO"] [ix]
        mycat["Flux00.0"] = secat["FLUX_AUTO"]   [ix]
        mycat["FErr00.0"] = secat["FLUXERR_AUTO"][ix]
        mycat["Flags"   ] = secat["FLAGS"]       [ix]
        mycat["Alpha"   ] = secat["ALPHA_J2000"] [ix]
        mycat["Delta"   ] = secat["DELTA_J2000"] [ix]
        for k, a in enumerate(apstr):
            mycat[f"Mag{a}" ] = secat["MAG_APER"]    [ix, k]
            mycat[f"Err{a}" ] = secat["MAGERR_APER"] [ix, k]
            mycat[f"Flux{a}"] = secat["FLUX_APER"]   [ix, k]
            mycat[f"FErr{a}"] = secat["FLUXERR_APER"][ix, k]

        hdr = fits.getheader(scif)
        hdr["IMNAXIS1"] = hdr["NAXIS1"]
        hdr["IMNAXIS2"] = hdr["NAXIS2"]
        hdr["APERS"] = ','.join(apstr)
        hdr["NAPER"] = len(apstr)
        for k, a in enumerate(aper):
            hdr[f"APER{k+1:1d}"] = a

        pri_hdu = fits.PrimaryHDU(header=hdr, data=None)
        cat_hdu = fits.BinTableHDU(data=mycat)
        new_fits = fits.HDUList([pri_hdu, cat_hdu])
        new_fits.writeto(catf, overwrite=True)

        with open(txtf, "w") as ff:
            ff.write(
                f"{'#Num':>4s}  {'X':>8s} {'Y':>8s}  "
                f"{'Elong':>5s} {'FWHM':>5s}  "
                f"{'Mag00.0':>7s} {'Err00.0':>7s}  " +
                "  ".join([f"Mag{a} Err{a}" for a in apstr]) + 
                f"  {'Flags':>16s}  {'Alpha':>10s} {'Delta':>10s}\n")
            for s in mycat:
                ff.write((
                    "{s[Num]:4d}  {s[X]:8.3f} {s[Y]:8.3f}  "
                    "{s[Elong]:5.2f} {s[FWHM]:5.2f}  "
                    "{s[Mag00.0]:7.3f} {s[Err00.0]:7.4f}  " + 
                    "  ".join([f"{{s[Mag{a}]:7.3f}} {{s[Err{a}]:7.4f}}" for a in apstr]) + 
                    "  {s[Flags]:16b}  {s[Alpha]:10.6f} {s[Delta]:+10.6f}\n")
                    .format(s=s))
        logf.debug(f"{ns} objects dump to {catf}")

        if conf.draw_phot:
            img = fits.getdata(scif)
            _, imgm, imgs = sigma_clipped_stats(img, sigma=3)
            ix = np.where(mycat["Err00.0"] < conf.draw_phot_err)

            fig, ax = plt.subplots(figsize=(8, 8))
            ax.imshow(img-imgm, vmax=3*imgs, vmin=-3*imgs,
                origin='lower', cmap='gray')
            ax.scatter(mycat["X"][ix], mycat["Y"][ix],
                s=10, c="none", marker="o", edgecolors="red")
            ax.set_title(f"{scif} {ns:d}")
            fig.savefig(pngf)

    logf.info(f"{nf:3d} files photometry done")
