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
from .u_conf import config
from .u_workmode import workmode
from .u_log import init_logger
from .u_utils import loadlist, rm_ix, zenum, meanclip, tqdm_bar


def phot(
        conf:config,
        red_dir:str,
        obj:str,
        band:str,
        # se_cmd:str="source-extractor",
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
    # if os.system(f"which {se_cmd} > /dev/null"):
    #     raise OSError(f"Wrong Source-Extractor Command: {se_cmd}")
    
    # check installation of se
    if os.system(f"which sex > /dev/null") == 0:
        se_cmd = "sex"
    elif os.system(f"which source-extractor > /dev/null") == 0:
        se_cmd = "source-extractor"
    elif os.system(f"which sextractor > /dev/null") == 0:
        se_cmd = "sextractor"
    else:
        raise OSError("No Source-Extractor Installed")

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
    mode.start_lazy()
    ix = []
    for i, (scif, catf) in zenum(bf_fits_list, cat_fits_list):
        if mode.missing(scif, "corrected image", None) or \
           mode.exists(catf, "catalog", None):
            ix.append(i)
    # remove missing file
    mode.end_lazy(logf)
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
    # remark any aperture providec by caller
    no_aper = not aper
    # apertures
    if not aper:
        aper = [5.0]
    aper = aper if isinstance(aper, (list, tuple)) else [aper]
    apstr  = [] if no_aper else [f"{a:04.1f}" for a in aper]
    # se command
    se_call = f"{se_cmd} -c {se_sex} {{}} " \
              f"-parameters_name {se_par} " \
              f"-CATALOG_NAME {{}} " \
              f"-PHOT_APERTURES {{}} 2> /dev/null"
    # se command for fwhm testing
    se_fwhm = f"{se_cmd} -c {conf.here}bright.sex {{}} " \
              f"-parameters_name {conf.here}bright.param " \
              f"-CATALOG_NAME {{}} "\
              f"2>/dev/null"

    # catalog datatypes
    mycatdt = [
        ("Num",      np.uint16 ),
        ("X",        np.float64),
        ("Y",        np.float64),
        ("Elong",    np.float32),
        ("FWHM",     np.float32),
        ("MagAUTO",  np.float32),
        ("ErrAUTO",  np.float32),
        ("FluxAUTO", np.float32),
        ("FErrAUTO", np.float32),
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
    bc = max((max(aper) * 2, 20))

    # load images and process
    pbar = tqdm_bar(nf, f"PHOT {obj} {band}")
    for i, (scif, sef, catf, txtf, pngf) in zenum(
            bf_fits_list, se_fits_list, cat_fits_list, cat_txt_list, cat_png_list):

        # call se for fwhm
        logf.debug(f"SE {i+1:03d}/{nf:03d}: "+se_fwhm.format(scif, sef))
        os.system(se_fwhm.format(scif, sef))
        # load and processor catalog
        secat = fits.getdata(sef, 2)
        ix = np.where(
            (bc < secat["X_IMAGE"]) & (secat["X_IMAGE"] < nx - bc) &
            (bc < secat["Y_IMAGE"]) & (secat["Y_IMAGE"] < ny - bc) &
            (secat["MAGERR_AUTO"] < 0.05)) [0]
        fwhm, _ = meanclip(secat["FWHM_IMAGE"][ix])

        # call se
        apcomma = ",".join([f"{a:.1f}" if a > 0 else f"{-a*fwhm:.1f}" for a in aper])
        logf.debug(f"SE {i+1:03d}/{nf:03d}: "+se_call.format(scif, sef, apcomma))
        os.system(se_call.format(scif, sef, apcomma))

        # load and processor catalog
        secat = fits.getdata(sef, 2)
        # remove stars at border
        ix = np.where(
            (bc < secat["X_IMAGE_DBL"]) & (secat["X_IMAGE_DBL"] < nx - bc) &
            (bc < secat["Y_IMAGE_DBL"]) & (secat["Y_IMAGE_DBL"] < ny - bc) )[0]
        secat = secat[ix]

        ns = len(ix)
        mycat = np.empty(ns, mycatdt)
        mycat["Num"     ] = secat["NUMBER"]
        mycat["X"       ] = secat["X_IMAGE_DBL"] -1
        mycat["Y"       ] = secat["Y_IMAGE_DBL"] -1
        mycat["Elong"   ] = secat["ELONGATION"]
        mycat["FWHM"    ] = secat["FWHM_IMAGE"]
        mycat["MagAUTO" ] = secat["MAG_AUTO"]
        mycat["ErrAUTO" ] = secat["MAGERR_AUTO"]
        mycat["FluxAUTO"] = secat["FLUX_AUTO"]
        mycat["FErrAUTO"] = secat["FLUXERR_AUTO"]
        mycat["Flags"   ] = secat["FLAGS"]
        mycat["Alpha"   ] = secat["ALPHA_J2000"]
        mycat["Delta"   ] = secat["DELTA_J2000"]
        for k, a in enumerate(apstr):
            mycat[f"Mag{a}" ] = secat["MAG_APER"]    [:, k]
            mycat[f"Err{a}" ] = secat["MAGERR_APER"] [:, k]
            mycat[f"Flux{a}"] = secat["FLUX_APER"]   [:, k]
            mycat[f"FErr{a}"] = secat["FLUXERR_APER"][:, k]

        hdr = fits.getheader(scif)
        hdr["IMNAXIS1"] = hdr["NAXIS1"]
        hdr["IMNAXIS2"] = hdr["NAXIS2"]
        hdr["APERS"] = ','.join(apstr)
        hdr["NAPER"] = len(apstr)
        hdr["FWHM"] = fwhm
        for k, a in enumerate(apstr):
            hdr[f"APER{k+1:1d}"] = float(a)

        pri_hdu = fits.PrimaryHDU(header=hdr, data=None)
        cat_hdu = fits.BinTableHDU(data=mycat)
        new_fits = fits.HDUList([pri_hdu, cat_hdu])
        new_fits.writeto(catf, overwrite=True)

        with open(txtf, "w") as ff:
            ff.write(
                f"{'#Num':>4s}  {'X':>8s} {'Y':>8s}  "
                f"{'Elong':>5s} {'FWHM':>5s}  "
                f"{'MagAUTO':>7s} {'ErrAUTO':>7s}  " +
                "  ".join([f"Mag{a} Err{a}" for a in apstr]) + 
                f"  {'Flags':>16s}  {'Alpha':>10s} {'Delta':>10s}\n")
            for s in mycat:
                ff.write((
                    "{s[Num]:4d}  {s[X]:8.3f} {s[Y]:8.3f}  "
                    "{s[Elong]:5.2f} {s[FWHM]:5.2f}  "
                    "{s[MagAUTO]:7.3f} {s[ErrAUTO]:7.4f}  " + 
                    "  ".join([f"{{s[Mag{a}]:7.3f}} {{s[Err{a}]:7.4f}}" for a in apstr]) + 
                    "  {s[Flags]:16b}  {s[Alpha]:10.6f} {s[Delta]:+10.6f}\n")
                    .format(s=s))
        logf.debug(f"{ns} objects (FWHM={fwhm:5.2f}) dump to {catf}")

        if conf.draw_phot:
            img = fits.getdata(scif)
            _, imgm, imgs = sigma_clipped_stats(img, sigma=3)
            ix = np.where(mycat["ErrAUTO"] < conf.draw_phot_err)

            fig, ax = plt.subplots(figsize=(8, 8))
            ax.imshow(img-imgm, vmax=3*imgs, vmin=-3*imgs,
                origin='lower', cmap='gray')
            ax.scatter(mycat["X"][ix], mycat["Y"][ix],
                s=10, c="none", marker="o", edgecolors="red")
            ax.set_title(f"{scif} {ns:d}")
            fig.savefig(pngf)
            plt.close()
        pbar.update(1)
    pbar.close()

    logf.info(f"{nf:3d} files photometry done")
