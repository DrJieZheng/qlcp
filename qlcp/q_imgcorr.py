# -*- coding: utf-8 -*-
"""
    v1 201901, Dr. Jie Zheng, Beijing & Xinglong, NAOC
    v2 202101, Dr. Jie Zheng & Dr./Prof. Linqiao Jiang
    v3 202201, Zheng & Jiang
    v4 202304, Upgrade, restructure, Zheng & Jiang
    Quick_Light_Curve_Pipeline
"""


import os
import re
# import numpy as np
import astropy.io.fits as fits
from astropy import time, coordinates as coord, units as u
from PyAstronomy import pyasl
from qastutil import lst, azalt, airmass
from .u_conf import config, workmode
from .u_log import init_logger
from .u_utils import loadlist, rm_ix, zenum, uttimestr, hdr_dt, ra2hms, dec2dms


def imgcorr(
        conf:config,
        raw_dir:str,
        red_dir:str,
        obj:str,
        band:str,
        use_bias:str=None,
        use_flat:str=None,
        alt_bias:str=None,
        alt_flat:str=None,
        alt_coord:tuple[str]=None,
        mode:workmode=workmode(),
):
    """
    image bias & flat correction
    Scientific data correciton by bias & flat, and header update
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param obj: object
    :param band: band
    :param use_bias: bias filename
    :param use_flat: flat filename
    :param alt_bias: alternative bias filename
    :param alt_flat: alternative flat filename
    :param alt_coord: coordination of object, if not provided in fits header
                      tuple of ra and dec, hh:mm:ss.s +dd:mm:ss
    :param mode: input files missing or output existence mode
    :returns: Nothing
    """
    logf = init_logger("imgcorr", f"{red_dir}/log/imgcorr.log", conf)
    os.makedirs(f"{red_dir}/{obj}_{band}", exist_ok=True)
    mode.reset_append(workmode.EXIST_SKIP)

    # list file, and load list
    listfile = f"{red_dir}/lst/{obj}_{band}.lst"
    if mode.missing(listfile, f"{obj} {band} list", logf):
        return
    raw_list = loadlist(listfile, base_path=raw_dir)
    bf_fits_list = loadlist(listfile, base_path=red_dir,
                        suffix="bf.fits", separate_folder=True)

    # bias file, prefer use_bias, then bias today, then alt_bias
    bias_fits = use_bias if use_bias else f"{red_dir}/bias.fits"
    if workmode(workmode.MISS_SKIP).missing(bias_fits, "master bias", logf):
        bias_fits = alt_bias
    # flat file, prefer use_flat, then flat today, then alt_flat
    flat_fits = use_flat if use_flat else f"{red_dir}/flat_{band}.fits"
    if workmode(workmode.MISS_SKIP).missing(flat_fits, f"master flat {band}", logf):
        flat_fits = alt_flat

    # check file missing
    if mode.missing(bias_fits, "master bias", logf):
        return
    if mode.missing(flat_fits, f"master flat {band}", logf):
        return
    ix = []
    for i, (rawf, bff) in zenum(raw_list, bf_fits_list):
        if mode.missing(rawf, "raw image", logf) or \
           mode.exists(bff, "corrected image", logf):
            ix.append(i)
    # remove missing file
    rm_ix(ix, raw_list, bf_fits_list)
    nf = len(raw_list)

    # check coord
    if alt_coord:
        ra_exam = re.match("^[0-9]{2}:[0-9]{2}:[0-9]{2}(\\.[0-9]{0,4})?$", alt_coord[0])
        dec_exam = re.match("^[+-]?[0-9]{2}:[0-9]{2}:[0-9]{2}(\\.[0-9]{0,4})?$", alt_coord[1])
        if ra_exam is None:
            raise ValueError(f"RA format wrong (hh:mm:ss.ssss) <{alt_coord[0]}>")
        if dec_exam is None:
            raise ValueError(f"Dec format wrong (+/-dd:mm:ss.ssss) <{alt_coord[1]}>")

    if nf == 0:
        logf.info(f"SKIP {obj} {band} Nothing")
        return

    ###############################################################################

    logf.debug(f"{nf} files for {obj} {band}")

    # load bias and flat
    logf.debug(f"Loading Bias: {bias_fits}")
    data_bias = fits.getdata(bias_fits)
    logf.debug(f"Loading Flat: {flat_fits}")
    data_flat = fits.getdata(flat_fits)

    # image shape
    ny, nx = data_bias.shape

    # observatory info, this is used in transfer JD to HJD
    site = coord.EarthLocation(lat=conf.site_lat * u.deg, lon=conf.site_lon * u.deg,
                               height=conf.site_ele * u.m)
    # object
    hdr = fits.getheader(raw_list[0])
    obj_ra = alt_coord[0] if alt_coord else hdr.get("RA", None)
    obj_dec = alt_coord[1] if alt_coord else hdr.get("DEC", None)
    logf.debug(f"Coord({obj_ra} {obj_dec}) from {'PARAM' if alt_coord else raw_list[0]}")
    if obj_ra is None or obj_dec is None:
        obj = None
    else:
        obj = coord.SkyCoord(obj_ra, obj_dec, unit=(u.hour, u.deg), frame="icrs")

    # load images and process
    for i, (rawf, bff) in zenum(raw_list, bf_fits_list):
        logf.debug(f"Loading {i+1:03d}/{nf:03d}: {rawf:40s}")

        # process data
        dat = (fits.getdata(rawf) - data_bias) / data_flat
        # border cut
        bc = conf.border_cut
        if bc > 0:
            dat = dat[bc:-bc, bc:-bc]

        # load and handle header
        hdr = fits.getheader(rawf)
        # force check the header
        _ = hdr.tostring()

        # get time
        obs_dt = hdr_dt(hdr)
        obs_jd = time.Time(obs_dt, format='isot', scale='utc', location=site)

        # center of the exposure
        expt = hdr.get("EXPTIME", 0.0)
        obs_jd += expt / 2 / 86400

        obs_mjd = obs_jd.mjd
        obs_lst = coord.Angle(lst(obs_mjd, site.lon.deg), u.hour)
        if obj:
            # jd + ra/dec to bjd, hjd
            ltt_bary = obs_jd.light_travel_time(obj)
            obs_bjd = (obs_jd.tdb + ltt_bary).jd
            obs_hjd = pyasl.helio_jd(obs_jd.jd - 2.4e6, obj.ra.deg, obj.dec.deg) + 2.4e6
            # az, alt, airmass
            obs_ha = obs_lst - obj.ra
            obs_az, obs_alt = azalt(site.lat.deg, obs_lst.hour, obj.ra.deg, obj.dec.deg)
            obs_am = airmass(site.lat.deg, obs_lst.hour, obj.ra.deg, obj.dec.deg)
        else:
            obs_bjd = obs_hjd = obs_jd
            obs_az = obs_alt = obs_am = 0.0
            obs_ha = coord.Angle(0, u.hour)

        # add ra, dec, lst, jd, mjd, bjd, hjd, az, alt,
        hdr.update({
            "RA": ra2hms(obj.ra),
            "DEC": dec2dms(obj.dec),
            "LST": ra2hms(obs_lst),
            "HA": ra2hms(obs_ha),
            "DATE-OBS": obs_dt,
            "TIME-OBS": obs_dt,
            "JD": obs_jd.jd,
            "MJD": obs_mjd,
            "BJD": obs_bjd,
            "HJD": obs_hjd,
            "AZ": obs_az,
            "ALT": obs_alt,
            "AIRMASS": obs_am,
            "SITEELEV": site.height.value,
            "SITELAT": site.lat.deg,
            "SITELONG": site.lon.deg,
        })
        # adjust crpix for border cut
        hdr.update({
            "CRPIX1": hdr.get("CRPIX1", nx // 2) - conf.border_cut,
            "CRPIX2": hdr.get("CRPIX2", ny // 2) - conf.border_cut,
        })

        # add process time to header
        hdr.update(BZERO=0)
        hdr.append(('PROCTIME', uttimestr()))

        # save new fits
        fits.writeto(bff, data=dat, header=hdr, overwrite=True)
        logf.debug(f"Writing {i+1:03d}/{nf:03d}: {bff:40s}")

    logf.info(f"{nf:3d} files corrected")
