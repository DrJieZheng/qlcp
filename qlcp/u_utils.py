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
import time
import pickle
import numpy as np
# from scipy import stats as st
# from astropy.time import Time


def loadlist(listfile:str, base_path:str="", suffix:str="", separate_folder:bool=False)->list[str]:
    """
    Load file list from list file, add base path and suffix to each filename
    :param listfile:
    :param base_path: add base path to filenames in list
    :param suffix: if filename not ends with fits/fit/gz, then .fits will be append
    :param separate_folder: make a separate folder for each file
    :return: a list of filename
    """

    def _get_ext_(fn):
        # an inner function to split base name and ext
        # if the last ext is gz, then the former section is also part of the ext
        # 4x: xx.fits.gz xx.tar.gz
        sp = os.path.splitext(fn)
        base = sp[0]
        ext = sp[1]
        if ext == ".gz":
            spsp = os.path.splitext(base)
            ext = spsp[1] + ext
            base = spsp[0]
        return os.path.basename(base), ext

    # the last name of the list, if sperate_folder set, this is the folder name
    lstname = os.path.basename(os.path.splitext(listfile)[0])

    # load original list
    flst = [f.strip() for f in open(listfile, "r").readlines()]

    # format the suffix
    if suffix is not None and suffix != "" and not suffix.startswith("."):
        suffix = "." + suffix

    ori_path = [os.path.dirname(f) + "/" for f in flst]  # with ending /
    base_name = [_get_ext_(f)[0] for f in flst]  # pure name
    ori_ext = [_get_ext_(f)[1] for f in flst]  # ext with .

    # new suffix if provided
    if suffix is None or suffix.strip() == "":
        new_ext = ori_ext.copy()
    elif suffix in ('.fit', '.fits', '.fit.gz', '.fits.gz', ):
        # if new and ori suffix are both aliases of fits, keep original
        # different to last version, middle-fix is now a part of new suffix
        new_ext = [f if f in '.fit .fits .fit.gz .fits.gz' else suffix for f in ori_ext]
    else:
        new_ext = [suffix] * len(ori_ext)

    # if base path provided, use it, else keep original path
    if base_path is not None and base_path != "":
        new_path = [base_path] * len(flst)
    else:
        new_path = ori_path.copy()

    # add last name if separate_folder is set
    if separate_folder:
        new_lst = [f"{p}/{lstname}/{f}{e}" for p, f, e in zip(new_path, base_name, new_ext)]
    else:
        new_lst = [f"{p}/{f}{e}" for p, f, e in zip(new_path, base_name, new_ext)]

    return new_lst


def rm_ix(ix:list[int], *arr):
    """Remove a[ix] from all a in arr, no return, all action on same array"""
    # reverse sort ix
    ix.sort(reverse=True)
    # remove array items
    for i in ix:
        for a in arr:
            del a[i]


def zenum(*arr):
    """combine of enumerate and zip"""
    return enumerate(zip(*arr))


def list_exist(files:list[str])->list[bool]:
    """
    Check existance of a list of files
    :param files: a list containing filenames
    :return: a list of True/False
    """
    e = [os.path.isfile(f) for f in files]
    return e


def localtimestr() -> str:
    """
    Generate a string of current local time
    :return:
    """
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime())


def uttimestr() -> str:
    """
    Generate a string of current ut time
    :return:
    """
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())


def str2mjd(s:str) -> float:
    """
    ISO datetime string to mjd
    :param s: yyyy-mm-ddThh:mm:ss
    :return: mjd
    """
    # 40587: 1970-01-01T00:00:00
    return time.mktime(time.strptime(s, "%Y-%m-%dT%H:%M:%S")) / 86400.0 + 40587.0


def hdr_dt(hdr):
    """Find the datetime string from header, guess where the keyword is"""
    # 220607 use regular-express to find obs date and time
    dt_str = hdr.get("DATE-OBS", "") + "T" + hdr.get("TIME-OBS", "")
    dd = re.search("([0-9]{4}-[01][0-9]-[0-3][0-9])", dt_str)
    if dd:
        dd = dd[0]
    else:
        dd = re.search("([0-3][0-9]/[01][0-9]-[0-9][0-9])", dt_str)
        if dd:
            dd = dd[0]
            dd = f"20{dd[6:8]}-{dd[3:5]}-{dd[0:2]}"
        else:
            dd = "2019-10-18"
    tt = re.search("([012][0-9]:[0-5][0-9]:[0-5][0-9](\\.[0-9]*)?)", dt_str)
    if tt:
        tt = tt[0]
    else:
        tt = "00:00:00.000"
    dt_str = f"{dd}T{tt}"
    return dt_str


def meanclip(dat, nsigma=3.0, func=np.nanmean)->(float,float):
    """
    Compute clipped median and sigma of dat
    :param dat: data, can be list, tuple or np array, 1-d or n-d
    :param nsigma: how many sigma used in clipping
    :param func: method used to evaluate the median/mean
    :return:
    """
    if len(dat) == 0:
        m, s = np.nan, np.nan
    else:
        m1 = np.nanmedian(dat)
        s1 = np.nanstd(dat)
        g1 = np.abs(dat - m1) < nsigma * s1
        m2 = np.nanmedian(dat[g1])
        s2 = np.nanstd(dat[g1])
        g2 = np.abs(dat - m2) < nsigma * s2
        m = func(dat[g2])
        s = np.nanstd(dat[g2])
    return m, s


def unmatched(nall, matched):
    """
    选出未匹配到的目标
    :param nall:
    :param matched:
    :return:
    """
    # tag = np.ones(nall, bool)
    # tag[matched] = False
    # return np.where(tag)[0]
    return np.array(list( set(range(nall)) - set(matched) ))


def subset(ixall, matched):
    """
    find unmatched indices
    :param ixall:
    :param matched:
    :return:
    """
    return np.array(list( set(ixall) - set(matched) ))


def ra2hms(x):
    """
    from coord.ra to hms format
    """
    xx = x.hms
    return f"{int(xx.h):02d}:{int(xx.m):02d}:{xx.s:05.2f}"


def dec2dms(x):
    """
    from coord.dec to signed dms format
    """
    xx = x.signed_dms
    return f"{('+' if xx.sign>0 else '-'):1s}{int(xx.d):02d}:{int(xx.m):02d}:{xx.s:04.1f}"


def pkl_dump(filename:str, *dat):
    """dump variables into pickle file"""
    with open(filename, "wb") as ff:
        pickle.dump(dat, ff)


def pkl_load(filename:str):
    """load var from pickle file"""
    with open(filename, "rb") as ff:
        dat = pickle.load(ff)
    return dat


def cat2txt(filename:str, cat:np.ndarray, fmt:dict={}):
    """dump catalog into txt file"""
    # get column names, scalar or array, and array size
    scol = [k for k in cat.dtype.names if len(cat.dtype[k].shape) == 0]
    acol = [k for k in cat.dtype.names if len(cat.dtype[k].shape) == 1]
    ncol = cat.dtype[acol[0]].shape[0] if acol else 0
    with open(filename, "w") as ff:
        # write header
        ff.write("#")
        for k in scol:
            ff.write(f" {k}")
        for i in range(ncol):
            for k in acol:
                ff.write(f" {k}_{i+1:02d}")
        ff.write("\n")

        # write data
        for s in cat:
            for k in scol:
                ff.write(f" {s[k]}")
            for i in range(ncol):
                for k in acol:
                    ff.write(f" {s[k][i]}")
            ff.write("\n")
