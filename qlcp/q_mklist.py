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
from glob import glob
from .u_conf import config
from .u_workmode import workmode
from .u_log import init_logger


def mklist(
        conf:config,
        raw_dir:str,
        red_dir:str,
        obj:str=None,
        band:str=None,
        mode:workmode=workmode(),
        save_list_file:bool=True,
):
    """
    list making function, return a dict of object-band
    :param conf: config object
    :param raw_dir: raw files dir
    :param red_dir: red files dir
    :param obj: object(s) to process
    :param band: band(s) to process
    :param mode: input files missing or output existence mode
    :param save_list_file: save list or just return object-band dict tree
    :returns: the band-object dict, key is band, and value is list of objects
    """
    os.makedirs(f"{red_dir}/lst/", exist_ok=True)
    logf = init_logger("makelist", f"{red_dir}/log/mklist.log", conf)

    ###############################################################################

    def _filename_parse_(fn):
        """
        parse filename into object, band and sn
        """
        o, b, n = "UNKNOWN", "X", 0
        # try all patterns, if any matches, use it
        for p in conf.patterns:
            x = re.search(p, fn)
            if x:
                o = x.groupdict().get("obj", "UNKNOWN")
                b = x.groupdict().get("band", "X")
                n = int(x.groupdict().get("sn", "0"))
                break
        return o, b, n

    # get files
    files = glob(f"{raw_dir}/*.fit") + glob(f"{raw_dir}/*.fits")
    files = [os.path.basename(f) for f in files]
    files.sort()
    logf.info(f"{len(files):3d} fits files found in {raw_dir}")

    # file name parse
    fileinfo = [_filename_parse_(f) for f in files]
    objs = ["bias" if "bias" in f[0].lower() else "flat" if "flat" in f[0].lower() else f[0]
            for f in fileinfo]
    bands = ["" if "bias" in f[0].lower() else f[1] for f in fileinfo]

    # handle specified obj and band
    if isinstance(obj, str):
        objall = [obj, 'bias', 'flat']
    elif isinstance(obj, (list, tuple)):
        objall = list(obj) + ['bias', 'flat']
    else:
        objall = None
    if isinstance(band,  (str, list, tuple)):
        bandall = list(band) + ['']
    else:
        bandall = None
    # dispatch files
    file4band = {}  # [band][obj] -> list of filenames
    list4band = {}  # [band] -> list of objects
    file4obj = {}  # [obj] -> list of filenames
    list4obj = {}  # [obj] -> list of bands
    for f, o, b in zip(files, objs, bands):
        # skip files not as specified
        if objall and o not in objall or bandall and b not in bandall:
            continue
        
        # for a new bands, create a sub dict
        if b not in file4band:
            file4band[b] = {}
            list4band[b] = []
        if o not in file4band[b]:
            file4band[b][o] = []
            list4band[b].append(o)
        file4band[b][o].append(f)

        # for a new object, create a sub dict
        if o not in file4obj:
            file4obj[o] = {}
            list4obj[o] = []
        if b not in file4obj[o]:
            file4obj[o][b] = []
            list4obj[o].append(b)
        file4obj[o][b].append(f)

    # # find bad bands, with only flat but no objects
    # bad_band = []
    # for b, ff in file4band.items():
    #     if len(ff) == 1 and "flat" in ff:
    #         bad_band.append(b)
    # # remove bad bands
    # for b in bad_band:
    #     del file4band[b]
    #     del list4band[b]
    #     if b in file4obj["flat"]:
    #         file4obj[o].remove(b)
    
    logf.info(f"{list4obj}")

    # updating lists
    if save_list_file:
        # output lists, lst contain all files, ldt contains file and adding time
        for b in file4band:
            for o in file4band[b]:
                bb = ('_' if b else '') + b
                f_lst = file4band[b][o]
                # normal list
                lst_fn = f"{red_dir}/lst/{o}{bb}.lst"
                with open(lst_fn, "w") as ff:
                    ff.write("\n".join(f_lst))
                logf.info(f"\t{len(f_lst):3d} ==> {lst_fn}")

    return list4obj

