# -*- coding: utf-8 -*-
"""
    v1 201901, Dr. Jie Zheng, Beijing & Xinglong, NAOC
    v2 202101, Dr. Jie Zheng & Dr./Prof. Linqiao Jiang
    v3 202201, Zheng & Jiang
    v4 202304, Upgrade, restructure, Zheng & Jiang
    Quick_Light_Curve_Pipeline
"""


import os
import logging
import configparser


class config:
    """
    A class of config
    In default mode, an instance of this class will be created and act as system config
    but user can set other parameters
    """

    def __init__(self, ini_file:list|tuple|str=None, extra_conf=None):
        """
        Init settings
        :param ini_file: external ini file
        :param extra_conf: extra config items, as dict
        """

        # path of the program
        self.here = os.path.realpath(os.path.dirname(__file__)) + "/"   # program path
        # default log level
        self.file_log = logging.DEBUG           # log level for file
        self.scr_log = logging.INFO             # log level for screen
        # observatory
        self.site_lon           = 117.57722     # longitude 117.34.38
        self.site_lat           = 40.395833     # latitude +40.23.45
        self.site_ele           = 960           # elevation above sea level
        self.site_tz            = 8             # timezone
        # flat level limit
        self.flat_limit_low     =  5000         # low limit for flat level
        self.flat_limit_high    = 50000         # high limit for flat level
        # image correction
        self.border_cut         = 0             # cut border pixels
        # draw phot result
        self.draw_phot          = False         # draw phot result or not
        self.draw_phot_err      = 0.05          # max error of stars to be drawn
        # offset max distance
        self.offset_max_dis     = 250           # max distance for offset
        # max matching distance
        self.match_max_dis      = 10.0          # max distance for object matching
        # star pick
        self.pick_err_max       = 0.10          # max error for pick stars
        self.pick_star_n        = 100           # max number of stars
        self.pick_bad_max       = 0.2           # factor of bad stars
        self.pick_var_n         = 5             # max number of variable stars
        self.pick_var_std       = 0.10          # std of the variance stars
        self.pick_var_rad       = 0.5           # radius of the variance stars
        self.pick_ref_n         = 10            # max snumber of reference stars
        self.pick_ref_std       = 0.02          # std of the reference stars
        self.pick_ref_dif       = 0.10          # max-min limit of the reference stars
        # wcs setting
        self.wcs_max_err        = 0.05          # mag-err limit for choose good stars
        self.wcs_max_n          = 1000          # brightest n stars will be used
        self.wcs_min_n          = 20            # brightest n stars will be used
        self.wcs_max_dis        = 10            # pixels for last matching
        # flux cali setting
        self.flux_max_err       = 0.025         # mag-err limit for choose calibrating stars
        self.flux_max_n         = 1000          # n of the brightest stars will be used
        self.flux_max_dis       = 10.0          # pixels for matching image and ref stars
        self.flux_chk_res       = True          # plot residual check image or not

        # filename pattern
        self.patterns = [
            # UYUMa-0003I.fit  bias-0001.fits
            "(?P<date>[0-9]{0,8})(?P<obj>[^-_]*)-(?P<sn>[0-9]{3,6})(?P<band>[a-zA-Z]{0,1}).fit(s{0,1})",
            # flat_R_003.fit TCrB_V_001.fits
            "(?P<date>[0-9]{0,8})(?P<obj>[^-_]*)_(?P<band>[a-zA-Z]{0,1})(_{0,1})(?P<sn>[0-9]{3,6}).fit(s{0,1})",
            # flat_R_1_R_001.fit // for auto flat only
            "(?P<date>[0-9]{0,8})(?P<obj>flat)_(?P<band>[a-zA-Z])_(?P<sn>[0-9]{1,2})_[a-zA-Z]_([0-9]{3,6}).fit(s{0,1})"
        ]

        # load external ini file
        if ini_file:
            # if not None or '', transfer to list
            if isinstance(ini_file, (list, tuple)):
                ini_file = (ini_file,)
            for ff in ini_file:
                # process ini files one by one
                if isinstance(ff, str) and os.path.isfile(ff):
                    self._load_ini_(ff)

        # extra config
        if isinstance(extra_conf, dict):
            self.__dict__.update(extra_conf)

    def _load_ini_(self, f):
        """Load data from ini file (f)"""

        # type convert function
        def totype(s, e):
            """Convert string s to the type of e, if failed, use e as the result"""
            try:
                if isinstance(e, int):
                    v = int(s)
                elif isinstance(e, float):
                    v = float(s)
                elif isinstance(e, bool):
                    if s.lower() in "yes true":
                        v = True
                    elif s.lower() in "no false":
                        v = False
                    else:
                        v = e
                elif isinstance(e, (list, tuple)):
                    v = [e.strip() for e in s.split(",")]
                else:
                    v = s
            except Exception(e):
                v = e
            return v

        # construct a parser
        cp = configparser.ConfigParser()
        # loading
        with open(f) as ff:
            cp.read_string("[DEFAULT]\n" + ff.read())
            # check existing keys only
            for k, x in self.__dict__.items():
                if cp.has_option("DEFAULT", k):
                    # if provided, transfer to correct type
                    v = totype(cp.get("DEFAULT", k), x)
                    self.__dict__[k] = v
