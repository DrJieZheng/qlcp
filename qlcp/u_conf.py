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
        self.file_log = logging.INFO     # log level for file
        self.scr_log = logging.ERROR     # log level for screen
        # observatory
        self.site_lon = 117.57722         # longitude 117.34.38
        self.site_lat = 40.395833         # latitude +40.23.45
        self.site_ele = 960               # elevation above sea level
        self.site_tz  = 8                 # timezone
        # flat level limit
        self.flat_limit_low  =  5000      # low limit for flat level
        self.flat_limit_high = 50000      # high limit for flat level
        # image correction
        self.border_cut = 10              # cut border pixels
        # draw phot result
        self.draw_phot = False            # draw phot result or not
        self.draw_phot_err = 0.05          # max error of stars to be drawn
        # offset max distance
        self.offset_max_dis = 250         # max distance for offset
        # max calibration distance
        self.cali_max_dis = 10.0          # max distance for calibration
        # wcs setting
        self.wcs_max_err = 0.05           # mag-err limit for choose good stars
        self.wcs_max_n = 1000             # brightest n stars will be used
        self.wcs_min_n = 20               # brightest n stars will be used
        self.wcs_max_dis = 10             # pixels for last matching
        # flux cali setting
        self.flux_max_err = 0.025         # mag-err limit for choose calibrating stars
        self.flux_max_n = 1000            # n of the brightest stars will be used
        self.flux_max_dis = 10.0          # pixels for matching image and ref stars
        self.flux_chk_res = True          # plot residual check image or not

        # filename pattern
        self.patterns = [
            # UYUMa-0003I.fit  bias-0001.fits
            "(?P<obj>[^-_]*)-(?P<sn>[0-9]{3,6})(?P<band>[a-zA-Z]{0,1}).fit(s{0,1})",
            # flat_R_003.fit TCrB_V_001.fits
            "(?P<obj>[^-_]*)_(?P<band>[a-zA-Z]{0,1})(_{0,1})(?P<sn>[0-9]{3,6}).fit(s{0,1})",
            # flat_R_1_R_001.fit // for auto flat only
            "(?P<obj>flat)_(?P<band>[a-zA-Z])_(?P<sn>[0-9]{1,2})_[a-zA-Z]_([0-9]{3,6}).fit(s{0,1})"
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


class workmode:
    """A class for work mode constants and checking"""

    # input  files, v exists, x missing
    # output files, v exists, x no exists
    # stat   files, s exists or missing
    # Skip with warning
    # eXception Over -skip             v x s
    EXIST_ERROR   :int = 0b00000001  #   X X
    EXIST_SKIP    :int = 0b00000010  #   - -
    EXIST_OVER    :int = 0b00000100  #   O O
    EXIST_APPEND  :int = 0b00001000  #   O O
    MISS_ERROR    :int = 0b00100000  # X   X
    MISS_SKIP     :int = 0b01000000  # O   -

    def __init__(self, mode:int=MISS_SKIP+EXIST_APPEND):
        """init object"""
        self.mode_r = self.mode = mode

    def missing(self, filename:str, filetype:str, logf:logging.Logger):
        """
        Check the file is missing or not, log or raise exception
        return: False - OK, True - no but continue, Exception - No and error
        """
        if not os.path.isfile(filename):
            if self.mode_r & workmode.MISS_ERROR:
                logf.error(f"Stop missing {filetype}: `{filename}`")
                raise FileNotFoundError(f"Missing {filetype}: `{filename}`")
            else:  # WMODE.MISS_SKIP or not set
                logf.warning(f"Skip missing {filetype}: `{filename}`")
                return True
        else:
            return False

    def exists(self, filename:str, filetype:str, logf:logging.Logger):
        """
        Check the file is existing or not, log or raise exception,
        append should be transfer to over or skip in program
        return: True - OK, False - exists but continue, Exception - exists and error
        """
        if os.path.isfile(filename):
            if self.mode_r & workmode.EXIST_ERROR:
                logf.error(f"Stop existing {filetype}: `{filename}`")
                raise FileExistsError(f"Existing {filetype}: `{filename}`")
            elif self.mode_r & workmode.EXIST_OVER:
                logf.info(f"Overwrite existing {filetype}: `{filename}`")
                return False
            else:  # WMODE.EXIST_SKIP or not set or APPEND
                logf.warning(f"Skip existing {filetype}: `{filename}`")
                return True
        else:
            return False

    def reset_append(self, mode_to:int):
        """Replace append mode bit with another"""
        if self.mode & workmode.EXIST_APPEND:
            self.mode_r = self.mode & ~ workmode.EXIST_APPEND | mode_to
        else:
            self.mode_r = self.mode
