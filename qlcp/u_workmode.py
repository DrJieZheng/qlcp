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
from .u_utils import fnbase


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
        # reset mode and original mode
        self.mode_r = self.mode = mode
        # postpone the report of warning files
        self.lazy_files = {}
    
    def __str__(self):
        """print work mode"""
        return f"mode: {self.mode:08b} moder: {self.mode_r:08b}"
    
    def start_lazy(self):
        """Start lazy mode, only report the warning files"""
        self.lazy_files.clear()

    def lazy_warning(self, reason:str, filename:str, logf:logging.Logger):
        """Add a warning file to lazy list"""
        # if logf is Logger, no lazy mode
        if logf:
            logf.warning(f"{reason}: {filename}")
        else:
            if reason not in self.lazy_files:
                self.lazy_files[reason] = []
            self.lazy_files[reason].append(filename)
    
    def end_lazy(self, logf:logging.Logger):
        """End lazy mode and report the warning files"""
        for reason in self.lazy_files:
            logf.warning(f"{reason} ({len(self.lazy_files[reason])}): " +
            ", ".join(fnbase(self.lazy_files[reason])))
        # remove the lazy files
        self.lazy_files.clear()

    def missing(self, filename:str, filetype:str, logf:logging.Logger=None):
        """
        Check the file is missing or not, log or raise exception
        return: False - OK, True - no but continue, Exception - No and error
        """
        if not os.path.isfile(filename):
            if self.mode_r & workmode.MISS_ERROR:
                logf.error(f"Stop missing {filetype}: `{filename}`")
                raise FileNotFoundError(f"Missing {filetype}: `{filename}`")
            else:  # WMODE.MISS_SKIP or not set
                # logf.debug(f"Skip missing {filetype}: `{filename}`")
                self.lazy_warning(f"Skip missing {filetype}", filename, logf)
                return True
        else:
            return False

    def exists(self, filename:str, filetype:str, logf:logging.Logger=None):
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
                # logf.debug(f"Overwrite existing {filetype}: `{filename}`")
                self.lazy_warning(f"Overwrite existing {filetype}", filename, logf)
                return False
            else:  # WMODE.EXIST_SKIP or not set or APPEND
                # logf.debug(f"Skip existing {filetype}: `{filename}`")
                self.lazy_warning(f"Skip existing {filetype}", filename, logf)
                return True
        else:
            return False

    def reset_append(self, mode_to:int):
        """Replace append mode bit with another"""
        if self.mode & workmode.EXIST_APPEND:
            self.mode_r = self.mode & ~ workmode.EXIST_APPEND | mode_to
        else:
            self.mode_r = self.mode
