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


def init_logger(logger_name:str, log_file:str, conf,):
    """
    Init a logger
    :param logger_name:
    :param log_file:
    :param conf:
    :return:
    """
    # logger
    logger = logging.getLogger(logger_name)
    # basic level
    logger.setLevel(logging.DEBUG)
    # create log dir
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    # handlers: stream and file
    log_sh = logging.StreamHandler()
    log_fh = logging.FileHandler(filename=log_file, mode="a")
    # levels
    log_sh.setLevel(conf.scr_log)
    log_fh.setLevel(conf.file_log)
    # formatter
    log_sh_fmt = logging.Formatter(
        fmt="%(levelname)-7s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    log_fh_fmt = logging.Formatter(
        fmt="%(asctime)s|%(name)s|%(levelname)-7s|%(filename)-8s:%(lineno)03d|%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",)
    # set formatters to handlers
    log_sh.setFormatter(log_sh_fmt)
    log_fh.setFormatter(log_fh_fmt)
    # add handlers to logger
    logger.handlers.clear()
    logger.addHandler(log_sh)
    logger.addHandler(log_fh)

    return logger
