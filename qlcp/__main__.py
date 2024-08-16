# -*- coding: utf-8 -*-
"""
    v1 201901, Dr. Jie Zheng, Beijing & Xinglong, NAOC
    v2 202101, Dr. Jie Zheng & Dr./Prof. Linqiao Jiang
    v3 202201, Zheng & Jiang
    v4 202304, Upgrade, restructure, Zheng & Jiang
    Quick_Light_Curve_Pipeline
"""


import sys
import argparse
import os
from .u_workmode import workmode


def main():
    """
    A cli tool to run the pipeline.
    """
    if len(sys.argv) == 1:
        print("""Quick Light Curve Pipeline
Usage: python -m qlcp command arguments
Commands:
    Get the x,y coordinates of stars by clicking on the image:
        py -m qlcp getxy <fits_file>
        py -m qlcp getxy -h
    Light Curve Process:
        py -m qlcp <action> <raw_dir> <output_dir> [<optons>]
        py -m qlcp -h
""")
    elif sys.argv[1].lower() == "getxy":
        parser = argparse.ArgumentParser(description="XY Picker")
        parser.add_argument("action", type=str, 
            help="getxy")
        parser.add_argument("fitsfile", type=str, 
            help="FITS file to be processed")
        parser.add_argument("--size", type=int, default=10, 
            help="Picking radius in pixels")
        args = parser.parse_args()

        from .getxy import getxy
        if not os.path.isfile(args.fitsfile):
            raise FileNotFoundError(f"Cannot find file {args.fitsfile}")
        getxy(args.fitsfile, size=args.size)

        return
    else:
        def str_or_int(value):
            """test str or int, for argparse"""
            try:
                return int(value)
            except ValueError:
                return value

        def pos_xy(coord_str):
            """converts a string 'x,y' to a float tuple (x, y)"""
            try:
                x, y = coord_str.split(',')
                return float(x.strip()), float(y.strip())
            except ValueError:
                raise argparse.ArgumentTypeError(f"Cannot convert '{coord_str}' to a x,y position.")

        # parse arguments
        parser = argparse.ArgumentParser(description="Quick Light Curve Pipeline")
        parser.add_argument("action", type=str, default="lbfiopcd", required=False,
            help="""Action to be performed, steps, each char stands for one step
    l = List
    b = Bias comb
    f = Flat comb
    i = Image correct
    o = Offset
    p = find stars and Photometry
    s = photometry with sep
    w = Wcs
    k = picK ref stars
    K = loading last picKing result
    c = Catalog
    d = Differential flux""")
        parser.add_argument("raw", type=str, 
            help="Raw data directory")
        parser.add_argument("red", type=str, 
            help="Reduced output directory")
        parser.add_argument("-o", "--obj", type=str, nargs="*", default=None,
            help="Object(s) to process")
        parser.add_argument("-b", "--band", type=str, nargs="*", default=None, 
            help="Band(s) to process")
        parser.add_argument("-i", "--base", type=str_or_int, default=None, 
            help="Base image index or filename")
        parser.add_argument("-a", "--aper", type=float, nargs="*", default=None, 
            help="Aperture(s)")
        parser.add_argument("-p", "--starxy", type=pos_xy, nargs="*", default=None,
            help="Target positions")
        parser.add_argument("-t", "--target", type=int, nargs="*", default=None,
            help="Index of target star")
        parser.add_argument("-r", "--ref", type=int, nargs="*", default=None,
            help="Index of reference star")
        parser.add_argument("-c", "--check", type=int, nargs="*", default=None,
            help="Index of check star")
        parser.add_argument("--use_bias", type=str, default=None,
            help="File name of use bias image")
        parser.add_argument("--use_flat", type=str, default=None,
            help="File name of use flat image")
        parser.add_argument("--alt_bias", type=str, default=None,
            help="File name of alternative bias image")
        parser.add_argument("--alt_flat", type=str, default=None,
            help="File name of alternative flat image")
        parser.add_argument("--alt_coord", type=str, default=None,
            help="alternative coordinate of the object")
        parser.add_argument("--log_screen", type=str.lower, default="info",
            choices=["error", "warning", "info", "debug"],
            help="Log level for screen output: error warning info debug")
        parser.add_argument("--log_file", type=str.lower, default="debug",
            choices=["error", "warning", "info", "debug"],
            help="Log level for file output: error warning info debug")
        parser.add_argument("--ini_file", type=str,  nargs="*", default=None,
            help="Configuration file(s)")
        parser.add_argument("--file_miss", type=str.lower, default="skip",
            choices=["error", "skip"],
            help="Action for missing files: error, skip")
        parser.add_argument("--file_exists", type=str.lower, default="append",
            choices=["error", "skip", "over", "append"],
            help="Action for existing files: error, skip, over, append")
        # configuration parameters
        parser.add_argument("--conf_site_lon"       , type=float, default=None, 
			help="longitude, default is 117.34.38")
        parser.add_argument("--conf_site_lat"       , type=float, default=None, 
			help="latitude, default is  +40.23.45")
        parser.add_argument("--conf_site_ele"       , type=float, default=None, 
			help="elevation above sea level, default is  940")
        parser.add_argument("--conf_site_tz"        , type=float, default=None, 
			help="timezone, default is  +8")
        parser.add_argument("--conf_flat_limit_low" , type=int  , default=None, 
			help="low limit for flat level")
        parser.add_argument("--conf_flat_limit_high", type=int  , default=None, 
			help="high limit for flat level")
        parser.add_argument("--conf_border_cut"     , type=int  , default=None, 
			help="cut border pixels")
        parser.add_argument("--conf_draw_phot"      , type=bool , default=None, 
			help="draw phot result or not")
        parser.add_argument("--conf_draw_phot_err"  , type=float, default=None, 
			help="max error of stars to be drawn")
        parser.add_argument("--conf_offset_max_dis" , type=int  , default=None, 
			help="max distance for offset")
        parser.add_argument("--conf_match_max_dis"  , type=float, default=None, 
			help="max distance for object matching")
        parser.add_argument("--conf_pick_err_max"   , type=float, default=None, 
			help="max error for pick stars")
        parser.add_argument("--conf_pick_star_n"    , type=int  , default=None, 
			help="max number of stars")
        parser.add_argument("--conf_pick_bad_img"   , type=float, default=None, 
			help="factor of bad image")
        parser.add_argument("--conf_pick_bad_max"   , type=float, default=None, 
			help="factor of bad stars")
        parser.add_argument("--conf_pick_var_n"     , type=int  , default=None, 
			help="max number of variable stars")
        parser.add_argument("--conf_pick_var_std"   , type=float, default=None, 
			help="std of the variance stars")
        parser.add_argument("--conf_pick_var_rad"   , type=float, default=None, 
			help="radius of the variance stars")
        parser.add_argument("--conf_pick_ref_n"     , type=int  , default=None, 
			help="max snumber of reference stars")
        parser.add_argument("--conf_pick_ref_std"   , type=float, default=None, 
			help="std of the reference stars")
        parser.add_argument("--conf_pick_ref_dif"   , type=float, default=None, 
			help="max_min limit of the reference stars")
        parser.add_argument("--conf_wcs_max_err"    , type=float, default=None, 
			help="mag_err limit for choose good stars")
        parser.add_argument("--conf_wcs_max_n"      , type=int  , default=None, 
			help="brightest n stars will be used")
        parser.add_argument("--conf_wcs_min_n"      , type=int  , default=None, 
			help="brightest n stars will be used")
        parser.add_argument("--conf_wcs_max_dis"    , type=int  , default=None, 
			help="pixels for last matching")
        parser.add_argument("--conf_flux_max_err"   , type=float, default=None, 
			help="mag_err limit for choose calibrating stars")
        parser.add_argument("--conf_flux_max_n"     , type=int  , default=None, 
			help="n of the brightest stars will be used")
        parser.add_argument("--conf_flux_max_dis"   , type=float, default=None, 
			help="pixels for matching image and ref stars")
        parser.add_argument("--conf_flux_chk_res"   , type=bool , default=None, 
			help="plot residual check image or not")
        
        args = parser.parse_args()
        # print(args)

        # process arguments to parameters for run
        kwarg = {}
        kwarg["raw_dir"] = args.raw
        kwarg["red_dir"] = args.red
        kwarg["steps"] = args.action
        kwarg["obj"] = args.obj
        kwarg["band"] = args.band
        kwarg["use_bias"] = args.use_bias
        kwarg["use_flat"] = args.use_flat
        kwarg["alt_bias"] = args.alt_bias
        kwarg["alt_flat"] = args.alt_flat
        kwarg["alt_coord"] = args.alt_coord
        kwarg["base_img"] = args.base
        kwarg["aper"] = args.aper
        kwarg["starxy"] = args.starxy
        kwarg["ind_tgt"] = args.target
        kwarg["ind_ref"] = args.ref
        kwarg["ind_chk"] = args.check
        fe = {
            "error":workmode.EXIST_ERROR,
            "skip":workmode.EXIST_SKIP,
            "over":workmode.EXIST_OVER,
            "append":workmode.EXIST_APPEND,
            }.get(args.file_exists, 0)
        fm = {
            "error":workmode.MISS_ERROR,
            "skip":workmode.MISS_SKIP,
            }.get(args.file_miss, 0)
        kwarg["mode"] = fe + fm
        kwarg["ini_file"] = args.ini_file
        # config parameters
        for k in args.__dict__:
            v = args.__dict__[k]
            if k.startswith("conf_") and v is not None:
                kwarg[k[5:]] = v

        from .j_run import run



if __name__ == "__main__":
    main()
