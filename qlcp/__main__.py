# -*- coding: utf-8 -*-
"""
    v1 201901, Dr. Jie Zheng, Beijing & Xinglong, NAOC
    v2 202101, Dr. Jie Zheng & Dr./Prof. Linqiao Jiang
    v3 202201, Zheng & Jiang
    v4 202304, Upgrade, restructure, Zheng & Jiang
    Quick_Light_Curve_Pipeline
"""


import sys


def main():
    """
    A cli tool to run the pipeline.
    """
    if len(sys.argv) == 1:
        print("""
            Usage: python -m qlcp command arguments
            Commands:
              getxy <fits_file>  Get the x,y coordinates of stars by clicking on the image.
              <action> <raw_dir> <output_dir> [<optons>]
        """)
    elif sys.argv[1].lower() == "getxy":
        if len(sys.argv) == 3:
            from .getxy import getxy
            getxy(sys.argv[2])
        else:
            print("Usage: python -m qlcp getxy <fits_file>")
        return
    else:
        def str_or_int(value):
            """test str or int, for argparse"""
            try:
                return int(value)
            except ValueError:
                return value
        
        import argparse

        def pos_xy(coord_str):
            """converts a string 'x,y' to a float tuple (x, y)"""
            try:
                x, y = coord_str.split(',')
                return float(x.strip()), float(y.strip())
            except ValueError:
                raise argparse.ArgumentTypeError(f"Cannot convert '{coord_str}' to a x,y position.")

        
        # parse arguments
        import argparse
        parser = argparse.ArgumentParser(description="Quick Light Curve Pipeline")
        parser.add_argument("action", type=str, 
            help="Action to be performed")
        parser.add_argument("raw", type=str, 
            help="Raw data directory")
        parser.add_argument("red", type=str, 
            help="Reduced output directory")
        parser.add_argument("-o", "--obj", type=str, nargs="*", default=None,
            help="Object(s) to process")
        parser.add_argument("-b", "--band", type=str, default=None, 
            help="Band(s) to process")
        parser.add_argument("-i", "--base", type=str_or_int, default=None, 
            help="Base image index or filename")
        parser.add_argument("-a", "--aper", type=float, nargs="*", default=None, 
            help="Aperture radius")
        parser.add_argument("-p", "--starxy", type=pos_xy, nargs="*", default=None,
            help="Target positions")
        parser.add_argument("-t", "--target", type=int, nargs="*", default=None,
            help="Index of target star")
        parser.add_argument("-r", "--ref", type=int, nargs="*", default=None,
            help="Index of reference star")
        parser.add_argument("-c", "--check", type=int, nargs="*", default=None,
            help="Index of check star")
        parser.add_argument("--use-bias", type=str, default=None,
            help="File name of use bias image")
        parser.add_argument("--use-flat", type=str, default=None,
            help="File name of use flat image")
        parser.add_argument("--alt-bias", type=str, default=None,
            help="File name of alternative bias image")
        parser.add_argument("--alt-flat", type=str, default=None,
            help="File name of alternative flat image")
        parser.add_argument("--alt-coord", type=str, default=None,
            help="alternative coordinate of the object")
        parser.add_argument("--log-screen", type=str.lower, default="info",
            choices=["error", "warning", "info", "debug"],
            help="Log level for screen output: error warning info debug")
        parser.add_argument("--log-file", type=str.lower, default="debug",
            choices=["error", "warning", "info", "debug"],
            help="Log level for file output: error warning info debug")
        parser.add_argument("--ini-file", type=str,  nargs="*", default=None,
            help="Configuration file(s)")
        parser.add_argument("--file-miss", type=str.lower, default="skip",
            choices=["error", "skip"],
            help="Action for missing files: error, skip")
        parser.add_argument("--file-exists", type=str.lower, default="append",
            choices=["error", "skip", "over", "append"],
            help="Action for existing files: error, skip, over, append")
        
        args = parser.parse_args()
        print(args)

        # from .j_run import run


if __name__ == "__main__":
    main()
