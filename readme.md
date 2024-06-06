# QLCP: Quick Light Curve Pipeline v2024.06

## Description:

QLCP is a quick, easy-to-use pipeline for processing light curves from images observed by any ground based telescope.

## Installation:

QLCP is available on PyPI and can be installed with pip:

```bash
pip install qlcp
```

## Usage:

```python
import qlcp
from qlcp import workmode
import logging

rawd = '/path/to/raw/data/'
redd = '/path/to/reduced/data/'

qlcp.run(
    None,                   # ini file, if none, use default
    rawd,                   # raw data directory
    redd,                   # reduced data directory
    steps='lbfoipwkcdg',    # steps to do
    aper=[3,5,7,9],         # apertures, up to 9 apertures
    starxy=[                # star positions on base image
        (927,1018),         # target
        (855,920),          # ref1
        (1255,861),         # chk
        (1107,1349),        # ref2
        (1161,1434),        # ref3
        (1220,1289),        # ref4
        (688,1050),         # ref5
    ],
    ind_tgt=[0,6],          # target index
    ind_ref=[1,3,4,5,6],    # reference index
    ind_chk=[2,3,4],        # checker index
    mode=workmode(workmode.EXIST_OVER+workmode.MISS_ERROR),
                            # how to handle overwrite and missing files
    offset_max_dis=100,     # max offset for offset
    se_cmd='sex',           # command to run SExtractor
    cali_max_dis=20.0,      # max distance for identify calibration stars
    scr_log=logging.DEBUG,  # logging level
    draw_phot=False,        # draw photometry result or not
)
```

