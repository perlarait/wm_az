#!/usr/bin/env python
# -*- coding: utf-8 -*-

# """WORMSSSSSS!"""

## __Import__
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import csv
import numpy as np
import scipy as sp
import shutil
from pylab import figure, axes, pie, title, show, savefig
# from bin import csvs
from bin import mats2
from bin import clnr
from bin import grph
from bin import rscrps

## enviromental variables
zplim = 2
znlim = -2
zxlim = 3

## __verify structure__

## initial data shuttling
os.mkdir('../final')
# os.mkdir('../final/plots')
shutil.move('../data', '../final/')
os.mkdir('../data')

## __Import data Recursively through folders__
# repX = raw("which replicate?")
bb = glob.glob('../final/data/*/*.csv')

## __Clean up__
##	_linearize_
## _platenames_
clnr.pltnmij(bb)
clnr.pltnmr(bb)
## injects rep names
clnr.repij(bb)
## _drops,changes_
clnr.clnrr(bb)
##  _well_let/num_
clnr.idxr(bb)
#  _annotate +/-ctrls_
clnr.annt(bb)

## __Calculations__
# calculates z-score
rscrps.zscrn(bb)
rscrps.zscrp(bb)

## cln again
clnr.clnr2(bb)

## _Creates master csv_
clnr.mrgr(bb)

#Remove CTRL Outliers (Automatic exclusion) (need to do recursively by plate)
mats2.ex('../final/all.csv', zxlim)

# re-calculates z-score
dd = glob.glob('../final/*_excl.csv')
rscrps.zscrn(dd)
rscrps.zscrp(dd)

## __Hit&Tox__
mats2.hits('../final/all_excl.csv', zplim)
mats2.tox('../final/all_excl.csv', znlim)
mats2.countr('../final/all_excl_all_hits.csv')
mats2.countr('../final/all_excl_all_tox.csv')
#mats2.autopl('../final/all_excl_all_hits.csv')
#mats2.autopl('../final/all_excl_all_tox.csv')
#mats2.autopl2('../final/all_excl_all_hits.csv')
#mats2.autopl2('../final/all_excl_all_tox.csv')

## Correlate Echo exceptions

## Plotting & Sum
a = '../final/all_excl.csv'
b = '../final/all.csv'

## Generate summaries
rscrps.summr(a)
rscrps.summr(b)
rscrps.summr_pl(a)
rscrps.summr_pl(b)
rscrps.summr_rp(a)
rscrps.summr_rp(b)

## __between plate analysis__
#	Interested in deviation between positive and neg controls
#	Interested in size of effect between hits and tox

## box plots of ctrls_
rscrps.ctrlbxa(b)
rscrps.ctrlbxz(a)
## raw area plots
rscrps.pltala(b)
## all zscore plots
rscrps.pltalz(a)
## plots mean lines
rscrps.pltalmpl(b)
rscrps.pltalmnl(b)

## exp zscore w/perplate & rep mean o/ctrls

## __Finishing__
a = 'hit limit = ' + str(zplim)
b = 'tox limit = ' + str(znlim)
c = 'excl limit = ' + str(zxlim)
f = open('meta.txt', 'w')
f.write(a + '\n' + b + '\n' + c)
