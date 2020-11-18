# Cloudy_Helper

Some useful [Cloudy](http://www.nublado.org/) (Ferland et al. 2013) outputs management scripts

## Overview

List section provides name of scripts and their applications. For
details of each script, please read introduction at beginning of each
script. Please make sure you have latest Python3 to run these scripts.

## List of Scripts

CoolingFunction -- Example code used to show analytic fit functions (cooling
function and cooling rate) in Wang et al. (2013).

CoolingSpliter.py -- Use to get H&He, metal, and ee cooling( defined
in Wang et al. 2013) from CLOUDY’s 'save cooling each' and 'save
overview' outputs.

DataGather.py -- Use to merge temporary files generate by 'grid'
command into a whole file.

MultiCloudy.py -- Use to run cloudy when input variables need to be
expressed as a function (such as log(T) + log(n_H) = 10). Run cloudy
in a multiprocessing way.
