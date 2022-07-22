import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Polygon
import sys

# Define root directory
ROOT_DIR = Path().absolute()

# Get bash inputs to define src, data, image, and output directories
SRC_DIR = Path(sys.argv[1])
DATA_DIR = Path(sys.argv[2])
IM_DIR = Path(sys.argv[3])
GEO_DIR = Path(sys.argv[4])
OUT_DIR = Path(sys.argv[5])
# ROOT_DIR = Path('/home/durbank/scratch/heximap/')
# SRC_DIR = Path('/home/durbank/Documents/Research/Glaciers/heximap/src/')
# DATA_DIR = Path('/home/durbank/Documents/Research/Glaciers/GSLR/data/')
# IM_DIR = DATA_DIR.joinpath('hexagon/declass-ii/imagery/')
# GEO_DIR = DATA_DIR.joinpath('refDEMs')
# OUT_DIR = Path('home/durbank/Documents/Research/Glaciers/GSLR/scratch/outputs/')

# Load custom pyHEX module
sys.path.append(SRC_DIR.as_posix())
import pyHEX

# Define (and if needed) create tmp processing directory
TMP_DIR = ROOT_DIR.joinpath('tmp')
TMP_DIR.mkdir(parents=True, exist_ok=True)

# Get full path and filename for reference DEM 
# ***(currently gets only first DEM in folder)***
GeoRefPath = list(GEO_DIR.glob('*.tif'))[0]

################## USER EDITS MAY BE REQUIRED HERE ##################
## Subset RGI glaciers to glaciers of interest

# Load all RGI glacier outlines/data
RGI_files = list(DATA_DIR.joinpath("RGI-data").glob("*/*.shp"))
gdfs = [gpd.read_file(path) for path in RGI_files]
RGI_all = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), crs=gdfs[0].crs)

# Subset RGI data to region of interest (in this case, Southern Andes and glaciers larger than 1 km^2)
RGI_set = RGI_all[RGI_all["O1Region"] == "17"]
RGI_set = RGI_set.query("Area > 1")
#####################################################################

# Load hexagon metadata (assumes only one shp file in directory)
if ('declass-ii' in IM_DIR.as_posix()) and (
    'declass-iii' not in IM_DIR.as_posix()):
    hex_file = list(
        DATA_DIR.joinpath("hexagon/declass-ii/metadata").glob('*.shp'))[0]
elif IM_DIR.as_posix().contains('declass-iii'):
    hex_file = list(
        DATA_DIR.joinpath("hexagon/declass-iii/metadata").glob('*.shp'))[0]
else:
    print('Incorrect path names in data directory. Could not import hexagon metadata')

# Import hexagon metadata
hex_meta = gpd.read_file(hex_file)

# Subset hexagon images to those that intersect with the selected glaciers
bbox = RGI_set.total_bounds
bpoly = Polygon(
    ((bbox[0],bbox[1]), (bbox[2],bbox[1]), 
    (bbox[2],bbox[3]), (bbox[0],bbox[3])))
set_idx = hex_meta.intersects(bpoly)
hex_subset = hex_meta.loc[set_idx]

# Generate hexagon pairs
hex_subset, hex_pairs = pyHEX.pairify(hex_subset)

################## USER EDITS MAY BE REQUIRED HERE ##################
## Filter/Subset Hexagon Images to those matching critera of interest

# Subset hexagon pairs to those that have >33% overlap and where both members have been downloaded previously
pairs_subset = hex_pairs.loc[
    (hex_pairs['fArea']>0.33) & (hex_pairs['Download Count']==2)]
#####################################################################

# Find glaciers within hexagon pairs
pairs_subset = pyHEX.find_glaciers(pairs_subset, RGI_set)

################## USER EDITS MAY BE REQUIRED HERE ##################
## The remainder of the functions in this script and the later heximap 
# functions are designed to process a single hexagon pair at a time.
# Much of the code is structured in such a way that updating the functions/
# scripts to serially process multiple Hexagon pairs should be relatively
#  straightforward, but at present (July 2022) the necessary modifications 
# have not yet been implemented.
# The lines in this USER EDIT box therefore further subset the paired hexagon data to a single pair (in this case, the pair used for Morgan's tutorial walkthrough of the original heximap) 

# Subset paired data to single hex pair
tut_idx = 38
pair_i = pairs_subset.loc[[tut_idx]]

# Get glaciers intersecting selected hex pair
glaciers_i = RGI_set.loc[pair_i['glacier_idx'].iloc[0]]
#####################################################################


# Save glaciers that intersect with selected Hexagon pair
GlacierPath = TMP_DIR.joinpath('hex-glaciers')
GlacierPath.mkdir(parents=True, exist_ok=True)
# Ensure glaciers are in WGS84 coordinates (necessary for heximap)
rgi_hex = glaciers_i.copy()
rgi_hex.to_crs(epsg=4326, inplace=True)
# Export glaciers as shapefile
rgi_hex.to_file(GlacierPath.joinpath('rgi-hex.shp'))

# Get bounds of glaciated regions within hex pair
ROIs = pyHEX.pairROIs(glaciers_i)
# ****Subset ROIs to decrease heximap processing time during development***
ROIs = ROIs.iloc[[0,2,3,7]]

# Save ROI data for tutorial images
pyHEX.SaveROIs(ROIs, SavePath=TMP_DIR)

# Save image metadata for tutorial images
pyHEX.hex_asmat(
    hex_subset, pair_i.loc[tut_idx,'Indices'], 
    exp_path=TMP_DIR.joinpath('metadata'))

# Export MATLAB params file for later use
pyHEX.param_asmat(
    root=ROOT_DIR, source=SRC_DIR, image=IM_DIR, 
    IM1_name=pair_i.IDs.iloc[0][0], IM2_name=pair_i.IDs.iloc[0][1], 
    georef=GeoRefPath, GlacierShp=GlacierPath, OutDir=TMP_DIR)
