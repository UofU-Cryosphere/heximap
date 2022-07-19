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
OUT_DIR = Path(sys.argv[4])
# ROOT_DIR = Path('/home/durbank/scratch/heximap/')
# SRC_DIR = Path('/home/durbank/Documents/Research/Glaciers/heximap/src/')
# DATA_DIR = Path('/home/durbank/Documents/Research/Glaciers/GSLR/data/')
# IM_DIR = Path('/home/durbank/Documents/Research/Glaciers/GSLR/data/hexagon/declass-ii/imagery/')
# OUT_DIR = Path('home/durbank/Documents/Research/Glaciers/GSLR/scratch/outputs/')

# Define and (if needed) create tmp processing directory
TMP_DIR = ROOT_DIR.joinpath('tmp')
TMP_DIR.mkdir(parents=True, exist_ok=True)

print("Here's the Python output print statements...")
print(SRC_DIR)
print(DATA_DIR)
print(IM_DIR)
print(OUT_DIR)
print(TMP_DIR)

# Load custom modules
sys.path.append(SRC_DIR.as_posix())
import pyHEX

# Load hexagon metadata
hex2 = gpd.read_file(DATA_DIR.joinpath(
    "hexagon/declass-ii/metadata/declassii_620eb22e9f76579.shp"))

# Load glacier outlines
RGI_files = list(DATA_DIR.joinpath("RGI-data").glob("*/*.shp"))
gdfs = [gpd.read_file(path) for path in RGI_files]
RGI = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), crs=gdfs[0].crs)

# Subset hex data to region of Southern Andes
RGI_test = RGI[RGI["O1Region"] == "17"]
RGI_test = RGI_test.query("Area > 1")
bbox = RGI_test.total_bounds
bpoly = Polygon(
    ((bbox[0],bbox[1]), (bbox[2],bbox[1]), 
    (bbox[2],bbox[3]), (bbox[0],bbox[3])))
set_idx = hex2.intersects(bpoly)
hex2_subset = hex2.loc[set_idx]

# Generate hexagon pairs
hex2, hex2_pairs = pyHEX.pairify(hex2_subset)
pairs2 = hex2_pairs.loc[
    (hex2_pairs['fArea']>0.33) & (hex2_pairs['Download Count']==2)]

# Find glaciers within hexagon pairs
pairs2 = pyHEX.find_glaciers(pairs2, RGI_test)

# Subset paired data to single hex pair from tutorial
pair_idx = 38
pair_i = pairs2.loc[[pair_idx]]
pair_glaciers = RGI_test.loc[pair_i['glacier_idx'].iloc[0]]

# Get bounds of glaciated regions within hex pair (subsetted for development)
ROIs = pyHEX.pairROIs(pair_glaciers)
ROI_set = ROIs.iloc[[0,2,3,7]]

# Save ROI data for tutorial images
pyHEX.SaveROIs(ROI_set, SavePath=TMP_DIR)

# Save image metadata for tutorial images
pyHEX.hex_asmat(
    hex2, pairs2.loc[pair_idx,'Indices'], 
    exp_path=TMP_DIR.joinpath('metadata'))

# Export MATLAB param file for later use
pyHEX.param_asmat(
    root=ROOT_DIR, source=SRC_DIR, image=IM_DIR, 
    IM1_name=pair_i.IDs.iloc[0][0], IM2_name=pair_i.IDs.iloc[0][1], 
    georef=DATA_DIR, OutDir=TMP_DIR)
