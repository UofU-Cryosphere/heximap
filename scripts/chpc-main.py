import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Polygon
import sys
# import matlab.engine

# Get bash inputs to define src, data, and output directories
SRC_DIR = Path(sys.argv[1])
DATA_DIR = Path(sys.argv[2])
OUT_DIR = Path(sys.argv[3])

# Define and (if needed) create scratch processing directory
SCRATCH_DIR = SRC_DIR.joinpath('scratch')
SCRATCH_DIR.mkdir(parents=True, exist_ok=True)

# Load custom modules
sys.path.append(SRC_DIR.joinpath('src').as_posix())
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
pyHEX.SaveROIs(ROI_set, SavePath=SCRATCH_DIR.joinpath('hexROIs.mat'))

# Save image metadata for tutorial images
pyHEX.save_asmat(
    hex2, pairs2.loc[pair_idx,'Indices'], exp_path=SCRATCH_DIR)




# # Call matlab script of heximap processes
# eng = matlab.engine.start_matlab()
# eng.triarea(nargout=0)