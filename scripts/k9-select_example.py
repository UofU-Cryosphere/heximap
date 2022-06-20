# Script to determine suitable KH-9 Hexagon stereo pairs for glacier DEM reconstruction

# %%
import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Polygon

# Define root and data directories
ROOT_DIR = Path().absolute().parents[0]
DATA_DIR = Path('/home/durbank/Documents/Research/Glaciers/GSLR/data')
print(DATA_DIR)

# Load custom module
import sys
sys.path.append(ROOT_DIR.joinpath('src').as_posix())
import pyHEX

# %% Data import

# hex metadata filenames
hex2 = gpd.read_file(DATA_DIR.joinpath(
    "declass-ii/declassii_620eb22e9f76579.shp"))
hex3 = gpd.read_file(DATA_DIR.joinpath(
    "declass-iii/declassiii_620ec4d67f3bf101.shp"))

# Load glacier outlines
RGI_files = list(DATA_DIR.joinpath("RGI-data").glob("*/*.shp"))
gdfs = [gpd.read_file(path) for path in RGI_files]
RGI = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), crs=gdfs[0].crs)

# %% tmp work for testing with tutorial data

# # Subset glaciers to regions contained with hex_subset
# bbox = hex_subset.total_bounds
# RGI_test = RGI[RGI["O1Region"] == "17"]
# RGI_test = RGI_test.loc[
#     (RGI_test['CenLon']>bbox[0]) & (RGI_test['CenLon']<bbox[2])]
# RGI_test = RGI_test.loc[
#     (RGI_test['CenLat']>bbox[1]) & (RGI_test['CenLat']<bbox[3])]
# RGI_test = RGI_test.query("Area > 1")


# Subset hex data to region of Southern Andes
RGI_test = RGI[RGI["O1Region"] == "17"]
RGI_test = RGI_test.query("Area > 1")
bbox = RGI_test.total_bounds
bpoly = Polygon(
    ((bbox[0],bbox[1]), (bbox[2],bbox[1]), 
    (bbox[2],bbox[3]), (bbox[0],bbox[3])))

# set_idx = (hex3['Center L_1'].astype(float)>RGI_test.total_bounds[1]) \
#     & (hex3['Center L_1'].astype(float)<RGI_test.total_bounds[2]) \
#     & (hex3['Center L_2'].astype(float)>RGI_test.total_bounds[0]) \
#     & (hex3['Center L_2'].astype(float)<RGI_test.total_bounds[2])
set_idx = hex2.intersects(bpoly)
hex2_subset = hex2.loc[set_idx]

set_idx = hex3.intersects(bpoly)
hex3_subset = hex3.loc[set_idx]

# %%
hex2, hex2_pairs = pyHEX.pairify(hex2_subset)
pairs2 = hex2_pairs.loc[
    (hex2_pairs['fArea']>0.1) & (hex2_pairs['Download Count']==2)]

# hex3, hex3_pairs = pyHEX.pairify(hex3_subset, hex_class=3)
# pairs3 = hex3_pairs.loc[
#     (hex3_pairs['fArea']>0.1) & (hex3_pairs['Download Count']==2)]

pairs2 = pyHEX.find_glaciers(pairs2, RGI_test)

# %% Save metadata for tutorial images

pair_idx = 38

pyHEX.save_asmat(
    hex2, pairs2.loc[pair_idx,'Indices'], 
    exp_path=ROOT_DIR.joinpath('data/tmp'))

# Get total bounds of glaciated regions within hex pair
glacier_poly = RGI_test.loc[
    pairs2.loc[pair_idx, 'glacier_idx']].unary_union.convex_hull

# %%[markdown]
# 
# 
# 
# 
# 
#  
# ## Older deprecated stuff that might still be useful
# 
# 
# 
# 
# 
# 
#  
# %% 

# # Subset data for development
# RGI_test = RGI[RGI["O1Region"] == "15"].query("Area > 1")
# RGI_test = RGI_test.iloc[-1000:].query("CenLon > 95")
# RGI_bnds = RGI_test.unary_union.convex_hull
# hex = gpd.GeoDataFrame(hex3[hex3.intersects(RGI_bnds)].convert_dtypes())

# # Reformat some hex data
# hex["Acquisitio"] = pd.to_datetime(hex["Acquisitio"])
# hex["Download A"] = hex["Download A"] == "Y"

# import time
# t_0 = time.perf_counter()
# glacier_dict = {}
# download_count = {}

# for idx in RGI_test.index:

#     glacier = RGI_test.loc[idx]
#     hex_i = hex[hex.contains(glacier.geometry)]
#     mission_grp = hex_i.groupby(["Mission", "Camera"])

#     grp_cnt = mission_grp.count()
#     miss_keep = grp_cnt[grp_cnt.iloc[:, 0] >= 2].index

#     for key in miss_keep:
#         miss_idx = mission_grp.groups.get(key)
#         miss_hex = hex_i.loc[miss_idx]

#         while miss_hex.shape[0] > 1:
        
#             hex_0 = miss_hex.iloc[0]
#             f_idx = abs(hex_0["Frame Numb"] - miss_hex["Frame Numb"]) <= 1
#             if sum(f_idx) < 2:
#                 miss_hex.drop(miss_hex.index[0], inplace=True)
#                 continue

#             hex_vals = list(miss_hex.loc[f_idx].index)

#             glacier_key = glacier["RGIId"]
#             if glacier_key in glacier_dict:
#                 glacier_dict[glacier_key].append(tuple(hex_vals))
#                 download_count[glacier_key].append(
#                     miss_hex.loc[hex_vals,'Download A'].sum())
#             else:
#                 glacier_dict[glacier_key] = [(), tuple(hex_vals)]
#                 glacier_dict[glacier_key] = glacier_dict[glacier_key][1:]

#                 download_count[glacier_key] = [
#                     (), miss_hex.loc[hex_vals,'Download A'].sum()]
#                 download_count[glacier_key] = download_count[glacier_key][1:]

#             miss_hex.drop(hex_vals[0], inplace=True)

# t_end = time.perf_counter()
# print(f"Run time: {t_end-t_0} s")

# %%

# # df with unique combinations of glacier and hex pair (repeats of both glaciers and hex pairs are possible)
# pair_df = pd.DataFrame(
#     {'RGIId':glacier_dict.keys(), 'hex_pairs':glacier_dict.values(), 
#     'Downloads':download_count.values()}).explode(['hex_pairs', 'Downloads'])

# # df with unqiue glacier entries giving all (unpaired) hex images with coverage
# AllHex_df = pair_df.groupby('RGIId').sum()


# # df showing number of glaciers covered by each unique hex pair
# GlacierCnt_df = pair_df.groupby('hex_pairs').count().sort_values(
#     'RGIId', ascending=False)
# GlacierCnt_df['RGI_cnt'] = GlacierCnt_df['RGIId']
# GlacierCnt_df.drop(columns=['RGIId', 'Downloads'], inplace=True)

# # For each hex pair, determine which hex images are already downloaded
# idx1, idx2 = map(list,zip(*[el for el in GlacierCnt_df.index]))
# GlacierCnt_df['IM1_purchased'] = hex['Download A'].loc[idx1].values
# GlacierCnt_df['IM2_purchased'] = hex['Download A'].loc[idx2].values


# # High priority glaciers (based on at least 1 member of pair already downloaded)
# priority_idx = AllHex_df.query('Downloads > 0').index
# HighPair_df = pair_df[pair_df['RGIId'].isin(priority_idx)].query(
#     'Downloads > 0')
