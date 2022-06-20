
# Module imports
import pandas as pd
import geopandas as gpd
import warnings
import numpy as np
import scipy.io


# Function to standardize Hexagon 3 data to Hexagon Declass 2 conventions
def standardize(hex_gdf, hex_class=2):

    # Get copy of gdf (this avoids issues with setting values on slices)
    hex_gdf = hex_gdf.copy()

    # Match naming conventions of Declass 3 to Declass 2
    if hex_class==3:
        # Rename columns to match hexagon declass 2 names
        hex_gdf.rename(columns={
            'Frame Numb':'Frame', 'Camera':'Camera Typ', 
            'Download A':'Down Load'}, inplace=True)
        
        # Relabel `Down Load` to match hexagon declass 2 convention
        hex_gdf['Down Load'].replace(
            to_replace=["Y", "N"], value=["Yes", "No"], inplace=True)

    # Reformat `Aquisitio` to datetime and `Down Load` to logical
    hex_gdf.loc[:,"Acquisitio"] = pd.to_datetime(hex_gdf["Acquisitio"])
    hex_gdf.loc[:,'Down Load'] = hex_gdf['Down Load'] == "Yes"

    # Get new groupings based on consecutive frame number sequences
    frame_diff = hex_gdf['Frame'].diff() - 1
    frame_diff.iloc[0] = 0.
    grp_vals = frame_diff.ne(0).cumsum()
    hex_gdf['Group'] = grp_vals
    # # Get new groupings based on ID
    # if hex_class==2:
    #     patt_drop = r'(.*-\d0*)'
    #     tmp = hex_gdf['Entity ID'].str.split(patt_drop)
    #     str_new = tmp.apply(lambda x: x[-1])
    #     hex_gdf['Group'] = str_new.str.extract(r'(^\d+)').astype('int')
    # elif hex_class==3:
    #     patt_drop = r'(.*-\d0*)'
    #     tmp = hex_gdf['Entity ID'].str.split(patt_drop)
    #     str_new = tmp.apply(lambda x: x[-1])
    #     hex_gdf['Group'] = str_new.str.extract(r'(^\d+)').astype('int')


    return hex_gdf



# Function to find pairs of hexagon imagery, given a df of hex data
def pairify(hex_gdf, hex_class=2):

    # Standardize data format
    hex = standardize(hex_gdf, hex_class=hex_class)


    # Group by mission and camera type and keep those with at least 2 members
    mission_grp = hex.groupby(["Mission", "Camera Typ", "Group"])
    grp_cnt = mission_grp.count()
    miss_keep = grp_cnt[grp_cnt.iloc[:, 0] >= 2].index

    # Iterate through mission/camera indices
    ls_pairs = []
    download_cnt = []
    for key in miss_keep:
        miss_idx = mission_grp.groups.get(key)
        miss_hex = hex.loc[miss_idx]

        # Match up hexagon images as pairs
        while miss_hex.shape[0] > 1:
            
            # Select first hexagon image
            hex_0 = miss_hex.iloc[0]

            # Determine which other images are within 1 frame of current
            f_val = abs(hex_0["Frame"] - miss_hex["Frame"])
            f_idx = f_val <= 1

            # Drop current image if no other images are within frame and
            #  continue to next image
            if sum(f_idx) < 2:
                miss_hex.drop(miss_hex.index[0], inplace=True)
                continue
            
            # Get indices of image pairs for current image
            hex_vals = list(miss_hex.loc[f_idx].index)
            ls_pairs.append(hex_vals)

            # Number of previous downloads in pair
            download_cnt.append(miss_hex.loc[hex_vals,'Down Load'].sum())

            miss_hex.drop(hex_vals[0], inplace=True)

    # Get Entity IDs for hex pairs
    ID_pairs = []
    for idx in ls_pairs:
        ID_pairs.append(hex.loc[idx, 'Entity ID'].to_list())

    # Group hex pairs into df
    pairs_df = pd.DataFrame(
        {'IDs':ID_pairs, 'Indices':ls_pairs, 'Download Count':download_cnt})

    # Get (approximate) region of overlap between hex pairs
    warnings.filterwarnings(
        action='ignore', 
        message='.*Geometry is in a geographic CRS.*')
    geoms = []
    geo_fArea = np.empty(len(ls_pairs))
    for i,idx in enumerate(ls_pairs):
        pair_i = hex.loc[idx]
        geom_i = pair_i.iloc[0:1].intersection(pair_i.iloc[1::], align=False)
        geoms.append(geom_i.geometry.iloc[0])
        geo_fArea[i] = (geom_i.area / pair_i.area.mean()).values
    warnings.filterwarnings(
        action='default', 
        message='.*Geometry is in a geographic CRS.*')
    pairs_df['fArea'] = geo_fArea
    pairs_gdf = gpd.GeoDataFrame(data=pairs_df, geometry=geoms, crs=hex.crs)

    return (hex, pairs_gdf)

# Function to save hexagon image pair metadata to matlab file
def save_asmat(hex_gdf, pair_idx, exp_path):
    hex = hex_gdf.loc[pair_idx]
    IM1 = hex.iloc[0]
    IM1.index = IM1.index.str.replace(' ', '_')
    IM2 = hex.iloc[1]
    IM2.index = IM2.index.str.replace(' ', '_')
    scipy.io.savemat(
        exp_path.joinpath((IM1['Entity_ID']+'_meta.mat')), IM1.to_dict())
    scipy.io.savemat(
        exp_path.joinpath((IM2['Entity_ID']+'_meta.mat')), IM2.to_dict())

# Function to find glaciers within overlap region of hexagon pairs
def find_glaciers(pair_gdf, glacier_gdf, buffer_sz=10000):

    # Reproject to Equal Earth crs (to ensure better area calculations)
    pair_gdf = pair_gdf.copy()
    pair_reproj = pair_gdf.to_crs(epsg=8857)
    glacier_reproj = glacier_gdf.copy().to_crs(epsg=8857)

    # Define ROI for each glacier (necessary because hex bounds are only 
    # accurate to within a few km)
    # Default ROI is 10 km radius around glacier center
    glacier_buff = glacier_reproj.centroid.buffer(buffer_sz)

    glacier_idx = []
    for idx, row in pair_reproj.iterrows():
        in_idx = glacier_buff.within(row['geometry'])
        glacier_idx.append(glacier_buff.index[in_idx].to_list())
    pair_gdf['glacier_idx'] = glacier_idx
    pair_gdf['Num_glaciers'] = pair_gdf['glacier_idx'].map(lambda x: len(x))

    return pair_gdf