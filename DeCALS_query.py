
from dl import authClient as ac, queryClient as qc
from dl.helpers.utils import convert
from getpass import getpass
import time
from astropy.io import fits
import numpy as np
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
import glob

slice_size=50000

cluster_centre=SkyCoord(str(354.41916666667 ), str(0.27666666666667), frame='icrs',unit=(u.deg,u.deg))

offset=0

sql_count_query=""" SELECT t.ra from ls_dr9.photo_z as pz, ls_dr9.tractor as 
t where Q3C_RADIAL_QUERY(ra,dec,355.9154,0.3308,2) AND t.ls_id = pz.ls_id"""

total_rows_str = qc.query(sql=sql_count_query, scalar=True)
result = total_rows_str.split('\n')[1:-1]
total_rows=len(result)

combined_result = []

# Iterate through the data in smaller slices
for offset in range(0, total_rows, slice_size):
    # Construct the SQL query with the current offset and limit (slice size)
    sql_query = f"""
    SELECT t.ra,t.dec,t.flux_w1,t.flux_w2,t.flux_w3,t.flux_w4,t.flux_ivar_w1, t.flux_ivar_w2, t.flux_ivar_w3,
     t.flux_ivar_w4, t.flux_g, t.flux_r,t.flux_z,  t.flux_ivar_g, t.flux_ivar_r,t.flux_ivar_z, t.mag_g, t.mag_r, t.mag_z, t.mag_w1, t.mag_w2,
     t.mag_w3,t.mag_w4,t.objid, t.type, t.snr_g, t.snr_r, t.snr_z, t.snr_w1, t.snr_w2, t.snr_w3, t.snr_w4, t.w1_w2, t.w2_w3, t.w3_w4, pz.ls_id, pz.objid, pz.z_phot_mean,
     pz.z_phot_median, pz.z_spec, pz.z_phot_l68, pz.z_phot_u68, pz.z_phot_l95, pz.z_phot_u95
          FROM ls_dr9.photo_z AS pz, ls_dr9.tractor AS t
    WHERE Q3C_RADIAL_QUERY(t.ra, t.dec, {cluster_centre.ra.value}, {cluster_centre.dec.value}, 2)
    AND t.ls_id = pz.ls_id
    LIMIT {slice_size} OFFSET {offset}
    """
    t1 = time.time()
    # Use qc.query() to execute the query and get the result
    result = qc.query(sql=sql_query)
    result1 = result.split('\n')[1:-1]
    t2 = time.time()
    print('Time',t2-t1, len(result1), slice_size, offset)
  
    # If the result is empty, it means we have processed all rows, so break the loop
    if not result:
        break
    else:
        keys = np.array(result.split('\n')[0].split(','))

        t2 = time.time()

        print('Query finished in ', t2-t1, 'seconds')
    # Add the current slice to the combined result list
        combined_result += result1

        data=[]

        for row_str in result1:
            values = row_str.split(',')
            data.append(dict(zip(keys, values)))

        table = Table(rows=data)

        output_filename = 'query_result_'+str(offset)+'.fits'

        table.write(output_filename, format='fits', overwrite=True)

        print('table '+ output_filename + ' has been written ')


