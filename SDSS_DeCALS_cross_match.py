from astropy.io import fits
from astropy.table import Table, vstack, join
from astropy.coordinates import SkyCoord,match_coordinates_sky
from astropy import units as units
import pandas as pd
import numpy as np

def match_coord(ra1, dec1, ra2, dec2, search_radius=1., nthneighbor=1, verbose=True):
        '''
        Match objects in (ra2, dec2) to (ra1, dec1).         Inputs: 
                RA and Dec of two catalogs;
                search_radius: in arcsec;
                (Optional) keep_all_pairs: if true, then all matched pairs are kept; otherwise, if more than
                one object in t2 is match to the same object in t1 (i.e. double match), only the closest pair
                is kept.        Outputs: 
                idx1, idx2: indices of matched objects in the two catalogs;
                d2d: distances (in arcsec);
                d_ra, d_dec: the differences (in arcsec) in RA and Dec; note that d_ra is the actual angular 
                separation;
        '''
        t1 = Table()
        t2 = Table()
        # protect the global variables from being changed by np.sort
        ra1, dec1, ra2, dec2 = map(np.copy, [ra1, dec1, ra2, dec2])
        t1['ra'] = ra1
        t2['ra'] = ra2
        t1['dec'] = dec1
        t2['dec'] = dec2
        t1['id'] = np.arange(len(t1))
        t2['id'] = np.arange(len(t2))
        # Matching catalogs
        sky1 = SkyCoord(ra1*units.degree,dec1*units.degree, frame='icrs')
        sky2 = SkyCoord(ra2*units.degree,dec2*units.degree, frame='icrs')
        idx, d2d, d3d = sky2.match_to_catalog_sky(sky1, nthneighbor=nthneighbor)
        # This finds a match for each object in t2. Not all objects in t1 catalog are included in the result.         # convert distances to numpy array in arcsec
        d2d = np.array(d2d.to(units.arcsec))
        matchlist = d2d<search_radius
        if np.sum(matchlist)==0:
                if verbose:
                        print('0 matches')
                return np.array([], dtype=int), np.array([], dtype=int), np.array([]), np.array([]), np.array([])
        t2['idx'] = idx
        t2['d2d'] = d2d
        t2 = t2[matchlist]
        init_count = np.sum(matchlist)
        #--------------------------------removing doubly matched objects--------------------------------
        # if more than one object in t2 is matched to the same object in t1, keep only the closest match
        t2.sort('idx')
        i = 0
        while i<=len(t2)-2:
                if t2['idx'][i]>=0 and t2['idx'][i]==t2['idx'][i+1]:
                        end = i+1
                        while end+1<=len(t2)-1 and t2['idx'][i]==t2['idx'][end+1]:
                                end = end+1
                        findmin = np.argmin(t2['d2d'][i:end+1])
                        for j in range(i,end+1):
                                if j!=i+findmin:
                                        t2['idx'][j]=-99
                        i = end+1
                else:
                        i = i+1     
        mask_match = t2['idx']>=0
        t2 = t2[mask_match]
        t2.sort('id')
        if verbose:
                print('Doubly matched objects = %d'%(init_count-len(t2)))
        # -----------------------------------------------------------------------------------------
        if verbose:
                print('Final matched objects = %d'%len(t2))
        # This rearranges t1 to match t2 by index.
        t1 = t1[t2['idx']]
        
        ##########################################
        return np.array(t1['id']), np.array(t2['id']), np.array(t2['d2d'])


def DeCALS_cross_match(input_catalog,catalog_DeCALS,sep,ra,dec):

    catalog1 = Table.read(input_catalog)
    catalog2 = Table.read(catalog_DeCALS)

    table = catalog2

    
    column_names= table.colnames

    for column in column_names:

        current_dtype = table[column].dtype

        # Only change the data type if it's not an integer or float
        if current_dtype.kind not in ['i', 'f', 'U']:

            try:
                table[column] = np.array(table[column], dtype='float64')
            except ValueError:
                # If the data type conversion fails due to non-numeric values, skip this column
                pass
     


    # Define the maximum separation distance for cross-matching in arcseconds
    index=match_coord(ra1=catalog1[ra], dec1=catalog1[dec], ra2=catalog2['ra'], dec2=catalog2['dec'], search_radius=sep)
    

    matched_catalog1 = catalog1[index[0]]
    matched_catalog1['id'] = range(len(index[0]))


    matched_catalog2 = catalog2[index[0]]
    matched_catalog2['id'] = range(len(index[0]))

    joined_table = join(matched_catalog1, matched_catalog2, keys='id')


    joined_table.write('Zwcl2341_DeCALS_SDSS_join.fits', format='fits', overwrite=True)

  

 

if __name__ == "__main__":

    input_catalog= 'zwcl2341_90_DR17.fits'
    catalog_DeCALS= 'Zwcl2341_DeCALS.fits'
    sep=2.
    ra='ra'
    dec='dec'

    DeCALS_cross_match(input_catalog,catalog_DeCALS,sep,ra,dec)
  



