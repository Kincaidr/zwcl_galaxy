#!/usr/bin/env python
from astropy import coordinates as coords
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
from astropy.nddata.utils import Cutout2D
import numpy as np
import requests
from io import BytesIO
from PIL import Image
from astropy.wcs import WCS
from astropy.stats import biweight_scale
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.pyplot import show, plot, draw
from multiprocessing import Process
from scipy.constants import pi
import scipy.constants as constants
import statistics as stats
import astropy.cosmology as cp
import sys
from astropy.io import fits
from matplotlib.patches import Circle
import pandas as pd
import csv
from scipy.stats import norm
from scipy.optimize import curve_fit
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck13, z_at_value
import matplotlib.cm as cm
import os
from getpass import getpass
import time
import numpy as np
from astropy.table import Table,vstack
import glob
from dl import authClient as ac, queryClient as qc
from dl.helpers.utils import convert
from getpass import getpass
from astropy.table import Table, Column, join
from astromatch import Catalogue
from astromatch import Match
from mocpy import MOC

def LR(folder_path,radio_cat,optical_cat,cluster_centre,cluster_name):

    
    def fix_col_radio(radio_cat):
    
        radio = fits.open(radio_cat)
        data=radio[1].data
        

        data_e_ra=data['E_RA']*3600
        data_e_dec=data['E_DEC']*3600

        data['E_RA']=data_e_ra
        data['E_DEC']=data_e_dec

        fits.writeto(radio_cat+'_fix.fits',data=data,overwrite=True)

        return(radio_cat+'_fix.fits')
    

    def fix_col_optical_HSC(optical_cat):

        optical = fits.open(optical_cat)
        data=optical[1].data
        t=Table(data)

        for column in t.keys():
            if  'isnull' in column:
                t.remove_columns(column)

        raerr = np.maximum.reduce([data['g_sdsscentroid_raerr'] , data['r_sdsscentroid_raerr'], data['i_sdsscentroid_raerr'], data['z_sdsscentroid_raerr'], data['y_sdsscentroid_raerr']])*1e-6
        decerr = np.maximum.reduce([data['g_sdsscentroid_decerr'], data['r_sdsscentroid_decerr'], data['i_sdsscentroid_decerr'], data['z_sdsscentroid_decerr'], data['y_sdsscentroid_decerr']])*1e-6

        ra_err=Column(name='ra_err',data= raerr)
        dec_err=Column(name='dec_err',data= decerr)

        t.add_column(ra_err)
        t.add_column(dec_err)
        # Replace NaN values with 1e-06
        for col_name in t.colnames:
            col_data = t[col_name]
            col_data = np.nan_to_num(col_data, nan=1e-06)
            t[col_name] = col_data


        t.write(optical_cat+'_fix.fits',overwrite=True)

        return(optical_cat+'_fix.fits')

    def fix_col_optical_DeCALS(optical_cat):

        def flux2mag(flux):
                return (22.5-2.5*np.log10(flux))

        t=Table.read(optical_cat)
        #import IPython;IPython.embed()

        mag_g=Column(name='mag_g', data=flux2mag(t['FLUX_G']))
        mag_r=Column(name='mag_r', data=flux2mag(t['FLUX_R']))
        mag_i=Column(name='mag_i', data=flux2mag(t['FLUX_I']))


        # t.add_column( mag_g)
        # t.add_column( mag_r)
        # t.add_column( mag_i)
        

        # mag_g=flux2mag(t['FLUX_G'])
        # mag_r=flux2mag(t['FLUX_R'])
        # mag_i=flux2mag(t['FLUX_I'])

        return(t)


    radio=radio_cat
    optical=optical_cat


    # try:
    #      optical=fix_col_optical_HSC(optical_cat)
    # except ValueError:
    #      optical=optical_cat

    #import IPython;IPython.embed()

    sdss_mags = [ 'mag_g','mag_r','mag_i']#,'mag_w1','mag_w2','mag_w3']
    


    search_radius_a = 2*u.deg
    search_radius_b = 1.5*u.deg

    moc_xxl = MOC.from_elliptical_cone(
        lon=cluster_centre.ra,
        lat=cluster_centre.dec,
        a=search_radius_a,
        b=search_radius_b,
        pa=0*u.deg,
        max_depth=14
    )



    rcat = Catalogue(radio,name='radio_cat',id_col='Source_id',coord_cols=['RA','DEC'],poserr_cols=None ,poserr_type='circle',area=moc_xxl) #['E_RA','E_DEC']

    ocat = Catalogue(optical,name='optical_cat',id_col='OBJID',coord_cols=['RA','DEC'],poserr_cols=None,poserr_type='circle',area=moc_xxl,mag_cols=sdss_mags) #['ra_err','dec_err']

    rcat.poserr.add_syserr(0.4*u.arcsec)
    ocat.poserr.add_syserr(0.4*u.arcsec)

    xm = Match(rcat, ocat)

    match_results_lr = xm.run(method='lr', radius=2.0*u.arcsec)


    # match_results=xm.run(method='nway', radius=3.0*u.arcsec, use_mags=True)
    # match_results=xm.run(method='lr', radius=3.0*u.arcsec, use_mags=True)
 
    #cutoff=0.5
    xm.set_best_matchs(cutoff=0.5)

    best_matches = xm.get_matchs(match_type='best') 

    print('best_matches',len(best_matches))

    #LR=lr_stats['cutoff'][np.nanargmax(lr_stats['CR'])]


    lr_stats = xm.stats()

    dist=best_matches['Separation_radio_cat_optical_cat']

    n, bins = np.histogram(dist, bins=30)

    bin_centers = (bins[1:] + bins[:-1]) / 2

    n,bins,_=plt.hist(dist, bins=bin_centers,alpha=0.7)

    x=np.linspace(0,2,1000)
    def gauss(x, mu, sigma, A):
     return A*np.exp(-(x-mu)**2/2/sigma**2)
    
    expected = ( .2, .7, 1)

    # params, cov = curve_fit(gauss, bins, n, expected)

    # plt.plot(x, gauss(dist, *params), color='red', lw=3, label='model')


    #import IPython;IPython.embed()





    # optical_id_best=np.where(best_matches['SRCID_optical_cat'].value)[0]

    # radio_id_best=np.where(best_matches['SRCID_radio_cat'].value)[0]

    # optical = fits.open(optical_cat)
    # data_optical=optical[1].data

    # t=Table(data_optical)[optical_id_best]

    # t.write(optical_cat.replace('_new.fits','_LR.fits'),overwrite=True)

    # radio = fits.open(radio_cat)
    # data_radio=optical[1].data

    # t=Table(data_radio)[radio_id_best]

    # t.write(radio_cat.replace('_radio_cat.fits','_radio_LR.fits'),overwrite=True)

    
    

    # dra, ddec = xm.offset('optical_cat','radio_cat', match_type='best')
    # x = np.linspace(-10, 10, num=100)

    # fig = plt.figure(figsize=(13, 5))

    # plt.subplot(121)
    # plt.hist(dra.value, bins='auto', density=True)
    # plt.plot(x, norm.pdf(x, np.mean(dra), np.std(dra)))
    # plt.xlabel('DRA / arcsec')

    # plt.subplot(122)
    # plt.hist(ddec.value, bins='auto', density=True)
    # plt.plot(x, norm.pdf(x, np.mean(ddec), np.std(ddec)))
    # plt.xlabel('DDEC / arcsec')


    fig = plt.figure(figsize=(13, 7))


    plt.plot(lr_stats['cutoff'], (1 - lr_stats['error_rate']) * 100, '--r', lw=3)

    plt.xlabel('Likelihood ratio (LR)',fontsize=20)
    plt.ylabel('Purity (%)',color='red',fontsize=20)
    plt.yticks(fontsize=12)
   

    ax2 = plt.twinx()
    # Using 'r' for a red line, you can change the color as needed
    plt.plot(lr_stats['cutoff'], lr_stats['completeness']*100, lw=3)
    ax2.set_ylabel('completeness (%)',color='blue',fontsize=20)
    plt.axvline(x=0.5, color='black', linestyle='--')
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.xlim(0,1)
    #plt.title(cluster_name,fontsize=17)

    plt.savefig(folder_path+cluster_name + '_LR_ratio.pdf')
    plt.show()



def DeCALS_query(folder_path,cluster_name,cluster_centre):

    
    slice_size=50000

    offset=0

    sql_count_query=f""" SELECT t.ra from ls_dr9.photo_z as pz, ls_dr9.tractor as 
    t where Q3C_RADIAL_QUERY(ra,dec,{cluster_centre.ra.value},{cluster_centre.dec.value},2) AND t.ls_id = pz.ls_id"""

    total_rows_str = qc.query(sql=sql_count_query, scalar=True)
    result = total_rows_str.split('\n')[1:-1]
    total_rows=len(result)
    print('total rows',total_rows)


    combined_result = []

    # Iterate through the data in smaller slices
    for offset in range(0, total_rows, slice_size):
        # Construct the SQL query with the current offset and limit (slice size)
        sql_query = f"""
        SELECT t.ra,t.dec,t.ra_ivar,t.dec_ivar,t.flux_w1,t.flux_w2,t.flux_w3,t.flux_w4,t.flux_ivar_w1, t.flux_ivar_w2, t.flux_ivar_w3,
        t.flux_ivar_w4, t.flux_g, p.apflux_g_2, t.flux_r, p.apflux_r_2, t.flux_z, p.apflux_z_2,  t.flux_ivar_g, t.flux_ivar_r,t.flux_ivar_z, t.mag_g, t.mag_r, t.mag_z, t.mag_w1, t.mag_w2,
        t.mag_w3,t.mag_w4,t.objid, t.type, t.snr_g, t.snr_r, t.snr_z, t.snr_w1, t.snr_w2, t.snr_w3, t.snr_w4, t.w1_w2, t.w2_w3, t.w3_w4, pz.ls_id, pz.objid, pz.z_phot_mean,
        pz.z_phot_median, pz.z_spec, pz.z_phot_l68, pz.z_phot_u68, pz.z_phot_l95, pz.z_phot_u95,pz.survey
            FROM ls_dr9.photo_z AS pz, ls_dr9.tractor AS t, ls_dr9.apflux AS p
        WHERE Q3C_RADIAL_QUERY(t.ra, t.dec, {cluster_centre.ra.value}, {cluster_centre.dec.value}, 2)
        AND t.ls_id = pz.ls_id AND pz.ls_id=p.ls_id
        LIMIT {slice_size} OFFSET {offset}
        """
        t1 = time.time()
        # Use qc.query() to execute the query and get the result
        result = qc.query(sql=sql_query)
        result1 = result.split('\n')[1:-1]
        t2 = time.time()
        print('Time',t2-t1, len(result1), slice_size, offset)
    
        # If the result is empty, it means we have processed all rows, so break the loop
        if not result1:
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

def stack_table(folder_path,cluster_name ):

    all_query_tables=glob.glob('query_result_*.fits')
    all_query_tables.sort()
    #table_list=[]
    stack = None
    for i, query_table in enumerate(all_query_tables):
        print(i, len(all_query_tables), query_table)
        table = Table.read(query_table, format='fits')
        if stack:
            stack = vstack([stack, table])
        else:
            stack = table
        stack.write(folder_path+cluster_name+'_stacked_table.fits', format='fits', overwrite=True)

    catalog = Table.read(folder_path+cluster_name+'_stacked_table.fits')

    for column in catalog.keys():
    # Check if the column dtype starts with 'U' (Unicode string)
        if catalog[column].dtype.kind in ['U', 'S']:
            try:
            # Convert the column to 'float64'
                catalog[column] = catalog[column].astype('float64')
            except ValueError:
                continue

        catalog.write(cluster_name+'_DeCALS.fits',format='fits',overwrite=True)



def HSC_query(folder_path,cluster_name,DeCALS_catalogue_fits):

    HSC=cluster_name
#   SELECT
#          f.object_id
#         , f.ra
#         , f.dec
#         , f.g_cmodel_flux
#         , f.g_cmodel_fluxerr
#         , f.r_cmodel_flux
#         , f.r_cmodel_fluxerr
#         , f.i_cmodel_flux
#         , f.i_cmodel_fluxerr
#         , f.z_cmodel_flux
#         , f.z_cmodel_fluxerr
#         , f.y_cmodel_flux
#         , f.y_cmodel_fluxerr
#         , f.g_cmodel_mag
#         , f.g_cmodel_magerr
#         , f.r_cmodel_mag
#         , f.r_cmodel_magerr
#         , f.i_cmodel_mag
#         , f.i_cmodel_magerr
#         , f.z_cmodel_mag
#         , f.z_cmodel_magerr
#         , f.y_cmodel_mag
#         , f.y_cmodel_magerr
#         , w.object_id
#         , w.photoz_median
#         , w.photoz_err68_min
#         , w.photoz_err68_max
#         , w.photoz_mean
#         , w.sfr
#         , w.sfr_err68_min
#         , w.sfr_err68_max
#         , w.stellar_mass
#         , w.stellar_mass_err68_min
#         , w.stellar_mass_err68_max
#         , p.g_apertureflux_15_mag
#         , p.g_apertureflux_15_magerr
#         , p.r_apertureflux_15_mag
#         , p.r_apertureflux_15_magerr
#         , p.i_apertureflux_15_mag
#         , p.i_apertureflux_15_magerr
#         , p.object_id
#     FROM
#             pdr3_wide.forced as f
#         JOIN
#             pdr3_wide.photoz_mizuki as w
#         ON
#             f.object_id = w.object_id
#         JOIN
#             pdr3_wide.forced3 as p
#         ON
#             f.object_id = p.object_id
#         WHERE conesearch(coord, 354.4191, 0.2766, 7200.0)



#   SELECT
#          f.object_id
#         , f.ra
#         , f.dec
#         , f.g_cmodel_flux
#         , f.g_cmodel_fluxerr
#         , f.r_cmodel_flux
#         , f.r_cmodel_fluxerr
#         , f.i_cmodel_flux
#         , f.i_cmodel_fluxerr
#         , f.z_cmodel_flux
#         , f.z_cmodel_fluxerr
#         , f.y_cmodel_flux
#         , f.y_cmodel_fluxerr
#         , f.g_cmodel_mag
#         , f.g_cmodel_magerr
#         , f.r_cmodel_mag
#         , f.r_cmodel_magerr
#         , f.i_cmodel_mag
#         , f.i_cmodel_magerr
#         , f.z_cmodel_mag
#         , f.z_cmodel_magerr
#         , f.y_cmodel_mag
#         , f.y_cmodel_magerr
#         , w.object_id
#         , w.photoz_median
#         , w.photoz_err68_min
#         , w.photoz_err68_max
#         , w.photoz_mean
#         , w.sfr
#         , w.sfr_err68_min
#         , w.sfr_err68_max
#         , w.stellar_mass
#         , w.stellar_mass_err68_min
#         , w.stellar_mass_err68_max
#     FROM
#             pdr3_wide.forced as f
#         JOIN
#             pdr3_wide.photoz_mizuki as w
#         ON
#             f.object_id = w.object_id
#         WHERE conesearch(coord, 354.4191, 0.2766, 7200.0)
        


def RA_DEC_seperation_optical(folder_path,cluster_name,DeCALS_catalogue_fits):


    
    DEC_sep=(galaxy_dec_SDSS-galaxy_dec_DeCALS)*3600

    RA_sep=(galaxy_ra_SDSS-galaxy_ra_DeCALS)*3600

    n, bins = np.histogram(RA_sep, bins=50)

    bin_DEC_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2

    n, bins = np.histogram(DEC_sep, bins=50)

    bin_RA_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2
    
    bin_RA_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2


    x = np.linspace(-3, 3, num=100)

    fig = plt.figure(figsize=(20, 7))
    plt.subplot(121)
    plt.hist(RA_sep-bin_RA_centre, bins='auto', density=True)
    plt.plot(x, norm.pdf(x, np.mean(RA_sep-bin_RA_centre), np.std(RA_sep-bin_RA_centre)))
    plt.xlabel('RA sep (arcsec)',fontsize=20)
    plt.title("RA seperation distance distribution \n between SDSS and DeCALS for " + cluster_name,fontsize=22)
    plt.subplot(122)
    
    plt.hist(DEC_sep-bin_DEC_centre, bins='auto', density=True)
    plt.plot(x, norm.pdf(x, np.mean(DEC_sep-bin_DEC_centre), np.std(DEC_sep-bin_DEC_centre)))
    plt.xlabel('DEC sep (arcsec)',fontsize=20)
    plt.title("DEC seperation distance distribution \n between SDSS and DeCALS for " + cluster_name,fontsize=22)
    plt.savefig(folder_path+cluster_name+'_RA_DEC_seperation_optical.pdf')
    plt.show()


def general_seperation(folder_path,cluster_name,optical_cat_fits,radio_cat_fits):
    
    optical_catalogue = fits.open(optical_cat_fits)
    optical_cat=optical_catalogue[1].data

    radio_catalogue = fits.open(radio_cat_fits)
    radio_cat=radio_catalogue[1].data

    ra=optical_cat['ra']
    dec=optical_cat['dec']

    ra1=radio_cat['RA']
    dec1=radio_cat['DEC']

    coo_optical = SkyCoord(ra*u.deg, dec*u.deg)
    coo_radio = SkyCoord(ra1*u.deg, dec1*u.deg)
    idx_sdss, d2d_sdss, d3d_sdss = coo_radio.match_to_catalog_sky(coo_optical)

    ras_sim = np.random.rand(len(coo_radio))*coo_optical.ra.ptp() + coo_optical.ra.min()
    decs_sim = np.random.rand(len(coo_radio))*coo_optical.dec.ptp() + coo_optical.dec.min()
    ras_sim, decs_sim

    coo_simulated = SkyCoord(ras_sim, decs_sim)  
    idx_sim, d2d_sim, d3d_sim = coo_simulated.match_to_catalog_sky(coo_optical)

    plt.hist(d2d_sim.arcsec, histtype='step', label='Simulated', linestyle='dashed',range=(0,10),bins=30)
    plt.hist(d2d_sdss.arcsec,  histtype='step', label='MeerKAT to HSC',range=(0,10),bins=30)
    plt.xlabel('separation [arcsec]',fontsize=15)
    plt.ylabel('Counts',fontsize=15)
    plt.title(cluster_name,fontsize=17)
    plt.legend(loc=0)
    plt.tight_layout()
    plt.savefig(folder_path+cluster_name+'_general_seperation.pdf')
    plt.show()




def RA_DEC_seperation_radio(folder_path,cluster_name,radio_DeCALS_catalogue):
    
    radio_DeCALS_catalogue_fits = fits.open(radio_DeCALS_catalogue)
    radio_DeCALS_cat=radio_DeCALS_catalogue_fits[1].data

    import IPython;IPython.embed()

    galaxy_ra_MeerKAT = radio_DeCALS_cat['RA_1']*3600
    galaxy_dec_MeerKAT = radio_DeCALS_cat['DEC_1']*3600
    galaxy_ra_DeCALS_radio =  radio_DeCALS_cat['ra_2']*3600
    galaxy_dec_DeCALS_radio =  radio_DeCALS_cat['dec_2']*3600

    DEC_sep=(galaxy_dec_MeerKAT-galaxy_dec_DeCALS_radio)

    RA_sep=(galaxy_ra_MeerKAT-galaxy_ra_DeCALS_radio)

    n, bins = np.histogram(RA_sep, bins=50)

    bin_DEC_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2

    n, bins = np.histogram(DEC_sep, bins=50)

    bin_RA_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2
    

    mean_RA = np.mean(RA_sep)
    std_RA = np.std(RA_sep)


    x_min_RA =np.min(RA_sep) - 3 * std_RA
    x_max_RA = np.max(RA_sep) + 3 * std_RA


    x_RA = np.linspace(x_min_RA, x_max_RA, 1000)


    mean_DEC = np.mean(DEC_sep)
    std_DEC = np.std(DEC_sep)

    x_min_DEC =np.min(DEC_sep) - 3 * std_DEC
    x_max_DEC = np.max(DEC_sep) + 3 * std_DEC

    x_DEC = np.linspace(x_min_DEC, x_max_DEC, 1000)
    

    fig = plt.figure(figsize=(20, 7))
    plt.subplot(121)
    plt.hist(RA_sep-bin_RA_centre, bins='auto', density=True)
    plt.plot(x_RA, norm.pdf(x_RA, np.mean(RA_sep-bin_RA_centre), np.std(RA_sep-bin_RA_centre)))
    plt.xlabel('RA sep (arcsec)',fontsize=20)
    plt.title("RA seperation distance distribution \n between DeCALS and MeerKAT for " + cluster_name,fontsize=22)

    plt.subplot(122)
    plt.hist(DEC_sep-bin_DEC_centre, bins='auto', density=True)
    plt.plot(x_DEC, norm.pdf(x_DEC, np.mean(DEC_sep-bin_DEC_centre), np.std(DEC_sep-bin_DEC_centre)))
    plt.xlabel('DEC sep (arcsec)',fontsize=20)
    plt.title("DEC seperation distance distribution \n between DeCALS and MeerKAT for " + cluster_name,fontsize=22)
    plt.savefig(folder_path+cluster_name+'_RA_DEC_seperation_radio.pdf')
    plt.show()


def radio_redshift_DeCLAS(folder_path,cluster_name,radio_DeCALS_catalogue):

    radio_DeCALS_catalogue_fits = fits.open(radio_DeCALS_catalogue)
    radio_DeCALS_cat=radio_DeCALS_catalogue_fits[1].data
    
    z=radio_DeCALS_cat['z_spec'] > 0

    z_spec=radio_DeCALS_cat['z_spec'][z]

    print('Length is ',len(z_spec))

    plt.title(cluster_name)
    plt.axvline(np.median(z_spec))
    n, bins = np.histogram(z_spec, bins=50)
    plt.hist(z_spec, bins=30, color='purple', alpha=0.5, label='redshift')
    plt.text(np.max(z_spec),np.argmax(n)+7,f'Mean z is {np.median(z_spec):.4f}', fontsize=15, ha='right', va='top')
    plt.text(np.max(z_spec),np.argmax(n)+12,f'Total galaxies is {len(z_spec):.0f}', fontsize=15, ha='right', va='top')
    plt.show()





def circle(folder_path,cluster_name,galaxy_ra_DeCALS,galaxy_dec_DeCALS,zsp,cluster_centre,zsp_min,zsp_max,search_radius):

    
    galaxy_ra=[]
    galaxy_dec_circle=[]
    redshift=[]

    #cluster 90arcmin radius


    for i in range(len(galaxy_ra_DeCALS)):
        if np.sqrt((galaxy_ra_DeCALS[i]-cluster_centre.ra.value)**2+(galaxy_dec_DeCALS[i]-cluster_centre.dec.value)**2)<= search_radius:
            galaxy_ra_circle.append(galaxy_ra_DeCALS[i])
            galaxy_dec_circle.append(galaxy_dec_DeCALS[i])
            redshift.append(zsp[i])
    redshift_circle=np.array(redshift)

    
    galaxy_ra_new=np.array(galaxy_ra_circle)[(np.array(redshift_circle) >= zsp_min) & (np.array(redshift_circle) <= zsp_max)]

    galaxy_dec_new =np.array(galaxy_dec_circle)[(np.array(redshift_circle) >= zsp_min) & (np.array(redshift_circle) <= zsp_max)]

    
    plt.scatter(galaxy_ra_new,galaxy_dec_new)

    plt.xlabel("RA in deg")
    plt.ylabel("DEC in deg")
    plt.title("1.5 deg SDSS DR17 galaxies of " + cluster_name)
    plt.savefig(folder_path + cluster_name + '_galaxy_distribution.pdf')
    

    return(galaxy_ra_circle,galaxy_dec_circle,redshift_circle)


def redshift_distribution(folder_path,cluster_name,galaxy_ra_circle,galaxy_dec_circle,redshift_circle,cluster_centre,zsp_min,zsp_max):



    zsp_min_1 = 0.24
    zsp_max_1 = 0.30

    redshift_1=np.array(redshift_circle)[(np.array(redshift_circle) >= zsp_min) & (np.array(redshift_circle)<= zsp_max)]
    redshift_2=np.array(redshift_circle)[(np.array(redshift_circle) >= zsp_min_1) & (np.array(redshift_circle)<= zsp_max_1)]


    fig = plt.figure(figsize=(12, 5))
    
    ax1 = fig.add_axes([0.1,0.3,0.8,0.6])
    ax2 = fig.add_axes([0.65, 0.75,0.25,0.15])

 
    ax1.hist(redshift_1,bins=100)
    ax2.hist(redshift_2,bins=50)
   
    ax1.set_xlabel('Redshift',fontsize=15)
    ax1.set_ylabel('Number of galaxies',fontsize=15)
    ax1.set_title("SDSS spectroscopic redshift distribution of galaxies within 1.5 deg of " + cluster_name)
    ax1.set_xlim(0.24,0.30)
    ax2.set_xlim(0,1)

    plt.savefig(folder_path + cluster_name + '_redshift_distribution.pdf')
    plt.show()



def velocity_dispersion(cluster_centre,R_200_Mpc,zsp_min,zsp_max,cluster_z_lit,DeCALS_catalog):
    
    DeCALS_cat=fits.open(DeCALS_catalog)[1].data

    galaxy_ra_SDSS=DeCALS_cat['ra']
    galaxy_dec_SDSS=DeCALS_cat['dec']
    #import IPython;IPython.embed()

    zsp=DeCALS_cat['z_spec']

    def R_200_to_deg(distance_Mpc,redshift):

        cosmo = FlatLambdaCDM(H0=70, Om0=0.27)

        def radians_to_degrees(radians):
            degrees = radians * (180 / np.pi)
            return degrees

        d_A = cosmo.angular_diameter_distance(z=redshift)

        theta_radian = distance_Mpc/d_A

        R_200=radians_to_degrees(theta_radian.value)

        return(R_200)

    redshift=np.array(zsp)
    
    R_200=R_200_to_deg(R_200_Mpc,cluster_z_lit)
    
    def histogram(sample,bins,color):
        n, bins, patches = plt.hist(sample, bins=bins, color=color, alpha=0.1, rwidth=0.85)
        cluster_z=bins[np.argmax(n)]
        return cluster_z


    
    ra_5Mpc_circle=[]
    dec_5Mpc_circle=[]
    zsp_R_200=[]
    for i in range(len(galaxy_ra_SDSS)):
        if np.sqrt((galaxy_ra_SDSS[i]-cluster_centre.ra.value)**2+(galaxy_dec_SDSS[i]-cluster_centre.dec.value)**2)<= R_200:
            ra_5Mpc_circle.append(galaxy_ra_SDSS[i])
            dec_5Mpc_circle.append(galaxy_dec_SDSS[i])
            zsp_R_200.append(redshift[i])
            
        
    zsp_R_200=np.array(zsp_R_200)[~np.isnan(zsp_R_200)]

    new_cluster_z=np.median(zsp_R_200)
    
    new_cluster_z=histogram(zsp_R_200,100,'#0504ab')

 
    print('new cluster z from histogram:' ,new_cluster_z )

    R_200_members=[]

    for i in range(len(zsp_R_200)):
        if (zsp_R_200[i] >= new_cluster_z - 0.015) and (zsp_R_200[i] <= new_cluster_z + 0.015):
            R_200_members.append(zsp_R_200[i]) 

    print(max(R_200_members),min(R_200_members),len(R_200_members))

    #Cluster members
    
    v_rec_vel=np.array(np.array(R_200_members)*constants.c)
    

    median_rec_vel=stats.median(v_rec_vel)
    
    v_rest_vel=np.array(v_rec_vel - median_rec_vel) /(1+median_rec_vel/constants.c)
            
            
    biscl_2 = biweight_scale(v_rest_vel)


    print(' rest velocity calculated from the biweight scale',(biscl_2/10**3))

    sigma_clip=np.array(v_rest_vel)[(v_rest_vel>-3*biscl_2)&(v_rest_vel<3*biscl_2)]


    sigma_cluster = biweight_scale(sigma_clip)/10**3
    
    print(' rest velocity calculated after clip',(sigma_cluster))

    sigma_cluster_z=(sigma_cluster*10**3)/constants.c

    new_cluster_z=np.median(zsp_R_200)

    print('new cluster z from median:' ,new_cluster_z )

    R_200=  R_200_to_deg(R_200_Mpc,new_cluster_z)#5 Mpc

    return(sigma_cluster_z,new_cluster_z,zsp_R_200,R_200)
 

def  cluster_members(folder_path,cluster_name,SDSS_cat,DeCALS_catalogue_fits,DEIMOS_catalogue_fits ,HeCS_catalogue_fits,R_200_Mpc,cluster_centre,sigma_cluster_z,cluster_z_lit):
    
    def R_200_to_deg(distance_Mpc,redshift):

        cosmo = FlatLambdaCDM(H0=70, Om0=0.27)

        def radians_to_degrees(radians):
            degrees = radians * (180 / np.pi)
            return degrees

        d_A = cosmo.angular_diameter_distance(z=redshift)

        theta_radian = distance_Mpc/d_A

        R_200=radians_to_degrees(theta_radian.value)

        return(R_200)
    
    def cluster_members_find(redshift,sigma_cluster_z,cluster_z_lit):


        cluster_member=(redshift>=cluster_z_lit-3*(sigma_cluster_z))&(redshift<=cluster_z_lit+3*(sigma_cluster_z)) 

        return(cluster_member)

    DeCALS_cat= fits.open(DeCALS_catalogue_fits)[1].data
  

    SDSS_cat= fits.open(SDSS_cat)[1].data
    

    galaxy_ra_SDSS = SDSS_cat['ra']

    galaxy_dec_SDSS = SDSS_cat['dec']

    galaxy_ra_DeCALS= DeCALS_cat['ra']

    galaxy_dec_DeCALS = DeCALS_cat['dec']

    # galaxy_ra_DeCALS = DeCALS_cat['ra'][(DeCALS_cat['survey'] != 'SDSS') & (DeCALS_cat['survey'] != 'BOSS') & (DeCALS_cat['survey'] != 'eBOSS-ELG')& (DeCALS_cat['survey'] != 'eBOSS-ELG') ]

    # galaxy_dec_DeCALS = DeCALS_cat['dec'][(DeCALS_cat['survey'] != 'SDSS') & (DeCALS_cat['survey'] != 'BOSS') & (DeCALS_cat['survey'] != 'eBOSS-ELG') & (DeCALS_cat['survey'] != 'eBOSS-ELG')]
    
    # redshift_DeCALS = DeCALS_cat['z_spec'][(DeCALS_cat['survey'] != 'SDSS') & (DeCALS_cat['survey'] != 'BOSS') & (DeCALS_cat['survey'] != 'eBOSS-ELG')& (DeCALS_cat['survey'] != 'eBOSS-ELG')]



    redshift_SDSS=SDSS_cat['zsp']

    mask=DeCALS_cat['z_spec'] > 0

    galaxy_ra_DeCALS= DeCALS_cat['ra'][mask]

    galaxy_dec_DeCALS = DeCALS_cat['dec'][mask]

    redshift_DeCALS=DeCALS_cat['z_spec'][mask]

   
  
    cluster_member_SDSS=cluster_members_find(redshift_SDSS,sigma_cluster_z,cluster_z_lit) 

    cluster_member_DeCALS=cluster_members_find(redshift_DeCALS,sigma_cluster_z,cluster_z_lit) 


    galaxy_members_SDSS=[np.array(galaxy_ra_SDSS)[cluster_member_SDSS] , np.array(galaxy_dec_SDSS)[cluster_member_SDSS] ]

    galaxy_members_DeCALS= [np.array(galaxy_ra_DeCALS)[cluster_member_DeCALS] , np.array(galaxy_dec_DeCALS)[cluster_member_DeCALS] ]


    print(' Number of cluster members SDSS',len(np.array(galaxy_ra_SDSS)[cluster_member_SDSS]) )

    print(' Number of cluster members DeCALS',len(np.array(galaxy_ra_DeCALS)[cluster_member_DeCALS]) )


    if cluster_name == 'Zwcl2341':

        DEIMOS_catalogue = fits.open(DEIMOS_catalogue_fits)
        DEIMOS_cat=DEIMOS_catalogue[1].data
        
        galaxy_ra_DEIMOS = DEIMOS_cat['RAJ2000']
        galaxy_dec_DEIMOS = DEIMOS_cat['DEJ2000']

        redshift_DEIMOS=DEIMOS_cat['zspec']


        cluster_member_DEIMOS=cluster_members_find(redshift_DEIMOS,sigma_cluster_z,cluster_z_lit)
        
        galaxy_members_DEIMOS=[np.array(galaxy_ra_DEIMOS)[cluster_member_DEIMOS] , np.array(galaxy_dec_DEIMOS)[cluster_member_DEIMOS] ]
                                    
        print(' Number of cluster members DEIMOS',len(np.array(galaxy_ra_DEIMOS)[cluster_member_DEIMOS]) )

        
        fig, ax = plt.subplots(figsize=(11, 7))
        pixscale = 0.5
        cutoutsize = 0.25

        size = int(cutoutsize*3600/pixscale)
        url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(cluster_centre.ra.degree)+'&dec='+str(cluster_centre.dec.degree)+'&layer=hsc-dr3&pixscale='+str(pixscale)+'&size='+str(size)
        response = requests.get(url)
        dra = (cutoutsize/2)/np.cos(cluster_centre.dec.degree*np.pi/180)
        ddec = cutoutsize/2
        if response.status_code==200:
            img = Image.open(BytesIO(response.content))
            ax.imshow(img,extent=[cluster_centre.ra.degree+dra,cluster_centre.ra.degree-dra,cluster_centre.dec.degree-ddec,cluster_centre.dec.degree+ddec])


        ax.scatter(galaxy_members_SDSS[0],galaxy_members_SDSS[1],color='green',alpha=0.7,marker="x",label='SDSS')
        ax.scatter(galaxy_members_DeCALS[0],galaxy_members_DeCALS[1],color='yellow',alpha=0.7,marker="x",label='DeCALS')
        ax.scatter(galaxy_members_DEIMOS[0],galaxy_members_DEIMOS[1],color='cyan',alpha=0.5,marker="x",label='DEIMOS')
        ax.set_xlim(cluster_centre.ra.degree-dra, cluster_centre.ra.degree+dra)
        ax.set_ylim(cluster_centre.dec.degree-ddec, cluster_centre.dec.degree+ddec)

        plt.legend()

        plt.xlabel("RA in deg",fontsize=15)
        plt.ylabel("DEC in deg",fontsize=15)
        plt.tight_layout()
        plt.title(cluster_name + r" Member galaxies at $-3 \sigma < z_{cluster} <  3 \sigma$",fontsize=17)
        circle=plt.Circle((cluster_centre.ra.value, cluster_centre.dec.value), R_200_to_deg(R_200_Mpc,cluster_z_lit), edgecolor= 'red',facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed') #, )
        plt.gca().add_patch(circle)
        ax.add_patch(circle)
        ax.text(0,0, "R_200", ha="left", va="top",fontsize=10)
        plt.savefig(folder_path+cluster_name+'_cluster_members.pdf')
        plt.show()
   


    if cluster_name == 'A2631':
    #cluster_member=np.array(redshift_circle)[(redshift_circle>=new_cluster_z-3*(sigma_cluster_z))&(redshift_circle<=new_cluster_z+3*(sigma_cluster_z))]  

        HeCS_catalogue = fits.open(HeCS_catalogue_fits)
        HeCS_cat=HeCS_catalogue[1].data



        def hms_to_degrees(h, m, s):
            return 15 * (h + m/60 + s/3600)  # Convert hours to degrees

        def dms_to_degrees(d, m, s):
            sign = 1 if d >= 0 else -1
            return sign * (abs(d) + m/60 + s/3600)
        
        ra=[]
        dec=[]

        for i in range(len(HeCS_cat['RAJ2000'])):
            h, m, s = np.float64(HeCS_cat['RAJ2000'][i].split())
            ra.append(hms_to_degrees(h, m, s))
            d,m,s = np.float64(HeCS_cat['DEJ2000'][i].split())
            dec.append(dms_to_degrees(d, m, s))

        galaxy_ra_HeCS=ra
        galaxy_dec_HeCS=dec

        # fits_path = '/home/kincaid/saraswati/A2631/A2631_galaxy_radio/image_DI_beam_2poly_beam_A2631.app.restored.fits'
        # hdul = fits.open(fits_path)
        # data = hdul[0].data
        # header=hdul[0].header
        # ls  = data[0,0,:,:]
        # del header['*3']
        # del header['*4']
        # header['WCSAXES']= 2
        # header['NAXIS']= 2
        #w = WCS(header)


        def plot_on_ax(ax, url):
            response = requests.get(url)
            if response.status_code==200:
                img = Image.open(BytesIO(response.content))
                ax.imshow(img,extend=[354,358,-1.5,1.5])
            ax.set_xticks([])
            ax.set_yticks([])
            optical_data = np.array(img)
            return ax,img
    


        cluster_centre=SkyCoord(str(354.41916666667 ), str(0.27666666666667), frame='icrs',unit=(u.deg,u.deg))
        position = SkyCoord(cluster_centre.ra, cluster_centre.dec,frame='fk5',equinox='J2000.0') 
    

        redshift_HeCS=HeCS_cat['cz']/(constants.c/10**3)

        cluster_member_HeCS=cluster_members_find(redshift_HeCS,sigma_cluster_z,cluster_z_lit)
        
        galaxy_members_HeCS=[np.array(galaxy_ra_HeCS)[cluster_member_HeCS] , np.array(galaxy_dec_HeCS)[cluster_member_HeCS] ]
        
        print(' Number of cluster members HeCS',len(np.array(galaxy_ra_HeCS)[cluster_member_HeCS]) )
        

        fig, ax = plt.subplots(figsize=(11, 7))
        pixscale = 0.5
        cutoutsize = 0.25

        size = int(cutoutsize*3600/pixscale)
        url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(cluster_centre.ra.degree)+'&dec='+str(cluster_centre.dec.degree)+'&layer=hsc-dr3&pixscale='+str(pixscale)+'&size='+str(size)
        response = requests.get(url)
        dra = (cutoutsize/2)/np.cos(cluster_centre.dec.degree*np.pi/180)
        ddec = cutoutsize/2
        if response.status_code==200:
            img = Image.open(BytesIO(response.content))
            ax.imshow(img,extent=[cluster_centre.ra.degree+dra,cluster_centre.ra.degree-dra,cluster_centre.dec.degree-ddec,cluster_centre.dec.degree+ddec])

        ax.scatter(galaxy_members_SDSS[0],galaxy_members_SDSS[1],color='green',alpha=0.5,marker="x",label='SDSS')
        ax.scatter(galaxy_members_DeCALS[0],galaxy_members_DeCALS[1],color='yellow',alpha=0.5,marker="x",label='DeCALS')
        ax.scatter(galaxy_members_HeCS[0],galaxy_members_HeCS[1],color='cyan',alpha=0.5,marker="x",label='HeCS')
        ax.set_xlim(cluster_centre.ra.degree-dra, cluster_centre.ra.degree+dra)
        ax.set_ylim(cluster_centre.dec.degree-ddec, cluster_centre.dec.degree+ddec)

        plt.legend()

        plt.xlabel("RA in deg",fontsize=15)
        plt.ylabel("DEC in deg",fontsize=15)
        plt.title(cluster_name + r" Member galaxies at $-3 \sigma < z_{cluster} <  3 \sigma$",fontsize=17)
        circle=plt.Circle((cluster_centre.ra.value, cluster_centre.dec.value), R_200_to_deg(R_200_Mpc,cluster_z_lit), edgecolor= 'red',facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed') #, )
        plt.gca().add_patch(circle)
        ax.add_patch(circle)
        ax.text(0,0, "R_200", ha="left", va="top",fontsize=10)
        plt.savefig(folder_path+cluster_name+'_cluster_members.pdf')
        plt.show()

        return()





def redshift_plots(folder_path,cluster_name,sigma_cluster_z,new_cluster_z, DeCALS_catalogue_fits):

    def clip(data, nsigma):
        '''iteratively removes data until all is within nsigma of the median, then returns the median and std'''
        lennewdata = 0
        lenolddata = data.size
        while lenolddata>lennewdata:
            lenolddata = data.size
            data = data[np.where((data<np.nanmedian(data)+nsigma*np.nanstd(data))&(data>np.nanmedian(data)-nsigma*np.nanstd(data)))]
            lennewdata = data.size
        return np.median(data), np.std(data)

    def ivar2errflux(ivar):
        df = np.sqrt(1/ivar) 
        return df

 
    z_lower=0.25
    z_high=0.30

    DeCALS_catalogue = fits.open(DeCALS_catalogue_fits)
    DeCALS_cat=DeCALS_catalogue[1].data

    mask = np.where((DeCALS_cat['z_spec'] >0 )) #
                    
    
    zsp = DeCALS_cat['z_spec'][mask] 
    zph = DeCALS_cat['z_phot_median'][mask] #z_phot_median photoz_median'



    photoz_err= (DeCALS_cat['z_phot_u68']-DeCALS_cat['z_phot_l68'])/2

    mask_filter= np.where(( photoz_err > 0) & (photoz_err/DeCALS_cat['z_phot_median'] * 100 < 10)& (DeCALS_cat['z_spec'] >0 ))


    # M_err_raw= np.array((np.log10(HSC_catalog['stellar_mass_err68_max'])-np.log10(HSC_catalog['stellar_mass_err68_min']))/2)
    # sfr_err_raw= np.array((np.log10(HSC_catalog['sfr_err68_max'])-np.log10(HSC_catalog['sfr_err68_min']))/2)



    #mask_filter= np.where(((ivar2errflux(DeCALS_cat['flux_ivar_g'])/DeCALS_cat['flux_g']) * 100 < 1) & ((ivar2errflux(DeCALS_cat['flux_ivar_r'])/DeCALS_cat['flux_r'])*100 < 1)& ((ivar2errflux(DeCALS_cat['flux_ivar_z'])/DeCALS_cat['flux_z'])*100 < 1) & (DeCALS_cat['z_spec'] >0 ))

    zsp_filter = DeCALS_cat['z_spec'][mask_filter] 
    zph_filter = DeCALS_cat['z_phot_median'][mask_filter] #

    def ivar2errflux(ivar):
        df = np.sqrt(1/ivar) 
        return df

    upper_bad= np.where((zph_filter > 0.25) & (zph_filter < 0.30) & (zsp_filter > new_cluster_z+5*sigma_cluster_z)) 

    sample= np.where((zph_filter > 0.25) & (zph_filter < 0.30) & (zsp_filter < new_cluster_z+5*sigma_cluster_z) & (zsp_filter > new_cluster_z-5*sigma_cluster_z))

    plt.scatter(zph_filter[upper_bad],zsp_filter[upper_bad])

    plt.scatter(zph_filter[sample],zsp_filter[sample])
    

    # flux_w1_data_bad = ivar2errflux(DeCALS_cat['flux_ivar_w1'][upper_bad])/(DeCALS_cat['flux_w1'][upper_bad] * 100)
    # flux_w2_data_bad = ivar2errflux(DeCALS_cat['flux_ivar_w2'][upper_bad])/(DeCALS_cat['flux_w2'][upper_bad] * 100)

    # flux_w3_data_bad = ivar2errflux(DeCALS_cat['flux_ivar_w3'][upper_bad])/(DeCALS_cat['flux_w3'][upper_bad] * 100)
    # flux_w4_data_bad = ivar2errflux(DeCALS_cat['flux_ivar_w4'][upper_bad])/(DeCALS_cat['flux_w4'][upper_bad] * 100)


    # flux_w1_sample = ivar2errflux(DeCALS_cat['flux_ivar_w1'][sample]) /(DeCALS_cat['flux_w1'][sample] * 100)                 
    # flux_w2_sample = ivar2errflux(DeCALS_cat['flux_ivar_w2'][sample] ) /(DeCALS_cat['flux_w2'][sample]  * 100)  
    # flux_w3_sample = ivar2errflux(DeCALS_cat['flux_ivar_w3'][sample]) /(DeCALS_cat['flux_w3'][sample] * 100)                 
    # flux_w4_sample = ivar2errflux(DeCALS_cat['flux_ivar_w4'][sample] ) /(DeCALS_cat['flux_w4'][sample]  * 100)                  

    
    # # plt.hist(flux_w1_data_bad, bins=30, color='blue', alpha=0.5, label='flux_w1 pollution')
    # # plt.hist(flux_w2_data_bad, bins=30, color='green', alpha=0.5, label='flux_w2 pollution')
    # plt.hist(flux_w3_data_bad, bins=30, color='purple', alpha=0.5, label='flux_w3 pollution')
    # plt.hist(flux_w4_data_bad, bins=30, color='yellow', alpha=0.5, label='flux_w4 pollution')

    # # plt.hist(flux_w1_sample, bins=30, color='red', alpha=0.5, label='flux_w1 sample',histtype='step')
    # # plt.hist(flux_w2_sample, bins=30, color='purple', alpha=0.5, label='flux_w2 sample',histtype='step')
    # plt.hist(flux_w3_sample, bins=30, color='black', alpha=0.5, label='flux_w3 sample',histtype='step')
    # plt.hist(flux_w4_sample, bins=30, color='cyan', alpha=0.5, label='flux_w4 sample',histtype='step')

    # plt.xlabel('Percentage error')
    # plt.ylabel('Frequency')
    # plt.legend()
    # plt.title('Histograms of errors in the individual DeCALS bands for sample sources and pollution sources')

    # # Show the plot
    # plt.show()



    print(len(zsp),'Length before filter')
    print(len(zsp_filter),'Length after filter')
    # target_region_filter = zsp_final[np.where((zsp_final>0.25)&(zsp_final<0.3)&(zph_final>0.2)&(zph_final<0.4))]
    
    # target_region_nofilter = zsp_final_nofilter[np.where((zsp_final_nofilter>0.25)&(zsp_final_nofilter<0.3)&(zph_final_nofilter>0.2)&(zph_final_nofilter<0.4))]
    
    mask_sigma=np.where((zph>z_lower)&(zph<z_high))


    zph_minus_zsp_before = np.array(zph-zsp)[mask_sigma]
    
    
    n, bins = np.histogram(zph_minus_zsp_before, bins=50)

    bin_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2

    #import IPython;IPython.embed()
    fig = plt.figure(figsize=(20, 7))

    mean = np.mean(zph_minus_zsp_before)
    std =np.std(zph_minus_zsp_before)


    x_min =np.min(zph_minus_zsp_before) - 3 * std
    x_max = np.max(zph_minus_zsp_before) + 3 * std


    x = np.linspace(x_min, x_max, 1000)
    FWHM=clip(zph_minus_zsp_before, 3)[1]*2.23
    plt.hist(zph_minus_zsp_before-bin_centre, bins=50, density=True,range=(-0.2,0.2)) 
    plt.plot(x, norm.pdf(x, clip(zph_minus_zsp_before, 3)[0]-bin_centre, clip(zph_minus_zsp_before, 3)[1]))
    plt.xlim(-0.4,0.4)
    plt.text(0.4, 25,f'FWHM is {FWHM:.4f}', fontsize=15, color='red', ha='right', va='top')
    plt.xlabel('zphot - zpec ',fontsize=20)
    plt.title(r"zphot $\sigma$ determination for " + cluster_name,fontsize=22)
    plt.savefig(folder_path+cluster_name+'_delta_photo_z_distribution.pdf')
    plt.show()
   


    fig = plt.figure(figsize=(20, 10))


    plt.subplot(221)
    plt.scatter(zph,zsp)
    plt.plot(zsp, zsp, 'k-', lw=2,label='one-to-one')
    plt.axhline(y=new_cluster_z-5*sigma_cluster_z, color='r', linestyle='-')
    plt.axhline(y=new_cluster_z+5*sigma_cluster_z, color='r', linestyle='-')
    plt.axvline(x=0, color='r', linestyle='-')
    plt.axvline(x=1, color='r', linestyle='-')
    plt.title(cluster_name + " photo z vs spec z for galaxies at 1 > z > 0")
    #plt.title('Outlier ratio is {}'.format(mag_1))
    plt.xlabel('photometric redshift')
    plt.ylabel('spectroscopic redshift')
    plt.xlim(0,1)
    plt.ylim(0,1)


    plt.subplot(222)
    plt.scatter(zph,zsp)
    plt.plot(zsp, zsp, 'k-', lw=2,label='one-to-one')
    plt.axhline(y=new_cluster_z-5*sigma_cluster_z, color='r', linestyle='-')
    plt.axhline(y=new_cluster_z+5*sigma_cluster_z, color='r', linestyle='-')
    plt.axvline(x=z_lower, color='r', linestyle='-')
    plt.axvline(x=z_high, color='r', linestyle='-')
    plt.title('Zoom in of photo z vs spec z for galaxies at  1 > z > 0 showing \n the redshift sample of interest' )
    #plt.title('Outlier ratio is {}'.format(mag_1))
    plt.xlabel('photometric redshift')
    plt.ylabel('spectroscopic redshift')
    plt.xlim(0.1,0.5)
    plt.ylim(0.2,0.35)

    plt.subplot(223)
    plt.scatter(zph_filter,zsp_filter)
    plt.plot(zsp_filter, zsp_filter, 'k-', lw=2,label='one-to-one')
    plt.axhline(y=new_cluster_z-5*sigma_cluster_z, color='r', linestyle='-')
    plt.axhline(y=new_cluster_z+5*sigma_cluster_z, color='r', linestyle='-')
    plt.axvline(x=0, color='r', linestyle='-')
    plt.axvline(x=1, color='r', linestyle='-')
    plt.title(cluster_name + " photo z vs spec z for galaxies at 1 > z > 0 \n after magnitude filter")
    #plt.title('Outlier ratio is {}'.format(mag_1))
    plt.xlabel('photometric redshift')
    plt.ylabel('spectroscopic redshift')
    plt.xlim(0,1)
    plt.ylim(0,1)


    plt.subplot(224)
    plt.scatter(zph_filter,zsp_filter)
    plt.plot(zsp_filter, zsp_filter, 'k-', lw=2,label='one-to-one')
    plt.axhline(y=new_cluster_z-5*sigma_cluster_z, color='r', linestyle='-')
    plt.axhline(y=new_cluster_z+5*sigma_cluster_z, color='r', linestyle='-')
    plt.axvline(x=z_lower, color='r', linestyle='-')
    plt.axvline(x=z_high, color='r', linestyle='-')
    plt.title(cluster_name + " Zoom in of photo z vs spec z for galaxies at  1 > z > 0 showing \n the redshift sample of interest after magnitude filter")
    #plt.title('Outlier ratio is {}'.format(mag_1))
    plt.xlabel('photometric redshift')
    plt.ylabel('spectroscopic redshift')
    plt.xlim(0.1,0.5)
    plt.ylim(0.2,0.35)
    plt.savefig(folder_path+cluster_name+'_photo_z_determination.pdf')
    plt.show()

    
def M_SFR_HSC(folder_path,cluster_name,HSC_catalog_fits,radio_HSC_catalog_fits,cluster_z_lit,sigma_cluster_z):


    HCS_catalogue = fits.open(HSC_catalog_fits)
    HSC_cat=HCS_catalogue[1].data

    radio_HCS_catalogue = fits.open(radio_HSC_catalog_fits)
    radio_HSC_cat=radio_HCS_catalogue[1].data


    h=0.7
    H0=h*100
    cosmo = cp.FlatLambdaCDM(H0=h*100, Om0=0.30)



    def redshift_to_age(z):

        t=(Planck13.age(z)).value

        return(t)


    def popetso_MS1(logM,t):

        a0=0.20
        a1=-0.034
        b0=-26.134
        b1=4.722
        b2=-0.1925

        MS_1= (a1*t+b1)*logM + b2*logM**2+(b0+a0*t)

        return(MS_1)
    

    def AGN_MS(logM,z):
        p1=20.97
        p2=2.51
        p3=1.411

        AGN_MS= p1 + p2*np.log10(1+z)+p3*(logM-10)

        return(AGN_MS)


    
    SFR_constant = 3.18e-22


    M_raw=np.array(HSC_cat['stellar_mass'])

    sfr_raw=np.array(HSC_cat['sfr'])



    M_raw_radio=np.array(radio_HSC_cat['stellar_mass'])

    sfr_raw_radio=np.array(radio_HSC_cat['sfr'])

   

    photoz_err= (HSC_cat['photoz_err68_max']-HSC_cat['photoz_err68_min'])/2


    photoz_err_radio= (radio_HSC_cat['photoz_err68_min']-radio_HSC_cat['photoz_err68_max'])/2


    M_err_raw= np.array((np.log10(HSC_cat['stellar_mass_err68_max'])-np.log10(HSC_cat['stellar_mass_err68_min']))/2)
    sfr_err_raw= np.array((np.log10(HSC_cat['sfr_err68_max'])-np.log10(HSC_cat['sfr_err68_min']))/2)





    M_err_raw_radio= np.array((np.log10(radio_HSC_cat['stellar_mass_err68_max'])-np.log10(radio_HSC_cat['stellar_mass_err68_min']))/2)
    sfr_err_raw_radio = np.array((np.log10(radio_HSC_cat['sfr_err68_max'])-np.log10(radio_HSC_cat['sfr_err68_min']))/2)

    thresh=50

    mask= np.where(((HSC_cat['photoz_median']  > (cluster_z_lit-cluster_z_lit*sigma_cluster_z)) & (HSC_cat['photoz_median']  < (cluster_z_lit+cluster_z_lit*sigma_cluster_z)))) #& ((M_err_raw / M_raw) * 100 < thresh) & ((sfr_err_raw / sfr_raw) * 100  < thresh))
    #mask_radio= np.where(( radio_HSC_cat['photoz_median'] > 0) & (photoz_err/radio_HSC_cat['photoz_median'] * 100 < 10)& (radio_HSC_cat['photoz_median'] >0 )&(photoz > 0.25) & (photoz < 0.3))


    stellar_photoz=HSC_cat['photoz_median']

    radio_photoz=radio_HSC_cat['photoz_median']

    mask_radio = np.where(np.isin(radio_photoz, stellar_photoz))
    
    print('min photoz',np.min(stellar_photoz))

    print('max photoz',np.max(stellar_photoz))

    print('Total length of photoz before',len(HSC_cat['photoz_median']))

    print('Total length of photoz radio before',len(radio_photoz))
 

    M=M_raw[mask]
    #M_err=M_err_raw[mask]

    M_radio=M_raw_radio[mask_radio]
    #M_err_radio=M_err_raw_radio[mask_radio]
  
    sfr=sfr_raw[mask]
    #sfr_err=sfr_err_raw[mask]

    sfr_radio=sfr_raw_radio[mask_radio]
    #sfr_err_radio=sfr_err_raw_radio[mask_radio]

    print('Total length of photoz after',len(stellar_photoz))

    print('Total length of photoz radio after',len(radio_photoz[mask_radio]))
    

    z = [cluster_z_lit]

    tages = redshift_to_age(z)

    

    mask1=(np.log10(sfr_radio) > -5)   &  (np.log10(sfr_radio) < -3)

    mask2=(np.log10(sfr_radio) > -3) &  (np.log10(sfr_radio) < -1)

    mask3=(np.log10(sfr_radio) > -1) &  (np.log10(sfr_radio) < 3)


    colors = cm.Spectral(np.linspace(2, 0, len(tages)))

    fig = plt.figure(figsize=(13, 5))

    for t_ind,  t in enumerate(tages):

        logM=np.linspace(8.5,12,100)

        plt.plot(logM,popetso_MS1(logM,t) ,label='MS eqn 10 from Popesso at $t$='+"{0:.1f}".format(t)+ ' Gyr ' ,color=colors[t_ind],linestyle='--')


    plt.scatter(np.log10(M),np.log10(sfr),label='Galaxy stellar MS',alpha=0.7,s=15,color='blue')


    plt.scatter(np.log10(M_radio)[mask1],np.log10(sfr_radio)[mask1],label='Main Sequence galaxies',marker='x',color='orange')

    plt.scatter(np.log10(M_radio)[mask2],np.log10(sfr_radio)[mask2],label='Intermediate SFR galaxies',marker='x',color='red')

    plt.scatter(np.log10(M_radio)[mask3],np.log10(sfr_radio)[mask3],label='Quiescent low SFR galaxies',marker='x',color='green')
   
    #plt.errorbar(np.log10(M),np.log10(sfr), xerr=np.log10(M_err), fmt='none', capsize=4,color='orange')
    #plt.errorbar(np.log10(M),np.log10(sfr), yerr=np.log10(sfr_err), fmt='none', capsize=4,color='black')
    plt.xlim(7,12)

    # plt.ylim(-2,3)
    plt.title(f"{cluster_name} mizuki derived Solar masses and SFR for member galaxies ", fontsize=14)
    plt.xlabel('$log M/M_\odot$ using HCS photometric redshift',fontsize=17)
    plt.ylabel('$log$ SFR',fontsize=17)
    plt.tight_layout()
    plt.legend()
    plt.savefig(folder_path+cluster_name+'_galaxy_stellar_MS.pdf')
    plt.show()
	

  



def radio_SFR_DeCALS(folder_path,cluster_name,FAST_cat,radio_DeCALS_catalogue,cluster_z_lit ):

    
    radio_DeCALS=Table.read(radio_DeCALS_catalogue)
    columns_to_keep = ['id' ,'z','l68_z', 'u68_z', 'ltau' , 'l68_ltau' , 'u68_ltau' ,  'metal', 'l68_metal', 'u68_metal' , 'lage' , 'l68_lage' , 'u68_lage'  , 'Av' , 'l68_Av' , 'u68_Av', 'lmass' ,'l68_lmass' ,'u68_lmass' , 'lsfr','l68_lsfr','u68_lsfr' ,  'lssfr' ,'l68_lssfr' , 'u68_lssfr' ,  'la2t' ,  'l68_la2t' , 'u68_la2t' ,  'chi2']
    
    FAST_catalog = pd.read_table(FAST_cat, sep="\s+", usecols=columns_to_keep)

    t2 = Table.from_pandas(FAST_catalog)

    
    T=join(t2, radio_DeCALS,keys_left='id',keys_right='Isl_id')

    
    flux_MeerKAT=T['Total_flux']

    flux_err_MeerKAT=T['E_Total_flux']

    #spec_indx=radio_DeCALS['spec_indx']

    def flux2mag(flux):
        return (22.5-2.5*np.log10(flux))

    mag_W1=flux2mag(T['FLUX_W1'])
    mag_W2=flux2mag(T['FLUX_W2'])


    lmass=T['lmass']

    u68_lmass=T['u68_lmass']

    l68_lmass=T['l68_lmass']



    FAST_id=T['id']

    z=T['z']
    
    def redshift_to_age(z):

        t=(Planck13.age(z)).value

        return(t)



    def AGN_MS(logM,z):
        p1=20.97
        p2=2.51
        p3=1.411

        AGN_MS= p1 + p2*np.log10(1+z)+p3*(logM-10)

        return(AGN_MS)

    def Wise_R90(W1,W2):
        a=np.array(0.650)
        b=np.array(0.153)
        y=np.array(13.86)
        
        if W2.any() > y:

            R90_W=a*np.exp(b*(W2-y)**2)
        else:
            R90_W=np.array(a)

        return(R90_W)


    def AGN_SF_ratio(logM,z):

        q1=-0.63
        q2=-0.05
        q3=0.61

        AGN_MS= q1 + q2*np.log10(1+z)+q3*(logM-10)

        return(AGN_MS)


    def radio_lum(flux,z,si,freq):


        S_14=(np.array(flux)/(freq**si)) * (1400e6**si)

        lumo_distance=np.array(cosmo.luminosity_distance(z))

        radio_lumo=(4*np.pi*(S_14*10**-26)*(lumo_distance*3.086e22)**2)*(1/(1+(z**(si+1))))

        return(radio_lumo)


    def popetso_MS1(logM,t):

        a0=0.20
        a1=-0.034
        b0=-26.134
        b1=4.722
        b2=-0.1925

        MS_1= (a1*t+b1)*logM + b2*logM**2+(b0+a0*t)

        return(MS_1)

    def popetso_MS2(M,t):

        a0=2.693
        a1=-0.186
        a2=10.85
        a3=-0.0729
        a4=0.99

        MS_2 = a0+a1*t - np.log10(1+(M/(10**(a2+a3*t)))**-a4)

        return(MS_2)
    

    h=0.7
    H0=h*100
    cosmo = cp.FlatLambdaCDM(H0=h*100, Om0=0.30)

    #import IPython;IPython.embed()
    radio_lumo=radio_lum(flux_MeerKAT,z,-0.7,1283e6)

    radio_lumo_err=radio_lum(flux_err_MeerKAT,z,-0.7,1283e6)


    index_radio_loud=radio_lumo < 10**23

    index_wise=(mag_W1-mag_W2) <Wise_R90(mag_W1,mag_W2)
    
    combined_mask =  np.logical_or(index_radio_loud,index_wise)


    SFR_constant=3.18e-22
    SFR = 3.18e-22*np.array(radio_lumo[combined_mask] )


    SFR_err=3.18e-22*np.array(radio_lumo_err[combined_mask])

    Mstar_SDSS=lmass[combined_mask]

    Mstar_SDSS_err=((u68_lmass-l68_lmass)/2)[combined_mask]


    #mask= np.where((Mstar_SDSS_err_x/Mstar_SDSS_x * 100 < 30) & (SFR_err/SFR * 100 < 30))

    # SFR=SFR[mask]
    # SFR_err= SFR_err[mask]


    # Mstar_SDSS=Mstar_SDSS_x[mask]


    # Mstar_SDSS_err=Mstar_SDSS_err_x[mask]



    print("Length after filter",len(Mstar_SDSS))

    #import IPython;IPython.embed()

    z = [cluster_z_lit]

    tages = redshift_to_age(z)


    colors = cm.Spectral(np.linspace(2, 0, len(tages)))
    fig = plt.figure(figsize=(13, 5))
    for t_ind,  t in enumerate(tages):
    
        logM=np.linspace(8.5,12,100)
        M=10**(logM)
        
        plt.plot(logM,popetso_MS1(logM,t),label='MS eqn 10 from Popesso at $t$='+"{0:.1f}".format(t)+ 'Gyr ' ,color=colors[t_ind],linestyle='--')
        
        plt.plot(logM,popetso_MS2(M,t),label='MS eqn 14 from Popesso at $t$='+"{0:.1f}".format(t)+ 'Gyr ' ,color=colors[t_ind])
        
        plt.plot(logM,np.log10(SFR_constant*10**(AGN_MS(logM,t))),label='AGN Main Sequence at $t$='+"{0:.1f}".format(t)+ 'Gyr ' ,color=colors[t_ind],linestyle=':')
        mid_x=np.mean(logM)
        mid_y=np.mean(np.log10(SFR_constant*10**(AGN_MS(logM,t))))

        plt.text(mid_x, mid_y+0.2, 'AGN MS', fontsize=12, color='blue', ha='center', va='center',rotation=18)

    plt.scatter(Mstar_SDSS[np.where(SFR)],np.log10(SFR),marker='o',color='b',label='Galaxy MS')

    plt.errorbar(Mstar_SDSS[np.where(SFR)],np.log10(SFR), xerr=np.abs(Mstar_SDSS_err[np.where(SFR)]), fmt='none', capsize=4,color='orange')
    plt.errorbar(Mstar_SDSS[np.where(SFR)],np.log10(SFR), yerr=np.abs(np.log10(SFR_err)), fmt='none', capsize=4,color='black')

 
    plt.xlim(8.5,12)
    # plt.ylim(-1,3)
    plt.title(f"{cluster_name}", fontsize=15)     
    plt.xlabel('$log M/M_\odot$ using DeCALS photometric redshift',fontsize=17)
    plt.ylabel('$log$ Radio SFR',fontsize=17)
    plt.tight_layout()
    plt.legend()
    plt.savefig(folder_path+cluster_name+'_galaxy_MS.pdf')
    plt.show()  


    
def red_sequence(folder_path,cluster_name,optical_catalog,cluster_centre,R_200):

    def flux2mag(flux):
        return (22.5-2.5*np.log10(flux))

    mag_filter=27
    optical_cat=fits.open(optical_catalog)[1].data

    mask= np.where((np.sqrt((optical_cat['ra']-cluster_centre.ra.value)**2+(optical_cat['dec']-cluster_centre.dec.value))**2<= 15/60) & (optical_cat['Z_PHOT_MEDIAN'] > 0.15) &  (optical_cat['Z_PHOT_MEDIAN'] < 0.35) & (flux2mag(optical_cat['APFLUX_G'][:,0]) < mag_filter) & (flux2mag(optical_cat['APFLUX_R'][:,0]) < mag_filter))

    optical_cat=optical_cat[mask]


    mag_g=flux2mag(optical_cat['APFLUX_G'][:,0])
    mag_r=flux2mag(optical_cat['APFLUX_R'][:,0])

    redshift=optical_cat['Z_PHOT_MEDIAN']

    #import IPython;IPython.embed()
    tot_mag_r=flux2mag(optical_cat['FLUX_R'])


    
    def ivar2errflux(ivar):
        df = np.sqrt(1/ivar) 
        return df

    
    galaxy_ra=np.array(optical_cat['ra'])
    galaxy_dec=np.array(optical_cat['dec'])
   
    

    scatter=plt.scatter(tot_mag_r,mag_g-mag_r,c=redshift, cmap='viridis')


    scatter=plt.scatter(tot_mag_r,mag_g-mag_r,c=redshift, cmap='viridis')
    
    plt.colorbar(scatter, label='redshift')


    plt.xlabel(r'$mag_r$')
    plt.ylabel(r'$mag_g - mag_r$')
    # plt.ylim(-2,3)
    plt.title(cluster_name+ ' ' +' DeCALS red sequence')
    plt.legend()

    plt.show()


    
def SDSS_DeCALS_Mstar(folder_path,cluster_name,SDSS_FAST_cat,DeCALS_FAST_cat):
                
    columns_to_keep = ['id' ,'z','l68_z', 'u68_z', 'ltau' , 'l68_ltau' , 'u68_ltau','metal', 'l68_metal', 'u68_metal' , 'lage' , 'l68_lage' , 'u68_lage'  , 'Av' , 'l68_Av' , 'u68_Av', 'lmass' ,'l68_lmass' ,'u68_lmass' , 'lsfr','l68_lsfr','u68_lsfr' ,  'lssfr' ,'l68_lssfr' , 'u68_lssfr' ,  'la2t' ,  'l68_la2t' , 'u68_la2t' ,  'chi2']

    SDSS_zwcl = pd.read_table(SDSS_FAST_cat, sep="\s+", usecols=columns_to_keep)

    DeCALS_zwcl = pd.read_table(DeCALS_FAST_cat, sep="\s+", usecols=columns_to_keep)

    # SDSS_A2631 = pd.read_table(A2631_SDSS_cat, sep="\s+", usecols=columns_to_keep)

    # DeCALS_A2631= pd.read_table(A2631_DeCALS_cat, sep="\s+", usecols=columns_to_keep)



    # Mstar_SDSS= np.concatenate((np.array(SDSS_zwcl['lmass']) , np.array(SDSS_A2631['lmass'])))

    # Mstar_DeCALS= np.concatenate((np.array(DeCALS_zwcl['lmass']) , np.array(DeCALS_A2631['lmass'])))

    # l68_Mstar_SDSS= np.concatenate((np.array(SDSS_zwcl['l68_lmass']) , np.array(SDSS_A2631['l68_lmass'])))

    # u68_Mstar_SDSS= np.concatenate((np.array(SDSS_zwcl['u68_lmass']) , np.array(SDSS_A2631['u68_lmass'])))

    # l68_Mstar_DeCALS=  np.concatenate((np.array(DeCALS_zwcl['l68_lmass']) , np.array(DeCALS_A2631['l68_lmass'])))

    # u68_Mstar_DeCALS=  np.concatenate((np.array(DeCALS_zwcl['u68_lmass']) , np.array(DeCALS_A2631['u68_lmass'])))


    Mstar_SDSS= np.array(SDSS_zwcl['lmass']) 

    Mstar_DeCALS= np.array(DeCALS_zwcl['lmass']) 

    l68_Mstar_SDSS= np.array(SDSS_zwcl['l68_lmass']) 

    u68_Mstar_SDSS= np.array(SDSS_zwcl['u68_lmass']) 

    l68_Mstar_DeCALS=  np.array(DeCALS_zwcl['l68_lmass']) 

    u68_Mstar_DeCALS=  np.array(DeCALS_zwcl['u68_lmass']) 


    # Mstar_SDSS= np.array(SDSS_A2631['lmass']) 

    # Mstar_DeCALS= np.array(DeCALS_A2631['lmass']) 

    # l68_Mstar_SDSS= np.array(SDSS_A2631['l68_lmass']) 

    # u68_Mstar_SDSS= np.array(SDSS_A2631['u68_lmass']) 

    # l68_Mstar_DeCALS=  np.array(DeCALS_A2631['l68_lmass']) 

    # u68_Mstar_DeCALS=  np.array(DeCALS_A2631['u68_lmass']) 



    error_upper_DeCALS=u68_Mstar_DeCALS-Mstar_DeCALS

    error_lower_DeCALS=Mstar_DeCALS-l68_Mstar_DeCALS



    error_upper_DeCALS_perc_err=(error_upper_DeCALS/Mstar_DeCALS)*100

    error_lower_DeCALS_perc_err=(error_lower_DeCALS/Mstar_DeCALS)*100

    filtered_indices = (error_upper_DeCALS_perc_err <= 10) & (error_lower_DeCALS_perc_err <= 10)

    filtered_error_upper_DeCALS= error_upper_DeCALS[filtered_indices]

    filtered_error_lower_DeCALS= error_lower_DeCALS[filtered_indices]

    error_DeCALS= np.array([filtered_error_upper_DeCALS,filtered_error_lower_DeCALS])



    error_upper_SDSS=u68_Mstar_SDSS-Mstar_SDSS

    error_lower_SDSS=Mstar_SDSS-l68_Mstar_SDSS



    error_upper_SDSS_perc_err=(error_upper_SDSS/Mstar_SDSS)*100

    error_lower_SDSS_perc_err=(error_lower_SDSS/Mstar_SDSS)*100

    filtered_indices = (error_upper_SDSS_perc_err <= 10) & (error_lower_SDSS_perc_err <= 10)

    filtered_error_upper_SDSS= error_upper_SDSS[filtered_indices]

    filtered_error_lower_SDSS= error_lower_SDSS[filtered_indices]

    error_SDSS= np.array([filtered_error_upper_SDSS,filtered_error_lower_SDSS])


    print("len SDSS",len(error_SDSS[0]))

    print("len DeCALS",len(error_DeCALS[0]))

    fig = plt.figure(figsize=(6, 4))
    plt.scatter(Mstar_SDSS[:len(error_SDSS[0])],Mstar_DeCALS[:len(error_SDSS[0])])
    plt.errorbar(Mstar_SDSS[:len(error_SDSS[0])],Mstar_DeCALS[:len(error_SDSS[0])], yerr=np.resize(error_DeCALS,error_SDSS.shape),xerr=error_SDSS,fmt='k.') 
    plt.title(r'FAST derived Solar masses at $ 0 < z < 0.5 $ for ' + cluster_name)
    plt.xlabel('log$M/M_\odot$ using DeCALS spectroscopic redshift')
    plt.ylabel('log$M/M_\odot$ using DeCALS photometric redshift')
    plt.plot(Mstar_SDSS,Mstar_SDSS)
    plt.savefig(folder_path+cluster_name+'_FAST_derivied_Mstar.pdf')

    plt.show()
        


