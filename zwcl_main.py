import sys
from zwcl2341_plots import DeCALS_query,circle, redshift_distribution,cluster_members,velocity_dispersion,redshift_plots,RA_DEC_seperation_optical,SDSS_DeCALS_Mstar,red_sequence,circle,RA_DEC_seperation_radio, radio_SFR_DeCALS,LR,general_seperation,M_SFR_HSC,radio_redshift_DeCLAS,stack_table
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import numpy as np
from astropy.io import fits

def zwcl_main( argv ):
    
    cluster_name='Zwcl2341'


    folder_path= cluster_name+'/plots/'

    output_filename_DeCALS = folder_path+cluster_name+'DeCALS.fits'
    # ---------------------------------------------------------
    #  SDSS and DeCALS combined catalogues for photo vs spec tests
    # ---------------------------------------------------------

    DeCALS_catalog = 'Zwcl2341_DeCALS_dr10_full.fits'

    radio_cat='Zwcl2341_srl.fits'

    SDSS_cat='zwcl2341_90_DR17.fits'

    HSC_DeCALS_catalogue_fits = 'Zwcl2341_HSC_DeCALS_join.fits'

    HSC_catalog_fits='Zwcl2341_HSC.fits'
    radio_HSC_catalog_fits='Zwcl2341_HSC_DeCALS_FAST.fits'
    radio_DeCALS_catalogue = 'Zwcl2341_DeCALS_radio_join.fits'
    DEIMOS_catalogue_fits='DEISMOS_table.fits'
    HeCS_catalogue_fits='HeCS_table.fits'

    # ---------------------------------------------------------
                 
    # ---------------------------------------------------------
    #  FAST catalogues for 0 < z< 0.5 to compare SDSS vs DeCALS Mstar
    # ---------------------------------------------------------
                                 
    zspec_FAST_Mstar_cat = 'observations_zwcl_zspec.fout'

    zphot_FAST_Mstar_cat = 'observations_zwcl_zphot_median.fout'


    # ---------------------------------------------------------
    #  FAST catalogues for MS
    # ---------------------------------------------------------

    FAST_cat = 'observations_zwcl_DeCALS_SFR.fout'

    columns_to_keep = ['id' ,'z','l68_z', 'u68_z', 'ltau' , 'l68_ltau' , 'u68_ltau' ,  'metal', 'l68_metal', 'u68_metal' , 'lage' , 'l68_lage' , 'u68_lage'  , 'Av' , 'l68_Av' , 'u68_Av', 'lmass' ,'l68_lmass' ,'u68_lmass' , 'lsfr','l68_lsfr','u68_lsfr' ,  'lssfr' ,'l68_lssfr' , 'u68_lssfr' ,  'la2t' ,  'l68_la2t' , 'u68_la2t' ,  'chi2']

    SDSS = pd.read_table(FAST_cat, sep="\s+", usecols=columns_to_keep)

    u68_lmass=SDSS['u68_lmass']

    l68_lmass=SDSS['l68_lmass']

    lmass=SDSS['lmass']

    FAST_id=SDSS['id']

    FAST_z=SDSS['z']

    # ---------------------------------------------------------------
    #  Radio and DeCALS combined catalogues 
    # ---------------------------------------------------------------



    # ---------------------------------------------------------------
    #  cluster general paramaters 
    # ---------------------------------------------------------------


    cluster_centre=SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))
    
    zsp_min=0.24
    zsp_max=0.30

    zsp_min_1 = 0
    zsp_max_1 = 1
    
    search_radius = 1.5

    cluster_z_lit=0.2694

    R_200_Mpc=  1.62 

  
    sigma_cluster_z =  0.0018


    # ---------------------------------------------------------

    plotting_scripts(folder_path,cluster_name,output_filename_DeCALS,cluster_centre,HSC_catalog_fits,radio_HSC_catalog_fits,DeCALS_catalog,HSC_DeCALS_catalogue_fits,DEIMOS_catalogue_fits,HeCS_catalogue_fits,
                     zsp_min,zsp_max,search_radius,R_200_Mpc,zspec_FAST_Mstar_cat ,
                     zphot_FAST_Mstar_cat,FAST_cat,radio_DeCALS_catalogue, radio_cat,sigma_cluster_z,cluster_z_lit,SDSS_cat)
     
def A2631_main( argv ):

    cluster_name='A2631'
    folder_path= cluster_name+'/plots/'

    output_filename_DeCALS = folder_path+cluster_name+'DeCALS.fits'

    # ---------------------------------------------------------
    #  SDSS and DeCALS combined catalogues for photo vs spec tests
    # ---------------------------------------------------------

    SDSS_cat='A2631_90_DR17.fits'

    DeCALS_catalog = 'A2631_DeCALS_dr10_full.fits'

    HSC_DeCALS_catalogue_fits = 'A2631_HSC_DeCALS_join.fits'

    radio_cat='A2631_srl.fits'

    HSC_catalog_fits='A2631_HSC.fits'

    radio_HSC_catalog_fits='A2631_HSC_DeCALS_FAST.fits'

    DEIMOS_catalogue_fits='DEISMOS_table.fits'

    HeCS_catalogue_fits='HeCS_table.fits'

    radio_DeCALS_catalogue = 'A2631_DeCALS_radio_join.fits'
    # ---------------------------------------------------------
    #  catalogues for 0 < z< 0.5 to compare SDSS vs DeCALS Mstar
    # ---------------------------------------------------------

    zspec_FAST_Mstar_cat = 'observations_A2631_zspec.fout' 

    zphot_FAST_Mstar_cat  ='observations_A2631_zphot_median.fout'

    FAST_cat='observations_A2631_DeCALS_SFR.fout'

    # ---------------------------------------------------------


    # ---------------------------------------------------------
    #  FAST catalogues for MS
    # ---------------------------------------------------------
    DeCALS_FAST_cat = 'observations_A2631_DeCALS_SFR.fout'

    columns_to_keep = ['id' ,'z','l68_z', 'u68_z', 'ltau' , 'l68_ltau' , 'u68_ltau' ,  'metal', 'l68_metal', 'u68_metal' , 'lage' , 'l68_lage' , 'u68_lage'  , 'Av' , 'l68_Av' , 'u68_Av', 'lmass' ,'l68_lmass' ,'u68_lmass' , 'lsfr','l68_lsfr','u68_lsfr' ,  'lssfr' ,'l68_lssfr' , 'u68_lssfr' ,  'la2t' ,  'l68_la2t' , 'u68_la2t' ,  'chi2']
    SDSS = pd.read_table(DeCALS_FAST_cat, sep="\s+", usecols=columns_to_keep)

    u68_lmass=SDSS['u68_lmass']

    l68_lmass=SDSS['l68_lmass']

    lmass=SDSS['lmass']

    FAST_id=SDSS['id']

    FAST_z=SDSS['z']



    # ---------------------------------------------------------------
    #  cluster general paramaters 
    # ---------------------------------------------------------------


    cluster_centre=SkyCoord(str(354.41916666667 ), str(0.27666666666667), frame='icrs',unit=(u.deg,u.deg))
    
    cluster_z_lit=0.2765
    sigma_cluster_z =  0.0028

    zsp_min=0.24
    zsp_max=0.30
    
    search_radius = 1.5

    R_200_Mpc=  1.53 



    plotting_scripts(folder_path,cluster_name,output_filename_DeCALS,cluster_centre,HSC_catalog_fits,radio_HSC_catalog_fits,DeCALS_catalog,HSC_DeCALS_catalogue_fits,DEIMOS_catalogue_fits,HeCS_catalogue_fits,
                     zsp_min,zsp_max,search_radius,R_200_Mpc,zspec_FAST_Mstar_cat ,
                     zphot_FAST_Mstar_cat,FAST_cat,radio_DeCALS_catalogue, radio_cat,sigma_cluster_z,cluster_z_lit,SDSS_cat)





def  plotting_scripts(folder_path,cluster_name,output_filename_DeCALS,cluster_centre,HSC_catalog_fits,radio_HSC_catalog_fits,DeCALS_catalog,HSC_DeCALS_catalogue_fits,DEIMOS_catalogue_fits,HeCS_catalogue_fits,
                     zsp_min,zsp_max,search_radius,R_200_Mpc,zspec_FAST_Mstar_cat ,
                     zphot_FAST_Mstar_cat,FAST_cat,radio_DeCALS_catalogue, radio_cat,sigma_cluster_z,cluster_z_lit,SDSS_cat):
    

    # DeCALS_query(folder_path,cluster_name,cluster_centre)
    # stack_table(folder_path,cluster_name )
    #radio_redshift_DeCLAS(folder_path,cluster_name,radio_DeCALS_catalogue)
    #general_seperation(folder_path,cluster_name,DeCALS_catalog,radio_cat)
    #RA_DEC_seperation_optical(folder_path,cluster_name,SDSS_cat,DeCALS_cat)

    #RA_DEC_seperation_radio(folder_path,cluster_name,radio_DeCALS_catalogue)
    
    # galaxy_ra_circle,galaxy_dec_circle,redshift_circle=circle(folder_path,cluster_name,galaxy_ra_SDSS,galaxy_dec_SDSS,zsp,cluster_centre,zsp_min,zsp_max,search_radius)

    # redshift_distribution(folder_path,cluster_name,galaxy_ra_circle,galaxy_dec_circle,redshift_circle,cluster_centre,zsp_min,zsp_max)
   
    #sigma_cluster_z,new_cluster_z,zsp_R_200,R_200 =velocity_dispersion(cluster_centre,R_200_Mpc,zsp_min,zsp_max,cluster_z_lit,DeCALS_catalog)
    
    #cluster_members(folder_path,cluster_name,SDSS_cat,DeCALS_catalog,DEIMOS_catalogue_fits ,HeCS_catalogue_fits,R_200_Mpc,cluster_centre,sigma_cluster_z,cluster_z_lit)

    
    redshift_plots(folder_path,cluster_name,sigma_cluster_z,cluster_z_lit, DeCALS_catalog)


    #red_sequence(folder_path,cluster_name,DeCALS_catalog,cluster_centre,R_200_Mpc)

    #M_SFR_HSC(folder_path,cluster_name,HSC_catalog_fits,radio_HSC_catalog_fits,cluster_z_lit,sigma_cluster_z)

    #SDSS_DeCALS_Mstar(folder_path,cluster_name,zspec_FAST_Mstar_cat,zphot_FAST_Mstar_cat)

    #radio_SFR_DeCALS(folder_path,cluster_name,FAST_cat,radio_DeCALS_catalogue,cluster_z_lit )
    
    #LR(folder_path,radio_cat,DeCALS_catalog,cluster_centre,cluster_name)
    
    

if __name__ == "__main__":
    zwcl_main(sys.argv)
    A2631_main(sys.argv)


