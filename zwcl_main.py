
import sys
from zwcl2341_plots import circle, redshift_distribution,cluster_members,zwcl_galaxy_distribution,velocity_dispersion,redshift_plots,radio_SFR_plots,RA_DEC_seperation,AGN_relation,SDSS_DeCALS,colour_colour,circle
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import numpy as np
from astropy.io import fits

def main( argv ):
    

    SDSS_catalogue_fits = fits.open(r'Zwcl2341_DeCALS_SDSS.fits')
    SDSS_catalogue=SDSS_catalogue_fits[1].data
                                 

    #SDSS_DeCALS_catalogue=pd.read_csv(r'SDSS_Spec_Photo_DeCALS.csv')
    galaxy_ra=SDSS_catalogue['ra_1']
    galaxy_dec=SDSS_catalogue['dec_1']
    zph=SDSS_catalogue['z_spec']
    zsp=SDSS_catalogue['zsp']

    u_mag=SDSS_catalogue['u']
    u_err=SDSS_catalogue['u_err']
    g_mag=SDSS_catalogue['g']
    g_err=SDSS_catalogue['g_err']
    r_mag=SDSS_catalogue['r']
    r_err=SDSS_catalogue['r_err']
    i_mag=SDSS_catalogue['i']
    i_err=SDSS_catalogue['i_err']
    z_mag=SDSS_catalogue['z']
    z_err=SDSS_catalogue['z_err']
    
    # u_err=SDSS_DeCALS_catalogue['u_err']
    # g_err=SDSS_DeCALS_catalogue['g_err']
    # r_err=SDSS_DeCALS_catalogue['r_err']
    # i_err=SDSS_DeCALS_catalogue['i_err']
    # flux_i=SDSS_DeCALS_catalogue['FLUX_I']
    # ivar_i=SDSS_DeCALS_catalogue['FLUX_IVAR_I']
    # flux_g=SDSS_DeCALS_catalogue['FLUX_G']
    # ivar_g=SDSS_DeCALS_catalogue['FLUX_IVAR_G']
    # flux_r=SDSS_DeCALS_catalogue['FLUX_R']
    # ivar_r=SDSS_DeCALS_catalogue['FLUX_IVAR_R']
    


    radio_SDSS_catalogue = fits.open(r'Zwcl2341_DeCALS_radio.fits')
    radio_SDSS_catalogue=radio_SDSS_catalogue[1].data



    cross_match_galaxy_ra_radio=SDSS_catalogue['ra_1']
    cross_match_galaxy_dec_radio=SDSS_catalogue['dec_1']

    cross_match_galaxy_ra_SDSS=SDSS_catalogue['ra_2']
    cross_match_galaxy_dec_SDSS=SDSS_catalogue['dec_2']




    
    cluster_centre=SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))
    
    zsp_min=0.24
    zsp_max=0.30
    
    search_radius = 1.5

    R_200= 0.103 #5 Mpc

    mag1_filter=20
    
    
    mag2_filter=20

    galaxy_ra_circle,galaxy_dec_circle,redshift_circle=circle(galaxy_ra,galaxy_dec,zsp,cluster_centre,zsp_min,zsp_max,search_radius)

    # #zph_R_200,zsp_R_200=zwcl_galaxy_distribution(galaxy_ra_circle,galaxy_dec_circle,redshift_circle,cluster_centre,R_200,zsp_min,zsp_max)
   
    sigma_cluster_z,new_cluster_z,zph_R_200,zsp_R_200 =velocity_dispersion(galaxy_ra,galaxy_dec,zsp,zph,cluster_centre,R_200,zsp_min,zsp_max)
    
    cluster_members(galaxy_ra_circle,galaxy_dec_circle,redshift_circle,sigma_cluster_z,new_cluster_z)

    
    mag_1,mag_2,zsp_final,zph_final,target_region_filter=redshift_plots(zsp,zph,sigma_cluster_z,new_cluster_z,z_mag,r_mag,mag1_filter,mag2_filter)
    
    #radio_SFR_plots(radio_flux,solar_mass,zsp_1)
    
    #RA_DEC_seperation(cross_match_galaxy_ra_radio,cross_match_galaxy_dec_radio,cross_match_galaxy_ra_SDSS,cross_match_galaxy_dec_SDSS)
    #AGN_relation(radio_flux,zsp_1,solar_mass)
    
    #colour_colour(galaxy_ra,galaxy_dec,cluster_centre,R_200,zph,zsp,z_mag,g_mag,mag_filter)
    #SDSS_DeCALS(ivar_r,flux_r,r_err,band='r')
     
if __name__ == "__main__":
    main(sys.argv)


