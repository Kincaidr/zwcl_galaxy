
from optparse import OptionParser
import sys
from zwcl2341_plots import zwcl_galaxy_distribution,velocity_dispersion,redshift_plots,radio_SFR_plots,RA_DEC_seperation,AGN_relation,SDSS_DeCALS,colour_colour
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import csv
import numpy as np

def main( argv ):
    
    SDSS_catalogue = pd.read_csv(r'Zwcl2341_Spec_Photo.csv')
    SDSS_DeCALS_catalogue=pd.read_csv(r'SDSS_Spec_Photo_DeCALS.csv')
    galaxy_ra=SDSS_catalogue['ra']
    galaxy_dec=SDSS_catalogue['dec']
    zsp=SDSS_catalogue['zsp']
    zph=SDSS_catalogue['zph']
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
    
    u_err=SDSS_DeCALS_catalogue['u_err']
    g_err=SDSS_DeCALS_catalogue['g_err']
    r_err=SDSS_DeCALS_catalogue['r_err']
    i_err=SDSS_DeCALS_catalogue['i_err']
    flux_i=SDSS_DeCALS_catalogue['FLUX_I']
    ivar_i=SDSS_DeCALS_catalogue['FLUX_IVAR_I']
    flux_g=SDSS_DeCALS_catalogue['FLUX_G']
    ivar_g=SDSS_DeCALS_catalogue['FLUX_IVAR_G']
    flux_r=SDSS_DeCALS_catalogue['FLUX_R']
    ivar_r=SDSS_DeCALS_catalogue['FLUX_IVAR_R']
    
    radio_SDSS_catalogue = pd.read_csv(r'Zwcl2341_Spec_Photo_Radio.csv')

    cross_match_galaxy_ra_radio=radio_SDSS_catalogue['ra_1']
    cross_match_galaxy_dec_radio=radio_SDSS_catalogue['dec_1']
    cross_match_galaxy_ra_SDSS=radio_SDSS_catalogue['RA_2']
    cross_match_galaxy_dec_SDSS=radio_SDSS_catalogue['DEC_2']

    radio_SDSS_SFR_catalogue = pd.read_csv(r'Zwcl2341_Spec_Photo_SFR_Radio.csv')

    radio_flux=radio_SDSS_SFR_catalogue['Total_flux']
    solar_mass=radio_SDSS_SFR_catalogue['logMass']
    zsp_1=radio_SDSS_SFR_catalogue['zsp']
    

    
    cluster_centre=SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))
    
    zsp_min=0.26
    zsp_max=0.30
    
    R_200= 0.103 #5 Mpc

    mag_filter=20

    zph_R_200,zsp_R_200=zwcl_galaxy_distribution(galaxy_ra,galaxy_dec,zsp,cluster_centre,R_200,zsp_min,zsp_max)
   
    #zph_R_200,zsp_R_200,sigma_cluster_z,new_cluster_z=velocity_dispersion(galaxy_ra,galaxy_dec,zsp,zph,cluster_centre,R_200,zsp_min,zsp_max)
    
    #mag_1,mag_2,zsp_final,zph_final,target_region_filter=redshift_plots(zsp,zph,sigma_cluster_z,new_cluster_z,z_mag,r_mag,mag_filter)
    
    #radio_SFR_plots(radio_flux,solar_mass,zsp_1)
    
    #RA_DEC_seperation(cross_match_galaxy_ra_radio,cross_match_galaxy_dec_radio,cross_match_galaxy_ra_SDSS,cross_match_galaxy_dec_SDSS)
    #AGN_relation(radio_flux,zsp_1,solar_mass)
    
    #colour_colour(galaxy_ra,galaxy_dec,cluster_centre,R_200,zph,zsp,z_mag,g_mag,mag_filter)
    #SDSS_DeCALS(ivar_r,flux_r,r_err,band='r')
     
if __name__ == "__main__":
    main(sys.argv)


