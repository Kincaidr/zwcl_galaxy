
from optparse import OptionParser
import sys
from zwcl2341_plots import zwcl_galaxy_distribution,velocity_dispersion,redshift_plots,radio_SFR_plots,RA_DEC_seperation,AGN_relation,colour_colour
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import csv
import numpy as np

def main( argv ):
    
    SDSS_catalogue = pd.read_csv(r'zwcl2341DR17UncleanViral_kincaid_1.csv')
   
    radio_SDSS_catalogue = pd.read_csv(r'SDSS_radio_combined.csv')
    
    galaxy_ra=SDSS_catalogue['ra']
    galaxy_dec=SDSS_catalogue['dec']
    cross_match_galaxy_ra_radio=radio_SDSS_catalogue['ra_1']
    cross_match_galaxy_dec_radio=radio_SDSS_catalogue['dec_1']
    cross_match_galaxy_ra_SDSS=radio_SDSS_catalogue['RA_2']
    cross_match_galaxy_dec_SDSS=radio_SDSS_catalogue['DEC_2']
    zsp=SDSS_catalogue['zsp']
    zph=SDSS_catalogue['zph']
    u_mag=SDSS_catalogue['u']
    g_mag=SDSS_catalogue['g']
    r_mag=SDSS_catalogue['r']
    i_mag=SDSS_catalogue['i']
    z_mag=SDSS_catalogue['z']

    
    
    radio_SDSS_SFR_catalogue = pd.read_csv(r'solar_mass.csv')

    radio_flux=radio_SDSS_SFR_catalogue['Total_flux']
    solar_mass=radio_SDSS_SFR_catalogue['logMass']
    zsp_1=radio_SDSS_SFR_catalogue['zsp']
    
    cluster_centre= SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))
    
    R_200= 0.103 #5 Mpc

    mag_filter=20

    #zwcl_galaxy_distribution(galaxy_ra,galaxy_dec,zsp,cluster_centre,R_200)
   
    #zph_R_200,zsp_R_200,sigma_cluster_z,new_cluster_z=velocity_dispersion(galaxy_ra,galaxy_dec,zsp,zph,cluster_centre,R_200)
    
    #mag_1,mag_2,zsp_final,zph_final,target_region_filter=redshift_plots(zsp,zph,sigma_cluster_z,new_cluster_z,z_mag,r_mag,magnitude_filter)
    
    #radio_SFR_plots(radio_flux,solar_mass,zsp_1)
    
    #RA_DEC_seperation(cross_match_galaxy_ra_radio,cross_match_galaxy_dec_radio,cross_match_galaxy_ra_SDSS,cross_match_galaxy_dec_SDSS)
    AGN_relation(radio_flux,zsp_1,solar_mass)
    
    #colour_colour(galaxy_ra,galaxy_dec,cluster_centre,R_200,zph,zsp,z_mag,r_mag,mag_filter)
    
     
if __name__ == "__main__":
    main(sys.argv)