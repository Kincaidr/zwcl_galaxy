#!/usr/bin/env python
from astroquery.sdss import SDSS
from astropy import coordinates as coords
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
import numpy as np
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
import astropy.cosmology as cp
from scipy.stats import norm

def zwcl_galaxy_distribution(galaxy_ra,galaxy_dec,zsp,cluster_centre,R_200):


    def histogram(sample,bins,color):
        n, bins, patches = plt.hist(sample, bins=bins, color=color, alpha=0.1, rwidth=0.85)
        cluster_z=bins[np.argmax(n)]
        return cluster_z

    # pos = coords.SkyCoord('23h43m39.700s +0d19m51.000s', frame='icrs')
    # xid = SDSS.query_region(pos, radius='1deg',spectro=True,data_release=17,photoobj_fields=['objid','ra','dec','u','g','r','i'],specobj_fields=['z','class'])
    #xid = SDSS.query_sql(query)

  
    #xid.colnames

  
    
    redshift=np.array(zsp)[(zsp >= 0.26) & (zsp <= 0.29)]
    galaxy_ra=np.array(galaxy_ra)[(zsp >= 0.26) & (zsp <= 0.29)]
    
    galaxy_dec=np.array(galaxy_dec)[(zsp >= 0.26) & (zsp <= 0.29)]
    

    cluster_z=histogram(redshift,500,'#0504aa')
    cluster_z_n=histogram(redshift,100,'#0504aa')
    cluster_z_median= np.median(redshift)

    print('cluster z 500 bin size:',cluster_z)
    print('cluster z 100 bin size:',cluster_z_n)
    print('cluster z median:',cluster_z_median)

    print('total number of galaxies in SDSS is',len(redshift))
    print('min redshift in SDSS is',min(redshift))
    print('max redshift in SDSS is',max(redshift))            
    print('RA and DEC of galaxies in redshift range:',len(galaxy_ra))

    
    ra_radian = np.array(galaxy_ra)*np.pi/180

    dec_radian = np.array(galaxy_dec)*np.pi/180


    h=0.7
    H0=h*100
    cosmo = cp.FlatLambdaCDM(H0=h*100, Om0=0.30)


    sc = SkyCoord(ra_radian, dec_radian, unit='rad', representation_type='unitspherical')
    cartesian=sc.cartesian


    x_coord=(np.cos(dec_radian) * np.cos(ra_radian))
    y_coord=(np.cos(dec_radian) * np.sin(ra_radian))
    z_coord=(np.sin(dec_radian))

 

    z_filter=np.array(redshift)

    cluster_centre_ra_radian=cluster_centre.ra.radian
    cluster_centre_dec_radian=cluster_centre.dec.radian

    cluster_centre=[cluster_centre.ra.value,cluster_centre.dec.value]


    sc_cluster_centre=  SkyCoord(cluster_centre_ra_radian, cluster_centre_dec_radian, unit='rad', representation_type='unitspherical')
    cartesian_cluster_centre=sc_cluster_centre.cartesian

    cluster_centre_x=cartesian_cluster_centre.x.value
    cluster_centre_y=cartesian_cluster_centre.y.value
    cluster_centre_z=cartesian_cluster_centre.z.value

    comoving_centre_x=cosmo.comoving_distance((cluster_z_median)*cluster_centre_x).value
    comoving_centre_y=cosmo.comoving_distance((cluster_z_median)*cluster_centre_y).value
    comoving_centre_z=cosmo.comoving_distance((cluster_z_median)*cluster_centre_z).value


    # reference to centre

    comoving_x=((cosmo.comoving_distance(z_filter)*x_coord).value)
    comoving_y=((cosmo.comoving_distance(z_filter)*y_coord).value)
    comoving_z =((cosmo.comoving_distance(z_filter)*z_coord).value)




    # R_200= 0.1388
        
    fig1 = plt.figure()
    fig2 = plt.figure()
    fig3 = plt.figure()
    ay1 = fig1.add_subplot(1,1,1, aspect='auto',projection='3d')
    ay2 = fig2.add_subplot(1,1,1, aspect='equal')
    ay3 = fig3.add_subplot(1,1,1, aspect='equal')

    # 2D physical coordinates

    ay1.scatter(comoving_x, comoving_y, comoving_z, color='black', alpha=0.5)
    ay1.set_xlabel("X in Mpc")
    ay1.set_ylabel("Y in Mpc")
    ay1.set_zlabel("Z in Mpc")
    ay1.set_title("3-D distribution of galaxies in a 1.5 deg radius centered on Zwcl 2341")

    ay2.scatter(comoving_y, comoving_z,color='black', alpha=0.5)
    ay2.set_xlabel("Y in Mpc")
    ay2.set_ylabel("Z in Mpc")


    circle=plt.Circle((comoving_centre_y, comoving_centre_z), R_200, edgecolor= 'blue',
    facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed') #, )
    ay2.add_patch(circle)
    ay2.text(0,0, "R_200", ha="left", va="top",fontsize=10)
    fig2.savefig("Y_Z_Mpc")
    print(comoving_centre_y,comoving_centre_z)
    R_200= 0.103

    # 2D angular coordinates

    ay3.scatter(galaxy_ra, galaxy_dec,color='black', alpha=0.5)

    ay3.set_xlim(min(galaxy_ra),max(galaxy_ra))
    ay3.set_ylim(min(galaxy_dec),max(galaxy_dec))
    ay3.set_xlabel("RA in deg")
    ay3.set_ylabel("DEC in deg")
    ay3.set_title("1.5 deg SDSS DR17 galaxies of Zwcl 2341 (Shishir catalogue)")

    circle=plt.Circle((cluster_centre[0], cluster_centre[1]), R_200, edgecolor= 'blue',
    facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed') #, )
    ay3.add_patch(circle)
    ay3.text(0,0, "R_200", ha="left", va="top",fontsize=10)
    fig3.savefig("Y_Z_Mpc")

    plt.show()
    return(galaxy_ra,galaxy_dec,redshift)


def velocity_dispersion(galaxy_ra,galaxy_dec,zsp,cluster_centre,R_200):
    
    
    redshift=np.array(zsp)[(zsp >= 0.26) & (zsp <= 0.29)]
    galaxy_ra=np.array(galaxy_ra)[(zsp >= 0.26) & (zsp <= 0.29)]
    
    galaxy_dec=np.array(galaxy_dec)[(zsp >= 0.26) & (zsp <= 0.29)]
    
    def histogram(sample,bins,color):
        n, bins, patches = plt.hist(sample, bins=bins, color=color, alpha=0.1, rwidth=0.85)
        cluster_z=bins[np.argmax(n)]
        return cluster_z
    
    cluster_centre_ra_radian=cluster_centre.ra.radian
    cluster_centre_dec_radian=cluster_centre.dec.radian

    cluster_centre=[cluster_centre.ra.value,cluster_centre.dec.value]
    
    ra_5Mpc_circle=[]
    dec_5Mpc_circle=[]
    z_5Mpc_circle=[]

    for i in range(len(galaxy_ra)):
        if np.sqrt((galaxy_ra[i]-cluster_centre[0])**2+(galaxy_dec[i]-cluster_centre[1])**2)<= R_200:
            ra_5Mpc_circle.append(galaxy_ra[i])
            dec_5Mpc_circle.append(galaxy_dec[i])
            z_5Mpc_circle.append(redshift[i])

    print('arcmin circle',len(z_5Mpc_circle))

    
    new_cluster_z=histogram(z_5Mpc_circle,50,'#0504ab')


    print('new cluster z:' ,new_cluster_z)  


    cluster_z_filter=[]

    for i in range(len(z_5Mpc_circle)):
        if (z_5Mpc_circle[i] >= new_cluster_z - 0.015) and (z_5Mpc_circle[i] <= new_cluster_z + 0.015):
            cluster_z_filter.append(z_5Mpc_circle[i]) 
    print(max(cluster_z_filter),min(cluster_z_filter),len(cluster_z_filter))

    #Cluster members
    
    v_rec_vel=np.array(np.array(cluster_z_filter)*constants.c)
    
        
    median_rec_vel=stats.median(v_rec_vel)
    
    v_rest_vel=np.array(v_rec_vel - median_rec_vel) /(1+median_rec_vel/constants.c)
            
            
    biscl_2 = biweight_scale(v_rest_vel)


    sigma_clip=np.array(v_rest_vel)[(v_rest_vel>-3*biscl_2)&(v_rest_vel<3*biscl_2)]

    print(' rest velocity calculated from the biweight scale',(biscl_2/10**3))

    sigma_cluster_z = biweight_scale(sigma_clip)/10**3

    
    print(' rest velocity calculated from the biweight scale sigma clip',sigma_cluster_z )

    plt.show() 
    
    return(sigma_cluster_z,new_cluster_z )
 



def redshift_plots(zsp,zph,sigma_cluster_z,new_cluster_z,r_mag,z_mag):


    zph_filter=np.array(zph)[(z_mag < 20) & (r_mag < 20) ]
    zsp_filter=np.array(zsp)[(z_mag < 20) & (r_mag < 20) ]
    
    mask_array=~np.logical_or(np.isnan(zsp_filter), np.isnan(zph_filter))


    zph_final=zph_filter[mask_array]
    zsp_final=zsp_filter[mask_array]
    print(len(zph_final),'length spectro and photo array')
    sigma_z= (sigma_cluster_z*10**3) /constants.c


    plt.scatter(zph_final,zsp_final)
    plt.plot(zsp_final, zsp_final, 'k-', lw=2,label='one-to-one')
    plt.axhline(y=0.25, color='r', linestyle='-')
    plt.axhline(y=0.30, color='r', linestyle='-')
    plt.axvline(x=0.20, color='r', linestyle='-')
    plt.axvline(x=0.40, color='r', linestyle='-')
    plt.title('Photometric redshift vs spectrosopic redshift for galaxies at 1 > z > 0, r_mag,z_mag < 20 ')
    plt.xlabel('photometric redshift')
    plt.ylabel('spectroscopic redshift')
    plt.xlim(0,1)
    plt.ylim(0,1)
    # plt.xlim(new_cluster_z-(30*(sigma_z)),new_cluster_z+(30*(sigma_z)))
    # plt.ylim(new_cluster_z-(30*(sigma_z)),new_cluster_z+(30*(sigma_z)))
    plt.show()


def radio_SFR_plots(radio_flux,solar_mass,zsp_1):
    
    h=0.7
    H0=h*100
    cosmo = cp.FlatLambdaCDM(H0=h*100, Om0=0.30)

    lumo_distance=np.array(cosmo.luminosity_distance(zsp_1))

    SFR_constant = 3.18e-22


    S_14=(radio_flux/(1283e6**-0.8)) * (1400e6**-0.8)


    radio_lumo=4*np.pi*(S_14*10**-26)*(lumo_distance*3.086e22)**2

    print('radio lumonisties',radio_lumo)

    SFR = SFR_constant*radio_lumo


    print('SFR',SFR)

    plt.scatter( solar_mass,SFR, color='black', alpha=0.5)
    plt.xlabel("SFR in M_solar/yea r")
    plt.ylabel("Log(Solar mass)")
    plt.ylim(0,800)
    plt.title("SFR vs solar mass")
    plt.show()


def RA_DEC_seperation(cross_match_galaxy_ra_radio,cross_match_galaxy_dec_radio,cross_match_galaxy_ra_SDSS,cross_match_galaxy_dec_SDSS):
    
    DEC_sep=(cross_match_galaxy_dec_SDSS-cross_match_galaxy_dec_radio)*3600

    RA_sep=(cross_match_galaxy_ra_SDSS-cross_match_galaxy_ra_radio)*3600

    n, bins = np.histogram(RA_sep, bins=50)

    bin_DEC_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2

    n, bins = np.histogram(DEC_sep, bins=50)

    bin_RA_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2


    x = np.linspace(-10, 10, num=100)

    fig = plt.figure(figsize=(13, 5))
    plt.subplot(121)
    plt.hist(RA_sep-bin_RA_centre, bins='auto', density=True)
    plt.plot(x, norm.pdf(x, np.mean(RA_sep-bin_RA_centre), np.std(RA_sep-bin_RA_centre)))
    plt.xlabel('RA sep (arcsec)')
    plt.title("RA seperation distance distribution between radio and optical")
    plt.subplot(122)
    
    plt.hist(DEC_sep-bin_DEC_centre, bins='auto', density=True)
    plt.plot(x, norm.pdf(x, np.mean(DEC_sep-bin_DEC_centre), np.std(DEC_sep-bin_DEC_centre)))
    plt.xlabel('DEC sep (arcsec)')
    plt.title("DEC seperation distance distribution between radio and optical")

    plt.show()