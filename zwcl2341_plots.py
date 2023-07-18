#!/usr/bin/env python
#from astroquery.sdss import SDSS
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
from scipy.stats import norm
from scipy.optimize import curve_fit
cosmo = FlatLambdaCDM(H0=70, Om0=0.27)

def circle(cluster_name,galaxy_ra,galaxy_dec,zsp,cluster_centre,zsp_min,zsp_max,search_radius):

    galaxy_ra_circle=[]
    galaxy_dec_circle=[]
    redshift=[]

    #cluster 90arcmin radius


    for i in range(len(galaxy_ra)):
        if np.sqrt((galaxy_ra[i]-cluster_centre.ra.value)**2+(galaxy_dec[i]-cluster_centre.dec.value)**2)<= search_radius:
            galaxy_ra_circle.append(galaxy_ra[i])
            galaxy_dec_circle.append(galaxy_dec[i])
            redshift.append(zsp[i])
    redshift_circle=np.array(redshift)

    
    galaxy_ra_new=np.array(galaxy_ra_circle)[(np.array(redshift_circle) >= zsp_min) & (np.array(redshift_circle) <= zsp_max)]

    galaxy_dec_new =np.array(galaxy_dec_circle)[(np.array(redshift_circle) >= zsp_min) & (np.array(redshift_circle) <= zsp_max)]

    
    plt.scatter(galaxy_ra_new,galaxy_dec_new)

    plt.xlabel("RA in deg")
    plt.ylabel("DEC in deg")
    plt.title("1.5 deg SDSS DR17 galaxies of " + cluster_name)

    

    return(galaxy_ra_circle,galaxy_dec_circle,redshift_circle)


def redshift_distribution(cluster_name,galaxy_ra_circle,galaxy_dec_circle,redshift_circle,cluster_centre,R_200,zsp_min,zsp_max):

    zsp_min_1=0
    zsp_max_1=1

    redshift_1=np.array(redshift_circle)[(np.array(redshift_circle) >= zsp_min) & (np.array(redshift_circle)<= zsp_max)]
    redshift_2=np.array(redshift_circle)[(np.array(redshift_circle) >= zsp_min_1) & (np.array(redshift_circle)<= zsp_max_1)]

    fig = plt.figure(figsize=(5, 3))
    
    ax1 = fig.add_axes([0.1,0.1,0.9,0.9])
    ax2 = fig.add_axes([0.70, 0.72,0.25,0.25])

    ax1.hist(redshift_1,bins=100)
    ax2.hist(redshift_2,bins=50)
   
    ax1.set_xlabel('Redshift',fontsize=15)
    ax1.set_ylabel('Number of galaxies',fontsize=15)
    ax1.set_title("SDSS spectroscopic redshift distribution of galaxies within 1.5 deg of " + cluster_name)
    plt.savefig(cluster_name +  '_redshift_distribution.png')
   

def velocity_dispersion(galaxy_ra,galaxy_dec,zsp,zph,cluster_centre,R_200,zsp_min,zsp_max):
    

    redshift=np.array(zsp)
    
    def histogram(sample,bins,color):
        n, bins, patches = plt.hist(sample, bins=bins, color=color, alpha=0.1, rwidth=0.85)
        cluster_z=bins[np.argmax(n)]
        return cluster_z
    
    cluster_centre_ra_radian=cluster_centre.ra.radian
    cluster_centre_dec_radian=cluster_centre.dec.radian

    cluster_centre=[cluster_centre.ra.value,cluster_centre.dec.value]
    
    ra_5Mpc_circle=[]
    dec_5Mpc_circle=[]
    zsp_R_200=[]
    zph_R_200=[]
    for i in range(len(galaxy_ra)):
        if np.sqrt((galaxy_ra[i]-cluster_centre[0])**2+(galaxy_dec[i]-cluster_centre[1])**2)<= R_200:
            ra_5Mpc_circle.append(galaxy_ra[i])
            dec_5Mpc_circle.append(galaxy_dec[i])
            zsp_R_200.append(redshift[i])
            zph_R_200.append(zph[i])
        

    print('arcmin circle',zsp_R_200)
    

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

    return(sigma_cluster_z,new_cluster_z,zph_R_200,zsp_R_200 )
 

def  cluster_members(cluster_name,galaxy_ra_circle,galaxy_dec_circle,redshift_circle,sigma_cluster_z,new_cluster_z):
    
    cluster_centre=SkyCoord(str(354.41916666667 ), str(0.27666666666667), frame='icrs',unit=(u.deg,u.deg))

    #cluster_member=np.array(redshift_circle)[(redshift_circle>=new_cluster_z-3*(sigma_cluster_z))&(redshift_circle<=new_cluster_z+3*(sigma_cluster_z))]  
        
    cluster_member_3=((redshift_circle>=new_cluster_z-3*(sigma_cluster_z))&(redshift_circle<=new_cluster_z+3*(sigma_cluster_z)) )

    cluster_member_5=((redshift_circle>=new_cluster_z-5*(sigma_cluster_z))&(redshift_circle<=new_cluster_z+5*(sigma_cluster_z)) )

    # foreground= ((redshift_circle<new_cluster_z-3*(sigma_cluster_z))& (redshift_circle != cluster_member))
    
    # background= ((redshift_circle>new_cluster_z+3*(sigma_cluster_z))& (redshift_circle != cluster_member))


                                       
    print(' Number of cluster members',len(cluster_member_3) )



    galaxy_ra_members_3=np.array(galaxy_ra_circle)[cluster_member_3]
    galaxy_dec_members_3=np.array(galaxy_dec_circle)[cluster_member_3]


    galaxy_ra_members_5=np.array(galaxy_ra_circle)[cluster_member_5]
    galaxy_dec_members_5=np.array(galaxy_dec_circle)[cluster_member_5]

    fig, ax = plt.subplots()
    plt.scatter(galaxy_ra_members_3,galaxy_dec_members_3,color='blue',alpha=0.5,marker="x",label=r'Member galaxies at $-3 \sigma < z_{cluster} <  3 \sigma$')

  
    plt.legend()

    plt.xlabel("RA in deg")
    plt.ylabel("DEC in deg")
    plt.title("1.5 deg SDSS DR17 galaxies of " + cluster_name)
    circle=plt.Circle((cluster_centre.ra.value, cluster_centre.dec.value), R_200, edgecolor= 'red',facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed') #, )
    plt.gca().add_patch(circle)

    ax.add_patch(circle)
    ax.text(0,0, "R_200", ha="left", va="top",fontsize=10)
    plt.savefig('A2631_cluster_members.png')
    plt.show()




    return(sigma_cluster_z,new_cluster_z)


def zwcl_galaxy_distribution(galaxy_ra_new,galaxy_dec_new,redshift_circle,cluster_centre,R_200,zsp_min,zsp_max):

    h=0.7
    H0=h*100
    cosmo = cp.FlatLambdaCDM(H0=h*100, Om0=0.30)

    RA_DEC=SkyCoord(galaxy_ra_new,galaxy_dec_new,frame='icrs',unit=(u.deg,u.deg))

    ra_radian=RA_DEC.ra.radian
    dec_radian=RA_DEC.dec.radian

    sc = SkyCoord(ra_radian, dec_radian, unit='rad', representation_type='unitspherical')
    cartesian=sc.cartesian


    x_coord=(np.cos(dec_radian) * np.cos(ra_radian))
    y_coord=(np.cos(dec_radian) * np.sin(ra_radian))
    z_coord=(np.sin(dec_radian))


    z_filter=np.array(redshift_circle)

    cluster_centre_ra_radian=cluster_centre.ra.radian
    cluster_centre_dec_radian=cluster_centre.dec.radian


    sc_cluster_centre=  SkyCoord(cluster_centre_ra_radian, cluster_centre_dec_radian, unit='rad', representation_type='unitspherical')
    cartesian_cluster_centre=sc_cluster_centre.cartesian

    cluster_centre_x=cartesian_cluster_centre.x.value
    cluster_centre_y=cartesian_cluster_centre.y.value
    cluster_centre_z=cartesian_cluster_centre.z.value



    fig1 = plt.figure()
    fig3 = plt.figure()
    
    ay1 = fig1.add_subplot(1,1,1, aspect='auto',projection='3d')
    
    ay3 = fig3.add_subplot(1,1,1, aspect='equal')

    # 2D physical coordinates

    comoving_x=((cosmo.comoving_distance(z_filter)*x_coord).value)
    comoving_y=((cosmo.comoving_distance(z_filter)*y_coord).value)
    comoving_z =((cosmo.comoving_distance(z_filter)*z_coord).value)

    omoving_centre_x=cosmo.comoving_distance((cluster_z_median_2)*cluster_centre_x).value
    comoving_centre_y=cosmo.comoving_distance((cluster_z_median_2)*cluster_centre_y).value
    comoving_centre_z=cosmo.comoving_distance((cluster_z_median_2)*cluster_centre_z).value
    
    ay1.scatter(comoving_x, comoving_y, comoving_z, color='black', alpha=0.5)
    ay1.set_xlabel("X in Mpc")
    ay1.set_ylabel("Y in Mpc")
    ay1.set_zlabel("Z in Mpc")
    ay1.set_title("3-D distribution of galaxies in a 1.5 deg radius centered on Zwcl 2341")
   
    # 2D angular coordinates

    #ay3.scatter(galaxy_ra_foreground, galaxy_dec_foreground,color='red', alpha=0.5,label='foreground galaxies')
    #ay3.scatter(galaxy_ra_background, galaxy_dec_background,color='black', alpha=0.5,label='background galaxies')
    ay3.scatter(galaxy_ra_members,galaxy_dec_members,color='blue',alpha=0.5,marker="x",label='A2631 cluster members')
    
    ay3.set_xlim(min(galaxy_ra_circle),max(galaxy_ra_circle))
    ay3.set_ylim(min(galaxy_dec_circle),max(galaxy_dec_circle))
    ay3.set_xlabel("RA in deg")
    ay3.set_ylabel("DEC in deg")
    ay3.set_title("1.5 deg SDSS DR17 galaxies of Zwcl 2341")
    ay3.legend()
    circle=plt.Circle((cluster_centre.ra.value, cluster_centre.dec.value), R_200, edgecolor= 'red',
    facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed') #, )
    ay3.add_patch(circle)
    ay3.text(0,0, "R_200", ha="left", va="top",fontsize=10)

    plt.show()
    
    return(sigma_cluster_z,new_cluster_z)




def redshift_plots(cluster_name,zsp,zph,sigma_cluster_z,new_cluster_z,mag_1,mag_2,mag1_filter,mag2_filter):


    zph_filter=np.array(zph)[(mag_1 < mag1_filter) & (mag_2 < mag2_filter) ]
    zsp_filter=np.array(zsp)[(mag_1 < mag1_filter) & (mag_2 < mag2_filter) ]

    zph_nofilter=np.array(zph)
    zsp_nofilter=np.array(zsp)

    mask_array_1=~np.logical_or(np.isnan(zsp_filter), np.isnan(zph_filter))
    mask_array_2=~np.logical_or(np.isnan(zsp_nofilter), np.isnan(zph_nofilter))


    zph_final=zph_filter[mask_array_1]
    zsp_final=zsp_filter[mask_array_1]

    zph_final_nofilter=zph_nofilter[mask_array_2]
    zsp_final_nofilter=zsp_nofilter[mask_array_2]
    
    target_region_filter = zsp_final[np.where((zsp_final>0.25)&(zsp_final<0.3)&(zph_final>0.2)&(zph_final<0.4))]
    
    target_region_nofilter = zsp_final_nofilter[np.where((zsp_final_nofilter>0.25)&(zsp_final_nofilter<0.3)&(zph_final_nofilter>0.2)&(zph_final_nofilter<0.4))]
    
    zsp_minus_zph_before = np.array(zsp_final_nofilter-zph_final_nofilter)[np.where((zph_final_nofilter>0.2)&(zph_final_nofilter<0.4))]
    
    zsp_minus_zph_after = np.array(zsp_final-zph_final)[np.where((zph_final>0.2)&(zph_final<0.4))]
    

    FD=(100-(len(target_region_nofilter)-len(target_region_filter))/(len(target_region_nofilter))*100)
    
    
    perc_outlier_before=np.absolute(np.mean(zsp_minus_zph_before))
    perc_outlier_after=np.absolute(np.mean(zsp_minus_zph_after))
    
    fraction_outlier=(perc_outlier_after/perc_outlier_before)*100
    #(zsp_minus_zph<0.05)&(zsp_minus_zph>-0.05)
    # mag_1=f'{mag_1}'.split('=')[0]
    # mag_2=f'{mag_2}'.split('=')[0]
   
    print(len(zph_final),'length spectro and photo array')
    

    fig = plt.figure(figsize=(15, 7))

    # plt.xlim(new_cluster_z-(30*(sigma_z)),new_cluster_z+(30*(sigma_z)))
    # plt.ylim(new_cluster_z-(30*(sigma_z)),new_cluster_z+(30*(sigma_z)))

    plt.subplot(221)
    plt.scatter(zph_final_nofilter,zsp_final_nofilter)
    plt.plot(zsp_final_nofilter, zsp_final_nofilter, 'k-', lw=2,label='one-to-one')
    plt.axhline(y=new_cluster_z-5*sigma_cluster_z, color='r', linestyle='-')
    plt.axhline(y=new_cluster_z+5*sigma_cluster_z, color='r', linestyle='-')
    plt.axvline(x=0.25, color='r', linestyle='-')
    plt.axvline(x=0.34, color='r', linestyle='-')
    plt.title(cluster_name + " photo z vs spec z for galaxies at 1 > z > 0")
    #plt.title('Outlier ratio is {}'.format(mag_1))
    plt.xlabel('photometric redshift')
    plt.ylabel('spectroscopic redshift')
    plt.xlim(0,1)
    plt.ylim(0,1)


    plt.subplot(222)
    plt.scatter(zph_final_nofilter,zsp_final_nofilter)
    plt.plot(zsp_final_nofilter, zsp_final_nofilter, 'k-', lw=2,label='one-to-one')
    plt.axhline(y=new_cluster_z-5*sigma_cluster_z, color='r', linestyle='-')
    plt.axhline(y=new_cluster_z+5*sigma_cluster_z, color='r', linestyle='-')
    plt.axvline(x=0.25, color='r', linestyle='-')
    plt.axvline(x=0.34, color='r', linestyle='-')
    plt.title('Zoom in of photo z vs spec z for galaxies at 1 > z > 0 showing the redshift sample of interest' )
    #plt.title('Outlier ratio is {}'.format(mag_1))
    plt.xlabel('photometric redshift')
    plt.ylabel('spectroscopic redshift')
    plt.xlim(0.1,0.5)
    plt.ylim(0.2,0.35)

    plt.subplot(223)
    plt.scatter(zph_final,zsp_final)
    plt.plot(zsp_final, zsp_final, 'k-', lw=2,label='one-to-one')
    plt.axhline(y=new_cluster_z-5*sigma_cluster_z, color='r', linestyle='-')
    plt.axhline(y=new_cluster_z+5*sigma_cluster_z, color='r', linestyle='-')
    plt.axvline(x=0.25, color='r', linestyle='-')
    plt.axvline(x=0.34, color='r', linestyle='-')
    plt.title(cluster_name + "photo z vs spec z for galaxies at 1 > z > 0 after magnitude filter")
    #plt.title('Outlier ratio is {}'.format(mag_1))
    plt.xlabel('photometric redshift')
    plt.ylabel('spectroscopic redshift')
    plt.xlim(0,1)
    plt.ylim(0,1)


    plt.subplot(224)
    plt.scatter(zph_final,zsp_final)
    plt.plot(zsp_final, zsp_final, 'k-', lw=2,label='one-to-one')
    plt.axhline(y=new_cluster_z-5*sigma_cluster_z, color='r', linestyle='-')
    plt.axhline(y=new_cluster_z+5*sigma_cluster_z, color='r', linestyle='-')
    plt.axvline(x=0.25, color='r', linestyle='-')
    plt.axvline(x=0.34, color='r', linestyle='-')
    plt.title('Zoom in of photo z vs spec z for galaxies at 1 > z > 0 showing the redshift sample of interest after magnitude filter')
    #plt.title('Outlier ratio is {}'.format(mag_1))
    plt.xlabel('photometric redshift')
    plt.ylabel('spectroscopic redshift')
    plt.xlim(0.1,0.5)
    plt.ylim(0.2,0.35)

    plt.show()


    

    fig = plt.figure(figsize=(15, 7))
    plt.subplot(121)
    x = np.linspace(min(zsp_minus_zph_before), max(zsp_minus_zph_before), num=100)

    plt.hist(zsp_minus_zph_before, bins=50)
    plt.plot(x, norm.pdf(x, np.mean(zsp_minus_zph_before), np.std(zsp_minus_zph_before)))
    plt.title('Distribution of (photo z - spec z) for galaxies at 1 > z > 0')
    plt.vlines(np.median(zsp_minus_zph_before), 0, 1)
    plt.subplot(122)
    x = np.linspace(min(zsp_minus_zph_after), max(zsp_minus_zph_after), num=100)
    
    #plt.text(.1, .99, "fraction_outlier = {:.3f}%".format(fraction_outlier), ha='center', va='baseline')
    plt.plot([], [], ' ', label="fraction_outlier = {:.3f}%".format(fraction_outlier))
    plt.hist(zsp_minus_zph_after, bins=50)
    plt.plot(x, norm.pdf(x, np.mean(zsp_minus_zph_after), np.std(zsp_minus_zph_after)))
    plt.title(f"Distribution of  (photo z - spec z) for galaxies at 1 > z > 0 with + 'mag_1 + < {mag1_filter}, + mag_2 + < {mag2_filter}")
    plt.vlines(np.median(zsp_minus_zph_after), 0, 1)
    plt.legend()
    plt.show()
    
    return(zsp_final,zph_final,mag_1,mag_2)
    
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
    plt.ylabel(r'SFR ($\mathcal{M}_\odot / year)$')
    plt.xlabel(r'Log($\mathcal{M}_\odot)$')
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
    
    bin_RA_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2


    x = np.linspace(-3, 3, num=100)

    fig = plt.figure(figsize=(13, 5))
    plt.subplot(121)
    plt.hist(RA_sep-bin_RA_centre, bins='auto', density=True)
    plt.plot(x, norm.pdf(x, np.mean(RA_sep-bin_RA_centre), np.std(RA_sep-bin_RA_centre)))
    plt.xlabel('RA sep (arcsec)')
    plt.title("RA seperation distance distribution between SDSS and DeCALS for Zwcl2341")
    plt.subplot(122)
    
    plt.hist(DEC_sep-bin_DEC_centre, bins='auto', density=True)
    plt.plot(x, norm.pdf(x, np.mean(DEC_sep-bin_DEC_centre), np.std(DEC_sep-bin_DEC_centre)))
    plt.xlabel('DEC sep (arcsec)')
    plt.title("DEC seperation distance distribution between SDSS and DeCALS for Zwcl2341")

    plt.show()

def AGN_relation(radio_flux,zsp_1,solar_mass):
    
    def func2(x,m,c):
        lsq=(m*x) + c
        return lsq

    h=0.7
    H0=h*100
    cosmo = cp.FlatLambdaCDM(H0=h*100, Om0=0.30)

    lumo_distance=np.array(cosmo.luminosity_distance(zsp_1))

    SFR_constant = 3.18e-22


    S_14=(radio_flux/(1283e6**-0.8)) * (1400e6**-0.8)


    radio_lumo=np.log10(4*np.pi*(S_14*10**-26)*(lumo_distance*3.086e22)**2)
    

    #AGN relation
    
    mass = np.array([9.5,10.25,10.75,11.5])
    lumo1 = np.mean(np.array([22.59,22.93,23.28]))
    lumo2 = np.mean(np.array([23,23.36]))
    lumo3 = np.mean(np.array([22.99,23.43,23.87,24.19]))
    lumo4 = np.mean(np.array([23.40,23.82,24.17]))


    lumo=np.array([lumo1,lumo2,lumo3,lumo4])

    
    coeff2, var2 = curve_fit(func2,mass,lumo)
    
    
    yfit = func2(solar_mass,coeff2[0],coeff2[1])

    plt.plot(solar_mass,yfit,'k-',linewidth=2.0,label='AGN trend for 0.1 < z < 0.7')

    plt.scatter( solar_mass,radio_lumo, color='blue', alpha=0.5,label='Zwcl galaxies at zsp > 1 ')
    
    #plt.scatter( mass,lumo, color='black', alpha=0.5)
    
    #plt.ylabel("Log (L_1.4 GHz)")
    
    plt.ylabel(r'Log($\mathcal{L}^{AGN}_{1.4GHz}[W / Hz^{-1}]$')
    plt.xlabel(r'Log($\mathcal{M} / \mathcal{M}_\odot)$')
    plt.xlim(np.min(solar_mass),np.max(solar_mass))
    plt.title("Low redshift AGN correlation overplotted over ZwCl galaxies")
    plt.legend()
    plt.show()
    
    
def colour_colour(galaxy_ra,galaxy_dec,cluster_centre,R_200,zph,zsp,mag_1,mag_2,mag_filter):
    
    galaxy_ra=np.array(galaxy_ra)
    galaxy_dec=np.array(galaxy_dec)
    zph=np.array(zph)
    zsp=np.array(zsp)
    
    
    ra_5Mpc_circle=[]
    dec_5Mpc_circle=[]
    zsp_R_200=[]
    zph_R_200=[]
    for i in range(len(galaxy_ra)):
        if np.sqrt((galaxy_ra[i]-cluster_centre.ra.value)**2+(galaxy_dec[i]-cluster_centre.dec.value)**2)<= R_200:
            zph_R_200.append(zph[i])
            zsp_R_200.append(zsp[i])
    
    mask_array_1=~np.logical_or(np.isnan(zsp_R_200), np.isnan(zph_R_200))
    
   
    zsp_R_200=np.array(zsp_R_200)[mask_array_1]
    
    zph_R_200=np.array(zph_R_200)
    
    mag_2=np.array(mag_2)[mag_2 < mag_filter]
    
    mag_1=np.array(mag_1)[mag_1 < mag_filter]
    
    colour_1_zph=np.array(mag_1)[np.where((zph_R_200 > 0.2) & (zph_R_200 < 0.40))]
    
    colour_2_zph=np.array(mag_2)[np.where((zph_R_200 > 0.2) & (zph_R_200 < 0.40))]
    
    colour_1_zsp=np.array(mag_1)[np.where((zsp_R_200 > 0.26) & (zsp_R_200 < 0.30))]
    
    colour_2_zsp=np.array(mag_2)[np.where((zsp_R_200 > 0.26) & (zsp_R_200 < 0.30))]
    
    
    
    if len(colour_1_zph) > len(colour_2_zph):
        colour_1_zph=np.array(colour_1_zph)[0:len(colour_2_zph)]
    else:
        colour_2_ph=np.array(colour_2_zph)[0:len(colour_1_zph)]
        
    

    if len(colour_1_zsp) > len(colour_2_zsp):
        colour_1_zsp=np.array(colour_1_zsp)[0:len(colour_2_zsp)]
    else:
        colour_2_zsp=np.array(colour_2_zsp)[0:len(colour_1_zsp)]
  

    
    plt.scatter( colour_1_zsp,np.abs(colour_1_zsp-colour_2_zsp), color='blue', alpha=0.5,label='Spectra from 0.26 <z < 0.30')
    
    plt.scatter( colour_1_zph,np.abs(colour_1_zph-colour_2_zph), color='black', alpha=0.3,label='Photometry from 0.20 <z < 0.40',marker="x")
    
    
    plt.ylabel(r'$mag_1-mag_2$')
    plt.title(r"Colour-Colour plot of ZwCl at $R_{200} \sim 5$ Mpc filtered by $~mag_1,mag_2 < 20$")
    plt.legend()
    plt.xlabel(r'$mag_1$')
    plt.legend()
    plt.show()
    
def SDSS_DeCALS(ivar,flux_1,mag1_err,band):
                
        def ivar2errmag(ivar, flux_1):
            df = np.sqrt(1/ivar) # ivar = 1 / (df^2)
            # dm  = (2.5 / ln(10)) * (df / f)
            return (2.5 / np.log(10))  * (df / flux_1)
        
        mag_err_decal=ivar2errmag(ivar,flux_1)
        
        SNR_decal=(1/mag_err_decal)*(2.5/np.log(10))
        SNR_SDSS=(1/mag1_err)*(2.5/np.log(10))
        
        plt.scatter(SNR_decal,SNR_SDSS)
        plt.plot(1,1)
        plt.xlabel(r'$SNR_{DeCALS}$')
        plt.ylabel(r'$SNR_{SDSS}$')
        plt.savefig("SNR SDSS vs DeCALS for band" + ' ' + str(band))
        plt.title("SNR SDSS vs DeCALS for band" + ' ' + str(band))
            
        plt.show()
        
        

def R_200_to_deg(distance_Mpc,redshift):
            
    def radians_to_degrees(radians):
        degrees = radians * (180 / math.pi)
        return degrees

    d_A = cosmo.angular_diameter_distance(z=redshift)

    theta_radian = distance_Mpc/d_A

    R_200=radians_to_degrees(theta_radian.value)

    return(R_200)
        
