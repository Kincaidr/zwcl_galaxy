import config
import socket
import os, sys
import glob
import requests
import matplotlib.pyplot as plt
import matplotlib as mpl
import RA_DEC_names as radecname
import pandas as pd
from astropy.io import fits
import aplpy
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
mpl.use('agg')
from skimage.transform import resize
from io import BytesIO
from PIL import Image
import panstarrs_utils as ps
from astroquery.gaia import GaiaClass
from astropy.coordinates import SkyCoord
import numpy as np
from gaiaxpy import calibrate

from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy import units as u
from scipy.ndimage import gaussian_filter

from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = -1


#if 'lesta' in socket.gethostname():
#    connhandler = tapconn.ConnectionHandler(host='proxy.unige.ch', port=3128, sslport=3128)
#    tap_plus_conn_handler.set_tunnel('gea.esac.esa.int')
#    Gaia  = GaiaClass(tap_plus_conn_handler=connhandler)

#else:
Gaia  = GaiaClass() #tap_plus_conn_handler=connhandler)


#u, p = clemon, !m73D!h5XE+V)mb9
#Gaia.login(user=config.username, password=config.password)


#table = XMatch.query(cat1=open('/tmp/pos_list.csv'),
#                     cat2='vizier:II/246/out',
#                     max_distance=1.5 * u.arcsec)

def get_main_table(ra, dec, radius):
    radius *= 1/3600
    query = f"SELECT * \
    FROM gaiadr3.gaia_source \
    WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS',"+str(ra)+","+str(dec)+","+str(radius)+"))=1;"

    job     = Gaia.launch_job_async(query) #, dump_to_file=True)
    results = job.get_results()
    N = len(results)
    return results, N


def get_stellar_density(ra, dec):
    _, count = get_main_table(ra, dec, radius=100)
    count /= (np.pi*(100/3600)**2)
    return int(count)
    

def qso_candidate_table(source_id):
    query = f"SELECT * \
    FROM gaiadr3.qso_candidates \
    WHERE gaiadr3.qso_candidates.source_id="+str(source_id)+""

    job     = Gaia.launch_job_async(query) #, dusmp_to_file=True)
    results = job.get_results()
    return results

def plot_on_ax(ax, url, title):
    response = requests.get(url)
    if response.status_code==200:
        img = Image.open(BytesIO(response.content))
        ax.imshow(img)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    optical_data = np.array(img)
    return ax,img


def panstarrs_on_ax(ax, ra, dec):
    #download which filters are available imaging
    title_colours = 'grz'
    filters = list(np.unique(list(ps.getimages(ra,dec,size=5,filters="grizy")['filter'])))

    #sort into blue to red wavelengths
    order = ['g', 'r', 'i', 'z', 'y']
    filter_positions = [order.index(filt) for filt in filters]
    filters = list(np.array(filters)[np.argsort(filter_positions)])

    #choose the best 3 (grz or whatever is available)
    if len(filters)>2:
        while len(filters)>3:
            if 'y' in filters:
                filters.remove('y')
            else:
                if 'i' in filters:
                    filters.remove('i')

        #download colour image and plot it
        title_colours = ''.join(filters)
        im = ps.getcolorim(ra, dec, size=60, output_size=None, filters=title_colours, format="png")
        ax.imshow(im)

    #remove x and y ticks even when data is not available
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Pan-STARRS '+title_colours)
    return ax

def gaia_on_ax(ax, ra, dec):
    data, N= get_main_table(ra, dec, radius=10)
    if N:
        ras = (ra-data['ra'])*np.cos(dec*np.pi/180.)*3600
        decs = (data['dec']-dec)*3600
        sizes = 23-data['phot_g_mean_mag']
        pmsigs = ((data['pmra']/data['pmra_error'])**2. + (data['pmdec']/data['pmdec_error'])**2.)**0.5
        varflag = data['phot_variable_flag']
        var = ['VAR' if flag=='VARIABLE' else 'N/A' for flag in varflag]
        s = 0
        for i in range(len(ras)):
            if data['in_qso_candidates'][i]:
                s+=1
                qso_data = qso_candidate_table(data['source_id'][i])
                z = qso_data['redshift_qsoc'].filled(np.nan)[0]
                if np.isfinite(z):
                    z = round(z, 3)
                zscore = qso_data['zscore_qsoc'].filled(np.nan)[0]
                if np.isfinite(zscore):
                    zscore = round(zscore, 3)
                qscore = data['classprob_dsc_combmod_quasar'].filled(np.nan)[i]
                if np.isfinite(qscore):
                    qscore = round(qscore, 3)
                ax.scatter(ras[i], decs[i], s=sizes[i]*25, facecolors='none', edgecolors='b')
                ax.text(s=str(s), x=ras[i]-0.25, y=decs[i]+0.25)
                ax.text(x=0.5, y=1.04-0.1*s, s=str(s)+': '+str(z)+' ('+str(zscore)+', '+str(qscore)+') '+var[i], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            else:
                ax.scatter(ras[i], decs[i], s=sizes[i]*25, facecolors='none', edgecolors='r')
            
            if pmsigs[i]:
                ax.arrow(ras[i], decs[i], dx=data['pmra'][i]/data['pmra_error'][i]/5, dy=data['pmdec'][i]/data['pmdec_error'][i]/5, width=0.03, fc='black', ec='black', color='black')
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Gaia DR3 z (zp, qp) var')


def gaia_specs_on_ax(ax, ra, dec):
    """ Return plot of Gaia DR3 spectrum, number of spectra
    Queries the Gaia DR3 table and returns only objects within "radius" arcseconds
    which have BP/RP spectra, and plots these on the provided ax
    """
    verbose = True
    radius = 5./3600
    query = f"SELECT * \
    FROM gaiadr3.gaia_source \
    WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS',"+str(ra)+","+str(dec)+","+str(radius)+"))=1 AND gaiadr3.gaia_source.has_xp_continuous='True';"

    if verbose:
        print('Querying the gaiadr3.gaia_source table for detections with has_xp_continuous=True')
    job     = Gaia.launch_job_async(query) #, dump_to_file=True)
    results = job.get_results()
    if len(results)==0:
        if verbose:
            print('No results found')
    else:
        if verbose:
            print(len(results), 'results found. Downloading XP_CONT spectra')
        retrieval_type = 'XP_CONTINUOUS' # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
        #retrieval_type = 'ALL'
        data_structure = 'INDIVIDUAL'    # Options are: 'INDIVIDUAL', 'COMBINED', 'RAW'
        data_release   = 'Gaia DR3'      # Options are: 'Gaia DR3' (default), 'Gaia DR2'

        #get the BP/RP spectra
        datalink = Gaia.load_data(ids=results['source_id'], data_release = data_release, retrieval_type=retrieval_type, data_structure = data_structure, verbose = verbose, output_file = None)
        dl_keys  = [inp for inp in datalink.keys()]
        dl_keys.sort()
        
        for key in dl_keys:
            datalink[key][0].to_table().write(key, format='votable', overwrite=True)
            calibrated_spectra, sampling = calibrate(key)
            os.remove(key)
            x, y, e = sampling, calibrated_spectra.flux[0], calibrated_spectra.flux_error[0]
            ax.plot(sampling, calibrated_spectra.flux[0])
            ax.fill_between(x, y - e, y + e, alpha=0.2)


        ymin, ymax = ax.get_ylim()[0], ax.get_ylim()[1]
        '''for final_line in final_lines:
            wavline = (1+zsrc[i])*linedict[final_line]
            if (wavline/10.>ax.get_xlim()[0]) & (wavline/10.<ax.get_xlim()[1]):
                ax.vlines(wavline/10., ymin=ymin, ymax=ymax, ls='--', color='k', alpha=0.3)
                ax.text((wavline+150.)/10., (ymax-ymin)*0.9 + ymin, final_line, fontsize=8, horizontalalignment='left', verticalalignment='top', color='k')'''
        ax.set_ylim(ymin, ymax)
    ax.set_xlabel('wavelength')
    ax.set_ylabel('flux')
    ax.set_title('Gaia BP/RP spectra')
    return ax


def gaia_lightcurves_on_ax(ax, ra, dec):
    """ Return plot of Gaia DR3 spectrum, number of spectra
    Queries the Gaia DR3 table and returns only objects within "radius" arcseconds
    which have BP/RP spectra, and plots these on the provided ax
    """
    verbose = True
    radius = 5./3600
    query = f"SELECT * \
    FROM gaiadr3.gaia_source \
    WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS',"+str(ra)+","+str(dec)+","+str(radius)+"))=1 AND gaiadr3.gaia_source.has_xp_continuous='True';"

    if verbose:
        print('Querying the gaiadr3.gaia_source table for detections with has_xp_continuous=True')
    job     = Gaia.launch_job_async(query) #, dump_to_file=True)
    results = job.get_results()
    if len(results)==0:
        if verbose:
            print('No results found')
    else:
        if verbose:
            print(len(results), 'results found. Downloading XP_CONT spectra')
        retrieval_type = 'EPOCH_PHOTOMETRY' # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
        #retrieval_type = 'ALL'
        data_structure = 'INDIVIDUAL'    # Options are: 'INDIVIDUAL', 'COMBINED', 'RAW'
        data_release   = 'Gaia DR3'      # Options are: 'Gaia DR3' (default), 'Gaia DR2'

        #get the BP/RP spectra
        datalink = Gaia.load_data(ids=results['source_id'], data_release = data_release, retrieval_type=retrieval_type, data_structure = data_structure, verbose = verbose, output_file = None)
        dl_keys  = [inp for inp in datalink.keys()]
        dl_keys.sort()
        
        for key in dl_keys:
            dl_out = datalink[key][0].to_table()
            phot_set = dl_out[dl_out['band']=='G']
            plt.plot(phot_set['time'], phot_set['mag'], 'x', color='black')

        ymin, ymax = ax.get_ylim()[0], ax.get_ylim()[1]
        '''for final_line in final_lines:
            wavline = (1+zsrc[i])*linedict[final_line]
            if (wavline/10.>ax.get_xlim()[0]) & (wavline/10.<ax.get_xlim()[1]):
                ax.vlines(wavline/10., ymin=ymin, ymax=ymax, ls='--', color='k', alpha=0.3)
                ax.text((wavline+150.)/10., (ymax-ymin)*0.9 + ymin, final_line, fontsize=8, horizontalalignment='left', verticalalignment='top', color='k')'''
        ax.set_ylim(ymin, ymax)
    ax.set_xlabel(f'JD date [days]')
    ax.set_ylabel(f'magnitude [mag]')
    ax.set_title('Gaia G lightcurve')
    return ax

def sdss_spec_on_ax(ax, ra, dec):
    pos = coords.SkyCoord(ra, dec,  unit='deg')

    #this fails sometimes, even when get_spectra works
    xid_tab = SDSS.query_region(pos, spectro=True)

    #seems to fail at some positions and returns a table with a htmlhead column
    if xid_tab:
        if 'htmlhead' in xid_tab.colnames:
            xid_tab = None

    #this is only DR16 though
    #xid_tab = Vizier.query_region(pos, radius=8*u.arcsec, catalog='V/154/sdss16', column_filters={'zsp': '!='})
    #if xid_tab:
    #    xid_tab = xid_tab[0]

    ymin, ymax = 0, 1
    xmin, xmax = 3800, 10000
    colours = ['red', 'blue', 'black', 'yellow', 'purple']
    if xid_tab:
        xid_tab = xid_tab[np.unique(xid_tab['ra'], True)[1]]
        xmins, xmaxs = [], []
        ymins, ymaxs = [], []
        for k in range(len(xid_tab)):
            xid = xid_tab[k]
            sp = SDSS.get_spectra(plate=xid['plate'], mjd=xid['mjd'], fiberID=xid['fiberID'])
            
            #sp = SDSS.get_spectra(plate=int(xid['Sp-ID'].split('-')[0]), mjd=int(xid['Sp-ID'].split('-')[1]), fiberID=int(xid['Sp-ID'].split('-')[2]), data_release=16)
            flux, model, logwav, ivar = sp[0][1].data['flux'], sp[0][1].data['model'], sp[0][1].data['loglam'], sp[0][1].data['ivar']
            smoothspec = gaussian_filter(flux, 5)
            #get ylims
            kappa = np.where((logwav>np.log10(xmin))&(logwav<np.log10(xmax))&(ivar>np.nanpercentile(ivar, 1)))[0]
            ymin, ymax = np.min(smoothspec[kappa]), np.max(smoothspec[kappa])

            #plot the data and best fit model from SDSS
            ax.plot(10.**logwav, smoothspec, color=colours[k], alpha=0.8)
            ax.plot(10.**logwav, model, color='black', ls='--', alpha=0.5)

            xmins.append(np.min(10.**logwav))
            xmaxs.append(np.max(10.**logwav))
            ymins.append(ymin)
            ymaxs.append(ymax)

        xmin, xmax = np.min(xmins), np.max(xmaxs)
        ymin, ymax = np.min(ymins), np.max(ymaxs)
    #ax.set_title('RA, DEC='+str(ra)+', '+str(dec))
    ax.set_xlim(np.max([3800, xmin]), np.min([10000, xmax]))
    ax.set_ylim(ymin, ymax)

    #if plot_lines:
    #    for line in lines:
    #        ax.vlines(x=line*(1+z), ymin=0, ymax=1, color='black', ls='--', transform=ax.get_xaxis_transform(), alpha=0.5)

    return ax



def rms_value(fits_image):

    image = fits.open(fits_image)
    prihdr = image[0].header
    data = image[0].data
    image_data = data#[0,0,:,:]
    image_data=np.nan_to_num(image_data)

    img_min, img_max = np.min(image_data),np.max(image_data)
    img_mean=np.mean(image_data);img_median=np.median(image_data);img_std=np.std(image_data)
    print (img_min,img_max,img_mean,img_median,img_std)

    sigclip = SigmaClip(sigma = 3.0)#, iters=5)
    image_data_sigclip = sigclip(image_data)
    image_min_clip, image_max_clip = np.min(image_data_sigclip),np.max(image_data_sigclip)
    image_mean_clip=np.mean(image_data_sigclip);image_median_clip=np.median(image_data_sigclip);image_std_clip=np.std(image_data_sigclip)
    print (image_min_clip,image_max_clip,image_mean_clip,image_median_clip,image_std_clip)


    image_mean, image_median, image_stddev = sigma_clipped_stats(image_data, sigma = 3.0)#, iters = 5)
    print (image_mean,image_median)
    print (image_stddev)
    rms=image_stddev 

    return(rms)





   
def plot_image(fits_image,ra,dec):

    rms=rms_value(fits_image)

    hdu = fits.open(fits_image)[0]
    data = hdu.data
    header = hdu.header
    new_data = data[0,0,:,:]

    del header['*3']
    del header['*4']
    header['WCSAXES']= 2
    header['NAXIS']= 2

    contour_levels = np.linspace(np.min(new_data), np.max(new_data), num=10)
    position = SkyCoord(ra*u.deg, dec*u.deg,frame='fk5',equinox='J2000.0') 

    
    w = WCS(header)
    #cutout = np.rot90(np.fliplr(Cutout2D(new_data,position=position,size=size,wcs=w,mode='partial').data))
    #cutout = Cutout2D(new_data,position=position,size=size,wcs=w,mode='partial').data

    return(new_data,contour_levels,position,w)


def sdss_spec_positions_on_ax(ax, ra, dec):
    pos = coords.SkyCoord(ra, dec,  unit='deg')

    xid_tab = SDSS.query_region(pos, spectro=True)
    if xid_tab:
        if 'htmlhead' in xid_tab.colnames:
            xid_tab = None


    #xid_tab = Vizier.query_region(pos, radius=8*u.arcsec, catalog='V/154/sdss16', column_filters={'zsp': '!='})
    #if xid_tab:
    #    xid_tab = xid_tab[0]

    colours = ['red', 'blue', 'black', 'yellow', 'purple']
    ax.text(x=-0.8, y=0.95, s='SDSS spectra', horizontalalignment='center', verticalalignment='top', transform=ax.transAxes, color='k', fontsize=15)
    if xid_tab:
        xid_tab = xid_tab[np.unique(xid_tab['ra'], True)[1]]
        for k in range(len(xid_tab)):
            xid = xid_tab[k]
            sp = SDSS.get_spectra(plate=xid['plate'], mjd=xid['mjd'], fiberID=xid['fiberID'])
            #sp = SDSS.get_spectra(plate=int(xid['Sp-ID'].split('-')[0]), mjd=int(xid['Sp-ID'].split('-')[1]), fiberID=int(xid['Sp-ID'].split('-')[2]), data_release=16)

            class_sdss = sp[0][2].data['CLASS'][0]
            if sp[0][2].data['SUBCLASS'][0]:
                class_sdss += '-'+sp[0][2].data['SUBCLASS'][0]
            z_sdss = round(sp[0][2].data['Z'][0], 3)

            x = (ra-xid_tab['ra'][k])*np.cos(dec*np.pi/180.)*3600
            y = (xid_tab['dec'][k]-dec)*3600

            #x = (ra-xid_tab['RA_ICRS'][k])*np.cos(dec*np.pi/180.)*3600
            #y = (xid_tab['DE_ICRS'][k]-dec)*3600
            print(colours[k])
            ax.scatter(x, y, s=[250, 200, 150][k], facecolors='none', edgecolors=colours[k])
            ax.text(x=-0.8, y=0.75-0.1*k, s=class_sdss+' (z='+str(z_sdss)+')', horizontalalignment='center', verticalalignment='top', transform=ax.transAxes, color=colours[k])
                
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.set_xticks([])
    ax.set_yticks([])

def make_full_cutout(ra, dec, savename):
    fig = plt.figure(figsize=(10, 10))
    
    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(ra)+'&dec='+str(dec)+'&layer=ls-dr10&pixscale=0.06'
    ax = plt.axes([0.05, 0.74, 0.2, 0.2])
    plot_on_ax(ax, url, 'DECALS DR10')

    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(ra)+'&dec='+str(dec)+'&layer=sdss&pixscale=0.06'
    ax = plt.axes([0.28, 0.74, 0.2, 0.2])
    plot_on_ax(ax, url, 'model')

    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(ra)+'&dec='+str(dec)+'&layer=unwise-neo7&pixscale=0.25&size=256'
    ax = plt.axes([0.51, 0.74, 0.2, 0.2])
    plot_on_ax(ax, url, 'WISE FOV1')

    ax = plt.axes([0.74, 0.74, 0.2, 0.2])
    panstarrs_on_ax(ax, ra, dec)


    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(ra)+'&dec='+str(dec)+'&layer=ls-dr10&pixscale=0.25&size=256'
    ax = plt.axes([0.05, 0.51, 0.2, 0.2])
    plot_on_ax(ax, url, 'DECAL FOV1')

    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(ra)+'&dec='+str(dec)+'&layer=ls-dr10&pixscale=0.25&size=600'
    ax = plt.axes([0.28, 0.51, 0.2, 0.2])
    plot_on_ax(ax, url, 'DECAL FOV2')

    ax = plt.axes([0.51, 0.51, 0.2, 0.2])
    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(ra)+'&dec='+str(dec)+'&layer=sdss&pixscale=0.06'
    plot_on_ax(ax, url, 'SDSS')

    ax = plt.axes([0.74, 0.51, 0.2, 0.2])
    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(ra)+'&dec='+str(dec)+'&layer=hsc-dr2&pixscale=0.06'
    plot_on_ax(ax, url, 'HSC')

    ax = plt.axes([0.05, 0.28, 0.2, 0.2])

    cutout = Cutout2D(new_data,position=position,size=u.Quantity((15.36,15.36), u.arcsec),wcs=w,mode='partial').data
    ax.imshow(cutout,cmap='gray',vmin=np.nanpercentile(cutout, 0.1),origin='lower', vmax=np.nanpercentile(cutout, 99))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('MeerKAT FOV1')

    ax = plt.axes([0.28, 0.28, 0.2, 0.2])
    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(ra)+'&dec='+str(dec)+'&layer=vlass1.2&pixscale=0.06&size=256'
    plot_on_ax(ax, url, 'VLASS FOV1')

    
    ax = plt.axes([0.51, 0.28, 0.2, 0.2])
    cutout = Cutout2D(new_data,position=position,size=u.Quantity((132,132), u.arcsec),wcs=w,mode='partial').data
    ax.imshow(cutout,cmap='gray',vmin=np.nanpercentile(cutout, 0.1),origin='lower', vmax=np.nanpercentile(cutout, 99))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('MeerKAT FOV3')



    ax = plt.axes([0.74, 0.28, 0.2, 0.2])
    url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra='+str(ra)+'&dec='+str(dec)+'&layer=vlass1.2&pixscale=0.06&size=2200'
    plot_on_ax(ax, url, 'VLASS FOV3')



    ax = plt.axes([0.05, 0.1, 0.41, 0.15])
    gaia_specs_on_ax(ax, ra, dec)

    ax = plt.axes([0.53, 0.1, 0.41, 0.15])
    gaia_lightcurves_on_ax(ax, ra, dec)


    ax = plt.axes([0.05, -0.1, 0.41, 0.15])
    sdss_spec_on_ax(ax, ra, dec)

    ax = plt.axes([0.72, -0.1, 0.15, 0.15])
    sdss_spec_positions_on_ax(ax, ra, dec)

    titlestring = radecname.Jname([ra], [dec])[0]+': '+str(ra)[:10]+', '+str(dec)[:10]
    titlestring += '\n'
    titlestring += 'local Gaia density = '+str(get_stellar_density(ra, dec))+' per sq. deg.'
    plt.suptitle(titlestring, y=1.01)

    plt.savefig(savename, bbox_inches='tight')
    plt.close()



if __name__ == '__main__':


    fits_image=sys.argv[1]

    cat=sys.argv[2]

    # correct_axes(fits_image)	
    cat = pd.read_table(cat, sep="\s+")

    ras=np.array(cat['RA'])
    decs=np.array(cat['DEC'])
    # ras = [354.5760, ]
    # decs = [-0.9362,]

    for i in range(50):

        ra, dec = ras[i], decs[i]

        new_data,contour_levels,position,w=plot_image(fits_image,ra,dec)
     
        print(ra,dec)
        name = radecname.Jname([ra],[dec])[0]
        #savename = './full_cutouts/'+str(i)+'_'+str(ra)+'_'+str(dec)+'_full_cutout.png'
        if not os.path.exists('./robert_cutouts/'+str(i)+'.png'):
            make_full_cutout(ra, dec, savename='./cutouts/'+name+'_'+str(i)+'.png')