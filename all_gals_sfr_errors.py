import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import os,sys,string, math
from astropy.io import ascii
from astropy.table import Table, Column, vstack
import glob
import re


def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


#fit_files=glob.glob('all_gals/*.fit')
#fit_files=glob.glob('all_gals_fors2/*.fit')
fit_files=glob.glob('*.fit')
#fit_files=glob.glob('all_gals_vimos/*.fit')
fit_files=sorted(fit_files, key=natural_key)
#sed_files=sorted(sed_files, key=natural_key)

print( len(fit_files))

index=[]
for i in fit_files:
    #x=i.split('/')[1] #for other galaxies
    y=x.split('.')[0]
    #y=i.split('.')[0]
    index=np.append(index,y)
    print(index)

sfr_sigmas = np.array([])
ssfr_sigmas = np.array([])
ms_sigmas = np.array([])
sfr = np.array([])
ssfr = np.array([])
ms = np.array([])


for j,file in enumerate(fit_files):
    print( j)
    #if (index[j]=='67')|(index[j]=='83')|(index[j]=='155')|(index[j]=='290')|\
       #(index[j]=='325')|(index[j]=='603')|(index[j]=='1095')|(index[j]=='1179'):
    if (index[j]=='118')|(index[j]=='196'):#|(index[j]=='200')|(index[j]=='188'):
    #if (index[j]=='10'):
    #if (index[j]=='133')|(index[j]=='131')|(index[j]=='112')|(index[j]=='110')|(index[j]=='108')|(index[j]=='105')|(index[j]=='102')|(index[j]=='97')|(index[j]=='75')|(index[j]=='56')|(index[j]=='55')|(index[j]=='51'):
    #if 0:
        print( 'Galaxy %s \n' %index[j])        
        sfr_coeff = [-99,-99,-99]
        ssfr_coeff = [-99,-99,-99]
        ms_coeff = [-99,-99,-99]
    #else:
    elif 1:
    #if 1:
        f=open(file)
        lines = f.readlines()
        fluxes = lines[2]
        u = np.float(fluxes.split()[0])
        lines_values = lines[10]
        sSFR = np.float(lines_values.split()[4])
        Mstar = np.float(lines_values.split()[5])
        SFR = np.float(lines_values.split()[15])    
        lines_sfr = lines[619:679]#.strip('\n')
        lines_ms = lines[209:269]#.strip('\n')
        lines_ssfr = lines[136:206]#.strip('\n')
        sfr_val = np.array([np.float(lines_sfr[i].split()[0]) for i in range(len(lines_sfr))])
        sfr_pdf = np.array([np.float(lines_sfr[i].split()[1]) for i in range(len(lines_sfr))])
        ssfr_val = np.array([np.float(lines_ssfr[i].split()[0]) for i in range(len(lines_ssfr))])
        ssfr_pdf = np.array([np.float(lines_ssfr[i].split()[1]) for i in range(len(lines_ssfr))])
        ms_val = np.array([np.float(lines_ms[i].split()[0]) for i in range(len(lines_ms))])
        ms_pdf = np.array([np.float(lines_ms[i].split()[1]) for i in range(len(lines_ms))])
        f.close()
        
        maxima_sfr = argrelextrema(sfr_pdf[sfr_val > 0.], np.greater)
        #print maxima_sfr
        #print sfr_val[sfr_val > 0.][maxima_sfr]
        #print sfr_pdf[sfr_val > 0.][maxima_sfr]
        
        max_pdf = sfr_pdf[sfr_val > 0.][maxima_sfr]
        
        condition = max_pdf[max_pdf>0.09]
        #print condition
        
        max_val = sfr_val[sfr_val > 0.][maxima_sfr]
        
        condition2 = max_val[np.where(max_pdf>0.09)]
        #print condition2            
            
        if 0:
            print('Galaxy %s \n' %index[j])
            
            #max_sfr=np.maximum(sfr_pdf)
            #max_sfr = sfr_pdf.index(max(sfr_pdf))
            #max_ssfr = ssfr_pdf.index(max(ssfr_pdf))
            #max_ms = ms_pdf.index(max(ms_pdf))
                    
            sfr_coeff,sfr_varmatrix = curve_fit(gauss,sfr_val,sfr_pdf,p0=[0.15,sfr_val[sfr_pdf.argmax()],0.05])
            sfr_gauss = gauss(sfr_val,*sfr_coeff)
            
            ms_coeff,ms_varmatrix = curve_fit(gauss,ms_val,ms_pdf,p0=[0.45,ms_val[ms_pdf.argmax()],0.5])
            ms_gauss = gauss(ms_val,*ms_coeff)
                
            ssfr_coeff,ssfr_varmatrix = curve_fit(gauss,ssfr_val,ssfr_pdf,p0=[0.15,ssfr_val[ssfr_pdf.argmax()],1])
            ssfr_gauss = gauss(ssfr_val,*ssfr_coeff)
        
            
            print('Sigma SFR : ', sfr_coeff[2],'\n','Sigma sSFR : ', ssfr_coeff[2],'\n','Sigma Ms : ', ms_coeff[2])
            print('---------SFR---------','\n')
            print('Gaussian fit :', sfr_coeff[1],'|','magphys value :', np.log10(SFR),'\n')
            print('---------sSFR---------','\n')
            print('Gaussian fit :', ssfr_coeff[1],'|','magphys value :', np.log10(sSFR),'\n')
            print('---------Mstar---------','\n')
            print('Gaussian fit :', ms_coeff[1],'|','magphys value :', np.log10(Mstar),'\n')
            print('\n','\n')
            
            #plt.subplot(3,1,1)
            #plt.plot(sfr_val,sfr_pdf)
            #plt.plot(sfr_val,sfr_gauss,color='k')
            #plt.xlabel('SFR')
        
            #plt.subplot(3,1,2)
            #plt.plot(ssfr_val,ssfr_pdf)
            #plt.plot(ssfr_val,ssfr_gauss,color='k')
            #plt.xlabel('sSFR')
        
            #plt.subplot(3,1,3)
            #plt.plot(ms_val,ms_pdf)
            #plt.plot(ms_val,ms_gauss,color='k')
            #plt.xlabel('Ms')
            #plt.tight_layout()
            #plt.show()           

        elif 0:
        #elif (int(index[j]) in bad_results):
            print('Galaxy %s \n' %index[j])
            print('in bad results')
            
            sfr_pdf[sfr_val > 1.5] = 0.
                    
            sfr_coeff,sfr_varmatrix = curve_fit(gauss,sfr_val,sfr_pdf,p0=[0.15,sfr_val[sfr_pdf.argmax()],0.05])
            sfr_gauss = gauss(sfr_val,*sfr_coeff)
            
            ms_coeff,ms_varmatrix = curve_fit(gauss,ms_val,ms_pdf,p0=[0.45,ms_val[ms_pdf.argmax()],0.5])
            ms_gauss = gauss(ms_val,*ms_coeff)
                
            ssfr_coeff,ssfr_varmatrix = curve_fit(gauss,ssfr_val,ssfr_pdf,p0=[0.15,ssfr_val[ssfr_pdf.argmax()],1])
            ssfr_gauss = gauss(ssfr_val,*ssfr_coeff)
        
            #sfr_coeff[1] = sfr_val[sfr_val > 0.][maxima_sfr[0]][0]
            
            print('Sigma SFR : ', sfr_coeff[2],'\n','Sigma sSFR : ', ssfr_coeff[2],'\n','Sigma Ms : ', ms_coeff[2])
            print('---------SFR---------','\n')
            print('Gaussian fit :', sfr_coeff[1],'|','magphys value :', np.log10(SFR),'\n')
            print('---------sSFR---------','\n')
            print('Gaussian fit :', ssfr_coeff[1],'|','magphys value :', np.log10(sSFR),'\n')
            print('---------Mstar---------','\n')
            print('Gaussian fit :', ms_coeff[1],'|','magphys value :', np.log10(Mstar),'\n')
            print('\n','\n')
            
            if (index[j]=='166')|(index[j]=='1035')|(index[j]=='1155')|(index[j]=='1551')|(index[j]=='2022'):
                sfr_coeff = [-99,-99,-99]
                ssfr_coeff = [-99,-99,-99]
                ms_coeff = [-99,-99,-99]
                            
            #plt.subplot(3,1,1)
            #plt.plot(sfr_val,sfr_pdf)
            #plt.plot(sfr_val,sfr_gauss,color='k')
            #plt.xlabel('SFR')
        
            #plt.subplot(3,1,2)
            #plt.plot(ssfr_val,ssfr_pdf)
            #plt.plot(ssfr_val,ssfr_gauss,color='k')
            #plt.xlabel('sSFR')
        
            #plt.subplot(3,1,3)
            #plt.plot(ms_val,ms_pdf)
            #plt.plot(ms_val,ms_gauss,color='k')
            #plt.xlabel('Ms')
            #plt.tight_layout()
            #plt.show()  


        #elif ((len(condition) >= 2)&(len(condition2) >= 2)):#&(np.abs(condition2[0]-condition2[-1])>=0.5):#&(u<1.e-30):
        #elif ((len(condition) >= 2)&(len(condition2) >= 2))|((index[j]=='128')|(index[j]=='130')):#&(np.abs(condition2[0]-condition2[-1])>=0.5):#&(u<1.e-30):
        if 0:
            #if (np.abs(condition2[0]-condition2[-1])>=0.5):
            if (np.abs(condition2[0]-condition2[-1])>=0.5)|(index[j]=='3')|(index[j]=='18'):
            #if (np.abs(condition2[0]-condition2[-1])>=0.5)|(index[j]=='81')|(index[j]=='125')|(index[j]=='128')|(index[j]=='130'):
                print('Galaxy %s \n' %index[j])
                print('double peak condition')
                
                sfr_pdf[sfr_val > 1.65] = 0.
                        
                sfr_coeff,sfr_varmatrix = curve_fit(gauss,sfr_val,sfr_pdf,p0=[0.15,sfr_val[sfr_pdf.argmax()],0.05])
                sfr_gauss = gauss(sfr_val,*sfr_coeff)
                
                ms_coeff,ms_varmatrix = curve_fit(gauss,ms_val,ms_pdf,p0=[0.45,ms_val[ms_pdf.argmax()],0.5])
                ms_gauss = gauss(ms_val,*ms_coeff)
                    
                ssfr_coeff,ssfr_varmatrix = curve_fit(gauss,ssfr_val,ssfr_pdf,p0=[0.15,ssfr_val[ssfr_pdf.argmax()],1])
                ssfr_gauss = gauss(ssfr_val,*ssfr_coeff)
            
                #sfr_coeff[1] = sfr_val[sfr_val > 0.][maxima_sfr[0]][0]
                
                print('Sigma SFR : ', sfr_coeff[2],'\n','Sigma sSFR : ', ssfr_coeff[2],'\n','Sigma Ms : ', ms_coeff[2])
                print('---------SFR---------','\n')
                print('Gaussian fit :', sfr_coeff[1],'|','magphys value :', np.log10(SFR),'\n')
                print('---------sSFR---------','\n')
                print('Gaussian fit :', ssfr_coeff[1],'|','magphys value :', np.log10(sSFR),'\n')
                print('---------Mstar---------','\n')
                print('Gaussian fit :', ms_coeff[1],'|','magphys value :', np.log10(Mstar),'\n')
                print('\n','\n')
                
                #sfr_coeff = [-99,-99,-99]
                #ssfr_coeff = [-99,-99,-99]
                #ms_coeff = [-99,-99,-99]
                            
                #plt.subplot(3,1,1)
                #plt.plot(sfr_val,sfr_pdf)
                #plt.plot(sfr_val,sfr_gauss,color='k')
                #plt.xlabel('SFR')
            
                #plt.subplot(3,1,2)
                #plt.plot(ssfr_val,ssfr_pdf)
                #plt.plot(ssfr_val,ssfr_gauss,color='k')
                #plt.xlabel('sSFR')
            
                #plt.subplot(3,1,3)
                #plt.plot(ms_val,ms_pdf)
                #plt.plot(ms_val,ms_gauss,color='k')
                #plt.xlabel('Ms')
                #plt.tight_layout()
                #plt.show()                    
        
            else:
                print( 'Galaxy %s \n' %index[j])
                print( 'no double peak problems')
                #max_sfr=np.maximum(sfr_pdf)
                print( sfr_pdf.argmax())
                print( ms_val[ms_pdf.argmax()])
                        
                sfr_coeff,sfr_varmatrix = curve_fit(gauss,sfr_val,sfr_pdf,p0=[0.15,sfr_val[sfr_pdf.argmax()],0.05])
                sfr_gauss = gauss(sfr_val,*sfr_coeff)
                
                ssfr_coeff,ssfr_varmatrix = curve_fit(gauss,ssfr_val,ssfr_pdf,p0=[0.5,ssfr_val[ssfr_pdf.argmax()],1])
                ssfr_gauss = gauss(ssfr_val,*ssfr_coeff)
                
                ms_coeff,ms_varmatrix = curve_fit(gauss,ms_val,ms_pdf,p0=[0.25,ms_val[ms_pdf.argmax()],0.5])
                ms_gauss = gauss(ms_val,*ms_coeff)
            
                
                print('Sigma SFR :  ',sfr_coeff[2],'\n','Sigma sSFR : ',ssfr_coeff[2],'\n','Sigma Ms :   ',ms_coeff[2])
                print( '---------SFR---------','\n')
                print( 'Gaussian fit :',sfr_coeff[1],'|','magphys value :', np.log10(SFR),'\n')
                print( '---------sSFR---------','\n')
                print( 'Gaussian fit :',ssfr_coeff[1],'|','magphys value :', np.log10(sSFR),'\n')
                print( '---------Mstar---------','\n')
                print( 'Gaussian fit :',ms_coeff[1],'|','magphys value :', np.log10(Mstar),'\n')
                print( '\n','\n')
                
                #plt.subplot(3,1,1)
                #plt.plot(sfr_val,sfr_pdf)
                #plt.plot(sfr_val,sfr_gauss,color='k')
                #plt.xlabel('SFR')
            
                #plt.subplot(3,1,2)
                #plt.plot(ssfr_val,ssfr_pdf)
                #plt.plot(ssfr_val,ssfr_gauss,color='k')
                #plt.xlabel('sSFR')
            
                #plt.subplot(3,1,3)
                #plt.plot(ms_val,ms_pdf)
                #plt.plot(ms_val,ms_gauss,color='k')
                #plt.xlabel('Ms')
                #plt.tight_layout()
                #plt.show()  


        elif ((max(sfr_pdf)<0.7)&(max(sfr_pdf)>0.001))|((max(ms_pdf)<0.85)&(max(ms_pdf)>0.1)):#(max(sfr_pdf)<0.75):
            print( 'Galaxy %s \n' %index[j])
            print( 'maxima problem condition')
            #max_sfr=np.maximum(sfr_pdf)
            print( sfr_pdf.argmax())
            print( ms_val[ms_pdf.argmax()])
                    
            sfr_coeff,sfr_varmatrix = curve_fit(gauss,sfr_val,sfr_pdf,p0=[0.15,sfr_val[sfr_pdf.argmax()],0.05])
            sfr_gauss = gauss(sfr_val,*sfr_coeff)
            
            ssfr_coeff,ssfr_varmatrix = curve_fit(gauss,ssfr_val,ssfr_pdf,p0=[0.5,ssfr_val[ssfr_pdf.argmax()],1])
            ssfr_gauss = gauss(ssfr_val,*ssfr_coeff)
            
            ms_coeff,ms_varmatrix = curve_fit(gauss,ms_val,ms_pdf,p0=[0.25,ms_val[ms_pdf.argmax()],0.5])
            ms_gauss = gauss(ms_val,*ms_coeff)
        
            #plt.subplot(3,1,1)
            #plt.plot(sfr_val,sfr_pdf)
            #plt.plot(sfr_val,sfr_gauss,color='k')
            #plt.xlabel('SFR')
        
            #plt.subplot(3,1,2)
            #plt.plot(ssfr_val,ssfr_pdf)
            #plt.plot(ssfr_val,ssfr_gauss,color='k')
            #plt.xlabel('sSFR')
        
            #plt.subplot(3,1,3)
            #plt.plot(ms_val,ms_pdf)
            #plt.plot(ms_val,ms_gauss,color='k')
            #plt.xlabel('Ms')
            #plt.tight_layout()
            #plt.show()  
            
            print('Sigma SFR :  ',sfr_coeff[2],'\n','Sigma sSFR : ',ssfr_coeff[2],'\n','Sigma Ms :   ',ms_coeff[2])
            print('---------SFR---------','\n')
            print( 'Gaussian fit :',sfr_coeff[1],'|','magphys value :', np.log10(SFR),'\n')
            print( '---------sSFR---------','\n')
            print( 'Gaussian fit :',ssfr_coeff[1],'|','magphys value :', np.log10(sSFR),'\n')
            print( '---------Mstar---------','\n')
            print( 'Gaussian fit :',ms_coeff[1],'|','magphys value :', np.log10(Mstar),'\n')
            print( '\n','\n')
        
        else :
            print( 'Galaxy %s \n' %index[j])
            print( 'Sigma SFR :  ',0,'\n','Sigma sSFR : ',0,'\n','Sigma Ms :   ',0)
            print( '---------SFR---------','\n')
            print( 'magphys value :', SFR,'\n')
            print( '---------Mstar---------','\n')
            print( 'magphys value :', Mstar,'\n')
            print( '---------sSFR---------','\n')
            print( 'magphys value :', sSFR,'\n')
            print( '\n','\n')
            
            sfr_coeff = [-99, SFR, -99]
            ssfr_coeff = [-99, sSFR, -99]
            ms_coeff = [-99, Mstar, -99]
            

    sfr_sigmas = np.append(sfr_sigmas, sfr_coeff[2])
    ssfr_sigmas = np.append(ssfr_sigmas, ssfr_coeff[2])
    ms_sigmas = np.append(ms_sigmas, ms_coeff[2])
    sfr = np.append(sfr, sfr_coeff[1])
    ssfr = np.append(ssfr, ssfr_coeff[1])
    ms = np.append(ms, ms_coeff[1])
    #sfr=np.append(sfr,np.log10(SFR))
    #ssfr=np.append(ssfr,np.log10(sSFR))
    #ms=np.append(ms,np.log10(Mstar))
    
    
ids=range(len(sfr))    
#sys.exit()
       
data = Table([index,sfr,sfr_sigmas,ssfr,ssfr_sigmas,ms,ms_sigmas],
             names=['#filename','sfr','sigma_sfr','ssfr','sigma_ssfr','Mstar','sigma_Mstars'])

#ascii.write(data,output='all_gals_cl1411_magtot_final.dat', overwrite=True)
#ascii.write(data,output='fors2_magtot_cl1411_final.dat', overwrite=True)
ascii.write(data,output='mmt_magtot_cl1411_final.dat', overwrite=True)
#ascii.write(data,output='vimos_magtot_cl1411_final.dat', overwrite=True)

