import pandas as pd
import numpy as np
ids = []
cat =pd.read_csv(r'SDSS_Spec_Photo_DeCALS.csv')


def flux2mag(flux):
    return (22.5-2.5*np.log10(flux))


# def ivar2errmag(ivar, flux):
#     df = np.sqrt(1/ivar) 
#     return (2.5 / np.log(10)) * (df / flux)


def ivar2errflux(ivar):
    df = np.sqrt(1/ivar) 
    return df

fluxes=['FLUX_G','FLUX_R', 'FLUX_I', 'FLUX_Z','FLUX_W1','FLUX_W2','FLUX_W3','FLUX_W4']
flux_uncertanties=['FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_I' ,'FLUX_IVAR_Z','FLUX_IVAR_W1','FLUX_IVAR_W2','FLUX_IVAR_W3','FLUX_IVAR_W4']

f = open("observations.dat", "w")

line = '# HEADER'
f.write(line)
f.write('\n')
for i in range(10):
    line  = '  '
    line += str(cat['objid_1'][i]).replace('1237663278462','')
    line += '       '
    line += str(cat['zph'][i])
    line += '       '
    for j in range(len(fluxes)):
            if cat[fluxes[j]][i] < 0:
                line += str(0)  
                line += '       '
                line += str(0)  
                line += '       '
            else:
                line += str(cat[fluxes[j]][i]*3.631*10**-6)
                line += '       '
                line += str(ivar2errflux(cat[flux_uncertanties[j]][i])*3.631*10**-6)
                line += '       '  
    print(line)
    f.write(line)
    line = line.rstrip()
    f.write('\n')