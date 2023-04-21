import pandas as pd
import numpy as np
ids = []
cat =pd.read_csv(r'SDSS_Spec_Photo_DeCALS.csv')


def flux2mag(flux):
    return (22.5-2.5*np.log10(flux))


def ivar2errmag(ivar, flux):
    df = np.sqrt(1/ivar) 
    return (2.5 / np.log(10)) * (df / flux)


f = open("observations_new.dat", "w")
for i in range(10):
    line  = '  '
    line += str(cat['objid_1'][i])
    line += '       '
    line += str(cat['zph'][i])
    line += '       '
    for col in (['FLUX_G','FLUX_IVAR_G' ,'FLUX_R', 'FLUX_IVAR_R','FLUX_I', 'FLUX_IVAR_I','FLUX_Z', 'FLUX_IVAR_Z','FLUX_W1','FLUX_IVAR_W1','FLUX_W2','FLUX_IVAR_W2','FLUX_W3','FLUX_IVAR_W3','FLUX_W4','FLUX_IVAR_W4']):
        if (i % 2)  == 0:
            print(i)
            if cat[col][i] < 0:
                line += str(0)  
                line += '       '
            else:
                line += str(flux2mag(cat[col][i]))
                line += '       '
        else:#elif (i % 3)  == 0 and (i != 0):
            print(i)
            if  (cat[col][i-1]) < 0 :
                line += str(0) 
                line += '       '
            else:
                line += str(ivar2errmag(cat[col][i],flux2mag(cat[col][i-1])))
                line += '       '                                               
    print(line)
    f.write(line)
    line = line.rstrip()
    f.write('\n')