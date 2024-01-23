import numpy as np
from astropy import units as u

unitL=u.erg/u.s/u.Hz
unitF=u.erg/u.s/u.Hz/(u.m**2)

def ABmagtoFlux(ABmag,error):
    """
    Convert AB magnitude to flux in erg/s/cm^2/Hz

    Parameters
    ----------
    ABmag : float
        AB magnitude
    error : float
        error on the AB magnitude
    """
    return np.array((10**((ABmag+48.60)/-2.5),-3.34407e-20*np.exp(-0.921034*ABmag)*error))*u.erg/(u.s*u.Hz*u.cm**2)

def WISEmagtoFlux(WISEmag,error,band):
    """
    The WISE dataproducts giva a magnitude with a zeropoint based on Vega.
    This can be converted using a simple offset.
    
    Parameters
    ----------
    WISEmag : float
        WISE magnitude
    error : float
        error on the WISE magnitude
    band : int
        WISE filterband (1,2,3 or 4)
    """
    if type(band)!=int or band==0:
        raise TypeError('Please give the WISE filterband as an integer')
    offset=[0,2.699,3.339,5.174,6.620]
    ABmag=WISEmag+offset[band]
    return ABmagtoFlux(ABmag,error)

def binning(bins,time,data,error,clean=True):
    """
    Bins the data in the given bins and returns the binned data

    input the dataset as seperate lists/1D-numpy arrays
    give the bins as a list of border values

    Parameters
    ----------
    bins : list
        list of border values of the bins
    time : list
        list of time values 
    data : list
        list of data values
    error : list
        list of error values
    clean : bool
        if True, the data will be cleaned from outliers before binning
    """
    binned_time=[]
    binned_data=[]
    binned_error=[]
    mask=np.invert(np.isnan(data)+np.isnan(error))
    time=time[mask]
    data=data[mask]
    error=error[mask]

    for i in range(len(bins)-1):
        mask=(time>=bins[i])&(time<=bins[i+1])
        if clean:
            median=np.mean(data[mask])
            std=np.std(data[mask])
            #std=np.max(error[mask])
            mask*=(data>median-3*std)&(data<median+3*std)

        if np.linalg.norm(mask,0)!=0:
            mean=wmean(data[mask],error[mask])
            #binned_time.append(wmean(time[mask],error[mask])[0])
            binned_time.append(np.mean(time[mask]))
            binned_data.append(mean[0])
            binned_error.append(mean[1])
            #if clean:
            #    mask2=(data[mask]>0.1*binned_data[-1])&(data[mask]<10*binned_data[-1])
            #    if np.linalg.norm(mask2,0)!=0:
            #        mean=wmean(data[mask][mask2],error[mask][mask2])
            #        binned_data[-1]=mean[0]
            #        binned_error[-1]=mean[1]
        elif i!=0 and i<len(bins)-2:
            binned_time.append(bins[i])
            mask=(time>=bins[i-1])&(time<=bins[i+2])
            if clean:
                median=np.mean(data)
                std=np.std(data)
                #std=np.max(error)
                mask*=(data>median-std)&(data<median+std)

            mean=wmean(data[mask],error[mask])
            binned_data.append(mean[0])
            binned_error.append(mean[1])
        else:
            binned_time.append(bins[i])
            binned_data.append(0)
            binned_error.append(0)



    return np.array((binned_time,binned_data,binned_error))

def FluxtoLum(flux,distance=1171*u.Mpc):
    """
    Convert flux to luminosity
    
    Parameters
    ----------
    flux : astropy.units.quantity.Quantity
        flux
    distance : astropy.units.quantity.Quantity
        distance
    """
    return (4*np.pi*distance**2*flux).to(unitL)

def LumtoFlux(lum,distance=1171*u.Mpc):
    """
    Convert luminosity to flux
    
    Parameters
    ----------
    lum : astropy.units.quantity.Quantity
        luminosity
    distance : astropy.units.quantity.Quantity
        distance
    """
    return (lum/(4*np.pi*distance**2)).to(unitF)

def bb_temp(flux,wl,Tmin=1e-1,Tmax=100000,Tstep=1):
    """
    Returns the temperature of a blackbody with the same relative flux in the 2 bands given as inputs
    """
    flux1,flux2=flux
    wl1,wl2=wl
    from astropy.modeling.physical_models import BlackBody
    x,y=0,0
    #Tmin=1e-1
    #while x<0.01 and y<0.01:
    #    Tmin*=10
    #    #Tmax*=10
    #    bb1=BlackBody(temperature=Tmin*u.K)
    #    x=bb1(wl1)/bb1(bb1.lambda_max)
    #    bb2=BlackBody(temperature=Tmin*u.K)
    #    y=bb2(wl2)/bb2(bb2.lambda_max)
    temp=np.arange(Tmin,Tmax,Tstep)
    BB=BlackBody(temperature=temp*u.K)
    index=np.argmin(np.abs((flux1/flux2)-(BB(wl1)/BB(wl2))))

    return temp[index]

def wmean(values,errors):
    """Calculate the weighted mean of a set of values with an uncertainty"""
    mask=np.invert(np.isnan(values)+np.isnan(errors))
    values=values[mask]
    errors=errors[mask]
    if np.linalg.norm(mask,0)==0.:
        return [0.,np.inf]
    if np.sum(np.divide(1,np.power(errors,2)))==0.:
        return [0.,np.inf]
    return np.average(values,weights=np.divide(1,np.power(errors,2))), np.sqrt(np.divide(1,np.sum(np.divide(1,np.power(errors,2)))))

