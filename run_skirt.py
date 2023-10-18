import os
import glob
import astropy.units as u
from astropy.io import fits
import numpy as np
from scipy.interpolate import splrep, splev
from scipy.integrate import quad
from tqdm import tqdm
import matplotlib.pyplot as plt

import sys
sys.path.append('./PTS')
import pts.simulation as sm
import pts.do
pts.do.initializePTS()

wl_w1=3.368*u.um
wl_w2=4.618*u.um

def get_lightcurve(L_data,T_data,t_data,output_t,shell,grain,FWHM,spaceBins,distance=False,geometries=['initShell.fits','Shell.fits'],skiname="initShell.ski",Si=False,prefix='v6',OUTFILES='results/',plot_SED=False,SKIRTpath=None):
    
    #set the simulation parameters
    skifile = sm.SkiFile(skiname)

    #px_scale=shell[0]/spaceBins #for the geometry .fits files

    #Setup the Pre- and Post-Peak Config Files
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/dynamicStateOptions/DynamicStateOptions/recipes/GrainSizeDustDestructionRecipe[@FWHM]','FWHM',"{0}".format(FWHM))
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/grid/Sphere1DSpatialGrid/meshRadial/LogMesh[@numBins]','numBins',str(spaceBins))
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/grid/Sphere1DSpatialGrid/meshRadial/LogMesh[@centralBinFraction]','centralBinFraction',str(shell[0]/shell[1]))
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/grid/Sphere1DSpatialGrid[@maxRadius]','maxRadius',str(shell[1]))
    
    
    if distance!=False:
        if distance==1171*u.Mpc:
            skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/cosmology/FlatUniverseCosmology[@redshift]','redshift','0.2343272')
        else:
            from astropy.cosmology import FlatLambdaCDM
            import astropy.cosmology.units as cu
            skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/cosmology/FlatUniverseCosmology[@redshift]','redshift',str(distance.to(cu.redshift,cu.redshift_distance(FlatLambdaCDM(70*u.km/u.s/u.Mpc, 0.3,0.7),kind='luminosity')).value))
    
    a=np.logspace(np.log10(grain[1]),np.log10(grain[0]),10)
    mass=[np.abs(quad(lambda x: x**(-1*grain[2]),a[0],10*a[0])[0])+np.abs(quad(lambda x: x**(-1*grain[2]),10*a[0],np.inf)[0])]
    for i in range(1,len(a)):
        mass.append(quad(lambda x: x**(-1*grain[2]),a[i],a[i-1])[0])
    
    if Si!=False:
        mass=np.concatenate((np.multiply((1-Si),mass),np.multiply(Si,mass)))
        a=np.concatenate((a,a))

    mass=np.multiply(mass,shell[3]/np.sum(mass))

    for i in range(len(mass)):
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/normalization/MassMaterialNormalization[@mass]'.format(i+1),'mass','{0} Msun'.format(mass[i]))
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/materialMix/FragmentDustMixDecorator/dustMix/ConfigurableDustMix/populations/GrainPopulation/sizeDistribution/SingleGrainSizeDistribution[@size]'.format(i+1),'size','{0} micron'.format(a[i]))
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/geometry/ShellGeometry[@exponent]'.format(i+1),'exponent',str(shell[2]))
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/geometry/ShellGeometry[@minRadius]'.format(i+1),'minRadius','{0} pc'.format(shell[0]))
            
    #Run SKIRT
    SED=[]
    #Note: the lightcurve file is based on relative times, D=0 corresponds to direct light.
    lightcurve=[]
    temp=[]
    t_peak=t_data[np.argmax(L_data[1])]
    radius=[shell[0]]*len(mass)
    #radius_t=[]
    static=False
    for t in tqdm(range(len(t_data)),desc="SKIRT Runs"):
        lightcurve_t,SED_t,radius_t,temp_t,wavelengths=runSKIRT(L_data[1,t],T_data[1,t],t_data[t],skifile,OUTFILES=OUTFILES,SKIRTpath=SKIRTpath)
        SED.append(SED_t)
        lightcurve.append(lightcurve_t)
        temp.append(temp_t)
        if static==False:
            if t_data[t]>=t_peak:
                static=True

            for i in range(len(radius_t)):
                if radius_t[i]<=shell[0]:
                    pass
                elif radius_t[i]<radius[i]:
                    pass
                elif radius[i]<radius_t[i] and radius_t[i]<shell[1]:
                    mass_fraction=quad(lambda r: r**(-1*shell[2]), radius_t[i], shell[1])[0]/quad(lambda r: r**(-1*shell[2]), radius[i], shell[1])[0]
                    if mass_fraction<1.:
                        mass=np.multiply(mass_fraction,mass)
                    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/geometry/ShellGeometry[@minRadius]'.format(i+1),'minRadius','{0} pc'.format(radius_t[i]))
                    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/normalization/MassMaterialNormalization[@mass]'.format(i+1),'mass','{0} Msun'.format(mass[i]))
                    radius[i]=radius_t[i]
                else:
                    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/geometry/ShellGeometry[@minRadius]'.format(i+1),'minRadius','{0} pc'.format(shell[1]-1e-5))
                    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/materialMix/FragmentDustMixDecorator[@initialDensityFraction]'.format(i+1),'initialDensityFraction','0')
                    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/normalization/MassMaterialNormalization[@mass]'.format(i+1),'mass','1e-50 Msun')
                    radius[i]=radius_t[i]
    
    #Timesteps of the SKIRT output as defined in the config files
    t_bin=1*u.day
    timesteps=[]
    for i in np.arange(len(lightcurve[0][0,1:])):
        timesteps.append(i*t_bin.value)

    #Interpolate between timesteps to get emmission for every day
    representation=np.zeros((len(lightcurve[0][0,:]),len(wavelengths),3,len(t_data)+4))
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))
    for wl in tqdm(range(len(wavelengths)),leave=True):
        for i in range(1,len(lightcurve[0][0,:])-1):#first column contains wavelengths and final contains overflow, therefore these are discarded
            y=[]
            time=[]
            for t in range(len(t_data)):
                if lightcurve[t][wl,i]!=0.:
                    y.append(lightcurve[t][wl,i])
                    time.append(t_data[t])
            if np.linalg.norm(y)==0.:
                pass
            if len(y)>3:
                a=splrep(time,y,k=1)
            elif len(y)>1:
                a=splrep(time,y,k=1)
            else:
                continue
            representation[i,wl,0,:len(a[0])]=a[0]
            representation[i,wl,1,:len(a[1])]=a[1]
            representation[i,wl,2,0]=a[2]
            if i in [1,2,5,10,50,100]:
                if wl==w1:
                    a=splev(np.linspace(t_data[0],t_data[-1],30),(representation[i,w1,0,:],representation[i,w1,1,:],int(representation[i,wl,2,0])))
    lc_total=np.zeros((len(output_t),len(wavelengths)))
    for i in tqdm(range(len(output_t)),desc="Compile Lightcurve",leave=False,position=2):
        for t in np.arange(0,int(np.floor(output_t[i]-t_data[0])),1):
            for wl in range(len(wavelengths)):
                if np.linalg.norm(representation[int(t/1)+1,wl,0,:])!=0:
                    lc_total[i,wl]+=splev(output_t[i]-t,(representation[int(t/1)+1,wl,0,:],representation[int(t/1)+1,wl,1,:],int(representation[int(t/1)+1,wl,2,0])))

    if plot_SED==True:
        colors=['black','grey','red','orange','green','lime','blue','navy','pink','purple','brown']*10

        plt.figure()
        for i in [1,5,10,100,500,1000]:
            plt.plot(wavelengths,lightcurve[1][:,i],label=str(i))
        plt.legend()
        plt.xlabel('Wavelength (micron)')
        plt.xscale('log')
        plt.yscale('log')
        #plt.ylim(1e-30,1)
        plt.savefig('plots/'+prefix+'SEDcompare_tbin.pdf')

        fig = plt.figure()
        ax = plt.subplot(111)
        for t in range(len(lightcurve)):
            ax.plot(range(len(lightcurve[0][0,1:]))*t_bin+(t_data[t]-t_data[0])*u.day,np.sum(lightcurve[t][:,1:],axis=0),label="MJD="+str(t_data[t]))
        ylim=ax.get_ylim()
        #ax.plot([2*(inner*u.pc/c.c).to(u.day).value]*2,ylim,linestyle='dashed')
        ax.set_yscale('log')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_xlabel('Time since emission')
        ax.set_ylabel('Flux (Jy)')
        ax.set_title('Lightcurves of SKIRT instances compared')
        fig.savefig('plots/'+prefix+'t_bin_compare.pdf')

        plt.figure()
        for i in range(len(t_data)):
            plt.plot(lightcurve[i][w1,:],label=t_data[i])
        plt.xlabel('Days after emission')
        plt.ylabel('Flux (Jy)')
        plt.yscale('log')
        plt.savefig('plots/'+prefix+'t_bin_compare_w1.pdf')

    return lc_total*u.Jy,wavelengths,temp,radius

def runSKIRT(L,T,MJD,skifile,packets=1e6,OUTFILES="",savefiles=None,SKIRTpath=None):
    """
    This function uses the SKIRT program (https://skirt.ugent.be/) to simulate a lightcurve for a variable blackbody source with a dust geometry around it. The source is defined by the blackbody temperature (T) and the luminosity in the Johnson V band. The source has to be defined in every timestep for which the lightcurve should be determined (note that datapoints before 1ltt has passed are a lower limit as light emitted earlier is not considered).
    alpha defines for quickly the relevance of a timestep's emission' falls off after ltt has passed
    """

    #initialize skirt
    if SKIRTpath==None:
        skirt = sm.Skirt(path="SKIRT/release/SKIRT/main/skirt")
    else:
        skirt = sm.Skirt(path=SKIRTpath)

    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/sourceSystem/SourceSystem/sources/PointSource/normalization/IntegratedLuminosityNormalization[@integratedLuminosity]','integratedLuminosity',"{0} erg/s".format(L))
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/sourceSystem/SourceSystem/sources/PointSource/sed/BlackBodySED[@temperature]','temperature',"{0:.2E} K".format(T))

    if os.path.isdir(OUTFILES+str(int(MJD))+"/")==False:
        os.mkdir(OUTFILES+str(int(MJD))+"/")
    else:
        files = glob.glob(OUTFILES+str(int(MJD))+'/*')
        for f in files:
            os.remove(f)

    skifile.saveTo(OUTFILES+str(int(MJD))+"/torus.ski")
    simulation=skirt.execute(OUTFILES+str(int(MJD))+"/torus.ski",outDirPath=OUTFILES+str(int(MJD))+"/", console='brief')

    SED=np.loadtxt(OUTFILES+str(int(MJD))+'/torus_instrument1_sed.dat')
    lightcurve=np.loadtxt(OUTFILES+str(int(MJD))+'/torus_instrument1_lc.dat')

    temp=[]
    radius=[]
    i=0
    while os.path.isfile(OUTFILES+str(int(MJD))+'/torus_medium-temperature_{0}_T_xy.fits'.format(i)):
        datafile=fits.open(OUTFILES+str(int(MJD))+'/torus_medium-temperature_{0}_T_xy.fits'.format(i))[0]
        temperature=datafile.data
        size=int(np.ceil(0.5*len(temperature[0,:])))
        if temperature[size,size]!=0.:
            radius.append(0)
            temp.append(temperature[size,size])
            i+=1
        elif np.linalg.norm(temperature)==0.:
            radius.append(size*datafile.header['CDELT1'])
            temp.append(0.)
            i+=1
        else:
            for r in range(size):
                if temperature[size,size+r]!=0.:
                    radius.append(r*datafile.header['CDELT1'])
                    temp.append(temperature[size,size+r])
                    i+=1
                    break
    
    wavelengths=SED[:,0]
    return lightcurve,SED,radius,temp,wavelengths
