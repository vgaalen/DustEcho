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

def get_lightcurve(L_data,T_data,t_data,output_t,shell,grain,FWHM,spaceBins,output_wavelengths=None,distance=False,skiname="initShell.ski",Si=False,prefix='v6',
                   OUTFILES='results/',plot_SED=False,SKIRTpath=None):
    """
    This function uses the SKIRT program (https://skirt.ugent.be/) to simulate a lightcurve for a variable blackbody source with a dust geometry around it. 
    The source is defined by the blackbody temperature (T) and the integrated luminosity. 
    Both these parameters have to be defined in every timestep for which the lightcurve should be determined though a SKIRT simulation. 
    The full lightcurve is determined afterwards through interpolation to all points in output_t.

    Parameters
    ----------
    L_data : array
        Array containing the integrated luminosity of the source in erg/s for every timestep.
    T_data : array
        Array containing the blackbody temperature of the source in K for every timestep.
    t_data : array
        Array containing the time of every timestep in MJD.
    output_t : array
        Array containing the times for which the lightcurve should be determined in MJD (through interpolation).
    shell : array
        Array containing the parameters of the dust geometry; [inner radius, outer radius, exponent, total mass].
    grain : array
        Array containing the parameters of the dust grain size distribution; [max grain size, min grain size, exponent].
    FWHM : float
        The full width at half maximum of the transients lightcurve in days.
    spaceBins : int
        The number of spatial bins used in the SKIRT simulation.
    output_wavelength : array
        Array describing which wavelengths should be in the output; [min wavelength, max wavelength, number of bins].
    distance : float, optional
        The distance to the transient in Mpc. If not given, the distance is set to 1171 Mpc.
    skiname : str, optional
        The name of the SKIRT config file that should be used. Default is 'initShell.ski'.
    Si : float, optional
        The fraction of the dust mass that is made up of silicate grains, the rest is made up of carbonaceous grains. Default is False, which means that only carbonaceous grains are used.
    prefix : str, optional
        The prefix of the output files. Default is 'v6'.
    OUTFILES : str, optional
        The path to the output directory. Default is 'results/'.
    plot_SED : bool, optional
        If True, the SED of the first timestep is plotted. Default is False.
    SKIRTpath : str, optional
        The path to the SKIRT executable. Default is None, which means that the default path is used.
    """

    # Set the simulation parameters
    skifile = sm.SkiFile(skiname)

    # Setup the Config Files
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/dynamicStateOptions/DynamicStateOptions/recipes/GrainSizeDustDestructionRecipe[@FWHM]','FWHM',"{0}".format(FWHM))
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/grid/Sphere1DSpatialGrid/meshRadial/LogMesh[@numBins]','numBins',str(spaceBins))
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/grid/Sphere1DSpatialGrid/meshRadial/LogMesh[@centralBinFraction]','centralBinFraction',str(shell[0]/shell[1]))
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/grid/Sphere1DSpatialGrid[@maxRadius]','maxRadius',str(shell[1]))

    # Setup the wavelength grid
    if output_wavelengths is not None:
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/SEDInstrument/wavelengthGrid/LinWavelengthGrid[@minWavelength]','minWavelength','{0} micron'.format(output_wavelengths[0]))
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/SEDInstrument/wavelengthGrid/LinWavelengthGrid[@maxWavelength]','maxWavelength','{0} micron'.format(output_wavelengths[1]))
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/SEDInstrument/wavelengthGrid/LinWavelengthGrid[@numWavelengths]','numWavelengths','{0}'.format(output_wavelengths[2]))
    
    if distance!=False:
        if distance==1171*u.Mpc:
            skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/cosmology/FlatUniverseCosmology[@redshift]','redshift','0.2343272')
        else:
            from astropy.cosmology import FlatLambdaCDM
            import astropy.cosmology.units as cu
            skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/cosmology/FlatUniverseCosmology[@redshift]','redshift',str(distance.to(cu.redshift,cu.redshift_distance(FlatLambdaCDM(70*u.km/u.s/u.Mpc, 0.3,0.7),kind='luminosity')).value))
    
    # Determine the mass for each dust grain size bin
    a=np.logspace(np.log10(grain[1]),np.log10(grain[0]),10)
    mass=[np.abs(quad(lambda x: x**(-1*grain[2]),a[0],10*a[0])[0])+np.abs(quad(lambda x: x**(-1*grain[2]),10*a[0],np.inf)[0])]
    for i in range(1,len(a)):
        mass.append(quad(lambda x: x**(-1*grain[2]),a[i],a[i-1])[0])
    if Si!=False:
        mass=np.concatenate((np.multiply((1-Si),mass),np.multiply(Si,mass)))
        a=np.concatenate((a,a))
    mass=np.multiply(mass,shell[3]/np.sum(mass))

    # Setup the dust grain size distribution
    for i in range(len(mass)):
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/normalization/MassMaterialNormalization[@mass]'.format(i+1),'mass','{0} Msun'.format(mass[i]))
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/materialMix/FragmentDustMixDecorator/dustMix/ConfigurableDustMix/populations/GrainPopulation/sizeDistribution/SingleGrainSizeDistribution[@size]'.format(i+1),'size','{0} micron'.format(a[i]))
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/geometry/ShellGeometry[@exponent]'.format(i+1),'exponent',str(shell[2]))
        skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/geometry/ShellGeometry[@minRadius]'.format(i+1),'minRadius','{0} pc'.format(shell[0]))
            
    #Run SKIRT
    SED=[]
    #Note: the lightcurve file is based on relative times, therefore the zeropoint corresponds to direct light.
    lightcurve=[]
    temp=[]
    t_peak=t_data[np.argmax(L_data[1])]
    radius=[shell[0]]*len(mass)
    static=False # when True this parameter prevents any further sublimation (or actually prevent non-physical regrowth of dust when the transient falls off)
    for t in tqdm(range(len(t_data)),desc="SKIRT Runs"):
        # Run SKIRT
        lightcurve_t,SED_t,radius_t,temp_t,wavelengths,simulation=runSKIRT(L_data[1,t],T_data[1,t],t_data[t],skifile,OUTFILES=OUTFILES,SKIRTpath=SKIRTpath)
        # Save the output
        if np.sum(lightcurve_t)==0.:
            return 0,0,0,0, simulation
        SED.append(SED_t)
        lightcurve.append(lightcurve_t)
        temp.append(temp_t)

        # Update the dust grain size distribution
        if static==False:
            # If the peak has passed, the dust should no longer sublimate
            if t_data[t]>=t_peak:
                static=True
            # If the peak has not passed, the dust grain size distribution should be updated
            for i in range(len(radius_t)):
                # Only update if the inner radius of a shell has increased (due to finite spatial resolution in SKIRT this is not always the case, even before the transient peak)
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
                # If the inner radius of a shell has increased beyond the outer radius, the shell should be removed (as true removal is challanging in this implementation, the shell is set to a very small radius and mass)
                else:
                    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/geometry/ShellGeometry[@minRadius]'.format(i+1),'minRadius','{0} pc'.format(shell[1]-1e-5))
                    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/materialMix/FragmentDustMixDecorator[@initialDensityFraction]'.format(i+1),'initialDensityFraction','0')
                    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/mediumSystem/MediumSystem/media/GeometricMedium[{0}]/normalization/MassMaterialNormalization[@mass]'.format(i+1),'mass','1e-50 Msun')
                    radius[i]=radius_t[i]
    
    # Create an array with the timesteps as written by SKIRT (changes to the SKIRT config file should be mirrored here)
    t_bin=1*u.day
    timesteps=[]
    for i in np.arange(len(lightcurve[0][0,1:])):
        timesteps.append(i*t_bin.value)

    # Interpolate between timesteps to get emmission for every day
    representation=np.zeros((len(lightcurve[0][0,:]),len(wavelengths),3,len(t_data)+4))
    #w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    #w2=np.argmin(np.abs(wavelengths-wl_w2.value))
    for wl in tqdm(range(len(wavelengths)),leave=True):
        for i in range(1,len(lightcurve[0][0,:])-1): # first column contains wavelengths and final contains overflow, therefore these are discarded
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
            #if i in [1,2,5,10,50,100]:
            #    if wl==w1:
            #        a=splev(np.linspace(t_data[0],t_data[-1],30),(representation[i,w1,0,:],representation[i,w1,1,:],int(representation[i,wl,2,0])))
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
            plt.plot(lightcurve[i][0,:],label=t_data[i])
        plt.xlabel('Days after emission')
        plt.ylabel('Flux (Jy)')
        plt.yscale('log')
        plt.savefig('plots/'+prefix+'t_bin_compare_w1.pdf')

    return lc_total*u.Jy,wavelengths,temp,radius, simulation

def runSKIRT(L,T,MJD,skifile,OUTFILES="",SKIRTpath=None):
    """
    This function uses the SKIRT program (https://skirt.ugent.be/) to simulate one timestep in a lightcurve for a variable blackbody source with a dust geometry around it. 
    The source is defined by the blackbody temperature (T) and the integrated luminosity. 

    Parameters
    ----------
    L : float
        The integrated luminosity of the source in erg/s.
    T : float
        The blackbody temperature of the source in K.
    MJD : float
        The absolute time of the requested lightcurve in MJD.
    skifile : str
        The SKIRT config file that should be used.
    OUTFILES : str, optional
        The path to the output directory. Default is ''.
    SKIRTpath : str, optional
        The path to the SKIRT executable. Default is None, which means that the default path is used.
    """

    # Initialize skirt
    if SKIRTpath==None:
        skirt = sm.Skirt(path="SKIRT/release/SKIRT/main/skirt")
    else:
        skirt = sm.Skirt(path=SKIRTpath)

    # Update the SKIRT config file with this timestep's parameters
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/sourceSystem/SourceSystem/sources/PointSource/normalization/IntegratedLuminosityNormalization[@integratedLuminosity]','integratedLuminosity',"{0} erg/s".format(L))
    skifile.setStringAttribute('//skirt-simulation-hierarchy/MonteCarloSimulation/sourceSystem/SourceSystem/sources/PointSource/sed/BlackBodySED[@temperature]','temperature',"{0:.2E} K".format(T))

    # Create the output directory and clear it if necessary
    if os.path.isdir(OUTFILES+str(int(MJD))+"/")==False:
        os.mkdir(OUTFILES+str(int(MJD))+"/")
    else:
        files = glob.glob(OUTFILES+str(int(MJD))+'/*')
        for f in files:
            os.remove(f)

    # Run SKIRT
    skifile.saveTo(OUTFILES+str(int(MJD))+"/run.ski")
    simulation = skirt.execute(OUTFILES+str(int(MJD))+"/run.ski",outDirPath=OUTFILES+str(int(MJD))+"/", console='brief')

    if os.path.isfile(OUTFILES+str(int(MJD))+'/run_instrument1_sed.dat')==False or os.path.isfile(OUTFILES+str(int(MJD))+'/run_instrument1_lc.dat')==False:
        print("ERROR: SKIRT did not produce output files")
        return 0,0,0,0,0, simulation
    # Load the output
    SED=np.loadtxt(OUTFILES+str(int(MJD))+'/run_instrument1_sed.dat')
    lightcurve=np.loadtxt(OUTFILES+str(int(MJD))+'/run_instrument1_lc.dat')

    # Find the inner radius of every dust shell and the temperature at that radius
    temp=[]
    radius=[]
    i=0
    while os.path.isfile(OUTFILES+str(int(MJD))+'/run_medium-temperature_{0}_T_xy.fits'.format(i)):
        datafile=fits.open(OUTFILES+str(int(MJD))+'/run_medium-temperature_{0}_T_xy.fits'.format(i))[0]
        temperature=datafile.data
        size=int(np.ceil(0.5*len(temperature[0,:])))
        # Handle the no sublimation case
        if temperature[size,size]!=0.:
            radius.append(0)
            temp.append(temperature[size,size])
            i+=1
        # Handle the full sublimation case
        elif np.linalg.norm(temperature)==0.:
            radius.append(size*datafile.header['CDELT1'])
            temp.append(0.)
            i+=1
        # All other cases
        else:
            for r in range(size):
                if temperature[size,size+r]!=0.:
                    radius.append(r*datafile.header['CDELT1'])
                    temp.append(temperature[size,size+r])
                    i+=1
                    break
    
    wavelengths=SED[:,0]
    return lightcurve,SED,radius,temp,wavelengths,simulation
