import matplotlib.pyplot as plt
import matplotlib.colors as mpl_c
import matplotlib.cm as cm
import numpy as np
import astropy.constants as c
import astropy.units as u
import os
import seaborn as sns

wl_w1=3.368*u.um
wl_w2=4.618*u.um

def plot(data,W1,W2,values=[],label='',save='test.pdf'):
    """
    Plot one or multiple simulated lightcurves over the WISE observations.
    If multiple lightcurves are given, they are plotted with a colorbar.
    The color of each curve is determined by the value in the values array.
    The meaning of the numbers in the values array should be given by the label, which is printed on the colorbar.

    Parameters
    ----------
    data : list of arrays
        list of arrays with MJDs in the first column and luminosities in the second
    W1 : array
        array with MJDs in the first column, luminosities in the second and errors in the third
    W2 : array
        array with MJDs in the first column, luminosities in the second and errors in the third
    label : string
        label for the colorbar
    save : string
        filename for the plot
    """
    fig=plt.figure()
    ax=fig.get_axes

    cmap1=mpl_c.Colormap('blues')
    cmap2=mpl_c.Colormap('reds')

    plt.errorbar(W1[0],(c.c/wl_w1).to(u.Hz).value*W1[1],yerr=np.abs((c.c/wl_w1).to(u.Hz).value*W1[2]),fmt='.',capsize=3,color='red',label='IR-WISE/W1')
    plt.errorbar(W2[0],(c.c/wl_w2).to(u.Hz).value*W2[1],yerr=np.abs((c.c/wl_w2).to(u.Hz).value*W2[2]),fmt='.',capsize=3,color='darkred',label='IR-WISE/W2')
    for i in range(len(data)):
        plt.plot(data[i][0],(c.c/wl_w1).to(u.Hz).value*data[i][1],color=cmap1[values[i]])
        plt.plot(data[i][0],(c.c/wl_w2).to(u.Hz).value*data[i][2],color=cmap2[values[i]])

    cbar1=fig.colorbar(cmap1,ax=ax)
    cbar2=fig.colorbar(cmap2,ax=ax)
    cbar2.set_label(label)

    ax.set_ylabel(r'$\nu L_{\nu} [erg/s]$')
    ax.set_xlabel('MJD')

    fig.savefig(save)

def coveringfactorPlot(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))
    for i in range(len(amax)):
        for j in range(len(Lbol)):
            cf=np.max(W1[1])/np.max(Grid1[j,i,:]) #covering factor
            fig=plt.figure()
            ax=plt.subplot(111)
            ax.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[j,i,:]*cf,label='W1 fit',color='slategrey')
            ax.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[j,i,:]*cf,label='W2 fit',color='indianred')
            ax.errorbar(W1[0],(c.c/wl_w1).to(u.Hz).value*W1[1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2],fmt='.',capsize=3,label='W1 data',color='slategrey')
            ax.errorbar(W2[0],(c.c/wl_w2).to(u.Hz).value*W2[1],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2],fmt='.',capsize=3,label='WISE data',color='indianred')
            ax.plot([L[0][np.argmax(L[1,:])]]*2,ax.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
            #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
            ax.set_yscale('log')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            ax.set_xlabel('MJD')
            #ax.set_ylabel('ν L_ν [erg/s]')
            ax.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
            #fig.suptitle('ASASSN15lh Dust Echo Model')
            ax.set_title(target+' Dust Echo Model')
            fig.savefig(PLOTFILES+'lightcurve-'+str(Lbol[j])+'-'+str(amax[i])+'.pdf')

def compare_amax(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))
    
    #cmap1=mpl_c.Colormap('blues')
    #cmap2=mpl_c.Colormap('reds')
    
    fig=plt.figure(figsize=[7.5, 4])
    ax=plt.subplot(111)
    
    color=np.log10(amax)+(2-np.min(np.log10(amax)))
    color/=np.max(color)
    
    verybest_loss=np.inf
    
    for i in range(len(amax)):
        best_loss=np.inf
        best=False
        for j in range(len(Lbol)):
            if np.linalg.norm(Grid1[j,i,:])==0.:
                cf=1
            else:
                cf=np.max(W1[1])/np.max(Grid1[j,i,:]) #covering factor
                if cf>1:
                    cf=1

            loss=0
            for t in range(len(W1[0])):
                loss+=(W1[1,t]-cf*Grid1[j,i,np.argmin(np.abs(timesteps-W1[0,t]))])**2
            for t in range(len(W2[0])):
                    loss+=(W2[1,t]-cf*Grid2[j,i,np.argmin(np.abs(timesteps-W2[0,t]))])**2
            if loss<best_loss:
                best_loss=loss
                best=j
                if loss<verybest_loss:
                    verybest_loss=loss
                    verybest=[i,j]
        
        j=best
        print('Lbol: ',Lbol[j])
        if np.linalg.norm(Grid1[j,i,:])==0.:
            cf=1
        else:
            cf=np.max(W1[1])/np.max(Grid1[j,i,:]) #covering factor
            if cf>1:
                cf=1

        ax.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[j,i,:]*cf,label=r'W1, '+np.format_float_scientific(amax[i], unique=False, exp_digits=1,precision=2),color=(0,0.5-0.5*color[i],0.5,color[i]))#color=cmap1(amax[i]))
        ax.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[j,i,:]*cf,label=r'W2, '+np.format_float_scientific(amax[i], unique=False, exp_digits=1,precision=2),color=(0.5,0.5-0.5*color[i],0.1,color[i]))#color=cmap2(amax[i]))
    
    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    ax.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
    ax.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='WISE data',color='indianred')
    ax.plot([L[0][np.argmax(L[1,:])]]*2,ax.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
    ax.set_yscale('log')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'Band, $a_{\mathrm{min}}$')
    ax.set_xlabel('MJD')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax.set_title(target+' Dust Echo Model')
    fig.savefig(PLOTFILES+'compare_amax.pdf')
    
    print(target+' - very best: ')
    print('Lbol: ',Lbol[verybest[1]])
    print('amax: ',amax[verybest[0]])

def compare_amax2(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))
    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    
    #cmap1=mpl_c.Colormap('blues')
    #cmap2=mpl_c.Colormap('reds')
    
    fig=plt.figure(figsize=[10, 4])
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)
    
    vmin=np.log10(np.min(amax))
    vmax=np.log10(np.max(amax))
    #-0.1*(vmax-vmin)
    #cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax+0.1*(vmax-vmin)),cmap='cubehelix')
    #cmap1=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True))
    #cmap2=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True))
    cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax),cmap=sns.color_palette("Spectral", as_cmap=True))
    verybest_loss=np.inf
    
    for j in range(len(amax)):
        best_loss=np.inf
        best=[0]*5
        for i in range(len(Lbol)):
            for k in range(len(grainpowerlaw)):
                for l in range(len(total_mass)):
                    for m in range(len(alpha)):
                        for n in range(len(Si)):

                            if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
                                cf=1
                            else:
                                cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
                                if cf>1:
                                    cf=1

                            loss=0
                            for t in range(len(W1[0][maskw1])):
                                if cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0][maskw1][t]))]!=0.:
                                    loss+=(W1[1][maskw1][t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0][maskw1][t]))])**2/(cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0][maskw1][t]))])
                                else:
                                    loss+=(W1[1][maskw1][t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0][maskw1][t]))])**2
                            for t in range(len(W2[0][maskw2])):
                                if cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0][maskw2][t]))]!=0.:
                                    loss+=(W2[1][maskw2][t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0][maskw2][t]))])**2/(cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0][maskw2][t]))])
                                else:
                                    loss+=(W2[1][maskw2][t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0][maskw2][t]))])**2
                            if loss<best_loss:
                                best_loss=loss
                                best=[i,k,l,m,n]
                                if loss<verybest_loss:
                                    verybest_loss=loss
                                    verybest=[i,j,k,l,m,n]

        [i,k,l,m,n]=best
        print('best ',amax[j],' - ',best)
        print('loss: ',loss)
        if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
            cf=1
        else:
            cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
            if cf>1:
                cf=1
        
        ax1.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[i,j,k,l,m,n,:]*cf,label=np.format_float_scientific(amax[j], unique=False, exp_digits=1,precision=1),color=cmap.to_rgba(np.log10(amax[j])))#color=cmap1(amax[i]))
        ax2.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[i,j,k,l,m,n,:]*cf,label=np.format_float_scientific(amax[j], unique=False, exp_digits=1,precision=1),color=cmap.to_rgba(np.log10(amax[j])))#color=cmap2(amax[i]))
                               
    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    ax1.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
    ax2.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='WISE data',color='indianred')
    ax1.plot([L[0][np.argmax(L[1,:])]]*2,ax1.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    ax2.plot([L[0][np.argmax(L[1,:])]]*2,ax2.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'$a_{\mathrm{max}}$')
    ax1.set_xlabel('MJD')
    ax2.set_xlabel('MJD')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax1.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax1.set_title('W1')
    ax2.set_title('W2')

    fig.subplots_adjust(right=0.85)

    fig.suptitle(target+' Dust Echo Model')
    fig.savefig(PLOTFILES+'compare_amax(1).pdf')
    
    print(target+' - very best: ')
    print('amax: ',amax[verybest[1]])
    print(verybest)
    print('loss: ',verybest_loss)



def compare_amax3(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))
    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])

    #cmap1=mpl_c.Colormap('blues')
    #cmap2=mpl_c.Colormap('reds')

    fig=plt.figure(figsize=[10, 4])
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)

    vmin=np.log10(np.min(amax))
    vmax=np.log10(np.max(amax))
    #-0.1*(vmax-vmin)
    #cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax+0.1*(vmax-vmin)),cmap='cubehelix')
    #cmap1=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True))
    #cmap2=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True))
    cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax),cmap=sns.color_palette("Spectral", as_cmap=True))
    best_loss=np.inf

    for j in range(len(amax)):
        for i in range(len(Lbol)):
            for k in range(len(grainpowerlaw)):
                for l in range(len(total_mass)):
                    for m in range(len(alpha)):
                        for n in range(len(Si)):

                            if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
                                cf=1
                            else:
                                cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
                                if cf>1:
                                    cf=1

                            loss=0
                            for t in range(len(W1[0][maskw1])):
                                if cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0][maskw1][t]))]!=0.:
                                    loss+=(W1[1][maskw1][t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0][maskw1][t]))])**2/(cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0][maskw1][t]))])
                                else:
                                    loss+=(W1[1][maskw1][t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0][maskw1][t]))])**2
                            for t in range(len(W2[0][maskw2])):
                                if cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0][maskw2][t]))]!=0.:
                                    loss+=(W2[1][maskw2][t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0][maskw2][t]))])**2/(cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0][maskw2][t]))])
                                else:
                                    loss+=(W2[1][maskw2][t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0][maskw2][t]))])**2
                            if loss<verybest_loss:
                                best_loss=loss
                                best=[i,j,k,l,m,n]

    for j in range(len(amax)):
        if np.linalg.norm(Grid1[best[0],j,best[2],best[3],best[4],best[5],:])==0.:
            cf=1
        else:
            cf=np.max(W1[1])/np.max(Grid1[best[0],j,best[2],best[3],best[4],best[5],:]) #covering factor
            if cf>1:
                cf=1

        ax1.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[best[0],j,best[2],best[3],best[4],best[5],:]*cf,label=np.format_float_scientific(amax[j], unique=False, exp_digits=1,precision=1),color=cmap.to_rgba(np.log10(amax[j])))#color=cmap1(amax[i]))
        ax2.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[best[0],j,best[2],best[3],best[4],best[5],:]*cf,label=np.format_float_scientific(amax[j], unique=False, exp_digits=1,precision=1),color=cmap.to_rgba(np.log10(amax[j])))#color=cmap2(amax[i]))

    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    ax1.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
    ax2.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='WISE data',color='indianred')
    ax1.plot([L[0][np.argmax(L[1,:])]]*2,ax1.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    ax2.plot([L[0][np.argmax(L[1,:])]]*2,ax2.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'$a_{\mathrm{max}}$')
    ax1.set_xlabel('MJD')
    ax2.set_xlabel('MJD')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax1.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax1.set_title('W1')
    ax2.set_title('W2')

    fig.subplots_adjust(right=0.85)

    fig.suptitle(target+' Dust Echo Model')
    fig.savefig(PLOTFILES+'compare_amax(1).pdf')

    print(target+' - very best: ')
    print('amax: ',amax[best[1]])
    print(verybest)
    print('loss: ',best_loss)



def compare_grainpowerlaw(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,gpl,total_mass,alpha,Si,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))

    #cmap1=mpl_c.Colormap('blues')
    #cmap2=mpl_c.Colormap('reds')

    fig=plt.figure(figsize=[10, 4])
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)

    vmin=np.log10(np.min(gpl))
    vmax=np.log10(np.max(gpl))
    #cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax+0.1*(vmax-vmin)),cmap='cubehelix')
    #cmap1=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True))
    #cmap2=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True))
    cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax),cmap=sns.color_palette("Spectral", as_cmap=True))
    verybest_loss=np.inf

    for k in range(len(gpl)):
        best_loss=np.inf
        best=[0]*5
        for i in range(len(Lbol)):
            for j in range(len(amax)):
                for l in range(len(total_mass)):
                    for m in range(len(alpha)):
                        for n in range(len(Si)):

                            if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
                                cf=1
                            else:
                                cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
                                if cf>1:
                                    cf=1

                            loss=0
                            for t in range(len(W1[0])):
                                loss+=(W1[1,t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0,t]))])**2
                            for t in range(len(W2[0])):
                                loss+=(W2[1,t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0,t]))])**2
                            if loss<best_loss:
                                best_loss=loss
                                best=[i,j,l,m,n]
                                if loss<verybest_loss:
                                    verybest_loss=loss
                                    verybest=[i,j,k,l,m,n]

        [i,j,l,m,n]=best
        print('best ',gpl[k],' - ',best)
        if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
            cf=1
        else:
            cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
            if cf>1:
                cf=1

        ax1.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[i,j,k,l,m,n,:]*cf,label=gpl[k],color=cmap.to_rgba(np.log10(gpl[k])))#color=cmap1(amax[i]))
        ax2.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[i,j,k,l,m,n,:]*cf,label=gpl[k],color=cmap.to_rgba(np.log10(gpl[k])))#color=cmap2(amax[i]))

    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    ax1.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
    ax2.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='WISE data',color='indianred')
    ax1.plot([L[0][np.argmax(L[1,:])]]*2,ax1.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    ax2.plot([L[0][np.argmax(L[1,:])]]*2,ax2.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'$\rho\propto a^{-\alpha}$')
    ax1.set_xlabel('MJD')
    ax2.set_xlabel('MJD')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax1.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax1.set_title('W1')
    ax2.set_title('W2')

    fig.subplots_adjust(right=0.85)

    fig.suptitle(target+' Dust Echo Model')
    fig.savefig(PLOTFILES+'compare_gpl.pdf')

    print(target+' - very best: ')
    print('gpl: ',gpl[verybest[2]])
    print(verybest)



def compare_grainpowerlaw2(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,gpl,total_mass,alpha,Si,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))

    #cmap1=mpl_c.Colormap('blues')
    #cmap2=mpl_c.Colormap('reds')

    fig=plt.figure(figsize=[10, 4])
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)

    vmin=np.log10(np.min(gpl))
    vmax=np.log10(np.max(gpl))
    #cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax+0.1*(vmax-vmin)),cmap='cubehelix')
    #cmap1=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True))
    #cmap2=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True))
    cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax),cmap=sns.color_palette("Spectral", as_cmap=True))
    best_loss=np.inf

    for k in range(len(gpl)):
        for i in range(len(Lbol)):
            for j in range(len(amax)):
                for l in range(len(total_mass)):
                    for m in range(len(alpha)):
                        for n in range(len(Si)):

                            if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
                                cf=1
                            else:
                                cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
                                if cf>1:
                                    cf=1

                            loss=0
                            for t in range(len(W1[0])):
                                loss+=(W1[1,t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0,t]))])**2
                            for t in range(len(W2[0])):
                                loss+=(W2[1,t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0,t]))])**2
                            if loss<best_loss:
                                best_loss=loss
                                best=[i,j,l,m,n]

    for k in range(len(gpl)):
        if np.linalg.norm(Grid1[best[0],best[1],k,best[3],best[4],best[5],:])==0.:
            cf=1
        else:
            cf=np.max(W1[1])/np.max(Grid1[best[0],best[1],k,best[3],best[4],best[5],:]) #covering factor
            if cf>1:
                cf=1

        ax1.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[best[0],best[1],k,best[3],best[4],best[5],:]*cf,label=gpl[k],color=cmap.to_rgba(np.log10(gpl[k])))#color=cmap1(amax[i]))
        ax2.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[best[0],best[1],k,best[3],best[4],best[5],:]*cf,label=gpl[k],color=cmap.to_rgba(np.log10(gpl[k])))#color=cmap2(amax[i]))

    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    ax1.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
    ax2.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='WISE data',color='indianred')
    ax1.plot([L[0][np.argmax(L[1,:])]]*2,ax1.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    ax2.plot([L[0][np.argmax(L[1,:])]]*2,ax2.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'$\rho\propto a^{-\alpha}$')
    ax1.set_xlabel('MJD')
    ax2.set_xlabel('MJD')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax1.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax1.set_title('W1')
    ax2.set_title('W2')

    fig.subplots_adjust(right=0.85)

    fig.suptitle(target+' Dust Echo Model')
    fig.savefig(PLOTFILES+'compare_gpl.pdf')

    print(target+' - very best: ')
    print('gpl: ',gpl[best[2]])
    print(best)



def compare_sillicate(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))

    #cmap1=mpl_c.Colormap('blues')
    #cmap2=mpl_c.Colormap('reds')

    fig=plt.figure(figsize=[10, 4])
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)

    vmin=np.min(Si)
    vmax=np.max(Si)
    #cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax+0.1*(vmax-vmin)),cmap='cubehelix')
    #cmap1=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True))
    #cmap2=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True))
    cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax),cmap=sns.color_palette("Spectral", as_cmap=True))
    verybest_loss=np.inf

    for n in range(len(Si)):
        best_loss=np.inf
        best=[0]*5
        for i in range(len(Lbol)):
            for j in range(len(amax)):
                for k in range(len(grainpowerlaw)):
                    for l in range(len(total_mass)):
                        for m in range(len(alpha)):

                            if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
                                cf=1
                            else:
                                cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
                                if cf>1:
                                    cf=1

                            loss=0
                            for t in range(len(W1[0])):
                                loss+=(W1[1,t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0,t]))])**2
                            for t in range(len(W2[0])):
                                loss+=(W2[1,t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0,t]))])**2
                            if loss<best_loss:
                                best_loss=loss
                                best=[i,j,k,l,m]
                                if loss<verybest_loss:
                                    verybest_loss=loss
                                    verybest=[i,j,k,l,m,n]

        [i,j,k,l,m]=best
        print('best ',Si[n],' - ',best)
        if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
            cf=1
        else:
            cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
            if cf>1:
                cf=1

        ax1.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[i,j,k,l,m,n,:]*cf,label=Si[n],color=cmap.to_rgba(Si[n]))#color=cmap1(amax[i]))
        ax2.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[i,j,k,l,m,n,:]*cf,label=Si[n],color=cmap.to_rgba(Si[n]))#color=cmap2(amax[i]))

    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    ax1.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
    ax2.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='WISE data',color='indianred')
    ax1.plot([L[0][np.argmax(L[1,:])]]*2,ax1.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    ax2.plot([L[0][np.argmax(L[1,:])]]*2,ax2.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'$Si fraction$')
    ax1.set_xlabel('MJD')
    ax2.set_xlabel('MJD')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax1.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax1.set_title('W1')
    ax2.set_title('W2')

    fig.subplots_adjust(right=0.85)

    fig.suptitle(target+' Dust Echo Model')
    fig.savefig(PLOTFILES+'compare_Si.pdf')

    print(target+' - very best: ')
    print('Si: ',Si[verybest[5]])
    print(verybest)



def compare_alpha(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))

    #cmap1=mpl_c.Colormap('blues')
    #cmap2=mpl_c.Colormap('reds')

    fig=plt.figure(figsize=[10, 4])
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)

    vmin=np.log10(np.min(alpha))
    vmax=np.log10(np.max(alpha))
    #cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax+0.1*(vmax-vmin)),cmap='cubehelix')
    #cmap1=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True))
    #cmap2=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True))
    cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax),cmap=sns.color_palette("Spectral", as_cmap=True))
    verybest_loss=np.inf

    for m in range(len(alpha)):
        best_loss=np.inf
        best=[0]*5
        for i in range(len(Lbol)):
            for j in range(len(amax)):
                for k in range(len(grainpowerlaw)):
                    for l in range(len(total_mass)):
                        for n in range(len(Si)):

                            if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
                                cf=1
                            else:
                                cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
                                if cf>1:
                                    cf=1

                            loss=0
                            for t in range(len(W1[0])):
                                loss+=(W1[1,t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0,t]))])**2
                            for t in range(len(W2[0])):
                                loss+=(W2[1,t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0,t]))])**2
                            if loss<best_loss:
                                best_loss=loss
                                best=[i,j,k,l,n]
                                if loss<verybest_loss:
                                    verybest_loss=loss
                                    verybest=[i,j,k,l,m,n]

        [i,j,k,l,n]=best
        print('best ',alpha[m],' - ',best)
        if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
            cf=1
        else:
            cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
            if cf>1:
                cf=1

        ax1.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[i,j,k,l,m,n,:]*cf,label=alpha[m],color=cmap.to_rgba(np.log10(alpha[m])))#color=cmap1(amax[i]))
        ax2.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[i,j,k,l,m,n,:]*cf,label=alpha[m],color=cmap.to_rgba(np.log10(alpha[m])))#color=cmap2(amax[i]))

    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    ax1.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
    ax2.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='WISE data',color='indianred')
    ax1.plot([L[0][np.argmax(L[1,:])]]*2,ax1.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    ax2.plot([L[0][np.argmax(L[1,:])]]*2,ax2.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'$\rho\propto r^{-\alpha}$')
    ax1.set_xlabel('MJD')
    ax2.set_xlabel('MJD')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax1.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax1.set_title('W1')
    ax2.set_title('W2')

    fig.subplots_adjust(right=0.85)

    fig.suptitle(target+' Dust Echo Model')
    fig.savefig(PLOTFILES+'compare_alpha.pdf')

    print(target+' - very best: ')
    print('alpha: ',alpha[verybest[4]])
    print(verybest)




def compare_mass(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))

    #cmap1=mpl_c.Colormap('blues')
    #cmap2=mpl_c.Colormap('reds')

    fig=plt.figure(figsize=[10, 4])
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)

    vmin=np.log10(np.min(total_mass))
    vmax=np.log10(np.max(total_mass))
    #cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax+0.1*(vmax-vmin)),cmap='cubehelix')
    #cmap1=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True))
    #cmap2=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True))
    cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax),cmap=sns.color_palette("Spectral", as_cmap=True))
    verybest_loss=np.inf

    for l in range(len(total_mass)):
        best_loss=np.inf
        best=[0]*5
        for i in range(len(Lbol)):
            for j in range(len(amax)):
                for k in range(len(grainpowerlaw)):
                    for m in range(len(alpha)):
                        for n in range(len(Si)):

                            if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
                                cf=1
                            else:
                                cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
                                if cf>1:
                                    cf=1

                            loss=0
                            for t in range(len(W1[0])):
                                loss+=(W1[1,t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0,t]))])**2
                            for t in range(len(W2[0])):
                                loss+=(W2[1,t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0,t]))])**2
                            if loss<best_loss:
                                best_loss=loss
                                best=[i,j,k,m,n]
                                if loss<verybest_loss:
                                    verybest_loss=loss
                                    verybest=[i,j,k,l,m,n]

        [i,j,k,m,n]=best
        print('best ',total_mass[l],' - ',best)
        if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
            cf=1
        else:
            cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
            if cf>1:
                cf=1

        ax1.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[i,j,k,l,m,n,:]*cf,label=total_mass[l],color=cmap.to_rgba(np.log10(total_mass[l])))#color=cmap1(amax[i]))
        ax2.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[i,j,k,l,m,n,:]*cf,label=total_mass[l],color=cmap.to_rgba(np.log10(total_mass[l])))#color=cmap2(amax[i]))

    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    ax1.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
    ax2.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='WISE data',color='indianred')
    ax1.plot([L[0][np.argmax(L[1,:])]]*2,ax1.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    ax2.plot([L[0][np.argmax(L[1,:])]]*2,ax2.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'$Total Mass$')
    ax1.set_xlabel('MJD')
    ax2.set_xlabel('MJD')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax1.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax1.set_title('W1')
    ax2.set_title('W2')

    fig.subplots_adjust(right=0.85)

    fig.suptitle(target+' Dust Echo Model')
    fig.savefig(PLOTFILES+'compare_mass.pdf')

    print(target+' - very best: ')
    print('total_mass: ',total_mass[verybest[3]])
    print(verybest)


def compare_Lbol(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES):
    w1=np.argmin(np.abs(wavelengths-wl_w1.value))
    w2=np.argmin(np.abs(wavelengths-wl_w2.value))

    #cmap1=mpl_c.Colormap('blues')
    #cmap2=mpl_c.Colormap('reds')

    fig=plt.figure(figsize=[10, 4])
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)

    vmin=np.log10(np.min(Lbol))
    vmax=np.log10(np.max(Lbol))
    #cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax+0.1*(vmax-vmin)),cmap='cubehelix')
    #cmap1=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True))
    #cmap2=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin-0.1*(vmax-vmin),vmax=vmax), cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True))
    cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=vmin,vmax=vmax),cmap=sns.color_palette("Spectral", as_cmap=True))
    verybest_loss=np.inf

    for i in range(len(Lbol)):
        best_loss=np.inf
        best=[0]*5
        for j in range(len(amax)):
            for k in range(len(grainpowerlaw)):
                for l in range(len(total_mass)):
                    for m in range(len(alpha)):
                        for n in range(len(Si)):

                            if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
                                cf=1
                            else:
                                cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
                                if cf>1:
                                    cf=1

                            loss=0
                            for t in range(len(W1[0])):
                                loss+=(W1[1,t]-cf*Grid1[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W1[0,t]))])**2
                            for t in range(len(W2[0])):
                                loss+=(W2[1,t]-cf*Grid2[i,j,k,l,m,n,np.argmin(np.abs(timesteps-W2[0,t]))])**2
                            if loss<best_loss:
                                best_loss=loss
                                best=[j,k,l,m,n]
                                if loss<verybest_loss:
                                    verybest_loss=loss
                                    verybest=[i,j,k,l,m,n]

        [j,k,l,m,n]=best
        print('best ',Lbol[i],' - ',best)
        if np.linalg.norm(Grid1[i,j,k,l,m,n,:])==0.:
            cf=1
        else:
            cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:]) #covering factor
            if cf>1:
                cf=1

        ax1.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[i,j,k,l,m,n,:]*cf,label=Lbol[i],color=cmap.to_rgba(np.log10(Lbol[i])))#color=cmap1(amax[i]))
        ax2.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[i,j,k,l,m,n,:]*cf,label=Lbol[i],color=cmap.to_rgba(np.log10(Lbol[i])))#color=cmap2(amax[i]))

    maskw1=(W1[1]>W1[2])
    maskw2=(W2[1]>W2[2])
    ax1.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
    ax2.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='WISE data',color='indianred')
    ax1.plot([L[0][np.argmax(L[1,:])]]*2,ax1.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    ax2.plot([L[0][np.argmax(L[1,:])]]*2,ax2.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
    #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'$L_{\mathrm{bol}}$')
    ax1.set_xlabel('MJD')
    ax2.set_xlabel('MJD')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax1.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax1.set_title('W1')
    ax2.set_title('W2')

    fig.subplots_adjust(right=0.85)

    fig.suptitle(target+' Dust Echo Model')
    fig.savefig(PLOTFILES+'compare_Lbol.pdf')

    print(target+' - very best: ')
    print('Lbol: ',Lbol[verybest[0]])
    print(verybest)

    return verybest





def SEDplot(luminosity,timesteps,wavelengths,peak,prefix):
    fig=plt.figure(figsize=[7.5, 4])
    ax=plt.subplot(111)
    #
    cmap=cm.ScalarMappable(norm=mpl_c.Normalize(vmin=(np.min(timesteps)-peak)*0.,vmax=(np.max(timesteps)-peak)*1.05),cmap='cubehelix')
    
    for t in np.linspace(3,len(timesteps)-1,num=10,dtype=int):
        ax.plot(wavelengths,luminosity[t,:],label='{0:.1f}'.format(timesteps[t]-peak),color=cmap.to_rgba(timesteps[t]-peak))
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),title='Time since peak [days]')
    ax.set_xlabel('Wavelength [micron]')
    #ax.set_ylabel('ν L_ν [erg/s]')
    ax.set_ylabel(r'$L_{\nu}$ [erg/s]')
    #cbar1=fig.colorbar(cmap1,ax=ax)
    #cbar2=fig.colorbar(cmap2,ax=ax)
    #cbar2.set_label(r'$a_{\mathrm{max}}$')
    #fig.suptitle('ASASSN15lh Dust Echo Model')
    ax.set_title(target+' Dust Echo Model SED')
    fig.savefig(PLOTFILES+'_SED.pdf')

if __name__=='__main__':
    target="ASASSN-14li"#"ASASSN-15lh"

    if target=='ASASSN-15lh':
        distance=1171*u.Mpc #for ASASSN15lh
        peak=57250 #for ASASSN15lh
        ylim=[2e41,3e45]
    elif target=='ASASSN-14li':
        distance=230*u.Mpc #for AT2019dsg (Stein+2021)
        peak=56906 #MJD for ASASSN-14li
        ylim=[6e39,3e43]
    elif target=='AT2019dsg':
        distance=89.5*u.Mpc #for ASASSN14li (Jiang+2016 https://iopscience.iop.org/article/10.3847/2041-8205/828/1/L14/pdf)
        peak=58605 #for AT2019dsg
        ylim=[5e41,5e43]
    else:
        raise ValueError("Target Distance unknown")
    
    day='Thesis'
    today='Thesis/Compile(9)'
    prefix='Thesis-Results/final'
    
    GRAINCUTOFF=True
    GPL=True #plot for different values for grainpowerlaw
    SI=True
    ALPHA=True
    MASS=True
    LBOL=True

    SED=False

    DATAFILES="data/"+target+'/'
    PLOTFILES="plots/"+target+'/'+str(today)+'/'
    if os.path.isdir(PLOTFILES)==False:
        os.makedirs(PLOTFILES[:-1])
    OUTFILES="results/"+target+'/'+str(today)+'/'
    if os.path.isdir(OUTFILES)==False:
        os.makedirs(OUTFILES[:-1])
    GRIDFILES="results/"+target+'/'+str(day)+'/'

    amax=10
    Lbol=1.279e45
    centralBin=0.01
    
    total_mass= 1 #optical_depth=100
    Si=0. #Silicates fraction (as opposed to graphites)
    
    if SED:
        #f = open(GRIDFILES+'GridSearch/'+str(day)+'_params.txt','r')
        #alpha = float(f.readline())
        #amin = float(f.readline())
        #x = np.fromstring(f.readline()[1:-1],sep=' ')
        #grainpowerlaw = float(f.readline())
        #x = np.fromstring(f.readline()[1:-1],sep=' ')
        #timesteps = np.fromstring(f.readline()[1:-1],sep=' ')
        #wavelengths = np.fromstring(f.readline()[1:-1],sep=' ')

        f = open(GRIDFILES+'GridSearch/'+str(day)+'_params.txt','r')
        alpha = np.fromstring(f.readline()[1:-2],sep=', ')[0]
        amin = np.fromstring(f.readline(),sep=' ')
        amax = np.fromstring(f.readline()[1:-2],sep=', ')[0]
        grainpowerlaw = np.fromstring(f.readline()[1:-2],sep=' ')[0]
        Si = np.fromstring(f.readline(),sep=' ')
        total_mass = np.fromstring(f.readline(),sep=' ')
        Lbol = np.fromstring(f.readline()[1:-2],sep=' ')[0]
        centralBin = np.fromstring(f.readline(),sep=' ')
        timesteps = np.fromstring(f.readline()[1:-2],sep=' ')
        wavelengths = np.fromstring(f.readline()[1:-2],sep=' ')

        #luminosity=np.loadtxt(GRIDFILES+'GridSearch/'+str(day)+'-'+np.format_float_scientific(Lbol, unique=False, exp_digits=2,precision=3)+'-'+np.format_float_scientific(amax, unique=False, exp_digits=2,precision=3)+'.txt')

        f=GRIDFILES+'GridSearch/'+day+'-'+np.format_float_scientific(Lbol, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(amax, unique=False, exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(amin, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(grainpowerlaw, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(total_mass, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(alpha, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(Si, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(centralBin, unique=False,exp_digits=2,precision=3) \
                        +'.txt'

        luminosity=np.loadtxt(f)
        SEDplot(luminosity,timesteps,wavelengths,peak,prefix)

    else:
        #folder=data+'('+str(total_mass)+'-'+str(Si)+'Si-'+str(grainpowerlaw)+')'
        L=np.genfromtxt(DATAFILES+target+'_lc_neutral.dat',skip_header=1).T[0:3]
        T=np.genfromtxt(DATAFILES+target+'_temp.dat',skip_header=1).T[0:3]
        W1=np.genfromtxt(DATAFILES+target+'_WISEproc.txt',skip_header=1)[:3]
        W2=np.genfromtxt(DATAFILES+target+'_WISEproc.txt',skip_header=1)[3:]

        f = open(GRIDFILES+'GridSearch/'+str(day)+'_params.txt','r')
        alpha = np.fromstring(f.readline()[1:-2],sep=', ')
        print(alpha)
        amin = np.fromstring(f.readline(),sep=' ')
        amax = np.fromstring(f.readline()[1:-2],sep=', ')
        grainpowerlaw = np.fromstring(f.readline()[1:-2],sep=' ')
        Si = np.fromstring(f.readline(),sep=' ')
        total_mass = np.fromstring(f.readline(),sep=' ')
        Lbol = np.fromstring(f.readline()[1:-2],sep=' ')
        centralBin = np.fromstring(f.readline(),sep=' ')
        timesteps = np.fromstring(f.readline()[1:-2],sep=' ')
        wavelengths = np.fromstring(f.readline()[1:-2],sep=' ')
        """
        if target=="ASASSN-14li":
            Lbol=[1.599e43,4.798e43,7.997e43,1.120e44,1.439e44,1.599e44]
            amax=['1.000e-02','5.000e-02','1.000e-01','5.000e-01','7.500e-01','1.000e+00','5.000e+00','7.500e+00','1.000e+01','5.000e+01']
            amin=['1.000e-04']
            grainpowerlaw=['1.500e+00','2.000e+00','2.500e+00','3.000e+00']
            total_mass=['1.000e+00']
            alpha=['1.000e+00','1.500e+00','2.000e+00','2.500e+00']
            Si=['0.000e+00','2.500e-01','5.000e-01','7.500e-01','9.999e-01']
            centralBin=['1.000e-02']
        """
        if False:
            pass
        else:
            import glob
            files=glob.glob(GRIDFILES+'GridSearch/Thesis-*')
            Lbol=[]
            amax=[]
            amin=[]
            grainpowerlaw=[]
            total_mass=[]
            alpha=[]
            Si=[]
            centralBin=[]
            name=len(GRIDFILES+'GridSearch/Thesis-')
            for f in files:
                if float(f[name:name+9]) not in Lbol:
                    Lbol.append(float(f[name:name+9]))
                if f[name+10:name+10+9] not in amax:
                    amax.append(f[name+10:name+10+9])
                if f[name+2*10:name+2*10+9] not in amin:
                    amin.append(f[name+2*10:name+2*10+9])
                if f[name+3*10:name+3*10+9] not in grainpowerlaw:
                    grainpowerlaw.append(f[name+3*10:name+3*10+9])
                if f[name+4*10:name+4*10+9] not in total_mass:
                    total_mass.append(f[name+4*10:name+4*10+9])
                if f[name+5*10:name+5*10+9] not in alpha:
                    alpha.append(f[name+5*10:name+5*10+9])
                if f[name+6*10:name+6*10+9] not in Si:
                    Si.append(f[name+6*10:name+6*10+9])
                if f[name+7*10:name+7*10+9] not in centralBin:
                    centralBin.append(f[name+7*10:name+7*10+9])


        w1=np.argmin(np.abs(wavelengths-wl_w1.value))
        w2=np.argmin(np.abs(wavelengths-wl_w2.value))
        Grid1=np.zeros((len(Lbol),len(amax),len(grainpowerlaw),len(total_mass),len(alpha),len(Si),len(timesteps)))
        Grid2=np.zeros((len(Lbol),len(amax),len(grainpowerlaw),len(total_mass),len(alpha),len(Si),len(timesteps)))

        #amax2=[]
        #for a in amax:
        #    if float(a)>=10 and float(a)<=10:
        #        amax2.append(a)
        # 
        #amax=amax2
        total_mass=['1.000e+00']
        print(centralBin)
        print(amin)
        
        Lbol=np.sort(np.array(Lbol))
        amax=np.sort(np.array(amax,dtype=np.float32))
        amin=np.sort(np.array(amin,dtype=np.float32))
        grainpowerlaw=np.sort(np.array(grainpowerlaw,dtype=np.float32))
        total_mass=np.sort(np.array(total_mass,dtype=np.float32))
        alpha=np.sort(np.array(alpha,dtype=np.float32))
        Si=np.sort(np.array(Si,dtype=np.float32))
        centralBin=np.sort(np.array(centralBin,dtype=np.float32))

        if target=='ASASSN-14li':
            Lbol=Lbol[(Lbol>1.5e43)&(Lbol<1.6e44)]
        
        counter=0
        for i in range(len(Lbol)):
            for j in range(len(amax)):
                for k in range(len(grainpowerlaw)):
                    for l in range(len(total_mass)):
                        for m in range(len(alpha)):
                            for n in range(len(Si)):
                                f=GRIDFILES+'GridSearch/'+day+'-'+np.format_float_scientific(Lbol[i], unique=False,exp_digits=2,precision=3) \
                                    +'-'+np.format_float_scientific(amax[j], unique=False, exp_digits=2,precision=3) \
                                    +'-'+np.format_float_scientific(amin[0], unique=False,exp_digits=2,precision=3) \
                                    +'-'+np.format_float_scientific(grainpowerlaw[k], unique=False,exp_digits=2,precision=3) \
                                    +'-'+np.format_float_scientific(total_mass[l], unique=False,exp_digits=2,precision=3) \
                                    +'-'+np.format_float_scientific(alpha[m], unique=False,exp_digits=2,precision=3) \
                                    +'-'+np.format_float_scientific(Si[n], unique=False,exp_digits=2,precision=3) \
                                    +'-'+np.format_float_scientific(centralBin[0], unique=False,exp_digits=2,precision=3) \
                                    +'.txt'
                                #f=GRIDFILES+'GridSearch/'+day+'-'+np.format_float_scientific(Lbol[i], unique=False,exp_digits=2,precision=3) \
                                #           +'-'+amax[j]+'-'+amin[0]+'-'+grainpowerlaw[k]+'-'+total_mass[l]+'-'+alpha[m]+'-'+Si[n]+'-'+centralBin[0]+'.txt'
                                
                                if os.path.isfile(f):
                                    luminosity=np.loadtxt(f)
                                else:
                                    luminosity=np.zeros((len(timesteps),len(wavelengths)))
                                    counter+=1
                                      
                                Grid1[i,j,k,l,m,n,:]=luminosity[:,w1]
                                Grid2[i,j,k,l,m,n:]=luminosity[:,w2]
                                        
        print(counter)
        print(np.sum(Grid1.shape))
        
        print(grainpowerlaw)
        print(grainpowerlaw.shape)
        print(Lbol)
        
        if GRAINCUTOFF:
            compare_amax2(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES)
                                                                        
        if GPL:
            compare_grainpowerlaw(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES)
                                                                            
        if SI:
            compare_sillicate(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES)
                                                                                
        if ALPHA:
            compare_alpha(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES)
                                                                                    
        if MASS:
            compare_mass(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES)
                                                                                        
        if LBOL:
            verybest=compare_Lbol(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,total_mass,alpha,Si,wavelengths,PLOTFILES)

        i,j,k,l,m,n=verybest
        fig=plt.figure(figsize=[10, 4])
        ax=plt.subplot(111)
        cf=np.max(W1[1])/np.max(Grid1[i,j,k,l,m,n,:])
        ax.plot(timesteps,(c.c/wl_w1).to(u.Hz).value*Grid1[i,j,k,l,m,n,:]*cf,label='W1 fit',color='slategrey')#color=cmap1(amax[i]))
        ax.plot(timesteps,(c.c/wl_w2).to(u.Hz).value*Grid2[i,j,k,l,m,n,:]*cf,label='W2 fit',color='indianred')#color=cmap2(amax[i]))

        maskw1=(W1[1]>W1[2])
        maskw2=(W2[1]>W2[2])
        ax.errorbar(W1[0][maskw1],(c.c/wl_w1).to(u.Hz).value*W1[1][maskw1],yerr=(c.c/wl_w1).to(u.Hz).value*W1[2][maskw1],fmt='.',capsize=3,label='W1 data',color='slategrey')
        ax.errorbar(W2[0][maskw2],(c.c/wl_w2).to(u.Hz).value*W2[1][maskw2],yerr=(c.c/wl_w2).to(u.Hz).value*W2[2][maskw2],fmt='.',capsize=3,label='W2 data',color='indianred')
        ax.plot([L[0][np.argmax(L[1,:])]]*2,ax.get_ylim(),linestyle='dashed',color='grey',label='TDE Peak')
        #ax.plot([t_data[np.argmax(L_abs[1,:])]+2*(inner*u.pc /c.c).to(u.day).value]*2,ax.get_ylim(),linestyle='dashed',color='lightgrey',label='Inner Radius LTT')
        ax.set_yscale('log')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),title=r'Band, $a_{\mathrm{min}}$')
        ax.set_xlabel('MJD')
        #ax.set_ylabel('ν L_ν [erg/s]')
        ax.set_ylabel(r'$\nu L_{\nu}$ [erg/s]')
        #cbar1=fig.colorbar(cmap1,ax=ax)
        #cbar2=fig.colorbar(cmap2,ax=ax)
        #cbar2.set_label(r'$a_{\mathrm{max}}$')
        #fig.suptitle('ASASSN15lh Dust Echo Model')
        ax.set_title(target+' Dust Echo Model')
        fig.savefig(PLOTFILES+'verybest.pdf')

    """
    if GRAINCUTOFF:
        #folder=data+'('+str(total_mass)+'-'+str(Si)+'Si-'+str(grainpowerlaw)+')'
        L=np.genfromtxt(DATAFILES+target+'_lc_neutral.dat',skip_header=1).T[0:3]
        T=np.genfromtxt(DATAFILES+target+'_temp.dat',skip_header=1).T[0:3]
        W1=np.genfromtxt(DATAFILES+target+'_WISEproc.txt',skip_header=1)[:3]
        W2=np.genfromtxt(DATAFILES+target+'_WISEproc.txt',skip_header=1)[3:]

        f = open(GRIDFILES+'GridSearch/'+str(day)+'_params.txt','r')
        alpha = np.fromstring(f.readline(),sep=' ')[1:-1][0]
        amin = np.fromstring(f.readline(),sep=' ')
        amax = np.fromstring(f.readline()[1:-1],sep=' ')
        grainpowerlaw = np.fromstring(f.readline()[1:-1],sep=' ')
        Si = np.fromstring(f.readline(),sep=' ')
        total_mass = np.fromstring(f.readline(),sep=' ')
        Lbol = np.fromstring(f.readline()[1:-1],sep=' ')
        centralBin = np.fromstring(f.readline(),sep=' ')
        timesteps = np.fromstring(f.readline()[1:-1],sep=' ')
        wavelengths = np.fromstring(f.readline()[1:-1],sep=' ')

        w1=np.argmin(np.abs(wavelengths-wl_w1.value))
        w2=np.argmin(np.abs(wavelengths-wl_w2.value))
        Grid1=np.zeros((len(Lbol),len(amax),len(timesteps)))
        Grid2=np.zeros((len(Lbol),len(amax),len(timesteps)))

        for i in range(len(Lbol)):
            for j in range(len(amax)):
                for k in [2]:#range(len(grainpowerlaw)):
                    f=GRIDFILES+'GridSearch/'+day+'-'+np.format_float_scientific(Lbol[i], unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(amax[j], unique=False, exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(amin, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(grainpowerlaw[k], unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(total_mass, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(alpha, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(Si, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(centralBin, unique=False,exp_digits=2,precision=3) \
                        +'.txt'

                    if os.path.isfile(f):
                        luminosity=np.loadtxt(f)
                    else:
                        luminosity=np.zeros((len(timesteps),len(wavelengths)))
                        print('x')
                    Grid1[i,j,:]=luminosity[:,w1]
                    Grid2[i,j,:]=luminosity[:,w2]

        compare_amax2(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,wavelengths,PLOTFILES)

    elif GPL:
        #folder=data+'('+str(total_mass)+'-'+str(Si)+'Si-'+str(grainpowerlaw)+')'
        L=np.genfromtxt(DATAFILES+target+'_lc_neutral.dat',skip_header=1).T[0:3]
        T=np.genfromtxt(DATAFILES+target+'_temp.dat',skip_header=1).T[0:3]
        W1=np.genfromtxt(DATAFILES+target+'_WISEproc.txt',skip_header=1)[:3]
        W2=np.genfromtxt(DATAFILES+target+'_WISEproc.txt',skip_header=1)[3:]

        f = open(GRIDFILES+'GridSearch/'+str(day)+'_params.txt','r')
        alpha = np.fromstring(f.readline()[1:-1],sep=' ')
        print(alpha)
        amin = np.fromstring(f.readline(),sep=' ')
        amax = np.fromstring(f.readline()[1:-1],sep=' ')
        grainpowerlaw = np.fromstring(f.readline()[1:-1],sep=' ')
        Si = np.fromstring(f.readline(),sep=' ')
        total_mass = np.fromstring(f.readline(),sep=' ')
        Lbol = np.fromstring(f.readline()[1:-1],sep=' ')
        centralBin = np.fromstring(f.readline(),sep=' ')
        timesteps = np.fromstring(f.readline()[1:-1],sep=' ')
        wavelengths = np.fromstring(f.readline()[1:-1],sep=' ')

        w1=np.argmin(np.abs(wavelengths-wl_w1.value))
        w2=np.argmin(np.abs(wavelengths-wl_w2.value))
        Grid1=np.zeros((len(Lbol),len(amax),len(grainpowerlaw),len(timesteps)))
        Grid2=np.zeros((len(Lbol),len(amax),len(grainpowerlaw),len(timesteps)))

        for i in range(len(Lbol)):
            for j in range(len(amax)):
                for k in range(len(grainpowerlaw)):
                    f=GRIDFILES+'GridSearch/'+day+'-'+np.format_float_scientific(Lbol[i], unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(amax[j], unique=False, exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(amin, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(grainpowerlaw[k], unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(total_mass, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(alpha, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(Si, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(centralBin, unique=False,exp_digits=2,precision=3) \
                        +'.txt'

                    if os.path.isfile(f):
                        luminosity=np.loadtxt(f)
                    else:
                        luminosity=np.zeros((len(timesteps),len(wavelengths)))
                        print('x')
                    Grid1[i,j,k,:]=luminosity[:,w1]
                    Grid2[i,j,k,:]=luminosity[:,w2]

        compare_grainpowerlaw(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,grainpowerlaw,wavelengths,PLOTFILES)

    elif SI:
        #folder=data+'('+str(total_mass)+'-'+str(Si)+'Si-'+str(grainpowerlaw)+')'
        L=np.genfromtxt(DATAFILES+target+'_lc_neutral.dat',skip_header=1).T[0:3]
        T=np.genfromtxt(DATAFILES+target+'_temp.dat',skip_header=1).T[0:3]
        W1=np.genfromtxt(DATAFILES+target+'_WISEproc.txt',skip_header=1)[:3]
        W2=np.genfromtxt(DATAFILES+target+'_WISEproc.txt',skip_header=1)[3:]

        f = open(GRIDFILES+'GridSearch/'+str(day)+'_params.txt','r')
        alpha = np.fromstring(f.readline(),sep=' ')
        amin = np.fromstring(f.readline(),sep=' ')
        amax = np.fromstring(f.readline()[1:-1],sep=' ')
        grainpowerlaw = np.fromstring(f.readline()[1:-1],sep=' ')
        Si = np.fromstring(f.readline(),sep=' ')
        total_mass = np.fromstring(f.readline(),sep=' ')
        Lbol = np.fromstring(f.readline()[1:-1],sep=' ')
        centralBin = np.fromstring(f.readline(),sep=' ')
        timesteps = np.fromstring(f.readline()[1:-1],sep=' ')
        wavelengths = np.fromstring(f.readline()[1:-1],sep=' ')

        w1=np.argmin(np.abs(wavelengths-wl_w1.value))
        w2=np.argmin(np.abs(wavelengths-wl_w2.value))
        Grid1=np.zeros((len(Lbol),len(amax),len(Si),len(timesteps)))
        Grid2=np.zeros((len(Lbol),len(amax),len(Si),len(timesteps)))

        for i in range(len(Lbol)):
            for j in range(len(amax)):
                for k in range(len(Si)):
                    f=GRIDFILES+'GridSearch/'+day+'-'+np.format_float_scientific(Lbol[i], unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(amax[j], unique=False, exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(amin, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(grainpowerlaw, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(total_mass, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(alpha, unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(Si[k], unique=False,exp_digits=2,precision=3) \
                        +'-'+np.format_float_scientific(centralBin, unique=False,exp_digits=2,precision=3) \
                        +'.txt'

                    if os.path.isfile(f):
                        luminosity=np.loadtxt(f)
                    else:
                        luminosity=np.zeros((len(timesteps),len(wavelengths)))
                        print('x')
                    Grid1[i,j,k,:]=luminosity[:,w1]
                    Grid2[i,j,k,:]=luminosity[:,w2]

        compare_sillicate(timesteps,Grid1,Grid2,W1,W2,Lbol,amax,Si,wavelengths,PLOTFILES)
    
"""
        
