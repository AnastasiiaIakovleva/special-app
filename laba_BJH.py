# -*- coding: utf-8 -*-

import argparse
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import BJH_functions as BJH

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('-o','--output', default="bjh", help="output image file")
    parser.add_argument('-i','--input', default="IsothermQ841.txt", help="input isotherm file")
    
    args = parser.parse_args()

    # Reading the input file
    inputFile = open(args.input)
    inputLines = inputFile.readlines()
    inputFile.close()
    # Initializing list with data
    pp0 = []
    adsorption = []
    
    values = []
    isIsotherm = False
    factorGas2Liq = 0.001547
    
    # Parsing lines    
    for line in inputLines:
        if "== Isotherm ==" in line:
        #if re.match('== Isotherm ==',line):
            isIsotherm = True
            continue
        if isIsotherm == True:
            if re.match('==',line):
                isIsotherm = False
                break
            if re.search('\d',line):
                values = list(map(float,line.split()))
                pp0.append(values[0])
                adsorption.append(values[1])
    # Finding adsorption and desorption branch
    max_point_idx = pp0.index(max(pp0))
    
    # Converting lists to numpy.ndarrays
    pp0d = np.array(pp0[max_point_idx:])
    des = np.array(adsorption[max_point_idx:])
    
    # Taking adsorption branch in reverse order
    values = pp0[:max_point_idx+1]
    values.reverse()
    pp0a = np.array(values)
    values = adsorption[:max_point_idx+1]
    values.reverse()
    ads = np.array(values)
    
    # BJH calculations
    
    # Vectorizing all functions
    
    BJH_thickness = np.vectorize(BJH.thickness)
    BJH_radiusCore = np.vectorize(BJH.radiusCore)
    BJH_radius2Core = np.vectorize(BJH.radius2Core)
    BJH_radiusPore = np.vectorize(BJH.radiusPore)
    BJH_rAverageMM = np.vectorize(BJH.rAverageMM)
    BJH_rAverageQCI = np.vectorize(BJH.rAverageQCI)
    
    # Core radius calculations via Kelvin's equation
    rka = BJH_radiusCore(pp0a)
    rkd = BJH_radiusCore(pp0d)
    
    # Statistical thickness calculations
    ta = BJH_thickness(pp0a)
    td = BJH_thickness(pp0d)
    
    # Pore radius calculations
    rpa = BJH_radiusPore(pp0a)
    rpd = BJH_radiusPore(pp0d)
    
    # Preparing arrays for each pressure interval
    rka_min = rka[1:]
    rka_max = rka[0:-1]
    ta_min = ta[1:]
    ta_max = ta[0:-1]
    rpa_min = rpa[1:]
    rpa_max = rpa[0:-1]

    rkd_min = rkd[1:]
    rkd_max = rkd[0:-1]
    td_min = td[1:]
    td_max = td[0:-1]
    rpd_min = rpd[1:]
    rpd_max = rpd[0:-1]     
    
    # Mean radius calculations
    rka_mean = BJH_rAverageQCI(rka_max,rka_min)
    rkd_mean = BJH_rAverageQCI(rkd_max,rkd_min)
    
    pp0a_mean = BJH_radius2Core(rka_mean)
    pp0d_mean = BJH_radius2Core(rkd_mean)
    
    rpa_mean = BJH_radiusPore(pp0a_mean)
    rpd_mean = BJH_radiusPore(pp0d_mean)
    
    # Calculation of c-factor for BJH procedure
    ca = rpa_mean/rka_mean
    cd = rpd_mean/rkd_mean
    
    # Square radius calculations for each interval
    R2a = (rpa_max/(rka_max+(ta_max-ta_min)/2.))**2
    R2d = (rpd_max/(rkd_max+(td_max-td_min)/2.))**2
    
    # Converting cc/g STP to mL(N2@77K)/g
    Vliqa = ads*factorGas2Liq
    Vliqd = des*factorGas2Liq
    
    Vliqa_min = Vliqa[1:]
    Vliqa_max = Vliqa[0:-1]
    
    Vliqd_min = Vliqd[1:]
    Vliqd_max = Vliqd[0:-1]
    
    # Preparing arrays for results
    Vpa = np.zeros(Vliqa_max.shape)
    Apa = np.zeros(Vliqa_max.shape)
    Vcuma = np.zeros(Vliqa_max.shape)
    Acuma = np.zeros(Vliqa_max.shape)
    
    Vpd = np.zeros(Vliqd_max.shape)
    Apd = np.zeros(Vliqd_max.shape)
    Vcumd = np.zeros(Vliqd_max.shape)
    Acumd = np.zeros(Vliqd_max.shape)
        
    # Main loop of calculations
    
    # Cumulative cAp 
    cApa = 0.0
    cApd = 0.0
    
    # Temporary variable
    dummy = 0.0
     
    for idx in np.arange(Vpa.shape[0]):
        # Pore volume calculation at a current step
        dummy = R2a[idx]*(Vliqa_max[idx]-Vliqa_min[idx]-(ta_max[idx]-ta_min[idx])*cApa/1000.)
        
        # Checking if it is greater than zero. If not, it's set zero sharp
        if dummy > 0:
            Vpa[idx]= dummy
        else:
            Vpa[idx]=0.0
        
        # Pore area calculation
        Apa[idx]= 2000.*Vpa[idx]/rpa_mean[idx]
        
        # Cumulative cAp calculation
        cApa += ca[idx]*Apa[idx]
    
    # reordering Rp for cumulative V and A plotting
    rpa_mean_reordered = np.flip(rpa_mean,0)    
    Vcuma = np.cumsum(np.flip(Vpa,0))
    Acuma = np.cumsum(np.flip(Apa,0))
    
       
    for idx in np.arange(Vpd.shape[0]):
        # Pore volume calculation at a current step
        dummy = R2d[idx]*(Vliqd_max[idx]-Vliqd_min[idx]-(td_max[idx]-td_min[idx])*cApd/1000.)
        
        # Checking if it is greater than zero. If not, it's set zero sharp
        if dummy > 0:
            Vpd[idx]= dummy
        else:
            Vpd[idx]=0.0
         
        # Pore area calculation
        Apd[idx]= 2000.*Vpd[idx]/rpd_mean[idx]
        
        # Cumulative cAp calculation
        cApd += cd[idx]*Apd[idx]
    
     # reordering Rp for cumulative V and A plotting
    rpd_mean_reordered = np.flip(rpd_mean,0)    
    Vcumd = np.cumsum(np.flip(Vpd,0))
    Acumd = np.cumsum(np.flip(Apd,0))
    
    #desorption branch
            
    # Plot dV/dlogD vs D
    Dd = 2*rpd_mean
    Dd_reordered = 2*rpd_mean_reordered
    dVlogd = 0.5*Vpd/(np.log10(rpd_max)-np.log10(rpd_min))
    
    fig1=plt.figure()
    ax1=fig1.add_subplot(1,1,1)
    ax1.set_xlabel(r'Диаметр, нм')
    ax1.set_xscale('log')
    ax1.set_ylabel(r'dV/dlogD, $\mathrm{см^3г^{-1}}$')
    ax1.plot(Dd, dVlogd, '-or',fillstyle = 'none')
    ax1.tick_params('x')
    ax1.set_xticks([0.5,1,2,3,5,10,20,30,50,100,200])
    ax1.set_xticklabels(['0.5','1','2','3','5','10','20','30','50','100','200'])
    ax1.set_xlim(0,200)
    ax1.set_ylim(ymin=0.0)  
    fig1.tight_layout()
    
    ax9=ax1.twinx()
    ax9.plot(Dd_reordered, Vcumd,'-sb',fillstyle = 'none')
    ax9.set_xlabel(r'Диаметр, нм')
    ax9.set_xlim(0,200)
    ax9.set_ylabel(r'V, $\mathrm{см^3г{-1}}$')
    ax9.set_ylim(ymin=0.0)  
    fig1.tight_layout()
    
    
    # Plot dV/dD vs D desorption branch 
    dVd=Vpd/(2*rpd_max-2*rpd_min)
    fig2=plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    ax2.plot(Dd, dVd, '-or',fillstyle = 'none')
    ax2.set_xlabel(r'Диаметр, нм')
    ax2.set_xlim(0,60)
    ax2.set_ylabel(r'dV/dD, $\mathrm{см^3г^{-1}нм{-1}}$')
    ax2.set_ylim(ymin=0.0) 
    
    #Plot V vs D
        
    ax3=ax2.twinx()
    ax3.plot(Dd_reordered, Vcumd,'-sb',fillstyle = 'none')
    ax3.set_xlabel(r'Диаметр, нм')
    ax3.set_xlim(0,60)
    ax3.set_ylabel(r'V, $\mathrm{см^3г{-1}}$')
    ax3.set_ylim(ymin=0.0)  
    fig2.tight_layout()
    
    # Produce output
    fig1.savefig(args.output+'1.pdf', dpi=200)
    fig2.savefig(args.output+'2.pdf', dpi=200)
    
    print ('Full volume ads: {:.4f}cc/g'.format(Vcuma[-1]))
    print ('Full surface area ads: {:.4f}m^2/г'.format(Acuma[-1]))
    print ('Full volume des: {:.4f}cc/g'.format(Vcumd[-1]))
    print ('Full surface area des: {:.4f}m^2/г'.format(Acumd[-1]))
    print ('Average diameter ads: {:.4f}nm'.format(4*Vcuma[-1]/Acuma[-1]*1000))
    print ('Average diameter des: {:.4f}nm'.format(4*Vcumd[-1]/Acumd[-1]*1000))
    
    #Plot dA/dD
    dAd=Apd/(2*rpd_max-2*rpd_min)
    fig6=plt.figure()
    ax6 = fig6.add_subplot(1,1,1)
    ax6.plot(Dd, dAd, '-or',fillstyle = 'none')
    ax6.set_xlabel(r'Диаметр, нм')
    ax6.set_xlim(0,60)
    ax6.set_ylabel(r'dA/dD, $\mathrm{м^2г^{-1}нм{-1}}$')
    ax6.set_ylim(ymin=0.0)
    
    #Plot A vs D
    ax7=ax6.twinx()    
    ax7.plot(Dd_reordered, Acumd,'-sb',fillstyle = 'none')
    ax7.set_xlabel(r'Диаметр, нм')
    ax7.set_xlim(0,60)
    ax7.set_ylabel(r'A, $\mathrm{м^2г^{-1}}$')
    ax7.set_ylim(ymin=0.0)  
    fig6.tight_layout()
    fig6.savefig(args.output+'6.pdf', dpi=200)
    
    #Plot dA/dlogD 
    dAlogd = Apd/(np.log10(rpd_max)-np.log10(rpd_min))
    fig8=plt.figure()
    ax8=fig8.add_subplot(1,1,1)
    ax8.set_xlabel(r'Диаметр, нм')
    ax8.set_xscale('log')
    ax8.set_ylabel(r'dA/dlogD, $\mathrm{м^2г^{-1}}$')
    ax8.plot(Dd, dAlogd, '-or',fillstyle = 'none')
    ax8.tick_params('x')
    ax8.set_xticks([0.5,1,2,3,5,10,20,30,50,100,200])
    ax8.set_xticklabels(['0.5','1','2','3','5','10','20','30','50','100','200'])
    ax8.set_xlim(0,200)
    ax8.set_ylim(ymin=0.0) 
    fig8.tight_layout()
    
    
    ax10=ax8.twinx()    
    ax10.plot(Dd_reordered, Acumd,'-sb',fillstyle = 'none')
    ax10.set_xlabel(r'Диаметр, нм')
    ax10.set_xlim(0,60)
    ax10.set_ylabel(r'A, $\mathrm{м^2г^{-1}}$')
    ax10.set_ylim(ymin=0.0)  
    fig8.tight_layout()
    
    fig8.savefig(args.output+'8.pdf', dpi=200)
    #Mode (desorption branch)
    print ('Des diameter moda: {:.4f}nm'.format(Dd[np.argmax(dVd)]))
    
    #adsorption branch
    
    # Plot dVa/dD vs D
    Da = 2*rpa_mean
    Da_reordered = 2*rpa_mean_reordered
    dVa=Vpa/(2*rpa_max-2*rpa_min)
    
    fig5=plt.figure()
    ax5 = fig5.add_subplot(1,1,1)
    ax5.plot(Da, dVa, '-or',fillstyle = 'none')
    ax5.set_xlabel(r'Диаметр, нм')
    ax5.set_xlim(0,100)
    ax5.set_ylabel(r'dVa/dD, $\mathrm{см^3г^{-1}нм{-1}}$')
    ax5.set_ylim(ymin=0.0) 
    print ('Ads diameter moda: {:.4f}nm'.format(Da[np.argmax(dVa)]))
    fig5.tight_layout()
    
    ax11=ax5.twinx()
    ax11.plot(Da_reordered, Vcuma,'-sb',fillstyle = 'none')
    ax11.set_xlabel(r'Диаметр, нм')
    ax11.set_xlim(0,100)
    ax11.set_ylabel(r'Va, $\mathrm{см^3г{-1}}$')
    ax11.set_ylim(ymin=0.0)  
    fig5.tight_layout()
    
    fig5.savefig(args.output+'5.pdf', dpi=200)
    
    # Plot dVa/dlogD vs D
    dVloga = 0.5*Vpa/(np.log10(rpa_max)-np.log10(rpa_min))
    
    fig12=plt.figure()
    ax12=fig12.add_subplot(1,1,1)
    ax12.set_xlabel(r'Диаметр, нм')
    ax12.set_xscale('log')
    ax12.set_ylabel(r'dVa/dlogD, $\mathrm{см^3г^{-1}}$')
    ax12.plot(Da, dVloga, '-or',fillstyle = 'none')
    ax12.tick_params('x')
    ax12.set_xticks([0.5,1,2,3,5,10,20,30,50,100,200])
    ax12.set_xticklabels(['0.5','1','2','3','5','10','20','30','50','100','200'])
    ax12.set_xlim(0,200)
    ax12.set_ylim(ymin=0.0)  
    fig12.tight_layout()
    
    ax13=ax12.twinx()
    ax13.plot(Da_reordered, Vcuma,'-sb',fillstyle = 'none')
    ax13.set_xlabel(r'Диаметр, нм')
    ax13.set_xlim(0,200)
    ax13.set_ylabel(r'Va, $\mathrm{см^3г{-1}}$')
    ax13.set_ylim(ymin=0.0)  
    fig12.tight_layout()
    
    fig12.savefig(args.output+'12.pdf', dpi=200)
    
    #Plot dAa/dD
    
    dAa=Apa/(2*rpa_max-2*rpa_min)
    fig14=plt.figure()
    ax14 = fig14.add_subplot(1,1,1)
    ax14.plot(Da, dAa, '-or',fillstyle = 'none')
    ax14.set_xlabel(r'Диаметр, нм')
    ax14.set_xlim(0,60)
    ax14.set_ylabel(r'dAa/dD, $\mathrm{м^2г^{-1}нм{-1}}$')
    ax14.set_ylim(ymin=0.0)
    
    #Plot A vs D
    ax15=ax14.twinx()    
    ax15.plot(Da_reordered, Acuma,'-sb',fillstyle = 'none')
    ax15.set_xlabel(r'Диаметр, нм')
    ax15.set_xlim(0,60)
    ax15.set_ylabel(r'Aa, $\mathrm{м^2г^{-1}}$')
    ax15.set_ylim(ymin=0.0)  
    fig14.tight_layout()
    
    fig14.savefig(args.output+'14.pdf', dpi=200)
    
    #Plot dAa/dlogD
    
    dAloga = Apa/(np.log10(rpa_max)-np.log10(rpa_min))
    fig16=plt.figure()
    ax16=fig16.add_subplot(1,1,1)
    ax16.set_xlabel(r'Диаметр, нм')
    ax16.set_xscale('log')
    ax16.set_ylabel(r'dAa/dlogD, $\mathrm{м^2г^{-1}}$')
    ax16.plot(Da, dAloga, '-or',fillstyle = 'none')
    ax16.tick_params('x')
    ax16.set_xticks([0.5,1,2,3,5,10,20,30,50,100,200])
    ax16.set_xticklabels(['0.5','1','2','3','5','10','20','30','50','100','200'])
    ax16.set_xlim(0,200)
    ax16.set_ylim(ymin=0.0) 
    fig16.tight_layout()
    
    ax17=ax16.twinx()    
    ax17.plot(Da_reordered, Acuma,'-sb',fillstyle = 'none')
    ax17.set_xlabel(r'Диаметр, нм')
    ax17.set_xlim(0,200)
    ax17.set_ylabel(r'Aa, $\mathrm{м^2г^{-1}}$')
    ax17.set_ylim(ymin=0.0)  
    fig16.tight_layout()
    
    fig16.savefig(args.output+'16.pdf', dpi=200)
    
    

    #Median (desorption branch)
    H,X1 = np.histogram(Dd[1:], bins=50)
    dx=X1[1]-X1[0]
    F1=np.cumsum(H)*dx
    F1=F1/np.max(F1)
    fig4=plt.figure()
    ax=fig4.add_subplot(111)
    ax.plot(X1[1:], F1, '-or',fillstyle = 'none') 
    print ('Des diameter mediana: {:.4f}nm'.format(X1[1:][np.abs(F1-0.5).argmin()]))
    
    #Median (adsorption branch)
    Da = 2*rpa_mean
    H,X2 = np.histogram(Da[1:], bins=50)
    dx=X2[1]-X2[0]
    F2=np.cumsum(H)*dx
    F2=F2/np.max(F2)
    fig5=plt.figure()
    ax=fig5.add_subplot(111)
    ax.plot(X2[1:], F2, '-or',fillstyle = 'none') 
    print ('Ads diameter mediana: {:.4f}nm'.format(X2[1:][np.abs(F2-0.5).argmin()]))
    
if __name__ == "__main__":
    main()
