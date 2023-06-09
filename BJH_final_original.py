# -*- coding: utf-8 -*-
"""
Created on Wed Jan 6 14:57:52 2021

@author: 79149
"""

import tkinter
from tkinter import filedialog as fd
import matplotlib.pyplot as plt
import numpy as np
import re
import math as m


def Pavg(x):
    return m.exp(-2*9.53/x)
def RC(x):
    return -9.53/m.log(x)
def TW(x):
    return 3.54*(-5/m.log(x))**(1/3)
def Davg(x,y):
    return 2*((x+y)*x*y/(x**2+y**2))
def CSAa(x,y):
    return m.pi*((x/2)+y)**2
def CSA(x,y):
    return m.pi*(y*(x+y))
def deltaTW(x,y,z):
    return 2*x/(m.sqrt((y*z*m.pi)**2 + 4*m.pi*x*y) + m.pi*y*z)
# def plotBJH(D1,dV,\
#             lblY1 = r'dV/dlogD, $\mathrm{см^3г^{-1}}$'):
    
#     fig,ax1 = plt.subplots()
    
#     ax1.set_xlabel(r'Диаметр, А')
#     ax1.set_xscale('log')
#     ax1.set_ylabel(lblY1)
#     ax1.plot(D1, dV,'-sb',fillstyle = 'none')
#     ax1.set_xticks([10,20,30,50,100,200,300,500,1000])
#     ax1.set_xticklabels(['10','20','30','50','100','200','300','500','1000'])
#     ax1.set_ylim(bottom=0.0)
      
    
#     plt.tight_layout()

def plotBJH(D1,dV,D2,V,\
                lblY2 = r'V(D), $\mathrm{см^3г^{-1}}$',\
                lblY1 = r'dV/dlogD, $\mathrm{см^3г^{-1}}$'):
    
    fig,ax1 = plt.subplots()
    
    ax1.set_xlabel(r'Диаметр, А')
    ax1.set_xscale('log')
    ax1.set_ylabel(lblY2)
    ax1.plot(D1, V,'-sb',fillstyle = 'none')
    ax1.set_xticks([10,20,30,50,100,200,300,500,1000])
    ax1.set_xticklabels(['10','20','30','50','100','200','300','500','1000'])
    ax1.set_ylim(bottom=0.0)
      
    ax2 = ax1.twinx()
    ax2.set_ylabel(lblY1)
    ax2.plot(D2, dV, '-or',fillstyle = 'none')
    ax2.set_ylim(bottom=0.0)
    
    
    plt.tight_layout()
    

def main():
    
    #finding the adsorption isotherm values ​​in the raw file
    
    root = tkinter.Tk()
    root.filename = fd.askopenfilename(initialdir='/Изотермы адсорбции', title='isotherms', filetypes=(('txt files ', '*.txt'),('all files', '*.*')))
    txt_file = root.filename
    root.mainloop()
    
    lines = []                                                                 
    
    with open(txt_file,'r') as text:
        lines = text.readlines()
    
    pr = []
    ads = []
    
    flag_isotherm = False    
    
    for k in lines:
                
                if re.search('== Isotherm ==',k):
                    flag_isotherm = True
                    continue
                
                if re.search('== Absolute Isotherm ==',k):
                    flag_isotherm = False
                    continue
                
                if flag_isotherm:
                    
                    if re.search('\d',k):
                      values = list(map(float,k.split()))
                      
                      if ~np.isnan(values[0]):
                          pr.append(values[0])
                          ads.append(values[1])
                          continue

    #splitting adsorption and desorption curves
    
    max_ads = max(ads)
    max_ads_index = ads.index(max_ads)
    
    adsorption = ads[:max_ads_index+1]
    pr_adsorption = pr[:max_ads_index+1]
    
    desorption = ads[max_ads_index:]
    
    pr_desorption = pr[max_ads_index:]
    
    
    
    #BJH calculations
    
    
    #Volume gas2liquid
    VI=[]
    for k in desorption:
        we = desorption.index(k)
        k = k*0.0015468
        VI.append(k) 
    
    #empty list for calcilations 
    Rci = []
    Rcj = []
    Daverage =[]
    LP = []
    
    #desorption branch
    for n in range(len(VI)):
        if n == 0:
            Vc_1 = VI[0] - VI[1]
            TW_1 = TW(pr_desorption[0])
            TW_2 = TW(pr_desorption[1])
            delta_Tw = TW_1-TW_2
            RC_1 = RC(pr_desorption[0])
            RC_2 = RC(pr_desorption[1])
            Davg_1 = Davg(RC_1,RC_2)
            Pravg_1 = Pavg(Davg_1)
            TWavg_1 = TW(Pravg_1)
            delta_Td = TWavg_1 - TW_2
            CSAa_1 = CSAa(Davg_1,delta_Td)
            LP_1 = Vc_1/CSAa_1
            LP.append(LP_1)
            RC_1 = RC_1 + delta_Tw
            Davg_1 = Davg_1 + 2*delta_Td
            Rci.append(RC_1)
            Rcj.append(RC_2)
            Daverage.append(Davg_1)
            
            
        if 0<n<(len(VI)-1):   
            delta_Tw = TW(pr_desorption[n])-TW(pr_desorption[n+1])
            Vdsumm = 0.0 
            for i in range(len(Daverage)):
                Vdsumm += LP[i]*(Daverage[i]+delta_Tw)
            v_d = m.pi*delta_Tw*Vdsumm
            
            if v_d > (VI[n] - VI[n+1]): #не происходит отккрытия новых поря
                lsumm = 0.0 #summary of length 
                ldsumm = 0.0 #summary of length multiply to D
                for g in range(len(Daverage)):
                    lsumm += LP[g]
                    ldsumm += LP[g]*Daverage[g]
                    
                delta_Tw = 2*(VI[n] - VI[n+1])/\
                    (m.sqrt((m.pi*ldsumm)**2 + 4*m.pi*(VI[n] - VI[n+1])*lsumm) + m.pi*ldsumm)
               
                for b in range(len(Daverage)):
                    Daverage[b] += 2*delta_Tw
                    Rci[b] += delta_Tw
                    Rcj[b] += delta_Tw 
            
            if v_d < (VI[n] - VI[n+1]):          #освобождение новых пор
                Vc_j = VI[n] - VI[n+1] - v_d
                RC_j = RC(pr_desorption[n])
                RC_k = RC(pr_desorption[n+1]) 
                Davg_j = Davg(RC_j,RC_k)
                TW_j = TW(pr_desorption[n])
                TW_k = TW(pr_desorption[n+1])
                delta_Tw = TW_j-TW_k
                Pravg_j = Pavg(Davg_j)
                TWavg_j = TW(Pravg_j)
                delta_Td = TWavg_j - TW_k
                CSAa_j = CSAa(Davg_j,delta_Td)
                LP_j = Vc_j/CSAa_j
                LP.append(LP_j)
                for b in range(len(Daverage)):
                    Daverage[b] += 2*delta_Tw
                    Rcj[b] += delta_Tw 
                
                Daverage.append(Davg_j+2*delta_Td)
                
                Rcj.append(RC_k)
                Rci.append(RC_j)
                for b in range(len(Rci)):
                    Rci[b] += delta_Tw
                    
    lsumm = 0.0 
    ldsumm = 0.0 
    for g in range(len(Daverage)):
        lsumm += LP[g]
        ldsumm += LP[g]*Daverage[g]
                   
    delta_Tw = 2*(VI[-1])/\
       (m.sqrt((m.pi*ldsumm)**2 + 4*m.pi*(VI[-1])*lsumm) + m.pi*ldsumm)
              
    for b in range(len(Daverage)):
        Daverage[b] += 2*delta_Tw
        Rci[b] += delta_Tw
        Rcj[b] += delta_Tw 
    
    Rci = np.array(Rci) 
    Rcj = np.array(Rcj)
    LP = np.array(LP)  
    Daverage = np.array(Daverage)   
    
    Di = 2*Rci
    Dj = 2*Rcj
    VPI = m.pi*LP*(Daverage/2)**2            
    SAI = m.pi*LP*Daverage*1.0e04
    
    
#    print('Di','Dj','Daverage','VPI','SAI')
 #   for k in np.arange(Di.shape[0]):
 #       print(Di[k], Dj[k], Daverage[k], VPI[k], SAI[k])
    
    print('По десорбции:')
    print('Удельный объем мезопор , см3/г: ', np.cumsum(VPI)[-1])
    print('Удельная поверхность мезопор , м2/г: ', np.cumsum(SAI)[-1])
    
    
#distribution

   

    dVdDi = VPI/(Di-Dj)
    
    
    
    dAdDi = SAI/(Di-Dj)
    
    
    dVdlogDi = VPI/np.log10(Di/Dj)
    
    
    dAdlogDi = SAI/np.log10(Di/Dj) 
    
    

    # plotBJH(Daverage,dVdDi,Daverage,np.flip(np.cumsum(VPI)),\
            
    #          lblY1 = r'dV/dD, $\mathrm{см^3г^{-1}A^{-1}}$',\
    #          lblY2 = r'V(D), $\mathrm{см^3г^{-1}}$')
        
    # print('Daverage', '      ', 'dVdDi')
    # for k in np.arange(Daverage.shape[0]):
    #            print(Daverage[k], '      ',dVdDi[k])
    
    # h = np.argmax(dVdDi)
    # print('Наиболее вероятный диаметр dV/dD по десорбционной ветке, A:', Daverage[h])
        
        
    # plotBJH(Daverage,dAdDi,Daverage,np.flip(np.cumsum(SAI)),\
            
    #          lblY1 = r'dA/dD, $\mathrm{м^2г^{-1}A^{-1}}$',\
    #          lblY2 = r'A(D), $\mathrm{см^3г^{-1}}$')
    
    # print('Daverage', '      ', 'dAdDi')    
    # for k in np.arange(Daverage.shape[0]):
    #           print(Daverage[k], '      ',dAdDi[k])
         
    # h = np.argmax(dAdDi)
    # print('Наиболее вероятный диаметр dA/dD по десорбционной ветке, A:', Daverage[h])
    
    # plotBJH(Daverage,dVdlogDi,Daverage,np.flip(np.cumsum(VPI)),\
            
    #          lblY1 = r'dV/dlogDi, $\mathrm{см^3г^{-1}}$',\
    #          lblY2 = r'V(D), $\mathrm{см^3г^{-1}}$')
        
    # print('Daverage', '      ', 'dVdlogDi')
    # for k in np.arange(Daverage.shape[0]):
    #         print(Daverage[k], '      ',dVdlogDi[k])
    
    # h = np.argmax(dVdlogDi)
    # print('Наиболее вероятный диаметр dV/dlogD по десорбционной ветке, A:', Daverage[h])    
    
    # plotBJH(Daverage,dAdlogDi,Daverage,np.flip(np.cumsum(SAI)),\
            
    #          lblY1 = r'dAdlogD, $\mathrm{см^3г^{-1}A^{-1}}$',\
    #          lblY2 = r'V(D), $\mathrm{см^3г^{-1}}$')
    
    # print('Daverage', '      ', 'dAdlogDi')
    # for k in np.arange(Daverage.shape[0]):
    #           print(Daverage[k], '      ',dAdlogDi[k]) 
        
    # h = np.argmax(dAdlogDi)
    # print('Наиболее вероятный диаметр dA/dlogD по десорбционной ветке, A:', Daverage[h])
        
  
    
  
    
  
    
  
    # adsorbtion branch
    

    VI=[]
    for k in adsorption:
        we = adsorption.index(k)
        k = k*0.0015468
        VI.append(k) 
    
    
    
    VI.reverse()
    pr_adsorption.reverse()
    
    
   
    Rci = []
    Rcj = []
    Daverage =[]
    LP = []
    
 
    for n in range(len(VI)):
        if n == 0:
            Vc_1 = VI[0] - VI[1]
            TW_1 = TW(pr_adsorption[0])
            TW_2 = TW(pr_adsorption[1])
            delta_Tw = TW_1-TW_2
            RC_1 = RC(pr_adsorption[0])
            RC_2 = RC(pr_adsorption[1])
            Davg_1 = Davg(RC_1,RC_2)
            Pravg_1 = Pavg(Davg_1)
            TWavg_1 = TW(Pravg_1)
            delta_Td = TWavg_1 - TW_2
            CSAa_1 = CSAa(Davg_1,delta_Td)
            LP_1 = Vc_1/CSAa_1
            LP.append(LP_1)
            RC_1 = RC_1 + delta_Tw
            Davg_1 = Davg_1 + 2*delta_Td
            Rci.append(RC_1)
            Rcj.append(RC_2)
            Daverage.append(Davg_1)
            
            
        if 0<n<(len(VI)-1):
            delta_Tw = TW(pr_adsorption[n])-TW(pr_adsorption[n+1])
            Vdsumm = 0.0 
            for i in range(len(Daverage)):
                Vdsumm += LP[i]*(Daverage[i]+delta_Tw)
            v_d = m.pi*delta_Tw*Vdsumm
            
            if v_d > (VI[n] - VI[n+1]):
                lsumm = 0.0 
                ldsumm = 0.0 
                for g in range(len(Daverage)):
                    lsumm += LP[g]
                    ldsumm += LP[g]*Daverage[g]
                    
                delta_Tw = 2*(VI[n] - VI[n+1])/\
                    (m.sqrt((m.pi*ldsumm)**2 + 4*m.pi*(VI[n] - VI[n+1])*lsumm) + m.pi*ldsumm)
               
                for b in range(len(Daverage)):
                    Daverage[b] += 2*delta_Tw
                    Rci[b] += delta_Tw
                    Rcj[b] += delta_Tw 
            
            if v_d < (VI[n] - VI[n+1]):
                Vc_j = VI[n] - VI[n+1] - v_d
                RC_j = RC(pr_adsorption[n])
                RC_k = RC(pr_adsorption[n+1]) 
                Davg_j = Davg(RC_j,RC_k)
                TW_j = TW(pr_adsorption[n])
                TW_k = TW(pr_adsorption[n+1])
                delta_Tw = TW_j-TW_k
                Pravg_j = Pavg(Davg_j)
                TWavg_j = TW(Pravg_j)
                delta_Td = TWavg_j - TW_k
                CSAa_j = CSAa(Davg_j,delta_Td)
                LP_j = Vc_j/CSAa_j
                LP.append(LP_j)
                for b in range(len(Daverage)):
                    Daverage[b] += 2*delta_Tw
                    Rcj[b] += delta_Tw 
                
                Daverage.append(Davg_j+2*delta_Td)
                
                Rcj.append(RC_k)
                Rci.append(RC_j)
                for b in range(len(Rci)):
                    Rci[b] += delta_Tw
                    
    lsumm = 0.0
    ldsumm = 0.0 
    for g in range(len(Daverage)):
        lsumm += LP[g]
        ldsumm += LP[g]*Daverage[g]
                   
    delta_Tw = 2*(VI[-1])/\
        (m.sqrt((m.pi*ldsumm)**2 + 4*m.pi*(VI[-1])*lsumm) + m.pi*ldsumm)
              
    for b in range(len(Daverage)):
        Daverage[b] += 2*delta_Tw
        Rci[b] += delta_Tw
        Rcj[b] += delta_Tw 
    
    Rci = np.array(Rci) 
    Rcj = np.array(Rcj)
    LP = np.array(LP)  
    Daverage = np.array(Daverage)   
    
    Di = 2*Rci
    Dj = 2*Rcj
    VPI = m.pi*LP*(Daverage/2)**2            
    SAI = m.pi*LP*Daverage*1.0e04  
    
    print('По адсорбции:')
    print('Удельный объем мезопор , см3/г: ', np.cumsum(VPI)[-1])
    print('Удельная поверхность мезопор , м2/г: ', np.cumsum(SAI)[-1])
    
    print('По адсорбции:')
    
 #   for k in np.arange(Di.shape[0]):
   #     print(Di[k], Dj[k], Daverage[k], VPI[k], SAI[k])
    
    
    
    

    dVdDi = VPI/(Di-Dj)
    
    dAdDi = SAI/(Di-Dj)
    dVdlogDi = VPI/np.log10(Di/Dj)
    dAdlogDi = SAI/np.log10(Di/Dj) 
    
    
    plotBJH(Daverage,dVdDi,Daverage,np.flip(np.cumsum(VPI)),\
            
              lblY1 = r'dV/dD, $\mathrm{см^3г^{-1}A^{-1}}$',\
              lblY2 = r'V(D), $\mathrm{см^3г^{-1}}$')
        
    print('Daverage', '      ', 'dVdDi')
    for k in np.arange(Daverage.shape[0]):
                print(Daverage[k], '      ',dVdDi[k])
    
    h = np.argmax(dVdDi)
    print('Наиболее вероятный диаметр dV/dD по десорбционной ветке, A:', Daverage[h])
        
        
    plotBJH(Daverage,dAdDi,Daverage,np.flip(np.cumsum(SAI)),\
            
              lblY1 = r'dA/dD, $\mathrm{м^2г^{-1}A^{-1}}$',\
              lblY2 = r'A(D), $\mathrm{см^3г^{-1}}$')
    
    print('Daverage', '      ', 'dAdDi')    
    for k in np.arange(Daverage.shape[0]):
              print(Daverage[k], '      ',dAdDi[k])
         
    h = np.argmax(dAdDi)
    print('Наиболее вероятный диаметр dA/dD по десорбционной ветке, A:', Daverage[h])
    
    plotBJH(Daverage,dVdlogDi,Daverage,np.flip(np.cumsum(VPI)),\
            
              lblY1 = r'dV/dlogDi, $\mathrm{см^3г^{-1}}$',\
              lblY2 = r'V(D), $\mathrm{см^3г^{-1}}$')
        
    print('Daverage', '      ', 'dVdlogDi')
    for k in np.arange(Daverage.shape[0]):
            print(Daverage[k], '      ',dVdlogDi[k])
    
    h = np.argmax(dVdlogDi)
    print('Наиболее вероятный диаметр dV/dlogD по десорбционной ветке, A:', Daverage[h])    
    
    plotBJH(Daverage,dAdlogDi,Daverage,np.flip(np.cumsum(SAI)),\
            
              lblY1 = r'dAdlogD, $\mathrm{см^3г^{-1}A^{-1}}$',\
              lblY2 = r'V(D), $\mathrm{см^3г^{-1}}$')
    
    print('Daverage', '      ', 'dAdlogDi')
    for k in np.arange(Daverage.shape[0]):
              print(Daverage[k], '      ',dAdlogDi[k]) 
        
    h = np.argmax(dAdlogDi)
    print('Наиболее вероятный диаметр dA/dlogD по десорбционной ветке, A:', Daverage[h])
        
    
   
   
    
    
if __name__ == "__main__":
    main()