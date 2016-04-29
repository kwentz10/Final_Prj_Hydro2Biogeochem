# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 09:21:55 2016

@author: Katherine
"""

# -*- coding: utf-8 -*-
"""
Alpine Tundra Model

Created on Sun Mar 27 15:51:07 2016

@author: Katherine
"""

from Parameters_with_Overland_Flow import *
import numpy as np

def Alpine_Tundra_Model(t,y,SWCSOIL,TSOIL,thetaS_mod, LAI, fVl, fVr, fVw, ZS, BD, Vw, R_SOC, delta_Vl, delta_Vr, delta_Vw, delta_LC, delta_SC, xi_Vl, xi_Vr, xi_Vw, xi_LC, xi_SC, KSc, fDOC, KL, KVm, KVl, KVr, KPp, KPw): 
    
    #Inital Conditions
    Vl=y[0] #leaf carbon
    Vw=y[1] #wood carbon
    Vr=y[2] #root carbon
    LC=y[3] #litter carbon
    SC=y[4] #soil carbon
    R_SOC=y[5] #recalcitrant carbon-->not used in model
    DOC=y[6] #dissolved organic carbon    
    NSavail=y[7] #dissolved inorganic nitrogen
    PSavail=y[8] #dissolved inorganic phosphorus
    PSpm=y[9] #phosphorus in parent materials    
    
    #Day of the Year (0-365; 0 is equal to 366)
    nday=int(t%DOY) #check to make sure this works
    
    #Do not run model for timestep if environmental conditions are null
    if SOLAR[nday]==-999:
        print (nday+1)
        return
    
    if TAIR[nday]==-999:
        print (nday+1)
        return
     
    #Set Control Environmental Conditions
    Ta=TAIR[nday]
    thetaS=SWCSOIL[nday]
    Ts=TSOIL[nday]
    #Scale moisture and temperature by modeled moisture in Overland Flow Model
    if nday+1>=GS0 and nday+1<=GS1:
        thetaS_model=thetaS_mod[int((GSL-1)-(GS1-(nday+1))+1)] #substitute thetaS from overland flow model
        SWCfrac_mod_obs=thetaS_model/thetaS #ratio between modeled and obseved soil moisture (given air temp)
        thetaS=thetaS_model #make soil moisture in the model the same as overland flow modeled soil moisture
        Ts=Ts/SWCfrac_mod_obs #scale soil temperature by soil moisture
        
    #Leaching
    if nday+1>=GS0 and nday+1<=GS1:
         Qlch=0.001252566*thetaS
    else: 
         Qlch=0.0
    
    #Porosity
    Poro=1.0-BD/PD
    
    #Maximum GPP
    UPAR=SOLAR[nday]*0.5
    fPAR=(1.0-np.exp(0.5*(-LAI)))
    GPP_max=UPAR*fPAR*epsilon
    
    #Soil Temperature Constraint on GPP
    GPP_Temp=(Tmax-Ta)/(Tmax-Topt)*((Ta/Topt)**(Topt/(Tmax-Topt)))   
    if Ta<0.0:
        GPP_Temp=0.0
    
    #Soil Moisture Constraint on GPP
    SA_eff=(thetaS-PAW)/(Poro-PAW) #effective soil saturation    
    
    if SA_eff<=0.0:
        SA_eff=0.0
        
    if SA_eff>=1.0:
        SA_eff=0.99
    
    GPP_SM=-4.0*SA_eff*SA_eff+4.0*SA_eff    
    
    if (nday+1)<GS0 or (nday+1)>GS1:  #set moisture limitation to zero during winter
        GPP_SM=1.0
     
    #Soil Nutrient Constraint on GPP
    GPP_T_SM=GPP_Temp*GPP_max*GPP_SM #constrained GPP with soil moisture and temperature but not nutrients          
    
    if GPP_T_SM>0.0:
        GPP_N_lim=thetaS*NSavail*ZS/(GPP_T_SM*(fVl/delta_Vl+fVw/delta_Vw+fVr/delta_Vr))
        GPP_P_lim=thetaS*PSavail*ZS/(GPP_T_SM*(fVl/xi_Vl+fVw/xi_Vw+fVr/xi_Vr))
        
        if (GPP_N_lim<GPP_P_lim):
            minimum=GPP_N_lim
        else:
            minimum=GPP_P_lim
        
        if minimum<1.0:
            GPP_nut=minimum
        else:
            GPP_nut=1.0
        
        if GPP_nut<0.0:
            GPP_nut=0.0
        
        #Calculate GPP
        GPP=GPP_T_SM*GPP_nut
    else:
        GPP=0.0
        GPP_N_lim=1.0 #not limiting if GPP_T_SM<0.0 because temperature and moisture are primary limiting factors
        GPP_P_lim=1.0 #not limiting if GPP_T_SM<0.0 because temperature and moisture are primary limiting factors
    
    #Soil Saturation
    SA=thetaS/Poro
        
    #Function of Temperature and Moisture
    F_TM=(0.642*0.642-(SA-0.642)*(SA-0.642))*1.514*np.exp(0.048*Ts)   
    if Ts<=0.0:
        F_TM=0.0
    
    #Temperature Function
    F_T=np.exp(0.0693*Ta)
    if Ta<0.0:
        F_T=0.0
    
    #Vegetation Carbon
    VC=Vl+Vw+Vr

    #Decomposition of Litter Carbon
    LC_decomp=KL*LC*F_TM
    
    #Autotrophic Respiration
    RVm=KVm*VC*F_T #vegetation maintenance respiration
    NPP_Vm=GPP-RVm #NPP taking into account only maintenance respiration, used to determine growth respiration
    
    if NPP_Vm>0:
        RVg=GR*NPP_Vm #vegetation growth respiration
    else: 
        RVg=0.0
    
    RA=RVm+RVg #total autotrophic respiration
    
    #Calculate NPP
    NPP=GPP-RA
    
    if NPP<0.0:
        NPP=0.0
        RA=GPP
        
    #Nitrogen & Phosphorus Precipitation/Deposition
    Nppt=N_PPT/DOY #converts g N m-2 yr-1 to g N m-2 day-1
    Pppt=P_PPT/DOY #converts g N m-2 yr-1 to g N m-2 day-1
    
    #Nitrogen & Phosphorus Fertilization
    if FTLZ_N=="false":
        N_fer=0
    elif FTLZ_N=="true":
        if t>DOY*(RUNY-1):
            if nday+1>=GS0 and nday+1<=GS1:
                N_fer=NFTLZ/GSL  
            else:
                N_fer=0
        else: 
            N_fer=0
            
    if FTLZ_P=="false":
        P_fer=0
    elif FTLZ_P=="true":
        if t>DOY*(RUNY-1):
            if nday+1>=GS0 and nday+1<=GS1:
                P_fer=PFTLZ/GSL  
            else:
                P_fer=0
        else:
            P_fer=0
    
    #Mass Balance of DOC
    DOCp=fDOC*BD*ZS*SC*KSc*F_TM #DOC production
    DOCRh=thetaS*ZS*DOC*KDOC*F_TM #DOC loss due to decomposition
    DOClch=Qlch*thetaS*DOC #DOC leaching
    
    #Mass Balance of Soil Inorganic Nitrogen
    NSnet=(1.0/delta_LC-fLc/delta_SC)*LC_decomp+fDIN*BD*ZS*SC*KSc*F_TM/delta_SC #net N mobilization/mineralization   
    NSupt=(fVl/delta_Vl+fVw/delta_Vw+fVr/delta_Vr)*NPP #N uptake by vegetation  
    NSlch=Qlch*thetaS*NSavail #N leaching
    NSdpt=Nppt #N deposition  
   
    #Mass Balance of Soil Inorganic Phosphorus
    PSnet=(1.0/xi_LC-fLc/xi_SC)*LC_decomp+fDIP*BD*ZS*SC*KSc*F_TM/xi_SC #net P mobilization/mineralization
    PSupt=(fVl/xi_Vl+fVw/xi_Vw+fVr/xi_Vr)*NPP #P uptake by vegetation
    PSppt=thetaS*ZS*KPp*PSavail*F_TM #P mineral precipitation
    PSwthr=KPw*PSpm*F_TM #P weathering
    PSlch=Qlch*thetaS*PSavail #P leaching
    PSdpt=Pppt #P deposition        
        
    #Distribution Coefficient of DOC (other distribution coefficients are in the parameters file, but need to establish BD before determining this one)
    dDOC=0.6*1.0/BD #(m3 g-1) linear distribution coefficent of DOC between sorbed and aqueous phase
    
    #Model ODEs
    f0=fVl*NPP-KVl*Vl #Vl
    f1=0.0 #Vw-->no change
    f2=fVr*NPP-KVr*Vr #Vr
    f3=-LC_decomp+KVl*Vl+KVr*Vr #LC
    f4=(fLc*LC_decomp-BD*ZS*SC*KSc*F_TM)/(BD*ZS) #SOC
    f5=(BD*ZS*SC*KSc*F_TM) #recalcitrant SOC-->DOESN"T DO ANYTHING IN MODEL
    f6=(DOCp-DOCRh-DOClch)/(thetaS*ZS + BD*ZS*dDOC) #DOC
    f7=(NSnet-NSupt-NSlch+NSdpt+N_fer)/(thetaS*ZS+BD*ZS*dN) #DIN
    f8=(PSnet-PSupt-PSppt-PSlch+PSdpt+PSwthr+P_fer)/(thetaS*ZS+BD*ZS*dP) #DIP
    f9=0.0 #PM-->no change
    
#    if t==0.0 or t>365.0:
#        print ("t=%.3f\t\tnday=%d\t\tN_change=%.6f\t\tNSavail=%.6f\t\tNSnet=%.6f\t\tNSupt=%.6f\t\tNSlch=%.6f" %(t,nday,(NSnet-NSupt-NSlch+NSdpt+N_fer)/(thetaS*ZS+BD*ZS*dN),NSavail,NSnet,NSupt,NSlch) ) 
#       
                  
    #Final Output of ODEs for Each Timestep--Use as y values for Next Model Run Timestep                
    return [f0,f1,f2,f3,f4,f5,f6,f7,f8,f9]
    
    
