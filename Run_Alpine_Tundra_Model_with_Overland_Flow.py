# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 09:20:47 2016

@author: Katherine
"""

# -*- coding: utf-8 -*-
"""
Run Alpine Tundra Model for 

Created on Mon Mar 28 12:45:15 2016

@author: Katherine
"""

from Alpine_Tundra_Model_with_Overland_Flow import Alpine_Tundra_Model
from scipy.integrate import ode
import numpy as np
from Parameters_with_Overland_Flow import *
import Overland_Flow_2D_Model_revise6 as of
import xlsxwriter as xlwt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import pyplot as plt


#INITIALIZE GRID VALUES THAT I WANT TO PLOT OVER TIME
NPP_grid_daily=[]
NPP_grid_total=[]

for r in range(of.num_cols*of.num_rows):
    File_Name='node_%d.xlsx'%(r+1)
    ISITE=of.node_sites["Sites_%d" %(r+1)]
    thetaS_mod=of.node_soilmois["SoilMois_%d" %(r+1)]

    #Define Type of Simulation 
    #if equal to 0 then not simulated; if equal to 1 then simulated
    #if all simulations are set to zero then control simulation is assumed
    optim=0 #optimization to estimate unknown parameters; if equal to 1 then control=1 and everything else =0
    printing=0 #print outputs to excel file

    #Site Specific Parameters
    if ISITE[0]==0.0: #dry meadow
        #parameters
        LAI=LAI_D
        fVl=fVl_D
        fVr=fVr_D
        fVw=fVw_D
        ZS=ZS_D
        BD=BD_D
        Vw=Vw_D
        R_SOC=R_SOC_D
        TSOIL=TSOIL_D
        SWCSOIL=SWCSOIL_D
        delta_Vl=delta_Vl_D
        delta_Vr=delta_Vr_D
        delta_Vw=delta_Vw_D
        delta_LC=delta_LC_D
        delta_SC=delta_SC_D
        xi_Vl=xi_Vl_D
        xi_Vr=xi_Vr_D
        xi_Vw=xi_Vw_D
        xi_LC=xi_LC_D
        xi_SC=xi_SC_D   
        #target values
        DOC_T=DOC_T_D
        DIN_T=DIN_T_D
        DIP_T=DIP_T_D
        VC_T=VC_T_D
        LC_T=LC_T_D
        SOC_T=SOC_T_D
        HA_T=HA_T_D  
        if optim==0.0: #optimization not happening
            #estimated parameters
            KSc=0.0000241857949278597
            fDOC = 0.00357666658099339
            KL = 0.0012178918824701
            KVm = 0.000813783228195743
            KVl = 0.00101553982856062
            KVr = KVl
            KPp = 0.00077331948120995
            KPw = 0.00000182876881156541
    elif ISITE[0]==1.0: #wet meadow
        #parameters
        LAI=LAI_W
        fVl=fVl_W
        fVr=fVr_W
        fVw=fVw_W
        ZS=ZS_W
        BD=BD_W
        Vw=Vw_W
        R_SOC=R_SOC_W
        TSOIL=TSOIL_W
        SWCSOIL=SWCSOIL_W
        delta_Vl=delta_Vl_W
        delta_Vr=delta_Vr_W
        delta_Vw=delta_Vw_W
        delta_LC=delta_LC_W
        delta_SC=delta_SC_W
        xi_Vl=xi_Vl_W
        xi_Vr=xi_Vr_W
        xi_Vw=xi_Vw_W
        xi_LC=xi_LC_W
        xi_SC=xi_SC_W
        #target values
        DOC_T=DOC_T_W
        DIN_T=DIN_T_W
        DIP_T=DIP_T_W
        VC_T=VC_T_W
        LC_T=LC_T_W
        SOC_T=SOC_T_W
        HA_T=HA_T_W
        if optim==0.0: #optimization not happenning
            #estimated parameters
            KSc=0.0000207629852717425
            fDOC=0.00474068205258691
            KL=0.00129254119148419
            KVm=0.00141500506767766
            KVl=0.00172376180149576
            KVr=KVl
            KPp=0.000572536431395864
            KPw=0.00000379613997296766
            


    #INITIAL CONDITIONS--------------------------------
    #Initial Conditions Vector for Model ODEs
    y0=[VC_T*fVl,Vw,VC_T*fVr,LC_T,SOC_T,R_SOC,DOC_T,DIN_T,DIP_T,PM_T]

    #Time Array for Model
    dt=1.0
    t0=0.0
    t=np.arange(t0,DOY*RUNY,dt) #(start, stop, timestep)
    
    #RUN MODEL-----------------------------------------
    soln=[]
    soln.append(np.array(y0))
    solver=ode(Alpine_Tundra_Model)
    solver.set_integrator('dopri5')
    solver.set_initial_value(y0).set_f_params(SWCSOIL,TSOIL,thetaS_mod, LAI, fVl, fVr, fVw, ZS, BD, Vw, R_SOC, delta_Vl, delta_Vr, delta_Vw, delta_LC, delta_SC, xi_Vl, xi_Vr, xi_Vw, xi_LC, xi_SC, KSc, fDOC, KL, KVm, KVl, KVr, KPp, KPw)
    
    for tt in t[1:]:
        soln.append(solver.integrate(tt))
        

    #All Years Solutions to ODEs
    soln_matrix=np.matrix(soln)
    ode_results=[soln_matrix[:,0].tolist(), soln_matrix[:,1].tolist(), soln_matrix[:,2].tolist(), soln_matrix[:,3].tolist(),soln_matrix[:,4].tolist(),soln_matrix[:,5].tolist(),soln_matrix[:,6].tolist(),soln_matrix[:,7].tolist(),soln_matrix[:,8].tolist(),soln_matrix[:,9].tolist() ]

    #Initialize Last Year Run Solutions to ODEs
    finalyr_ode_0=[]
    finalyr_ode_1=[]
    finalyr_ode_2=[]
    finalyr_ode_3=[]
    finalyr_ode_4=[]
    finalyr_ode_5=[]
    finalyr_ode_6=[]
    finalyr_ode_7=[]
    finalyr_ode_8=[]
    finalyr_ode_9=[]

    #Determine Index of When Final Year Begins and When Final Year Ends
    finalyear_begi=(len(ode_results[0])/RUNY)*(RUNY-1)-1
    finalyear_endi=finalyear_begi+DOY

    #Make Lists of Last Year Run Solutions to ODEs
    for i in np.arange(finalyear_begi,finalyear_endi+1,1):
        i=int(i)
        finalyr_ode_0+=ode_results[0][i]
        finalyr_ode_1+=ode_results[1][i]
        finalyr_ode_2+=ode_results[2][i]
        finalyr_ode_3+=ode_results[3][i]
        finalyr_ode_4+=ode_results[4][i]
        finalyr_ode_5+=ode_results[5][i]
        finalyr_ode_6+=ode_results[6][i]
        finalyr_ode_7+=ode_results[7][i]
        finalyr_ode_8+=ode_results[8][i]
        finalyr_ode_9+=ode_results[9][i]
    
    #Initialize Modeled Values Which Will Approximate Target Values at Steady State
    TotalDOC=0.0
    TotalDIN=0.0
    TotalDIP=0.0
    TotalSOC=0.0
    TotalVC=0.0
    TotalNPP=0.0
    TotalLC=0.0
    TotalGPP=0.0
    TotalRA=0.0
    TotalRH=0.0  

#####################################################            
    NPP_daily=[] 
    NPP_total=[0]
#####################################################

    #Empty Values for Printing to Excel File    
    if printing==1:
        y1=[]
        y2=[]
        y3=[]
        y4=[]
        y5=[]
        y6=[]
        y7=[]
        y8=[]
        y9=[]
        y10=[]
        y11=[]
        y12=[]
        y13=[]
        y14=[]
        y15=[]
        y16=[]
        y17=[]
        y18=[]
        y19=[]   
        
        #Calculate Constant Values Using Solved ODEs for Every Day During the Final Model Run Year  
    for tstep in range(1,len(finalyr_ode_0)):  #fix tstep because you cant have a negative tstep or else it will take last value of array
    
        Vl=finalyr_ode_0[tstep-1]
        Vw=finalyr_ode_1[tstep-1]
        Vr=finalyr_ode_2[tstep-1]
        LC=finalyr_ode_3[tstep-1]
        SC=finalyr_ode_4[tstep-1]
        DOC=finalyr_ode_6[tstep-1]
        NSavail=finalyr_ode_7[tstep-1]
        PSavail=finalyr_ode_8[tstep-1]
        PSpm=finalyr_ode_9[tstep-1] 
    
        #Do not run model for timestep if environmental conditions are null
        if SOLAR[tstep-1]==-999: #in this case, tstep-1 is for current iteration
            continue
    
        if TAIR[tstep-1]==-999: #in this case, tstep-1 is for current iteration
            continue
        
        #Set Control Environmental Conditions
        Ta=TAIR[tstep-1] #in this case, tstep-1 is for current iteration
        thetaS=SWCSOIL[tstep-1] #in this case, tstep-1 is for current iteration
        Ts=TSOIL[tstep-1] #in this case, tstep-1 is for current iteration
        
        #Scale moisture and temperature by modeled moisture in Overland Flow Model
        if tstep>=GS0 and tstep<=GS1:
            print "Observed Soil Moisture=",thetaS
            print "Observed Soil Temperature=",Ts
            thetaS_model=thetaS_mod[int((GSL-1)-(GS1-(tstep))+1)] #substitute thetaS from overland flow model
            print "Modeled Soil Moisture=", thetaS_model            
            SWCfrac_mod_obs=thetaS_model/thetaS #ratio between modeled and obseved soil moisture (given air temp)
            print "Ratio Between Modeled/Observed Soil Moisture=", SWCfrac_mod_obs            
            thetaS=thetaS_model #make soil moisture in the model the same as overland flow modeled soil moisture
            Ts=Ts/SWCfrac_mod_obs #scale soil temperature by soil moisture
            print "Modeled Soil Temperature=", Ts
        
        #Leaching
        if tstep>=GS0 and tstep<=GS1:
            Qlch=0.001252566*thetaS
        else: 
            Qlch=0.0
    
        #Porosity
        Poro=1.0-BD/PD
    
        #Maximum GPP
        UPAR=SOLAR[tstep-1]*0.5 #in this case, tstep-1 is for current iteration
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
    
        if (tstep)<GS0 or (tstep)>GS1:  #set moisture limitation to zero during winter
            GPP_SM=1.0
     
        #Soil Nutrient Constraint on GPP
        GPP_T_SM=GPP_Temp*UPAR*fPAR*epsilon*GPP_SM #constrained GPP with soil moisture and temperature but not nutrients          
    
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
                if tstep>=GS0 and tstep<=GS1:
                    N_fer=NFTLZ/GSL  
                else:
                    N_fer=0
            else: 
                N_fer=0
            
        if FTLZ_P=="false":
            P_fer=0
        elif FTLZ_P=="true":
            if t>DOY*(RUNY-1):
                if tstep>=GS0 and tstep<=GS1:
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
   
        #Heterotrophic Respiration    
        RH=(1.0-fLc)*LC_decomp+BD*ZS*SC*KSc*F_TM*(1-fDOC)+thetaS*ZS*DOC*KDOC*F_TM #only calculating heterotrophic respiration for the last year in run
    
        #N and P Limitations    
        if GPP_N_lim>1.0:
            GPP_N_lim=1.0
        if GPP_P_lim>1.0:
            GPP_P_lim=1.0
    
        #Accumulate Values     
        if (tstep)>=GS0 and (tstep)<=GS1: #averages and totals during the growing season
            TotalDOC+=finalyr_ode_6[tstep]
            TotalDIN+=finalyr_ode_7[tstep]
            TotalDIP+=finalyr_ode_8[tstep]
            TotalSOC+=finalyr_ode_4[tstep]
            TotalVC+=finalyr_ode_0[tstep]+finalyr_ode_1[tstep]+finalyr_ode_2[tstep]
            TotalNPP+=NPP #(from last ODE iteration)
            TotalGPP+=GPP #(from last ODE iteration)
            TotalLC+=finalyr_ode_3[tstep]
            TotalRH+=RH #(from last ODE iteration)
            TotalRA+=RA #(from last ODE iteration)
            
#####################################################            
            NPP_daily+=[NPP] 
            NPP_total+=[NPP_total[abs(int(GS0-tstep))]+NPP]
#####################################################
                
        #Print Limitations, Respiration, Carbon, and C,N,P Leaching For Every Day of Last Year Run 
        if printing==1:
            y1+=[tstep]
            y2+=[1.0-GPP_Temp]
            y3+=[1.0-GPP_SM]
            y4+=[1.0-GPP_N_lim]
            y5+=[1.0-GPP_P_lim]
            y6+=[GPP_max]
            y7+=[NPP]
            y8+=[GPP]
            y9+=[RA]
            y10+=[RH]
            y11+=[RA+RH]
            y12+=[Vl]
            y13+=[Vr]
            y14+=[VC]
            y15+=[LC]
            y16+=[SC]
            y17+=[DOClch]
            y18+=[NSlch]
            y19+=[PSlch]

    if printing==1:
    
        AvgDOC=TotalDOC/GSL
        AvgDIN=TotalDIN/GSL
        AvgDIP=TotalDIP/GSL
        AvgSOC=TotalSOC/GSL
        AvgTotalVC=TotalVC/GSL
        AvgTotalLC=TotalLC/GSL
    
        wb=xlwt.Workbook(File_Name)
    
        ws1=wb.add_worksheet('Daily Outputs')
        daily_names=['DOY','GPP_Temp','GPP_SM','GPP_N','GPP_P','GPP_Max','NPP','GPP','Ra', 'Rh',
                     'Ecosystem_Resp','Vl','Vr','VC','LC','SOC','DOC_leaching','N_leaching','P_leaching']
        daily_values=[y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19]
        
        for i in range(len(daily_names)):
            ws1.write(0,i,daily_names[i])
    
        for xi in range(len(daily_values)):
            for i in range(len(daily_values[xi])):
                ws1.write(1+i,xi,daily_values[xi][i])
    

        ws2=wb.add_worksheet('Average Outputs')
        avg_names=['AvgTotalVC','AvgTotalLC','AvgSOC','AvgDOC','AvgDIN','AvgDIP','TotalRH','TotalRA',
                    'TotalGPP','TotalNPP','TotalVC','TotalLC','SOC','DOC','DIN','DIP']
        avg_values=[AvgTotalVC,AvgTotalLC, AvgSOC, AvgDOC, AvgDIN, AvgDIP, TotalRH, TotalRA, TotalGPP, TotalNPP, TotalVC,
                 TotalLC, TotalSOC, TotalDOC, TotalDIN, TotalDIP]
    
        for i in range(len(avg_names)):
            ws2.write(0,i,avg_names[i]) 
        
        for i in range(len(avg_values)):
            ws2.write(1,i,avg_values[i])         
        
        wb.close()  
    
    #add values for each node 
    NPP_grid_daily+=[NPP_daily]
    NPP_grid_total+=[NPP_total[1:]]


#Make NPP Matrix
NPP_matrix_daily=np.matrix(NPP_grid_daily)
NPP_list_daily=NPP_matrix_daily.tolist()
Daily_Grid_NPP={}
for d in range(int(GSL)):
    Daily_Grid_NPP["NPP_{0}".format(d+1)]=[]
for d in range(int(GSL)):
    for i in range(len(NPP_list_daily)):
        Daily_Grid_NPP["NPP_%i" %(d+1)]+=[NPP_list_daily[i][d]]   

NPP_matrix_total=np.matrix(NPP_grid_total)
NPP_list_total=NPP_matrix_total.tolist()
Total_Grid_NPP={}
for d in range(int(GSL)):
    Total_Grid_NPP["NPP_{0}".format(d+1)]=[]
for d in range(int(GSL)):
    for i in range(len(NPP_list_total)):
        Total_Grid_NPP["NPP_%i" %(d+1)]+=[NPP_list_total[i][d]] 
    
#Plot NPP
fb3=plt.figure(3,figsize=(24,12)) #figure blueprint 
for p in range(len(Daily_Grid_NPP)):
    
    plt.close()  
    fb3=plt.figure(3,figsize=(24,12)) #figure blueprint
    
    NPP_raster_daily=of.mg.node_vector_to_raster(np.array(Daily_Grid_NPP["NPP_%d" %(p+1)]))
    NPP_raster_total=of.mg.node_vector_to_raster(np.array(Total_Grid_NPP["NPP_%d" %(p+1)]))
   
    figure3_ax1 = plt.subplot(1,2,1)
    figure3_ax1.set_title('Daily Net Primary Productivity',fontsize=20)
    figure3_ax1.set_xlabel('Distance (m)',fontsize=18)
    figure3_ax1.set_ylabel('Distance (m)',fontsize=18)
    ax_l3=figure3_ax1.imshow(NPP_raster_daily, extent=[0,of.num_cols*of.dx,0,of.num_rows*of.dx])
    
    divider=make_axes_locatable(figure3_ax1)
    cax=divider.append_axes("right",size="5%", pad=0.05)    
    
    ax_u3=plt.colorbar(mappable=ax_l3, cax=cax) #add in ticks=v for stable min/max of colorbar
    ax_u3.set_label('g C m-2 d-1',fontsize=18)
    
    figure3_ax2=plt.subplot(1,2,2)
    figure3_ax2.set_title('Total Net Primary Productivity', fontsize=20)
    figure3_ax2.set_xlabel('Distance (m)',fontsize=18)
    figure3_ax2.set_ylabel('Distance (m)',fontsize=18)
    ax_l4=figure3_ax2.imshow(NPP_raster_total, extent=[0,of.num_cols*of.dx,0,of.num_rows*of.dx])
    
    divider2=make_axes_locatable(figure3_ax2)
    cax2=divider2.append_axes("right",size="5%", pad=0.05) 
    
    ax_u4=plt.colorbar(mappable=ax_l4, cax=cax2) #add in ticks=v for stable min/max of colorbar
    ax_u4.set_label('g C m-2',fontsize=18)
    
    fb3.savefig('NPP_A_%d.png'%(p+1), bbox_inches='tight')  
    plt.pause(0.005)
    

        
