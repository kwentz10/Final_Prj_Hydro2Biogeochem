# -*- coding: utf-8 -*-
"""
Created on Wed Mar  13 2006

@author: Katherine
"""

from __future__ import division
import numpy as np
from landlab import RasterModelGrid
from matplotlib import pyplot as plt
from matplotlib import gridspec 
import sys
sys.stdout=open('print','w')
from mpl_toolkits.axes_grid1 import make_axes_locatable
#__________________________INITIALIZE______________________________#

##--Initialize Days that Receive Rainfall (used in run loop)

days=31 #length of growing season

##--Define Constants 

n_roughness=0.045 #roughness of bed surface 
K=0.0001 #hydraluic conductivity (ms-1)
rainrate_mean=5 #mean rain rate (mmhr-1)  
rainrate_max=15/1000/3600 #max rain rate (mmhr-1)
stormlength_mean=5*3600 #mean length of storm --12 hours--in (s)
stormlength_max=12*3600
por=0.5 #1.0-(800000.0 /2640000.0)
soil_height=1.5 #thickness of soil (m)


#Determine Infiltration/Soil Moisture Relationship with Linear Line
inf_max=200/1000/3600 #maximum infiltration rate (ms-1) when soil water content is minimal
inf_min=0/1000/3600 #minimum infiltration rate (ms-1) when soil water content is at saturation
swc_max=por #approximate maximum soil water content (effective soil saturation) (m3/m3)--from data; maximum soil water content=minimum infiltration rate
swc_min=0.08 #approximate minimum soil water content (m3/m3)--from data; minimum soil water content=maximum infiltration rate
m_inf_swc=(inf_max-inf_min)/(swc_min-swc_max) #slope of line: inf=m*swc+b (assuming linear relationship between infiltration rate and soil moisture)
b_inf_swc=inf_max-(swc_min*m_inf_swc) #y intercept of line ^^

##--Randomize Days that Contain Storm, Storm Intensity, and Storm Length

rand_rain=int(np.ceil(0.6*days)) #randomize number of growing season days that receive rain
index_rain=list(set(np.ceil(np.random.rand(rand_rain)*days))) #index of the days that receive rain (1 through # of days)
days_rain=len(index_rain)
rainrate=np.random.exponential(rainrate_mean, days_rain)/1000/3600 #rain rate in (ms-1)
for ii in range(len(rainrate)):
    if rainrate[ii]>rainrate_max:
        rainrate[ii]=rainrate_max
stormlength=np.random.exponential(stormlength_mean, days_rain) #length of storm in (s)
for ii in range(len(stormlength)):
    if rainrate[ii]>stormlength_max:
        rainrate[ii]=stormlength_max

##--Initialize Day Length (used in run loop)

max_hr=24 #hours in a day (hr)
t_min=0 #start of day (s)
t_max=3600*max_hr #length of day (s) 
dt=0.3 #model run timestep (s)
timearray=np.arange(t_min, t_max+dt, dt)
num_runs=len(timearray) #number of runs per day

##--Grid Parameters

dx=10 #spacing between grid cells (1 m^2 grid cell)
num_rows=12 #10 grid rows
num_cols=12 #10 grid columns

##--2D Grid Dimensions

mg=RasterModelGrid(num_cols,num_rows,dx) #make 2D grid

##--Initialize Overland Flow Model Data Fields

zs=mg.add_zeros('node','land_surface_elevation') #bedrock elevation (m)
hh=mg.add_zeros('node','water_thickness') #thickness of water (m)
zw=mg.add_empty('node','water_elevation') #water+bedrock elevation (m)
w_slope=mg.add_zeros('link','water_slope') #slope of water surface; make all links have zero slope
wsub_slope=mg.add_zeros('link','soil_moisture_slope')
q=mg.add_zeros('link','water_flux') #flux of water (ms-1)
sat_zone_thickness=mg.add_zeros('link','saturated_zone_thickness_at_link') #soil water thickness at links
thickness_node=mg.add_zeros('node','saturated_zone_thickness_at_node') #soil water thickness at nodes
subwater_height_node=mg.add_zeros('node','subwater_height_at_node')  #total height of soil water (thickness+bedrock)

##--Initialize Randomized Elevation

for c in range(1,num_cols):  #lowest row
    zs[c]=zs[c-1]+((np.random.rand(1)-0.5)*0.8) #from -30 to +30 m
for r in range(1,num_rows): #all other rows
    n=r*num_cols
    zs[n]=zs[n-num_cols]+((np.random.rand(1)-0.5)*0.8) #from -30 to +30 m
    for c in range(1,num_cols):
        n=(r*num_cols)+c
        zs[n]=((zs[n-num_cols]+zs[n-1])/2.0)+((np.random.rand(1)-0.5)*0.8) #from -30 to +30 m

##--Initialize Water+Bedrock Elevation Grid

zw[:]=zs+hh
bedrock_elev=zs-soil_height #1 m lower than zs--land surface

##--Initialize Site and Hydraulic Conductivity Data Field

site_grid=mg.add_empty('node','meadow_site') #type of alpine meadow; 0==dry meadow 1==wet meadow
sites=[]

for ii in range(len(zs)):
    if zs[ii]<np.mean(zs):
        site_grid[ii]=1
    else:
        site_grid[ii]=0
        
sites+=[site_grid] #put site_grid in a list (in order to write to dictionary)

##--Initialize Infiltration Data Field

inf_grid_timestep=mg.add_zeros('node','infiltration every second')

##--Initialize Soil Moisture Data Field 

soil_moisture_timestep=mg.add_empty('node','soil_moisture_timestep') #initial soil moisture grid (list adds soil moisture grid for every day timestep to put into biogeochem model)
for ii in range(len(soil_moisture_timestep)):
    if site_grid[ii]==0:
        soil_moisture_timestep[ii]= 0.1835 #(volumetric water content m3/m3) number obtained from real data (soil moisture prior to onset of growing season for dry meadow)
    elif site_grid[ii]==1:
        soil_moisture_timestep[ii]=0.2384 #(volumetric water content m3/m3) number obtained from real data (soil moisture prior to onset of growing season for wet meadow)

soil_moisture=[soil_moisture_timestep]

##--Boundary Conditions (ENWS)

mg.set_closed_boundaries_at_grid_edges(False, False, False, False) #water can flow out of all sides

##--Initialize Plotting

#Plot Elevation
elevation_raster=mg.node_vector_to_raster(zs)
fb1=plt.figure(1,figsize=(12,12)) #figure blueprint
figure1 = plt.subplot()
figure1.set_title('Elevation',fontsize=32)
figure1.set_xlabel('Distance (m)',fontsize=24)
figure1.set_ylabel('Distance (m)',fontsize=24)
l1=figure1.imshow(elevation_raster, extent=[0,num_cols*dx,0,num_rows*dx])

divider3=make_axes_locatable(figure1)
cax3=divider3.append_axes("right",size="5%", pad=0.05)   

u1=plt.colorbar(mappable=l1, cax=cax3) #add in ticks=v for stable min/max of colorbar
u1.set_label('Meters',fontsize=24)
fb1.savefig('Elevation_A.png', bbox_inches='tight') 

#Plot for Rainfall, Water Thickness, and Soil Moisture
fb2=plt.figure(2,figsize=(12,10))
rain_plt=[] #for rain plot
rain_time=[] #for rain plot


#__________________________RUN______________________________#

##--Run Through Growing Season Days (the index 0 initial values are declared above; I loop through 1-days+1 in order to account for all growing season days)

for d in range(1,days+1):#days+1
       
    ##--Run Through Each Second of Each Growing Season Day and Model Overland Flow 
    
    for t in range(num_runs): #num_runs 
        
        #Set Up infiltration Grid to Use in Soil Moisture Timestep to Determine Infiltration Rate (reset to zero after each dt)        
        inf_grid_timestep[:]=0.0
        
        #Calculate Overland Water Flux
        #Parameters for Manning Eq.
        w_slope=mg.calculate_gradients_at_links(zw)
        h_edge=mg.map_max_of_link_nodes_to_link('water_thickness')
        #Overland Flow
        q=-np.sign(w_slope)*(1.0/n_roughness)*(h_edge**(5.0/3.0))*((abs(w_slope))**(1.0/2.0)) #Manning Eq. 
        dqdx=mg.calculate_flux_divergence_at_nodes(q[mg.active_links])
                
        #Rain During Storm
        if d in index_rain:
            i=index_rain.index(d)
            if timearray[t]<=stormlength[i]: #cannot go beyond 24 hours; even if the random function produced a value greater than 24 hours for stormlength
                rain=rainrate[i]
            else:
                rain=0   
        else:
            rain=0
        
        #Infiltration Rate as a Function of Soil Moisture
        inf=m_inf_swc*soil_moisture_timestep+b_inf_swc 
        
        #Calculate hh
        dhdt=-dqdx+rain-inf 
        dh=dhdt[mg.core_nodes]*dt #change in water thickness only for core nodes
        hh[mg.core_nodes]+=dh #water thickness only for core nodes; h is additive
        hh[:]=np.maximum(hh,0) #water thickness cannot be negative; make sure bottom limit of z does not go below zs     
        
        #Calculate Infiltration at Each Time Step (additive over the day)
        for ii in range(len(hh)):
            if hh[ii]>=inf[ii]*dt: #if water thickness is equal to greater than maximum infiltration rate*timestep (m)
                inf_grid_timestep[ii]=inf[ii]*dt #to use in soil_moisture_timestep
            elif hh[ii]>0.0 and hh[ii]<inf[ii]*dt: #else if water thickness is greater than zero but less than maximum infiltration rate*timestep (m)
                inf_grid_timestep[ii]=hh[ii] #to use in soil_moisture_timestep
            else:
                inf_grid_timestep[ii]=0.0 #to use in soil_moisture_timestep
        
        #Additive Soil Moisture Time Step for Each dt of Each Day    
        soil_moisture_timestep=soil_moisture_timestep+(inf_grid_timestep*por/soil_height) #convert infiltration thickness to swc
        
        #Calculate Subsurface Flow
        #Parameters for Darcy's Law
        thickness_node[:]=soil_moisture_timestep*soil_height/por    
        subwater_height_node[:]=thickness_node+bedrock_elev
        sat_zone_thickness[:]=mg.map_value_at_max_node_to_link('subwater_height_at_node','saturated_zone_thickness_at_node')
        wsub_slope=mg.calculate_gradients_at_links((soil_moisture_timestep*soil_height/por)+bedrock_elev)
        #Subsurface Flow        
        q_sub=(-wsub_slope)*K*sat_zone_thickness #Darcy's Law
        dqsubdx=mg.calculate_flux_divergence_at_nodes(q_sub[mg.active_links]) #change in soil water flow over space
        dswhdt=dqsubdx[mg.core_nodes] #change in soil water thickness over time
        dswh=dswhdt*dt #change in soil water thickness
        
        #Recalculate Soil Moisture
        soil_moisture_timestep[mg.core_nodes]+=dswc*por/soil_height #convert soil water thickness to SWC
        soil_moisture_timestep[:]=np.maximum(soil_moisture_timestep,0) #swc cannot be negative
        
        for iii in range(len(soil_moisture_timestep)): #make sure soil moisture thickness cannot go beyond saturation
            if soil_moisture_timestep[iii]>por:
                hh[iii]+=soil_height*(soil_moisture_timestep[iii]-por)
                soil_moisture_timestep[iii]=por #make sure it has zero infiltration via a negative infiltration rat
        
        #Update Water+Bedrock Elevation
        zw[:]=hh+zs
        
    ##--Soil Moisture Grid 
    
    #Decay of Soil Moisture (by 15% each day)--run through list of soil_moisture and have nodes in each list index grid decay 
    for x in range(len(soil_moisture_timestep)):
        soil_moisture_timestep[x]=soil_moisture_timestep[x]*0.9
    
    soil_moisture+=[soil_moisture_timestep] #volumetric soil moisture; add inf_grid to soil_moisture (xm3/1m3 of soil); inf_grid=infiltrated water (m) per 1m by 1m grid cell--so inf_grid*1*1 equals inf_grid
    
    ##--Plotting
    
    plt.close()  
    fb2=plt.figure(2,figsize=(12,10))
    gs=gridspec.GridSpec(2,6)
    plt.ion()
    
    #Plot Rainfall 
    if d in index_rain:
        i=index_rain.index(d)
        rain=rainrate[i]
        storm=stormlength[i]
    else:
        rain=0
        storm=0
    
    if storm>24*3600:
        storm=24*3600
        
    rain_plt+=[rain]
    rain_plt2=[]
    
    for xx in range(len(rain_plt)):
        rain_plt2+=[rain_plt[xx]*1000*3600]
    rain_time+=[d]
    
    figure2_ax0=plt.subplot(gs[1,:])      
    figure2_ax0.plot(rain_time,rain_plt2,'-o')
    figure2_ax0.set_xlim(0,days+1)
    figure2_ax0.set_ylim(-5,rainrate_mean*4)
    figure2_ax0.set_xlabel('Growing Season: Days',fontsize=18)
    figure2_ax0.set_ylabel('Rain Intensity (mm/hr)',fontsize=18)
    figure2_ax0.set_title('Rainfall Over the Growing Season',fontsize=24)
    figure2_ax0.text(d,(rain*1000*3600)+2.5,'%.2f hr' % (storm/3600.0))
    
    #Plot Water Thickness
    water_thickness_raster=mg.node_vector_to_raster(hh)
    water_thickness_raster_mm=water_thickness_raster*1000
    figure2_ax1=plt.subplot(gs[0,0:3])
    figure2_ax1.set_title('Water Thickness',fontsize=24)
    figure2_ax1.set_xlabel('Distance (m)',fontsize=14)
    figure2_ax1.set_ylabel('Distance (m)',fontsize=14)
    ax_l1=figure2_ax1.imshow(water_thickness_raster_mm, extent=[0,num_cols*dx,0,num_rows*dx])
    ax_u1=plt.colorbar(mappable=ax_l1) #add in ticks=v for stable min/max of colorbar
    ax_u1.set_label('Millimeters',fontsize=14)
    
    #Plot Soil Moisture
    soil_moisture_raster=mg.node_vector_to_raster(soil_moisture[d]) 
    figure2_ax2=plt.subplot(gs[0,3:])
    figure2_ax2.set_title('Soil Moisture',fontsize=24)
    figure2_ax2.set_xlabel('Distance (m)',fontsize=14)
    figure2_ax2.set_ylabel('Distance (m)',fontsize=14)
    ax_l2=figure2_ax2.imshow(soil_moisture_raster, extent=[0,num_cols*dx,0,num_rows*dx])
    ax_u2=plt.colorbar(mappable=ax_l2) #add in ticks=v for stable min/max of colorbar
    ax_u2.set_label('Volumetric Water Content (m3/m3)',fontsize=14)
    
    #Finalize Plot
    plt.tight_layout()
    plt.savefig('OverlandFlow_A_%d.png'%d)
    plt.pause(.0005)
    

#__________________________FINALIZE______________________________#

##--Make Dictionary of Grid Node Soil Moisture; SoilMois_NODE:[day0,day1,day2...]

node_soilmois={}
for n in range(num_rows*num_cols):
    node_soilmois["SoilMois_{0}".format(n+1)]=[]
for n in range(num_rows*num_cols):
    for i in range(len(soil_moisture)):
        node_soilmois["SoilMois_%d" %(n+1)]+=[soil_moisture[i][n]]

##--Make Dictionary of Grid Node Sites; Site_NODE:[day0]--sites do not change over time

node_sites={}
for n in range(num_rows*num_cols):
    node_sites["Sites_{0}".format(n+1)]=[]
for n in range(num_rows*num_cols):
    for i in range(len(sites)):
        node_sites["Sites_%d" %(n+1)]+=[sites[i][n]]




        

