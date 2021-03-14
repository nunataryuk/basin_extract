'''
Created on 6 Oct 2020

@author: thomasgumbricht
'''

# imports

from __future__ import division

import os

import sys

import numpy as np
from scipy.spatial.distance import cdist
from operator import itemgetter
from params.be_params import Params
from ds_manage.datasource import DS
from scipy.signal._arraytools import even_ext

#import shapely.wkt

class ProcStage1BasinOutlets(DS):
    '''
    classdocs
    '''
    
    def __init__(self, params):
        '''
        Constructor
        '''
        
        # Initiate the data source (DS) manager
        DS.__init__(self)
        
        self.params = params 

        # Set the path to the new basins and data to same as the source path
        self.basinfp = self.datafp = self.params.fp
    
        if self.params.verbose:
            inforstr = '        Stage 1: Scripting up GRASS to find basin outlets'
            print (inforstr)
            
        self.stage1datafp = os.path.join(self.datafp, 'stage1')
        self.stage1scriptfp = os.path.join(self.stage1datafp, 'script')
        if not os.path.exists(self.stage1scriptfp):
            os.makedirs(self.stage1scriptfp)
    
        self._Grassscript()
        
        cmd = '\n# To run the output script you must have GRASS GIS open in a Terminal window session.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += 'chmod 755 %(fpn)s\n' %{'fpn': self.GRASSshFPN}
        cmd += '# Then execute the command from your GRASS terminal session:\n'
        cmd += '%(fpn)s\n\n'%{'fpn': self.GRASSshFPN}
        
        print (cmd)
                         
    def _Grassscript(self):    
        '''
        '''
        GRASSshFN = '%(s)s_grass_find_basin_outlets_stage1.sh' %{'s':self.params.region}

        self.GRASSshFPN = os.path.join(self.stage1scriptfp, GRASSshFN)
        
        cmd = '# Basin delineation: GRASS watershed analysis\n'
        cmd += '# Created by the Python package basin_extract\n\n'
                
        cmd += '# To run this script you must have GRASS GIS open in a Terminal window session.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += '# chmod 755 %(fpn)s\n' %{'fpn': self.GRASSshFPN}
        cmd += '# Then execute the command from your GRASS terminal session:\n'
        cmd += '# GRASS 7.x.3x ("region"):~ > %(fpn)s\n\n'%{'fpn': self.GRASSshFPN}

        cmd += '# Create destination path\n'
        
        cmd += 'mkdir -p %s\n\n' % (self.stage1scriptfp)
        
        cmd += '# Set region from the original DEM, if the layer is not at DEM@PERMANENT you have to manually update the layer reference.\n'
        cmd += 'g.region raster=%(dem)s\n\n' %{'dem': 'DEM@PERMANENT'}
                
        if self.params.grassDEM == 'hydro_fill_dem' and os.path.exists(self.params.hydroFillDEMFPN_s0):
            
            cmd += '# Import hydrologically corrected DEM from stage 0\n'
            
            cmd += 'r.in.gdal input=%(src)s output=%(dst)s\n\n' %{'src':self.params.hydroFillDEMFPN_s0, 'dst':'hydro_fill_dem'}
                        
        cmd += '# Single flow direction (SFD) watershed analysis\n'
        
        cmd += 'r.watershed -as elevation=%(dem)s accumulation=SFD_upstream drainage=SFD_drainage threshold=%(th)d --overwrite\n\n' %{'dem':self.params.grassDEM , 'th':self.params.basinCellThreshold}
        
        cmd += '# Mulitple flow directions (MFD) watershed analysis\n'
        
        cmd += 'r.watershed -a elevation=%(dem)s accumulation=MFD_upstream drainage=MFD_drainage threshold=%(th)d --overwrite\n\n' %{'dem':self.params.grassDEM , 'th':self.params.basinCellThreshold}
        
        cmd += '# Export SFD and MFD with color ramps by removing the "#" sign.\n\n'
        
        cmd += '# Convert SFD to 10 * natural log to get a Byte range\n'
        
        cmd += '# r.mapcalc "SFD_ln_upstream = 10*log(SFD_upstream)" --overwrite\n\n'
        
        cmd += '# Set color ramp\n'
        
        cmd += '# r.colors map=SFD_ln_upstream color=ryb\n\n'
        
        cmd += '# Export as geotiff \n'

        cmd += '# r.out.gdal -f input=SFD_ln_upstream format=GTiff type=Byte output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.upstreamLnSFDfpn_s1}
        
        cmd += '# Convert MFD to 10 * natural log to get a Byte range\n'
        
        cmd += '# r.mapcalc "MFD_ln_upstream = 10*log(MFD_upstream)" --overwrite\n\n'
        
        cmd += '# Set color ramp\n'
        
        cmd += '# r.colors map=MFD_ln_upstream color=ryb\n\n'
        
        cmd += '# Export as geotiff \n'
        
        cmd += '# r.out.gdal -f input=MFD_ln_upstream format=GTiff type=Byte output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.upstreamLnMFDfpn_s1}
        
        cmd += '# Identify terminal water body \n' 
        
        cmd += 'r.mapcalc "drain_terminal = if(isnull(%(dem)s), 1, null())" --overwrite\n\n'  %{'dem':self.params.grassDEM}
        
        cmd += '# Get the shoreline \n'
        
        cmd += 'r.buffer input=drain_terminal output=shoreline distances=%(dist)f units=meters --overwrite\n\n' %{'dist':self.params.adjacentDist}
        
        cmd += '# Extract accumulated SFD drainage for the shoreline \n'
        
        #cmd += 'r.mapcalc "shoreline_SFD_flowacc = if((shoreline > 1), SFD_upstream, null())" --overwrite\n\n'
        
        cmd += 'r.mapcalc "basin_SFD_outlets = if((shoreline > 1), if ((SFD_upstream >= %(bct)d),SFD_upstream,null() ) , null())" --overwrite\n\n' %{'bct':self.params.basinCellThreshold}

        cmd += '# Extract accumulated MFD drainage for the shoreline \n'

        #cmd += 'r.mapcalc "shoreline_MFD_flowacc = if((shoreline > 1), MFD_upstream, null())" --overwrite\n\n'
        
        cmd += 'r.mapcalc "basin_MFD_outlets = if((shoreline > 1), if ((MFD_upstream >= %(bct)d),MFD_upstream,null() ) , null())" --overwrite\n\n' %{'bct':self.params.basinCellThreshold}
        '''
        cmd += '# Threshold minimum drainage area \n'
        cmd += '# You must have created a reclass text file in advance, for example llke this \n'
        cmd += '### 2000 thru 99999999 = 1 ###\n'
        cmd += '### * = NULL ###\n\n'
        
        cmd += '# Reclass SFD shoreline to cells with upstream areas larger than the given threshold\n'
        
        cmd += 'r.reclass input=shoreline_SFD_flowacc output=basin_SFD_outlets rules="%(path)s" --overwrite\n\n' %{'path':self.params.flowaccReclassPath}
        
        cmd += '# Reclass MFD shoreline to cells with upstream areas larger than the given threshold\n'
        
        cmd += 'r.reclass input=shoreline_MFD_flowacc output=basin_MFD_outlets rules="%(path)s" --overwrite\n\n' %{'path':self.params.flowaccReclassPath}
        '''
        
        cmd += '# Combine SFD and MFD accumulated shoreline cells\n'
        
        cmd += 'r.mapcalc "basin_shoreline_outlets = if ((isnull(basin_SFD_outlets)), basin_MFD_outlets, basin_SFD_outlets)" --overwrite\n\n'
        '''
        cmd += '# Reclass the combined shoreline to cells with upstream areas larger than the given threshold\n'

        cmd += 'r.reclass input=shoreline_flowacc output=basin_shoreline_outlets rules="%(path)s" --overwrite\n\n' %{'path':self.params.flowaccReclassPath}
        '''
        cmd += '# Vectorize the combined shoreline basin outlet candidates\n'
        
        cmd += 'r.to.vect input=basin_shoreline_outlets output=basin_outlet_pt type=point --overwrite\n\n'
        
        cmd += '# For the standard approach the separate SFD and MFD solutions are not required,\n'
        cmd += '# to use them in the further analysis remove the "#" sign.\n\n'
        
        cmd += '# Vectorize SFD shoreline basin outlet candidates by removing the "#"\n'
        
        cmd += '# r.to.vect input=basin_SFD_outlets output=basin_SFD_outlets_pt type=point  --overwrite\n\n'
        
        cmd += '# Vectorize MFD shoreline basin outlet candidates by removing the "#"\n'
        
        cmd += '# r.to.vect input=basin_MFD_outlets output=basin_MFD_outlets_pt type=point --overwrite\n\n'
        
        cmd += '# Identify basin river mouths from cells adjacent to the terminal drainage \n'
        cmd += '# The default is to set the motuhs to cells with an elevation <= 0.\n'
        cmd += '# Change the value to include more upstream areas or if your terminal drainage is at another elevation.\n\n'
        
        cmd += 'r.mapcalc "lowlevel_DEM = if((%(dem)s <= 0), 0, null())"\n\n' %{'dem': self.params.grassDEM}
        
        cmd += '# Identify all "lowlevel" cells connected to a basin outlet candidate.\n'
        cmd += '# This cost grow analysis should join mouths belonging to the same basin but separated at the shoreline .\n'
        
        cmd += 'r.cost -n input=lowlevel_DEM output=lowlevel_outlet_costgrow start_points=basin_outlet_pt max_cost=0 --overwrite\n\n'
        
        cmd += '# Clump all "lowlevel" cells associated with a basin outlet candidate (generates a unique id for each basin).\n'
        
        cmd += 'r.clump -d input=lowlevel_outlet_costgrow  output=lowlevel_outlet_clumps --overwrite\n\n'

        cmd += '# Export the basin outlet clumps by removing the "#" sign.\n'
        
        cmd += '# r.out.gdal -f input=lowlevel_outlet_clumps format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.lowlevelOutletClumpsfpn_s1}
        
        cmd += '# Identify beachfront river mouths\n'
        cmd += '# in order to identify separate mouths even if hydrologically connected inside of the shoreline .\n'

        cmd += 'r.mapcalc "shoreline_DEM = if((shoreline > 1 && %(dem)s <= 0), 0, null())" --overwrite\n\n' %{'dem': self.params.grassDEM}
        
        cmd += '# This cost grow analysis joins the cells of separate mouth outlets along the shoreline.\n'
 
        cmd += 'r.cost -n input=shoreline_DEM output=shoreline_outlet_costgrow start_points=basin_outlet_pt max_cost=0 --overwrite\n\n'
        
        cmd += '# Clump all shoreline mouth cells (generates a unique id for each mouth).\n'
        
        cmd += 'r.clump -d input=shoreline_outlet_costgrow output=shoreline_outlet_clumps --overwrite\n\n'
        
        cmd += '# Export the mouth outlet clumps by removing the "#" sign.\n'

        cmd += '# r.out.gdal -f input=shoreline_outlet_clumps format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.shorelineOutletClumpsfpn_s1}
        
        ### From this point the separate SFD and MFD data sources are not needed for a combined solutions. The coding is however included for convenience. ###
        
        cmd += '# Vectorize shoreline mouths (all) outlet raste as points\n'
        
        cmd += 'r.to.vect input=shoreline_outlet_clumps output=basin_all_outlets_pt type=point --overwrite\n\n'
        
        cmd += '# Add columns to vector databases\n'
        
        cmd += 'v.db.addcolumn map=basin_all_outlets_pt columns="upstream DOUBLE PRECISION"\n\n'
        
        cmd += 'v.db.addcolumn map=basin_all_outlets_pt columns="sfdup DOUBLE PRECISION"\n\n'
        
        cmd += 'v.db.addcolumn map=basin_all_outlets_pt columns="mfdup DOUBLE PRECISION"\n\n'
        
        cmd += 'v.db.addcolumn map=basin_all_outlets_pt columns="elevation INT"\n\n'
        
        cmd += 'v.db.addcolumn map=basin_all_outlets_pt columns="mouth_id INT"\n\n'
        
        cmd += 'v.db.addcolumn map=basin_all_outlets_pt columns="basin_id INT"\n\n'
        
        cmd += '# or with multiple columns added in one command'
        
        cmd += '# v.db.addcolumn map=basin_SFD_outlets_pt columns="upstream DOUBLE PRECISION, elevation INT, mouth_id INT, basin_id INT"\n\n'
        
        cmd += '# v.db.addcolumn map=basin_MFD_outlets_pt columns="upstream DOUBLE PRECISION, elevation INT, mouth_id INT, basin_id INT"\n\n'
        
        cmd += '# Add data to vector databases\n\n'
        
        cmd += 'v.what.rast map=basin_all_outlets_pt raster=SFD_upstream column=sfdup\n\n'
        
        cmd += 'v.what.rast map=basin_all_outlets_pt raster=MFD_upstream column=mfdup\n\n'
        
        cmd += 'v.what.rast map=basin_all_outlets_pt raster=%(dem)s column=elevation\n\n' %{'dem': self.params.grassDEM}
        
        cmd += 'v.what.rast map=basin_all_outlets_pt raster=shoreline_outlet_clumps column=mouth_id\n\n'
        
        cmd += 'v.what.rast map=basin_all_outlets_pt raster=lowlevel_outlet_clumps column=basin_id\n\n'
        
        cmd += 'db.execute sql="UPDATE basin_all_outlets_pt SET upstream=mfdup"\n\n'
        
        cmd += 'db.execute sql="UPDATE basin_all_outlets_pt SET upstream=sfdup WHERE sfdup>mfdup"\n\n'
                
        cmd += '# v.what.rast map=basin_SFD_outlets_pt raster=SFD_upstream column=upstream\n\n'
        
        cmd += '# v.what.rast map=basin_SFD_outlets_pt raster=%(dem)s column=elevation\n\n' %{'dem': self.params.grassDEM}
        
        cmd += '# v.what.rast map=basin_SFD_outlets_pt raster=shoreline_outlet_clumps column=mouth_id\n\n'
        
        cmd += '# v.what.rast map=basin_SFD_outlets_pt raster=lowlevel_outlet_clumps column=basin_id\n\n'
        
        cmd += '# v.what.rast map=basin_MFD_outlets_pt raster=MFD_upstream column=upstream\n\n'
        
        cmd += '# v.what.rast map=basin_MFD_outlets_pt raster=%(dem)s column=elevation\n\n' %{'dem': self.params.grassDEM}   
        
        cmd += '# v.what.rast map=basin_MFD_outlets_pt raster=shoreline_outlet_clumps column=mouth_id\n\n'
        
        cmd += '# v.what.rast map=basin_MFD_outlets_pt raster=lowlevel_outlet_clumps column=basin_id\n\n'
        
        cmd += '# Upload the x and y coordinates to the vector db'
        
        cmd += 'v.to.db map=basin_all_outlets_pt option=coor columns=xcoord,ycoord\n\n'
        
        cmd += '# v.to.db map=basin_SFD_outlets_pt option=coor columns=xcoord,ycoord\n\n'
        
        cmd += '# v.to.db map=basin_MFD_outlets_pt option=coor columns=xcoord,ycoord\n\n'
        
        cmd += '# Export basin outlet candidates - "all outlets p"t required in ESRI Shape format for stage 2 processing\n'
        
        cmd += 'v.out.ogr type=point input=basin_all_outlets_pt format=ESRI_Shapefile output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.allOutletsfpn_s1}
        
        cmd += '# v.out.ogr type=point input=basin_SFD_outlets_pt format=ESRI_Shapefile output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.SFDOutletsfpn_s1}
        
        cmd += '# v.out.ogr type=point input=basin_MFD_outlets_pt format=ESRI_Shapefile output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.MFDOutletsfpn_s1}
        
        cmd += '# Build a virtual wall across the basin mouths that will be used to force the flow out\n'
        cmd += '# of the basin mouths to pass a single cell and facilitate hydrological balance calculations.\n\n'
        
        cmd += '# Buffer one cell from the identified mouths.\n\n'

        cmd += 'r.buffer input=shoreline_outlet_clumps output=mouth_buffer distances=%(dist)f units=meters --overwrite\n\n' %{'dist':self.params.adjacentDist}
        
        cmd += '# Reclass the buffer to get cells only on the "seaside" of the shoreline.\n'
        cmd += '# The wall will be too thick where it is not exactly horisontal or vertical.\n'
        
        cmd += 'r.mapcalc "thickwall = if((drain_terminal == 1 && mouth_buffer > 1), 1, null())" --overwrite\n\n'
        
        cmd += '# Export the thickwall by removing the "#" sign.\n'
        
        cmd += '# r.out.gdal -f format=GTiff type=Byte input=thickwall output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.thickwallfpn_s1} 
        
        cmd += '# Shrink the thickwall to 1 cell width.\n'
        cmd += '# Remove the thickwall from the terminal drainage.\n'
        
        cmd += 'r.mapcalc "drain_terminal_remain = if((drain_terminal == 1),  if ((isnull(thickwall)),1, null() ),null() ) " --overwrite\n\n' 

        cmd += '# Buffer from the remaining terminal drainage nd in over land, now including thickwall.\n'
        
        cmd += 'r.buffer input=drain_terminal_remain output=mouthshoreline distances=%(dist)f units=meters --overwrite\n\n' %{'dist':self.params.adjacentDist}
        
        cmd += '# Export by removing the "#" sign.\n'
        
        cmd += '# r.out.gdal -f format=GTiff type=Byte input=mouthshoreline output=output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.mouthshorelinefpn_s1}  
        
        cmd += '# Combine mouthshoreline and shorewall to get the final shorewall, set shorewall elevation to 9999.\n'
        
        cmd += 'r.mapcalc "shorewall = if((mouthshoreline > 1 && thickwall == 1), 9999, null())" --overwrite\n\n'
        
        cmd += '# Export the shorewall by removing the "#" sign.\n'
        
        cmd += '# r.out.gdal -f format=GTiff type=Byte input=shorewall output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.shorewallfpn_s1} 
        
        cmd += '# Vectorize the shorewall.\n'
        
        cmd += 'r.to.vect input=shorewall output=shorewall_pt type=point --overwrite\n\n'
        
        cmd += '# Export the shorewall vector - required in ESRI Shape format for stage 2 processing.\n'
        
        cmd += 'v.out.ogr type=point input=shorewall_pt format=ESRI_Shapefile output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.shorewallptfpn_s1} 
             
        cmd += '# Fill up the DEM with the holes created from the difference between thickwall and shorewall\n'
        
        cmd += 'r.mapcalc "fillholeDEM = if ((thickwall==1 && isnull(shorewall)),9999, null() ) " --overwrite\n\n'
        
        cmd += '# Export the fillholeDEM by removing the "#" sign.\n'

        cmd += '# r.out.gdal format=GTiff input=fillholeDEM output=%(fpn)s --overwrite\n\n' %{'fpn': self.params.shorefillDEMfpn_s1}        
        
        

        GRASSshF = open(self.GRASSshFPN,'w')
        
        GRASSshF.write(cmd)
        
        GRASSshF.close()