'''
Created on 6 Oct 2020

@author: thomasgumbricht
'''

# imports

from __future__ import division
import os
import sys
#from osgeo import ogr, osr

from params.be_params_v02 import Params
from ds_manage.datasource import DS
from basin_extract import basin_extract_stage0 as s0
from basin_extract import basin_extract_stage1 as s1
from basin_extract import basin_extract_stage2 as s2
from basin_extract import basin_extract_stage4 as s4
#import shapely.wkt

class ListSrcDstDS(DS):
    '''
    classdocs
    '''
    
    def __init__(self, params):
        '''
        Constructor
        '''
        
        # Set all the input parameters
        self.params = params 
        
        # Set the script path
        self.scriptfp = os.path.join(self.params.srcfp)
        
        if not os.path.exists(self.scriptfp):
            
            infostr = " There are no scripts for stage %(stage)d, should have been \
                 at %(f)s" %{'stage':self.params.stage,'f': self.scriptfp}
                 
            print (infostr)
            
        # Set the path to the new basins and data to same as the source path
        self.basinfp = self.datafp = self.params.srcfp
        
        # Go to the stage for this run
        if self.params.stage == 0:
            
            print("    Files produced at stage 0 (DEM fixing)")
            
            if os.path.exists(self.hydoFillPtDEMFPN_s0):
                infostr = '        DEM with single cell filled pits DONE:\n        %s' %(self.hydoFillPtDEMFPN_s0)
            else:
                infostr = '        DEM with single cell filled pits MISSING'
            print (infostr)
            
            if os.path.exists(self.hydoFillAreaDEMFPN_s0):
                infostr = '        DEM with area filled pits DONE:\n        %s' %(self.hydoFillAreaDEMFPN_s0)
            else:
                infostr = '        DEM with area filled pits MISSING'
            print (infostr)
            
            if os.path.exists(self.fillDEMtiles_s0):
                infostr = '        Polygon file with tiling system for DEM fixing DONE:\n        %s' %(self.fillDEMtiles_s0)
            else:
                infostr = '        Polygon file with tiling system for DEM fixing MISSING'
            print (infostr)

        if self.params.stage == 2:
            
            pass
                     
                
        if self.params.stage == 4:
            
            pass

class BasinExtract(DS):
    '''
    classdocs
    '''
    
    def __init__(self, params):
        '''
        Constructor
        '''
        
        # Initiate the data source (DS) manager
        DS.__init__(self)
        
        # Set all the input parameters
        self.params = params 
        
        # Set the script path
        self.scriptfp = os.path.join(self.params.fp)
        
        if not os.path.exists(self.scriptfp):
            
            os.makedirs(self.scriptfp)
            
        # Set the path to the new basins and data to same as the source path
        self.basinfp = self.datafp = self.params.fp
        
        # Go to the stage for this run
        if self.params.stage == 0:
            s0.ProcStage0(self.params)
            
        if self.params.stage == 1:
            s1.ProcStage1BasinOutlets(self.params)
        
        if self.params.stage == 2:
            
            # The target file was deleted if <overwrite> set to True
            if not os.path.exists(self.params.BasinOutletFpn_s2):
                
                self.ProcStage2()
                
            else:

                self.ProcStage2()


                                
        if self.params.stage == 4:
            
            self.ProcStage4()
                         
    def ProcStage2(self): 
        '''
        '''
        self.params.stage2datafp = os.path.join(self.datafp, 'stage2')
        
        self.params.stage2scriptfp = os.path.join(self.params.stage2datafp, 'script')
        
        if not os.path.exists(self.params.stage2scriptfp):
            
            os.makedirs(self.params.stage2scriptfp)
         
        s2.CreateGRASSprocess(self.params)
                        
    def ProcStage4(self): 
        '''
        '''

        self.stage4datafp = os.path.join(self.datafp, 'stage4')
        self.stage4scriptfp = os.path.join(self.stage4datafp, 'script')
        if not os.path.exists(self.stage4scriptfp):
            os.makedirs(self.stage4scriptfp)
        
        if not os.path.exists(self.params.Outletsfpn_s1):
            exitstr ='EXITING: the output file from stage 1 %s is missing.\n' %(self.params.Outletsfpn_s1)

            print (exitstr)
            SNULLE
        
        # Start stage 4
        s4.CleanBasinPolys(self.params)  
            
    
                     
if __name__ == "__main__":
    
    xmlFPN = '/Volumes/GRASS2020/GRASSproject/SRTM/region/basin/amazonia_drain_filldem1cell/0/stage0/script/grass_drainage_outlets_stage0.xml'
    xmlFPN = '/Volumes/GRASS2020/GRASSproject/SRTM/region/basin/amazonia_test_20201220/xml/grass_drainage_outlets_stage0.xml'
    #xmlFPN = '/Volumes/GRASS2020/GRASSproject/SRTM/region/basin/amazonia_test_20201220/xml/grass_drainage_outlets_stage1.xml'
    #xmlFPN = '/Volumes/GRASS2020/GRASSproject/SRTM/region/basin/amazonia_test_20201220/xml/grass_drainage_outlets_stage2.xml'
    #xmlFPN = '/Volumes/GRASS2020/GRASSproject/SRTM/region/basin/amazonia_test_20201220/xml/grass_drainage_outlets_stage4.xml'
    
    #ADD ROUTINE FOR FILLING OF THE NO DATA FROM STAGE 0
    
    
    #xmlFPN = '/Volumes/GRASS2020/GRASSproject/ESA-DUE/region/basin/nordichydro_v01/xml/grass_drainage_outlets_stage0.xml'
    #xmlFPN = '/Volumes/GRASS2020/GRASSproject/ESA-DUE/region/basin/nordichydro_v01/xml/grass_drainage_outlets_stage1.xml'
    #xmlFPN = '/Volumes/GRASS2020/GRASSproject/ESA-DUE/region/basin/nordichydro_v01/xml/grass_drainage_outlets_stage2.xml'
    xmlFPN = '/Volumes/GRASS2020/GRASSproject/ESA-DUE/region/basin/nordichydro_v01/xml/grass_drainage_outlets_stage4.xml'
        
    params = Params(xmlFPN)
    
    if params.filecheck:
        
        # list the source and destiation files
        
        ListSrcDstDS(params)
    
    else:
        # Start main process
        
        BasinExtract(params)
    