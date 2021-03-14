'''
Created on 18 Nov 2020

@author: thomasgumbricht
'''

import os
from ds_manage.datasource import DS

#import grass.script as grass

class ProcStage0(DS):
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
        
        # Set the parameter that determines whether to fill only single cells or areas (-f flag in GRASS package)
        if not self.params.fillDirCells:
            self.params.vt = False
        elif self.params.fillDirCells == 1:
            self.params.vt = 'pt'
        else:
            self.params.vt = 'area'
         
        # Set the path to the new basins and data to same as the source path
        self.basinfp = self.datafp = self.params.fp
    
        if self.params.verbose:
            inforstr = '        Stage 0: Patch up DEM'
            print (inforstr)
            
        self.stage0datafp = os.path.join(self.datafp, 'stage0')
        self.stage0scriptfp = os.path.join(self.stage0datafp, 'script')
        if not os.path.exists(self.stage0scriptfp):
            os.makedirs(self.stage0scriptfp)
    
        self._TerminalDrainage()
        
        self._FillDirTiling()
        
        self._PatchSuperimpose()
        
        cmd = '# To run the output scripts you must have GRASS GIS open in a Terminal window session.\n\n'
        cmd += '# Script for removing land locked no-data regions.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += 'chmod 755 %(fpn)s\n' %{'fpn': self.GRASS0shFPN}
        cmd += '# Then execute the command form your GRASS terminal session:\n'
        cmd += '%(fpn)s\n\n'%{'fpn': self.GRASS0shFPN}
        
        print (cmd)
        
        cmd = '# Script for tiled hydrolgoical DEM correction.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += 'chmod 755 %(fpn)s\n' %{'fpn': self.GRASS1shFPN}
        cmd += '# Then execute the command form your GRASS terminal session:\n'
        cmd += '%(fpn)s\n\n'%{'fpn': self.GRASS1shFPN}
        
        print (cmd)
        
        cmd = '# Two alternative scripts for superimposing hydrological corrections.\n\n'
        cmd += '# Alternativ 1: GRASS processing (slower).\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += 'chmod 755 %(fpn)s\n' %{'fpn': self.GRASS2shFPN}
        cmd += '# Then execute the command form your GRASS terminal session:\n'
        cmd += '%(fpn)s\n\n'%{'fpn': self.GRASS2shFPN}
        
        print (cmd)
        
        cmd = '# Alternativ 2: GDAL and ogr2ogr processing (faster, requires that GDAL is installed).\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += 'chmod 755 %(fpn)s\n' %{'fpn': self.GRASS3shFPN}
        cmd += '# Then execute the command form your GRASS terminal session:\n'
        cmd += '%(fpn)s\n\n'%{'fpn': self.GRASS3shFPN}
        
        print (cmd)
    
    def _PatchSuperimpose(self):
        '''
        '''
        if not self.params.fillDirCells and not self.params.invertedFillDirCells:
            
            infostr = '    No filling of pits (fillDirCell == 0)'
            
            print (infostr)
            
            return
        
        OGRshFN = '%(s)s_ogr_superimpose-dem-filldir-%(vt)s_stage0-step2.sh' %{'s':self.params.region, 'vt':self.params.vt}
        self.OGRshFPN  = OGRshFPN = os.path.join(self.stage0scriptfp, OGRshFN)
        self.OGRshF = open(OGRshFPN,'w')
        
        cmd = '# Basin delineation: Superimpose vector tiles with dem fill values using GDAL \n'
        cmd += '# Created by the Python package basin_extract\n\n'
                       
        cmd += '# Run this script from any terminal window session.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += '# chmod 755 %(fpn)s\n' %{'fpn': self.OGRshFPN}
        cmd += '# Then execute the command form your terminal:\n'
        cmd += '# $ %(fpn)s\n\n'%{'fpn': self.OGRshFPN}
        
        cmd += '# Copy the DEM to use as baseline for burning the DEM fill vectors\n'
        
        cmd += 'cp %(src)s %(dst)s\n\n' %{'src':self.params.inlandCompDEMFPN_s0, 'dst':self.params.hydroFillDEMFPN_s0}
        
        if self.params.fillDirCells:
            
            if self.params.fillDirCells == 1:
                
                fillid = 'inland_fill_pt'
                
                vOutFpn = self.params.hydroFillPtFPN_s0
                
                searchstr = 'hydro-fill-pt-*.shp'
                
            elif self.params.fillDirCells > 1:
                
                fillid = 'inland_fill_area'
                
                vOutFpn = self.params.hydroFillAreaFPN_s0
                
                searchstr = 'hydro-fill-area-*.shp'
            
            cmd += '# Shell script for patching upp all filled tiles\n'
            
            cmd += 'cd %(fp)s\n\n' %{'fp':self.stage0datafp} 
            
            cmd += 'files=(%(search)s)\n\n' %{'search':searchstr}
            
            cmd += 'first=${files[0]}\n\n'
    
            cmd += 'files=("${files[@]:1}")\n\n'
            
            cmd += 'ogr2ogr -skipfailures %(dstFPN)s "$first"\n\n' %{'dstFPN': vOutFpn}
    
            cmd += 'for file in "${files[@]}"; do\n'
            
            cmd += '    echo "$file"\n'
    
            cmd += '    ogr2ogr -append -skipfailures %(dstFPN)s "$file"\n\n' %{'dstFPN': vOutFpn}
    
            cmd += 'done\n'
    
            cmd += '# Rasterize the patched vectors with filled pits\n'
            
            if len(self.params.fillSQL) <= 1: # Whithout SQL
                   
                cmd += 'GDAL_rasterize -a filldem %(src)s %(dst)s\n\n'  %{'src':vOutFpn, 'dst':self.params.hydroFillDEMFPN_s0}           
            
            else: #With SQL
                
                cmd += 'GDAL_rasterize -a filldem -where \'%(sql)s\' %(src)s %(dst)s\n\n'  %{'src':self.params.hydroFillPtFPN_s0, 
                                        'dst':vOutFpn, 'sql':self.params.fillSQL}     
            
        if self.params.invertedFillDirCells:
            
            if self.params.invertedFillDirCells == 1:
                
                fillid = 'inverted_fill_pt'
                
                vOutFpn = self.params.invertedFillPtFPN_s0
                
                searchstr = 'inverted-fill-pt-*.shp'
                
            elif self.params.fillDirCells > 1:
                
                fillid = 'inverted_fill_area'
                
                vOutFpn = self.params.invertedFillAreaFPN_s0
                
                searchstr = 'inverted-fill-area-*.shp'
            
            cmd += '# flatten the peaks, it is recommended that you set a query to restrict the filling\n'
            cmd += '# thresholds below \'nbupmax>250 AND nbupmax<50000 AND nbupq3<250\' set to press down single peaks near streams and small\n'
            cmd += '# rivers, (nbmax>250 AND nbmax<50000) restricted to positions that have water only in one adjacent cell (nbup3<250).\n'
            cmd += '# The fill value is set to the 1st quantile of the DEM of the 3x3 cell neighborhood (nbdemq1).\n'
            
            
            cmd += 'files=(%(search)s)\n\n' %{'search':searchstr}
            
            cmd += 'first=${files[0]}\n\n'
    
            cmd += 'files=("${files[@]:1}")\n\n'
            
            cmd += 'ogr2ogr -skipfailures %(dstFPN)s "$first"\n\n' %{'dstFPN': vOutFpn}
    
            cmd += 'for file in "${files[@]}"; do\n'
            
            cmd += '    echo "$file"\n'
    
            cmd += '    ogr2ogr -append -skipfailures %(dstFPN)s "$file"\n\n' %{'dstFPN': vOutFpn}
    
            cmd += 'done\n'
    
            cmd += '# Rasterize the patched vectors with flattened peaks\n'
            
            if len(self.params.invertedFillSQL) <= 1:
                cmd += 'GDAL_rasterize -a %(fillval)s %(src)s %(dst)s\n' %{'src':vOutFpn, 
                                        'dst':self.params.hydroFillDEMFPN_s0, 'fillval':self.params.invertedFillValue}
            else:
                cmd += 'GDAL_rasterize -a %(fillval)s -where \'%(sql)s\' %(src)s %(dst)s\n' %{'src':vOutFpn, 
                                        'dst':self.params.hydroFillDEMFPN_s0, 'fillval':self.params.invertedFillValue, 'sql':self.params.invertedFillSQL}           
     
            self.OGRshF.write(cmd)
                        
            self.OGRshF.close()
        
        # Set the names of the script and text files
        GRASS2shFN = '%(s)s_grass_patch-tiles-dem-filldir-%(vt)s_stage0-step2.sh' %{'s':self.params.region, 'vt':self.params.vt}
        
        self.GRASS2shFPN = os.path.join(self.stage0scriptfp, GRASS2shFN)

        GRASS2shF = open(self.GRASS2shFPN,'w')
        
        cmd = '# Basin delineation: Patch up and superimpose filled cells from r.fill.dir using GRASS\n'
        cmd += '# Created by the Python package basin_extract\n\n'
                        
        cmd += '# To run this script you must have GRASS GIS open in a Terminal window session.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += '# chmod 755 %(fpn)s\n' %{'fpn': self.GRASS2shFPN}
        cmd += '# Then execute the command from your GRASS terminal session:\n'
        cmd += '# GRASS 7.x.x ("region"):~ > %(fpn)s\n\n'%{'fpn': self.GRASS2shFPN}
        
        GRASS2shF.write(cmd)
        
        cmd = '# Create a shell variable of all the input vectors\n'
        
        if self.params.fillDirCells:
        
            if self.params.fillDirCells == 1:
                
                fillid = 'inland_fill_pt'
                
                vOutFpn = self.params.hydroFillPtFPN_s0
                
            elif self.params.fillDirCells > 1:
                
                fillid = 'inland_fill_area'
                
                vOutFpn = self.params.hydroFillAreaFPN_s0
                
            cmd += 'MAPS=$(g.list type=vector separator=comma pat="%(id)s_*")\n\n' %{'id': fillid}
 
            cmd += '# Patch together all the input vectors\n'
            
            cmd += 'v.patch -e input=$MAPS output=%(id)s\n\n' %{'id': fillid}
            
            cmd += '# Rasterize the patched vectors to the value of column filldem\n'
            
            if len(self.params.fillSQL) > 1:
                           
                cmd += 'v.to.rast input=%(id)s where="%(sql)s" output=%(id)s use=attr attribute_column=filldem --overwrite\n\n' %{'id': fillid, 'sql':self.params.fillSQL}
                
            else:
                
                cmd += 'v.to.rast input=%(id)s output=%(id)s use=attr attribute_column=filldem --overwrite\n\n' %{'id': fillid}
            
            cmd += '# Superimpose the rasterized filldir DEM values \n'
                        
            cmd += 'r.mapcalc "hydrofill_DEM  = if( isnull(%(id)s),inland_comp_DEM, %(id)s)" --overwrite\n\n' %{'id': fillid} 
           
            cmd += '# Export the patched vector as an ESRI shape file by removing the "#"\n'
                        
            cmd += '# v.out.ogr input=%(id)s type=point format=ESRI_Shapefile output=%(out)s \n\n' %{'id': fillid, 'out':vOutFpn}
          
        if self.params.invertedFillDirCells:
            
            cmd += '# Create a shell variable of all the inverted flattening vectors\n'
            
            if self.params.invertedFillDirCells == 1:
                
                fillid = 'inverted_fill_pt'
                
                vOutFpn = self.params.invertedFillPtFPN_s0
                
            elif self.params.invertedFillDirCells > 1:
                
                fillid = 'inverted_fill_area'
                
                vOutFpn = self.params.invertedFillAreaFPN_s0
            
            cmd += 'MAPS=$(g.list type=vector separator=comma pat="%(id)s_*")\n' %{'id': fillid} 
        
            cmd += '# Patch together all the input vectors\n'
                
            cmd += 'v.patch -e input=$MAPS output=%(id)s\n'  %{'id': fillid} 
               
            cmd += '# Rasterize the patched vectors to the value of column filldem\n'
            
            if len(self.params.invertedFillSQL) > 1:
                           
                cmd += 'v.to.rast input=%(id)s where="%(sql)s" output=%(id)s use=attr attribute_column=%(val)s --overwrite\n\n' %{'sql':self.params.invertedFillSQL,
                                                        'val':self.params.invertedFillValue,'id': fillid}
                
            else:
                
                cmd += 'v.to.rast input=%(id)s output=%(id)s use=attr attribute_column=%(val)s --overwrite\n\n' %{'sql':self.params.invertedFillValue, 'id': fillid}

            cmd += '# Superimpose the rasterized filldir DEM values \n'
                 
            if  self.params.fillDirCells:
                
                cmd += 'r.mapcalc "hydroflat_DEM  = if( isnull(%(id)s),hydrofill_DEM, %(id)s)" --overwrite\n\n'  %{'id': fillid}
            
                cmd += 'r.mapcalc "hydrofill_DEM  =  hydroflat_DEM" --overwrite\n\n'
                
            else:    
                        
                cmd += 'r.mapcalc "hydrofill_DEM  = if( isnull(%(id)s),inland_comp_DEM, %(id)s)" --overwrite\n\n'  %{'id': fillid} 
    
            cmd += '# Export the patched vector as an ESRI shape file\n'
                            
            cmd += '# v.out.ogr input=%(id)s type=point format=ESRI_Shapefile output=%(out)s\n\n' %{'id': fillid, 'out':vOutFpn}
                                
        cmd += '# Export the hydrologically corrected DEM as a geotiff\n'
        
        cmd += 'r.out.gdal -f input=hydrofill_DEM format=GTiff output=%(out)s  --overwrite\n\n'  %{'out': self.params.hydroFillDEMFPN_s0}

        GRASS2shF.write(cmd)
                    
    def _FillDirTiling(self):
        '''
        '''
          
        # Set the script file names
        GRASS1shFN = '%(s)s_grass_dem-filldir-%(vt)s_stage0-step1.sh' %{'s':self.params.region, 'vt':self.params.vt}
        #GDALshFN = '%(s)s_GDAL_patch-tiles-dem-filldir-%(vt)s_stage0-step2.sh' %{'s':self.params.region,'vt':self.params.vt}
        
        self.GRASS1shFPN = GRASS1shFPN = os.path.join(self.stage0scriptfp, GRASS1shFN) 
        #self.GDALshFPN = GDALshFPN = os.path.join(self.stage0scriptfp, GDALshFN)
    
        self.GRASS1shF = open(GRASS1shFPN,'w')
        #self.GDALshF = open(GDALshFPN,'w')
        
        cmd = '# Basin delineation: Hydrological DEM fillup (r.fill.dir) with tiled DEM\n'
        cmd += '# Created by the Python package basin_extract\n\n'
        
        cmd += '# The tiling speeds up the process hundredfold compared to using the entire DEM.\n'
        cmd += '# Just make sure you have the parameter "tileCellOverlap" set to capture the pit' 
        cmd += '# sizes you want to fill (i.e. somewhat larger than the "FillDirCell" parameter). \n\n'
                
        cmd += '# To run this script you must have GRASS GIS open in a Terminal window session.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += '# chmod 755 %(fpn)s\n' %{'fpn': self.GRASS1shFPN}
        cmd += '# Then execute the command form your GRASS terminal session:\n'
        cmd += '# GRASS 7.x.3x ("region"):~ > %(fpn)s\n\n'%{'fpn': self.GRASS1shFPN}
        
        cmd += '# The script uses the filled DEM from stage0-step1 "inland_comp_DEM" as input.\n'
        cmd += '# To use another DEM as input manually edit the 1st r.mapcalc before the looping.\n\n'
        
        cmd += '# Copy the DEM to use as baseline for burning the DEM fill vectors\n'
        cmd += 'cp %(src)s %(dst)s\n\n' %{'src':self.params.inlandCompDEMFPN_s0, 'dst':self.params.hydroFillDEMFPN_s0}
        
        cmd += '# Set the region to the full map extent\n'
        cmd += 'g.region raster=%(dem)s\n\n' %{'dem': self.params.grassDEM}
        cmd += '# There is a bug in r.fill.dir and you must reclass all no data to a very low elevation (-1000)\n'
        cmd += 'r.mapcalc "rfilldir_DEM = if(isnull(inland_comp_DEM), -1000, inland_comp_DEM )" --overwrite\n\n'
                        
        if self.params.invertedFillDirCells:
            cmd += '#r.watershed to get upstream areas for the inverted height suppresion\n'
            cmd += 'r.watershed -a elevation=inland_comp_DEM accumulation=MFD_upstream_filldir threshold=%(th)d --overwrite\n\n' %{'dem':self.params.grassDEM , 'th':self.params.basinCellThreshold}

        self.GRASS1shF.write(cmd)

        tilesize = self.params.tileCellSize
        overlap = self.params.tileCellOverlap
        
        # featureD is a dict that hold the id and geometry of each polygon
        featureD = {}
        
        for c in range(0,int(self.params.cols/tilesize)+1):
            rW = max(self.params.west,self.params.west+c*tilesize*self.params.ewres-overlap*self.params.ewres)
            rE = min(self.params.east,tilesize*self.params.ewres+self.params.west+c*tilesize*self.params.ewres+overlap*self.params.ewres)
            for r in range(0,int(self.params.rows/tilesize)+1):
                rS = max(self.params.south,self.params.south+r*tilesize*self.params.nsres-overlap*self.params.nsres)
                rN = min(self.params.north,tilesize*self.params.nsres+self.params.south+r*tilesize*self.params.nsres+overlap*self.params.nsres)
                
                # Reset region
                
                cmd = 'g.region -a n=%(n)f s=%(s)f e=%(e)f w=%(w)f\n' %{'n':rN, 's':rS, 'e':rE, 'w':rW}
                
                cmd += 'data="$(r.stats -p input=inland_comp_DEM null_value=\'null\')"\n'
                cmd += 'if echo "$data" | grep -q "null\s100"; then\n'
                cmd += '    echo "Null tile - skipping"\n'
                cmd += 'else\n'
                cmd += '    echo "Valid tile"\n'
                
                self.GRASS1shF.write(cmd)
                
                tile_id = '%(c)d_%(r)d' %{'c':c,'r':r}
                
                featureD[tile_id] = {'id':tile_id, 'xmin':rW, 'xmax':rE, 'ymin':rS, 'ymax':rN}
                
                if self.params.fillDirCells:
                
                    if self.params.fillDirCells == 1:
                        
                        self._FillDirPoint(tile_id)
                        
                    else:
                        
                        self._FillDirArea(tile_id)
                      
                cmd = 'fi\n'
                
                self.GRASS1shF.write(cmd)
                
                if self.params.invertedFillDirCells:
                    
                    if self.params.invertedFillDirCells == 1:
                    
                        self._InvertedFillDirPoint(tile_id)
                    
                    else:
                    
                        NOTYET
                        
                '''
                # Remove all raster layers - takes too much space
                
                cmd = 'g.remove -f type=raster pattern="*_%(tid)s"\n' %{'tid':tile_id}
                
                GRASSshF.write(cmd)
                
                print (cmd)
                '''

        self.GRASS1shF.write('\n')
                
        cmd = 'g.region raster=%(dem)s\n\n' %{'dem': self.params.grassDEM}
        
        self.GRASS1shF.write(cmd)
               
        self.GRASS1shF.close()
        
        # Save the tiles used for the filldir 
        self.WriteStage0DS(self.params.fillDEMtiles_s0, None, featureD)
           
    def _FillDirPoint(self, tile_id):
        '''
        '''
  
        ptshpfn = 'inland-fill-pt-%(tid)s.shp' %{'tid':tile_id}
                
        ptshpfpn = os.path.join(self.stage0datafp, ptshpfn)
            
        # For single cells only with the -f flag    
        cmd = '    r.fill.dir -f input=rfilldir_DEM output=hydro_cellfill_DEM_%(tid)s direction=hydro_cellfill_draindir_%(tid)s areas=hydro_cellfill_problems_%(tid)s --overwrite\n' %{'tid':tile_id}
        
        # Check differences
        cmd += '    r.mapcalc "DEM_cellfill_diff_%(tid)s = hydro_cellfill_DEM_%(tid)s - rfilldir_DEM" --overwrite\n' %{'tid':tile_id}
    
        # get only diff areas
        cmd += '    r.mapcalc "inland_fill_cell_%(tid)s = if(DEM_cellfill_diff_%(tid)s != 0, 1, null() )" --overwrite\n' %{'tid':tile_id}
        
        # Convert to vector (area and pt)
        cmd += '    r.to.vect input=inland_fill_cell_%(tid)s output=inland_fill_pt_%(tid)s type=point --overwrite\n' %{'tid':tile_id}
        
        # Extract DEM fill value for points
        cmd += '    v.db.addcolumn map=inland_fill_pt_%(tid)s columns="filldem DOUBLE PRECISION" \n' %{'tid':tile_id}

        # Extract data from hydro_cellfill_DEM                 
        cmd += '    v.what.rast map=inland_fill_pt_%(tid)s column=filldem raster=hydro_cellfill_DEM_%(tid)s\n' %{'tid':tile_id}

        # Export always
        cmd += '    v.out.ogr input=inland_fill_pt_%(tid)s type=point format=ESRI_Shapefile output=%(pthshpfpn)s --overwrite\n\n' %{'tid':tile_id,'pthshpfpn':ptshpfpn}
                    
        self.GRASS1shF.write(cmd)
                      
    def _FillDirArea(self, tile_id): 
        '''
        '''
        areashpfn = 'inland-fill-area-%(tid)s.shp' %{'tid':tile_id}
                
        areashpfpn = os.path.join(self.stage0datafp, areashpfn) 
        
        # For all fill areas  without the _-f_ flag            
        cmd = 'r.fill.dir input=rfilldir_DEM output=hydro_areafill_DEM_%(tid)s direction=hydro_areafill_draindir_%(tid)s areas=hydro_areafill_problems_%(tid)s --overwrite\n' %{'tid':tile_id}
    
        # Check differences     
        cmd += 'r.mapcalc "DEM_areafill_diff_%(tid)s = hydro_areafill_DEM_%(tid)s - rfilldir_DEM" --overwrite\n' %{'tid':tile_id}
        
        # get only diff areas 
        cmd += 'r.mapcalc "inland_fill_area_%(tid)s = if(DEM_areafill_diff_%(tid)s != 0, 1, null() )" --overwrite\n' %{'tid':tile_id}
         
        # Convert to vector (area)
        cmd += 'r.to.vect input=inland_fill_area_%(tid)s output=inland_fill_area_%(tid)s type=area --overwrite\n' %{'tid':tile_id}
        
        # Get area of fill
        cmd += 'v.to.db map=inland_fill_area_%(tid)s type=centroid option=area columns=area_km2 units=kilometers\n' %{'tid':tile_id}

        cmd += 'v.db.addcolumn map=inland_fill_area_%(tid)s columns="filldem DOUBLE PRECISION" \n' %{'tid':tile_id}

        # Extract data from                    
        cmd += 'v.what.rast type=centroid map=inland_fill_area_%(tid)s column=filldem raster=hydro_areafill_DEM_%(tid)s\n' %{'tid':tile_id}

        # Export always, this layer to be used in GDAL_rasterize - simplest way to correct data  
        cmd += 'v.out.ogr input=inland_fill_area_%(tid)s type=area format=ESRI_Shapefile output=%(areahshpfpn)s --overwrite\n' %{'tid':tile_id,'areahshpfpn':areashpfpn}
         
        if tile_id == "0_0":
            
            cmd += '    ogr2ogr -skipfailures %(dstFPN)s %(srcFPN)s\n' %{'dstFPN': self.params.hydroFillAreaFPN_s0, 'srcFPN':areashpfpn}
        
        else:
            
            cmd += '    ogr2ogr -append  -skipfailures %(dstFPN)s %(srcFPN)s\n' %{'dstFPN': self.params.hydroFillAreaFPN_s0, 'srcFPN':areashpfpn}
            
        self.GDALshF.write(cmd)
        
    def _InvertedFillDirPoint(self, tile_id):
        '''
        '''
  
        ptshpfn = 'inverted-fill-pt-%(tid)s.shp' %{'tid':tile_id}
                
        ptshpfpn = os.path.join(self.stage0datafp, ptshpfn)
 
        
        cmd = '    r.mapcalc "invertedDEM_%(tid)s = 10000-rfilldir_DEM" --overwrite\n' %{'tid':tile_id}
                
        cmd += '    r.fill.dir -f input=invertedDEM_%(tid)s output=inverted_cellfill_DEM_%(tid)s direction=inverted_cellfill_draindir_%(tid)s areas=inverted_cellfill_problems_%(tid)s --overwrite\n' %{'tid':tile_id}
    
        # Check differences
        cmd += '    r.mapcalc "inverted_cellfill_diff_%(tid)s = inverted_cellfill_DEM_%(tid)s - invertedDEM_%(tid)s" --overwrite\n' %{'tid':tile_id}
        
        # get only diff areas
        cmd += '    r.mapcalc "inverted_fill_cell_%(tid)s = if(inverted_cellfill_diff_%(tid)s != 0, 1, null() )" --overwrite\n' %{'tid':tile_id}
        
        # Convert to vector (area and pt)
        cmd += '    r.to.vect input=inverted_fill_cell_%(tid)s output=inverted_fill_pt_%(tid)s type=point --overwrite\n' %{'tid':tile_id}

        #Extract DEM fill value for points
        cmd += '    v.db.addcolumn map=inverted_fill_pt_%(tid)s\
         columns="invdem DOUBLE PRECISION, filldem DOUBLE PRECISION, nbupmax DOUBLE PRECISION, nbupq3 DOUBLE PRECISION, nbdemq1 DOUBLE PRECISION, nbdemmin DOUBLE PRECISION" \n' %{'tid':tile_id}

        # Extract data from hydro_cellfill_DEM             
        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=invdem raster=inverted_cellfill_DEM_%(tid)s\n' %{'tid':tile_id}
        
        cmd += '    v.db.update map=inverted_fill_pt_%(tid)s column=filldem query_column="10000-(invdem)"\n' %{'tid':tile_id}

        cmd += '    r.neighbors input=inland_comp_DEM selection=inverted_fill_cell_%(tid)s\
         output=dem_neighbor_quart1_%(tid)s,dem_neighbor_min_%(tid)s\
          method=quart1,minimum --overwrite\n' %{'tid':tile_id}
          
        cmd += '    r.neighbors input=MFD_upstream_filldir selection=inverted_fill_cell_%(tid)s\
         output=inverted_neighbor_max_%(tid)s,inverted_neighbor_quart3_%(tid)s\
          method=maximum,quart3 --overwrite\n' %{'tid':tile_id}
        
        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=nbupmax raster=inverted_neighbor_max_%(tid)s\n' %{'tid':tile_id}
                
        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=nbupq3 raster=inverted_neighbor_quart3_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=nbdemq1 raster=dem_neighbor_quart1_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.what.rast map=inverted_fill_pt_%(tid)s column=nbdemmin raster=dem_neighbor_min_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.out.ogr input=inverted_fill_pt_%(tid)s type=point format=ESRI_Shapefile output=%(pthshpfpn)s --overwrite\n\n' %{'tid':tile_id,'pthshpfpn':ptshpfpn}
     
        #cmd = 'GDAL_rasterize -a filldem %(src_datasource)s %(dst)s\n'  %{'src_datasource':ptshpfpn, 'dst':self.params.invertedFillDEMFPN_s0}           
    
        #self.GDALshF.write(cmd)
        
        self.GRASS1shF.write(cmd)
        
    def _InvertedFillDirArea(self, tile_id):
        '''
        '''
  
        areashpfn = 'inverted-fill-area-%(tid)s.shp' %{'tid':tile_id}
                
        areashpfpn = os.path.join(self.stage0datafp, areashpfn)
 
        
        cmd = '    r.mapcalc "invertedDEM_%(tid)s = 10000-rfilldir_DEM" --overwrite\n' %{'tid':tile_id}
                
        cmd += '    r.fill.dir input=invertedDEM_%(tid)s output=inverted_areafill_DEM_%(tid)s direction=inverted_areafill_draindir_%(tid)s areas=inverted_areafill_problems_%(tid)s --overwrite\n' %{'tid':tile_id}
    
        # Check differences
        cmd += '    r.mapcalc "inverted_areafill_diff_%(tid)s = inverted_areafill_DEM_%(tid)s - invertedDEM_%(tid)s" --overwrite\n' %{'tid':tile_id}
        
        # get only diff areas
        cmd += '    r.mapcalc "inverted_fill_area_%(tid)s = if(inverted_areafill_diff_%(tid)s != 0, 1, null() )" --overwrite\n' %{'tid':tile_id}
        
        # Convert to vector (area and pt)
        cmd += '    r.to.vect input=inverted_fill_area_%(tid)s output=inverted_fill_area_%(tid)s type=area --overwrite\n' %{'tid':tile_id}

        #Extract DEM fill value for points
        cmd += '    v.db.addcolumn map=inverted_fill_area_%(tid)s\
         columns="invdem DOUBLE PRECISION, filldem DOUBLE PRECISION, nbupmax DOUBLE PRECISION, nbupq3 DOUBLE PRECISION, nbdemq1 DOUBLE PRECISION, nbdemmin DOUBLE PRECISION" \n' %{'tid':tile_id}

        # Extract data from hydro_cellfill_DEM             
        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=invdem raster=inverted_areafill_DEM_%(tid)s\n' %{'tid':tile_id}
        
        cmd += '    v.db.update map=inverted_fill_area_%(tid)s column=filldem query_column="10000-(invdem)"\n' %{'tid':tile_id}

        cmd += '    r.neighbors input=inland_comp_DEM \
         output=dem_neighbor_quart1_%(tid)s,dem_neighbor_min_%(tid)s\
          method=quart1,minimum --overwrite\n' %{'tid':tile_id}
          
        cmd += '    r.neighbors input=MFD_upstream_filldir \
         output=inverted_neighbor_max_%(tid)s,inverted_neighbor_quart3_%(tid)s\
          method=maximum,quart3 --overwrite\n' %{'tid':tile_id}
        
        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=nbupmax raster=inverted_neighbor_max_%(tid)s\n' %{'tid':tile_id}
                
        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=nbupq3 raster=inverted_neighbor_quart3_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=nbdemq1 raster=dem_neighbor_quart1_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.what.rast map=inverted_fill_area_%(tid)s column=nbdemmin raster=dem_neighbor_min_%(tid)s\n' %{'tid':tile_id}

        cmd += '    v.out.ogr input=inverted_fill_area_%(tid)s type=point format=ESRI_Shapefile output=%(areahshpfpn)s --overwrite\n\n' %{'tid':tile_id,'areahshpfpn':areashpfpn}
     
        if tile_id == "0_0":
            
            cmd += '    ogr2ogr -skipfailures %(dstFPN)s %(srcFPN)s\n\n' %{'dstFPN': self.params.invertedFillAreaFPN_s0, 'srcFPN':areashpfpn}
            
        else:
            
            cmd += '    ogr2ogr  -append -skipfailures %(dstFPN)s %(srcFPN)s\n\n' %{'dstFPN': self.params.invertedFillAreaFPN_s0, 'srcFPN':areashpfpn}

        self.GRASS1shF.write(cmd)
                                 
    def _TerminalDrainage(self):
        ''' GRASS shell script for terminal drainage "no data" fixing
        '''
        
        GRASS0shFN = '%(s)s_grass_fix-terminal-drainage-nodata_stage0-step0.sh' %{'s':self.params.region}

        self.GRASS0shFPN = os.path.join(self.stage0scriptfp, GRASS0shFN)

        GRASSshF = open(self.GRASS0shFPN,'w')
        
        cmd = '# Basin delineation: Assign elevation to land locked void cells\n'
        cmd += '# Created by the Python package basin_extract\n\n'
                
        cmd += '# To run this script you must have GRASS GIS open in a Terminal window session.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += '# chmod 755 %(fpn)s\n' %{'fpn': self.GRASS0shFPN}
        cmd += '# Then execute the command form your GRASS terminal session:\n'
        cmd += '# GRASS 7.x.3x ("region"):~ > %(fpn)s\n\n'%{'fpn': self.GRASS0shFPN}
        
        GRASSshF.write(cmd)
        
        GRASSshF.write('# Create destination path\n')
        
        cmd = 'mkdir -p %s\n\n' % (self.stage0scriptfp)
        
        GRASSshF.write(cmd)
        

        if not os.path.exists(self.params.inputDEMFPN_s0):
            
            GRASSshF.write('# Mapcalc to convert input DEM to DEM without negative values\n')
            
            cmd = 'r.mapcalc "inputDEMonlyp = if(%(dem)s<0, 0, %(dem)s)" --overwrite\n\n' %{'dem':self.params.grassDEM}
            
            cmd += '# Set color ramp to elevation\n'

            cmd += 'r.colors map=inputDEMonlyp color=elevation\n\n' %{'dem':self.params.grassDEM}
        
            cmd += '# Export as geotiff \n'
            
            cmd += 'r.out.gdal format=GTiff  input=inputDEMonlyp \
            output=%(fpn)s\n\n' %{'fpn':self.params.inputDEMFPN_s0}
            
            GRASSshF.write(cmd)
        
        GRASSshF.write('# Create a map of the terminal drainage\n')
        
        cmd = 'r.mapcalc "drain_terminal = if(isnull(\'%(dem)s\'), 1, null())" --overwrite\n\n' %{'dem':self.params.grassDEM}

        GRASSshF.write(cmd)
        
        GRASSshF.write('# Clump agglomerated cells of the terminal drainage\n')
        
        cmd = 'r.clump input=drain_terminal output=terminal_clumps --overwrite\n\n'
        
        GRASSshF.write(cmd)
        
        GRASSshF.write('# Vectorize the clumps\n')
        
        cmd = 'r.to.vect input=terminal_clumps output=terminal_clumps type=area --overwrite\n\n'

        GRASSshF.write(cmd)
        
        GRASSshF.write('# Add the area of each clump to the clump polygons attribute table\n')
        
        cmd = 'v.to.db map=terminal_clumps type=centroid option=area columns=area_km2 units=kilometers\n\n'
        
        GRASSshF.write(cmd)
        
        GRASSshF.write('# Export all clumps as ESRI shape files by removing the "#"\n')
        
        cmd = '# v.out.ogr input=terminal_clumps type=area format=ESRI_Shapefile \
        output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.terminalClumpsFPN_s0}
        
        GRASSshF.write(cmd)
        
        GRASSshF.write('# Rasterize all clumps below a threshhold size.\n')
        GRASSshF.write('# These raster cells will be changed to land areas with initial elevation = 0\n')
        
        cmd = 'v.to.rast input=terminal_clumps type=area where="area_km2< %(th)f" \
        output=drain_terminal_small use=val value=0 --overwrite\n\n' %{'th':self.params.terminalClumpAreakm2Threshold}
        
        GRASSshF.write(cmd)
        
        GRASSshF.write('# Export all small (land) clumps as ESRI shape files by removing the "#"\n')
        
        cmd = '# v.out.ogr input=drain_terminal_small type=area format=ESRI_Shapefile \
        output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.terminalClumpsSmallFPN_s0}
        
        GRASSshF.write(cmd)
        
        GRASSshF.write('# Rasterize all clumps above a threshhold size.\n')
        
        cmd = 'v.to.rast input=terminal_clumps type=area where="area_km2 >= %(th)f" \
        output=drain_terminal_large use=val value=0 --overwrite\n\n' %{'th':self.params.terminalClumpAreakm2Threshold}
        
        GRASSshF.write(cmd)
        
        GRASSshF.write('# Export all large (terminal drainage) clumps as ESRI shape files by removing the "#"\n')
        
        cmd = '# v.out.ogr input=drain_terminal_large type=area format=ESRI_Shapefile \
        output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.terminalClumpsLargeFPN_s0}
        
        GRASSshF.write(cmd)
        
        GRASSshF.write('# Mapcalc a new composit filled DEM over all terrestrial regions\n')

        cmd = 'r.mapcalc "inland_comp_DEM = if(( isnull(%(dem)s)), drain_terminal_small,%(dem)s )" \
        --overwrite\n\n' %{'dem':self.params.grassDEM}
        
        GRASSshF.write(cmd)
        
        GRASSshF.write('# Export the terrestrial composite DEM as a geotiff\n')
        
        cmd = 'r.out.gdal format=GTiff  input=inland_comp_DEM \
        output=%(fpn)s\n\n' %{'fpn':self.params.inlandCompDEMFPN_s0}
        
        GRASSshF.write(cmd)
        
        GRASSshF.close()