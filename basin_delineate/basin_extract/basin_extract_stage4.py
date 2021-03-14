'''
Created on 18 Nov 2020

@author: thomasgumbricht
'''

from ds_manage.datasource import DS
from osgeo import ogr, osr
from os import path
from sys import exit
#import shapely
from shapely.geometry.polygon import LinearRing
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.polygon import Polygon
from shapely.ops import unary_union
from shapely.wkb import loads

class CleanBasinPolys(DS):
    '''
    Remove incomplete basins and merge subbasins
    '''
    
    def __init__(self, params):
        '''
        Expects a parameter object
        '''
        
        # Initiate the data source (DS) manager
        DS.__init__(self)
        self.params = params
              
        self.attributeD = {}
        
        if self.params.verbose:
            inforstr = '        Stage 4: Remove incomplete basins and merge subbasins'
            print (inforstr)
             
        if not path.exists(self.params.BasinAreasFpn_s2):
            exitstr = 'EXITING: File with preliminary basin polygon file missing, %s' %(self.params.BasinAreasFpn_s2)
            exit(exitstr)
             
        #NOTE THATAT THERE ARE TWO POINT POINt FILES FROM STAGE 1 AND STAGE 
        
        if not path.exists(self.params.BasinOutletFpn_s2):
            exitstr = 'EXITING: File with preliminary basin outlet points missing, %s' %(self.params.BasinOutletFpn_s2)
            exit(exitstr)
            
        if self.params.verbose:    
            inforstr = '        Reading preliminary basin polygons from %s' %(self.params.BasinAreasFpn_s2)
            print (inforstr)

        ds = self.OpenDs(self.params.BasinAreasFpn_s2)
        
        # Get the source layer
        self.srcLayer = ds.GetLayer()
        
        # Get the spatial reference
        spatialRef = self.srcLayer.GetSpatialRef()

        if self.params.verbose:
            featureCount = self.srcLayer.GetFeatureCount()
            print ("        Number of preliminary basins: %d" % (featureCount))
            
            print ("        Building database")
        
        # Get the layer boundary
        self._GetLayerBoundary()
        
        # Loop 1: build a dict of all mouths and basins
        self._BuildMouthBasinDict()
        
        if self.params.verbose:
            
            print ("        Identifying polygons touching DEM perimeter")
        
        # Loop 2: remove all src basins touching the boundary
        self._FindBndTouchSrcFeats()
        
        if self.params.verbose:
            
            print ("        Identifying polygons draining to the same outlet")
            
        # Loop 3: join mouth polygons belonging to the same mouth outlets (Polygones -> Multipolygons)
        self._UnionMouths()
        
        if self.params.verbose:
            
            print ("        Identifying multi-mouth basins draining to more than one outlet")

        # Loop 4: join basins with multiple mouths
        self._UnionBasins()
        
        if self.params.verbose:
            
            print ("        Recalculating attributes for altered polygons")
        
        # Loop 5: Recalculate area of united polygons
        self._AreaBasins()
        
        if self.params.verbose:
            
            print ("        Writing results to new data sources")
            printstr = "            as: %s" % (self.params.basin_areasFPN_s4) 
            print (printstr)
            print ("        Omitted polygons (basins) saved")
            printstr = "            as: %s" % (self.params.basin_omitted_areasFPN_s4) 
            print (printstr)
        
        # Save the basins
        self.WriteStage4DS(self.params.basin_areasFPN_s4, spatialRef, self.featD, False)
        
        # Save the omitted basins
        self.WriteStage4DS(self.params.basin_omitted_areasFPN_s4, spatialRef, self.featD, True)

        if self.params.verbose:
            
            print ("        Comparing preliminary outlets with updated basins")
            
        # Compare outlet points from stage 2 and basins from stage 4
        self._FindOutletsWithoutBasins()
 
    def _GetLayerBoundary(self):
        ''' Get the boundary of the layer
        '''
        
        # Convert the layer extent to corner coordinates
        minX, maxX, minY, maxY = self.srcLayer.GetExtent()
        
        self.boundary = LinearRing([(minX, maxY), (maxX, maxY), (maxX, minY), (minX, minY)])

    def _BuildMouthBasinDict(self):
        ''' Build a dictionary of the input basins and mouths
        '''
        
        self.featD = {}
        
        self.basinD = {}
        
        self.mouthD = {}
        
        for x, feature in enumerate(self.srcLayer):
            
            mouth_id = feature.GetField("mouth_id")
            
            basin_id = feature.GetField("basin_id")
            
            area_km2 = feature.GetField("area_km2")
            
            x_mouth = feature.GetField("x_mouth")
            
            y_mouth = feature.GetField("y_mouth")

            self.featD[x] =  {'mouth_id':mouth_id, 'basin_id':basin_id, 
                              'area_km2':area_km2, 'x_mouth':x_mouth, 'y_mouth': y_mouth, 
                              'bndtouch': False, 'skip':False, 'newgeom': False}
            
            # Build a dict of all basins
            if basin_id in self.basinD:
                
                self.basinD[basin_id].append({'mouth_id':mouth_id, 'feat_id':x})
            
            else:
                
                self.basinD[basin_id] = [{'mouth_id':mouth_id, 'feat_id':x}]
                
            # Build a dict of all mouths
            if mouth_id in self.mouthD:
                
                self.mouthD[mouth_id].append({'mouth_id':mouth_id, 'feat_id':x})
            
            else:
                
                self.mouthD[mouth_id] = [{'mouth_id':mouth_id, 'feat_id':x}]
        
        self.srcLayer.ResetReading()  # reset the read position to the start
            
    def _FindBndTouchSrcFeats(self):
        ''' Identify basin polygons that touch the layer boundary
        '''
        
        intersectL = []
        
        for x, feature in enumerate(self.srcLayer):

            geom = loads(feature.GetGeometryRef().ExportToWkb())
            
            if geom.intersects(self.boundary):
                
                intersectL.append(x)
                
                self.featD[x]['bndtouch'] = True
                
                # This poly should be skipped when writing the new basin vectors
                self.featD[x]['skip'] = True
                
                if self.params.verbose > 1:
                    print ('            Touching boundary: %s %s' %(x, self.featD[x]))

        self.srcLayer.ResetReading()  # reset the read position to the start
        
        # Loop over the boundary touching polygons to see if there are any other linked polygons
        
        for i in intersectL:
            
            basin_id = self.featD[i]['basin_id']
            
            for j in self.basinD[basin_id]:
                
                x = j['feat_id'] 

                self.featD[x]['bndtouch'] = True
                
    def _UnionMouths(self):
        ''' Union polygons sharing the same mouth
        '''

        for mouth_id in self.mouthD:
 
            featL = []
            
            if len(self.mouthD[mouth_id]) > 1:

                # Get the largest polygon
                area = 0
                
                for j in self.mouthD[mouth_id]: 
                    
                    featL.append( j['feat_id'] )
                    
                    if self.featD[ j['feat_id']]['area_km2' ] > area:
                        
                        area = self.featD[ j['feat_id']]['area_km2' ]
                        
                        largestpoly = j['feat_id']
                
                self._UnionPolygons(featL,largestpoly)
 
    def _UnionBasins(self):
        ''' Union polygons sharing the same basin
        '''
        for basin_id in self.basinD:  
            
            featL =[]
            
            if len(self.basinD[basin_id]) > 1:
                
                # Get the largest polygon
                area = 0
                
                for j in self.basinD[basin_id]:
                    
                    featL.append( j['feat_id'] )
                    
                    if self.featD[ j['feat_id']]['area_km2' ] > area:
                        
                        area = self.featD[ j['feat_id']]['area_km2' ]
                        
                        largestpoly = j['feat_id']

                self._UnionPolygons(featL,largestpoly)

                       
    def _UnionPolygons(self, featL, largest):
        ''' Union polygons
        '''
        # If any of the polygons to be united is touching the boundary, all are infected (to be skipped)
        
        boundaryinfected = False
        
        if self.params.verbose > 1:
            print ('            Union (target): %s' %(self.featD[largest]))
        
        for f in featL:
            
            if self.featD[f]['bndtouch']:
                
                boundaryinfected = True
                
                break
        
        if boundaryinfected:
            
            for f in featL:
                
                self.featD[f]['skip'] = True
                
            # No need to union the polygons, they will all be omitted
            return
        
        i = featL.index(largest)
        
        # pop the largest basin from feature list
        featL.pop(i)
                
        # Start by merging polygons that share at least the side of one(1) cell
        n = -1
        
        while True:
            
            n += 1

            doneL = []

            for i in range( n,len(featL) ):
                
                f = featL[i]
                
                if f in doneL:
                    
                    continue
                
                if self.params.verbose > 1:
                    print ('                adding: %s' %(self.featD[f]))
                
                if i == 0: # Only happen in the first loop
                    
                    if self.featD[largest]['newgeom']: # if the geom is already united
                        
                        geomA = self.featD[largest]['newgeom']
                        
                    else: # pristine geom from file
                        
                        featureA = self.srcLayer.GetFeature(largest)
                    
                        geomA = loads(featureA.GetGeometryRef().ExportToWkb())
                
                if self.featD[f]['newgeom']:
  
                    geomB = self.featD[f]['newgeom']
                    
                else:
    
                    featureB = self.srcLayer.GetFeature(f)
                
                    geomB = loads(featureB.GetGeometryRef().ExportToWkb())
                        
                iSectGeom = self._BndIntersection(geomA,geomB)
                
                if iSectGeom == None: # None = not even touching
                    
                    continue 
                
                if iSectGeom.length > 0: 
                    
                    unionL = [geomA,geomB]
                    
                    doneL.append(f)
                    
                    # This geom is now united with largestand should be omitted when writing to file
                
                    self.featD[f]['skip'] = True
                    
                    if unary_union(unionL).is_valid:

                        geomA = unary_union(unionL)
                    
                    else:
                        infostr = '    Invalid union for mouths %s and %s, with intesect length of %s' %(largest,f,iSectGeom.length)
                        print (infostr )


            if n >= len(featL):
                break
                   
        # In the second loop geoms that are touching only at a point are united
        if len(featL) > len(doneL): # Featureremainin that are just touching at single point
            
            for i in range(len(featL) ):
                
                f = featL[i]
                
                if f in doneL:
                    
                    continue
                
                # This geom should be skipped when writing the new basin vectors
                self.featD[f]['skip'] = True
                
                featureB = self.srcLayer.GetFeature(f)
                
                geomB = loads(featureB.GetGeometryRef().ExportToWkb())

                unionL = [geomA,geomB]
                
                if unary_union(unionL).is_valid:
                    
                    geomA = unary_union(unionL)
                
                else:
                    
                    infostr = '    Invalid union for mouths %s and %s, with no intersect length (%s)' %(largest,f,iSectGeom.length)
                    print (infostr )
                    
                    if geomA.geom_type == 'MultiPolygon':
                        SNULLE
                    elif geomA.geom_type == 'Polygon':
                        
                        geomA = MultiPolygon( unionL )
                        
                    else:
                        
                        raise IOError('Shape is not a polygon.')
                    
        self.featD[largest]['newgeom'] = geomA
        
       
    def _BndIntersection(self, geomA,geomB):
        '''
        '''
        usebnd = False
        
        if geomA.intersects(geomB.boundary):
            
            if not geomA.is_valid:
                
                usebnd = True

            if not geomB.is_valid:
                
                usebnd = True
                                
            if usebnd:
                
                iSectGeom = geomA.boundary.intersection(geomB.boundary)

                return iSectGeom
            
            else:
                
                iSectGeom = geomA.intersection(geomB)

                return iSectGeom
        
        else:
            return None
         
    def _AreaBasins(self): 
        ''' Recalculate area of united basins
        '''
        
        for f in self.featD:
            
            if self.featD[f]['skip']:
                continue
            
            if self.featD[f]['newgeom']:
                self.featD[f]['area_km2'] = self.featD[f]['newgeom'].area/1000000
               
            # All joined basins smaller than the threshold to be skipped
            if self.featD[f]['area_km2'] < self.params.basinAreaThreshold:
                
                self.featD[f]['skip'] = True
                
    def _FindOutletsWithoutBasins(self):
        ''' Identify outlets from stage 2 that do not belong to a basin
        '''
        
        if self.params.verbose:    
            inforstr = '            Reading preliminary basin outlet points from %s' %(self.params.BasinOutletFpn_s2)
            print (inforstr)
            
        # Get the outlet points      
        dspt = self.OpenDs(self.params.BasinOutletFpn_s2)
        
        self.outletLayer = dspt.GetLayer()
        
        spatialRef = self.outletLayer.GetSpatialRef()

        if self.params.verbose:
            featureCount = self.outletLayer.GetFeatureCount()
            print ("            Number of outlets: %d" % (featureCount))
           
        # Get the basin layer created above
        dspoly = self.OpenDs(self.params.basin_areasFPN_s4)
        
        self.basinLayer = dspoly.GetLayer()
        
        if spatialRef == None:
            spatialRef = self.basinLayer.GetSpatialRef()

        if self.params.verbose:
            featureCount = self.basinLayer.GetFeatureCount()
            print ("            Number of basins: %d" % (featureCount))
            
        # Get the omitted basin layer created above
        dsomitted = self.OpenDs(self.params.basin_omitted_areasFPN_s4)
        
        self.omittedLayer = dsomitted.GetLayer()
        
        # Create a dict of the outlet points
        ptD = {}
        
        x = -1
            
        for ptfeature in self.outletLayer:
            
            x += 1
            
            ptgeom = loads(ptfeature.GetGeometryRef().ExportToWkb())
            
            mouth_id = ptfeature.GetField("mouth_id")
            
            basin_id = ptfeature.GetField("basin_id")
            
            upstream = ptfeature.GetField("upstream")
            
            xcoord = ptfeature.GetField("xcoord")
            
            ycoord = ptfeature.GetField("ycoord")
            
            x_orig = ptfeature.GetField("x_orig")
            
            y_orig = ptfeature.GetField("y_orig")
            
            ptD[x] = {'geom':ptgeom, 'upstream':upstream, 'mouth_id':mouth_id, 'basin_id':basin_id,
                      'x':xcoord, 'y':ycoord, 'x_orig':x_orig, 'y_orig':y_orig}
            
        #reset
        self.outletLayer.ResetReading()  # reset the read position to the start
        
        insideL = []
        duplicateL = []
        duplFeatD = {} 
        
        m = -1
        
        for polyfeature in self.basinLayer:
             
            m += 1
             
            polygeom = loads(polyfeature.GetGeometryRef().ExportToWkb())
            
            x = -1
            
            for ptfeature in self.outletLayer:
                
                x += 1
                
                ptgeom = loads(ptfeature.GetGeometryRef().ExportToWkb())
                
                if ptgeom.within(polygeom):
                    
                    if x in insideL:
                        
                        duplicateL.append(x)
                        
                    else:
                        
                        insideL.append(x)
            
            self.outletLayer.ResetReading()  # reset the read position to the start
                
        for pt in insideL:
            if pt in duplicateL:
                
                duplFeatD[pt] = ptD.pop(pt)
                
            else:
                
                ptD.pop(pt)

        # Loop over the omitted polygons and check if the remaining points belong to these
        
        for omitfeature in self.basinLayer:
             
             
            polygeom = loads(polyfeature.GetGeometryRef().ExportToWkb())
            
            x = -1
            
            for ptfeature in self.outletLayer:
                
                x += 1
                
                ptgeom = loads(ptfeature.GetGeometryRef().ExportToWkb())
                
                if ptgeom.within(polygeom):
                    
                    if x in insideL:
                        
                        duplicateL.append(x)
                        
                    else:
                        
                        insideL.append(x)
        
        if self.params.verbose:
            
            print ("            Number of duplicates found: %d" % (len(duplFeatD))) 
            printstr = "            Saved as: %s" % (self.params.basin_duplicate_s2_ptFPN_s4) 
            print (printstr)
            
            print ("            Number of omitted outlets found: %d" % (len(ptD)))
            printstr = "            Saved as: %s" % (self.params.basin_remain_s2_ptFPN_s4) 
            print (printstr) 
            
                  
        self.WriteStage2DS(self.params.basin_remain_s2_ptFPN_s4, spatialRef, ptD, True)
        
        self.WriteStage2DS(self.params.basin_duplicate_s2_ptFPN_s4, spatialRef, duplFeatD, True)