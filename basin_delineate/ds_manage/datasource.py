'''
Created on 31 Oct 2020

@author: thomasgumbricht
'''

from osgeo import ogr, osr, gdal
import os
import sys

class DS():
    '''
    classdocs
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        driverName = "ESRI Shapefile" # e.g.: GeoJSON, ESRI Shapefile
        self.driver = ogr.GetDriverByName(driverName)
            
    def DelDs(self,fpn):
        ''' Delete data source file
        '''

        if os.path.exists(fpn):
            self.driver.DeleteDataSource(fpn) 
            
    def OpenDs(self, fpn): 
        ''' Open existing data source
        '''
        if self.params.verbose:
            print ('        opening data source: %s' % (fpn))
        ds = self.driver.Open(fpn, 0) # 0 = read only
            
        if ds is None:
            exitstr =  ('Could not open %s' % (fpn))
            sys.exit(exitstr)
    
        return ds
        layer = ds.GetLayer()
        

        
        '''
        self.spatialRef = layer.GetSpatialRef()
        

        featureCount = layer.GetFeatureCount()
        print ("Number of features in %s: %d" % (os.path.basename(fpn),featureCount))
        '''
        
        
        return layer
        
    def CopyDs(self, fpn):
        ''' Copy an existing data source
        '''   
                
        self.OpenDs(fpn)
                          
        # Copy all features to finalOutletL
        
        finalOutletL = [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY(), feature.GetField("upstream"), 
                         feature.GetField("xcoord"), feature.GetField("ycoord"), feature.GetField("basin")) for feature in self.layer]
    
        return finalOutletL
      
    def GetDsTransform(self,fpn):
        '''
        '''
        raster = gdal.Open(fpn)
        gt =raster.GetGeoTransform()
        cols = raster.RasterXSize
        rows = raster.RasterYSize
        
        raster = None
        return gt,cols,rows
        '''
        print gt
        (258012.37107330866, 2.11668210080698, 0.0, 163176.6385398821, 0.0, -2.1168501270110074)
        pixelSizeX = gt[1]
        pixelSizeY =-gt[5]
        print pixelSizeX
        2.11668210080698
        print pixelSizeY
        2.1168501270110074 
        ''' 
      
    def CreateDs(self, fpn):
        ''' Create an empty DS
        ''' 
        driverName = "ESRI Shapefile" # e.g.: GeoJSON, ESRI Shapefile

        self.driver = ogr.GetDriverByName(driverName)
            
        self.dstds = self.driver.CreateDataSource(fpn)
        
    def CreateFields(self,mouth):
        '''
        '''
        self.dstlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        self.dstlayer.CreateField(ogr.FieldDefn('upstream', ogr.OFTReal))
        self.dstlayer.CreateField(ogr.FieldDefn('xcoord', ogr.OFTReal))
        self.dstlayer.CreateField(ogr.FieldDefn('ycoord', ogr.OFTReal))
        if mouth:
            self.dstlayer.CreateField(ogr.FieldDefn('basin_id', ogr.OFTInteger))
            self.dstlayer.CreateField(ogr.FieldDefn('mouth_id', ogr.OFTInteger))
            self.dstlayer.CreateField(ogr.FieldDefn('x_orig', ogr.OFTReal))
            self.dstlayer.CreateField(ogr.FieldDefn('y_orig', ogr.OFTReal))
            
    def CreateFeature(self, defn, mouth, i, f, feature):
        '''
        '''
        # Create a new feature (attribute and geometry)
        feat = ogr.Feature(defn)
        feat.SetField('id', i)
        feat.SetField('upstream',feature['upstream'])
        feat.SetField('xcoord',feature['x'])
        feat.SetField('ycoord',feature['y'])
        if mouth:
            feat.SetField('basin_id',feature['basin_id'])
            feat.SetField('mouth_id',feature['mouth_id'])
            feat.SetField('x_orig',feature['x_orig'])
            feat.SetField('y_orig',feature['x_orig'])
            
        return (feat)
    
    def WriteStage0DS(self, fpn, spatialRef, featureD):
        ''' Write data source to file
        '''
        
        self.CreateDs(fpn)
        
        if spatialRef == None:
            spatialRef = osr.SpatialReference()
            spatialRef.ImportFromProj4(self.params.proj4CRS)
        
        # Create the destination dataset layer
        self.dstlayer = self.dstds.CreateLayer('tiles', spatialRef, geom_type=ogr.wkbPolygon)

        self.dstlayer.CreateField(ogr.FieldDefn('id', ogr.OFTString))
        
        
        defn = self.dstlayer.GetLayerDefn()
                
        for f in featureD:
            
            
            feat = ogr.Feature(defn)
            
            feat.SetField('id', featureD[f]['id'])
            
            # Create ring
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(featureD[f]['xmin'], featureD[f]['ymax'])
            ring.AddPoint(featureD[f]['xmax'], featureD[f]['ymax'])
            ring.AddPoint(featureD[f]['xmax'], featureD[f]['ymin'])
            ring.AddPoint(featureD[f]['xmin'], featureD[f]['ymin'])
            ring.AddPoint(featureD[f]['xmin'], featureD[f]['ymax'])

            # Create polygon
            polygeom = ogr.Geometry(ogr.wkbPolygon)
            polygeom.AddGeometry(ring)
            
            #polygeom = ogr.CreateGeometryFromWkb(featureD[f]['newgeom'].wkb)
                   
            feat.SetGeometry(polygeom)

            self.dstlayer.CreateFeature(feat)
            
            feat = None  # destroy feat

    def WriteStage2DS(self, fpn, spatialRef, featureD, mouth=False):
        ''' Write data source to file
        '''
        
        self.CreateDs(fpn)
        
        if spatialRef == None:
            spatialRef = osr.SpatialReference()
            spatialRef.ImportFromProj4(self.params.proj4CRS)
        
        # Create the destination dataset layer
        self.dstlayer = self.dstds.CreateLayer('outlet', spatialRef, geom_type=ogr.wkbPoint)

        
        self.dstlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        self.dstlayer.CreateField(ogr.FieldDefn('upstream', ogr.OFTReal))
        self.dstlayer.CreateField(ogr.FieldDefn('xcoord', ogr.OFTReal))
        self.dstlayer.CreateField(ogr.FieldDefn('ycoord', ogr.OFTReal))
        if mouth:
            self.dstlayer.CreateField(ogr.FieldDefn('basin_id', ogr.OFTInteger))
            self.dstlayer.CreateField(ogr.FieldDefn('mouth_id', ogr.OFTInteger))
            self.dstlayer.CreateField(ogr.FieldDefn('x_orig', ogr.OFTReal))
            self.dstlayer.CreateField(ogr.FieldDefn('y_orig', ogr.OFTReal))
           
        defn = self.dstlayer.GetLayerDefn()
        
        i = -1
        for f in featureD:
            i +=1
            
            #feat = self._CreateFeature(defn, mouth, i, f, featureD[f])
            feat = ogr.Feature(defn)
            feat.SetField('id', i)
            feat.SetField('upstream',featureD[f]['upstream'])
            feat.SetField('xcoord',featureD[f]['x'])
            feat.SetField('ycoord',featureD[f]['y'])
            if mouth:
                feat.SetField('basin_id',featureD[f]['basin_id'])
                feat.SetField('mouth_id',featureD[f]['mouth_id'])
                feat.SetField('x_orig',featureD[f]['x_orig'])
                feat.SetField('y_orig',featureD[f]['x_orig'])
                
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(featureD[f]['x'], featureD[f]['y'])
            
            feat.SetGeometry(point)
            
            self.dstlayer.CreateFeature(feat)
            
            feat = None  # destroy feat
            
        '''
        i = -1
        for f in featureD:
            i +=1
            # Create a new feature (attribute and geometry)
            feat = ogr.Feature(defn)
            feat.SetField('id', i)
            feat.SetField('upstream',featureD[f]['upstream'])
            feat.SetField('xcoord',featureD[f]['x'])
            feat.SetField('ycoord',featureD[f]['y'])
            if mouth:
                feat.SetField('basin_id',featureD[f]['basin_id'])
                feat.SetField('mouth_id',featureD[f]['mouth_id'])
                feat.SetField('x_orig',featureD[f]['x_orig'])
                feat.SetField('y_orig',featureD[f]['x_orig'])
            '''
                
    def WriteStage4DS(self, fpn, spatialRef, featureD, omitted=False, mouth=False):
        ''' Write data source to file
        '''
        
        self.CreateDs(fpn)
        
        if spatialRef == None:
            spatialRef = osr.SpatialReference()
            spatialRef.ImportFromProj4(self.params.proj4CRS)
        
        # Create the destination dataset layer
        self.dstlayer = self.dstds.CreateLayer('basin', spatialRef, geom_type=ogr.wkbPolygon)

        self.dstlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        
        self.dstlayer.CreateField(ogr.FieldDefn('area_km2', ogr.OFTReal))
        
        self.dstlayer.CreateField(ogr.FieldDefn('basin_id', ogr.OFTInteger))
        
        if omitted:
            self.dstlayer.CreateField(ogr.FieldDefn('bndtouch', ogr.OFTInteger))
        
        if mouth:
            
            self.dstlayer.CreateField(ogr.FieldDefn('mouth_id', ogr.OFTInteger))
            
            self.dstlayer.CreateField(ogr.FieldDefn('x_mouth', ogr.OFTReal))
            
            self.dstlayer.CreateField(ogr.FieldDefn('y_mouth', ogr.OFTReal))
           
        defn = self.dstlayer.GetLayerDefn()
        
        i = -1
        
        for f in featureD:
            
            i +=1

            if omitted and not featureD[f]['skip']:
            
                continue
            
            elif not omitted and featureD[f]['skip']:
            
                continue
            
            feat = ogr.Feature(defn)
            
            feat.SetField('id', i)
            
            feat.SetField('area_km2',featureD[f]['area_km2'])
            
            feat.SetField('basin_id',featureD[f]['basin_id'])
        
            if omitted:
                    
                if featureD[f]['bndtouch']:
                    
                    feat.SetField('bndtouch',1)
                    
                else:
                    
                    feat.SetField('bndtouch',0)
            
            if mouth:
                
                feat.SetField('mouth_id',featureD[f]['mouth_id'])
                
                feat.SetField('x_mouth',featureD[f]['x_mouth'])
                
                feat.SetField('y_mouth',featureD[f]['y_mouth'])
                
            if featureD[f]['newgeom']:
                   
                polygeom = ogr.CreateGeometryFromWkb(featureD[f]['newgeom'].wkb)
                   
            else:
                
                srcFeature = self.srcLayer.GetFeature(f)
            
                polygeom = srcFeature.GetGeometryRef()
                
            feat.SetGeometry(polygeom)

            self.dstlayer.CreateFeature(feat)
            
            feat = None  # destroy feat
