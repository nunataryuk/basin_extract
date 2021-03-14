'''
Created on 1 Dec 2020

@author: thomasgumbricht
'''


from os import path

#import shapely
from shapely.geometry.polygon import LinearRing
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.polygon import Polygon
from shapely.ops import unary_union
from shapely.validation import explain_validity
from shapely.wkb import loads
from osgeo import ogr, osr

def PolyInspect(srcDS):
    '''
    '''
    if not path.exists(srcDS):
        exitstr = 'The polygonfile with all preliminary basins does not exist'
        exit(exitstr)

    print ('Opening input %s' %(srcDS))
    
    driverName = "ESRI Shapefile" # e.g.: GeoJSON, ESRI Shapefile
    driver = ogr.GetDriverByName(driverName)


    ds = driver.Open(srcDS, 0) # 0 = read only
    
    # Get the source layer
    srcLayer = ds.GetLayer()
    
    # Loop over the features in srcLayer
    n = 0
    for feature in srcLayer:
        n+=1
        print ('feature',n, feature.GetField("area_km2") )
        
        # Get the geometry
        
        geomA = loads(feature.GetGeometryRef().ExportToWkb())
        
        if not geomA.is_valid:
            print ('geomA is not valid')
            if geomA.boundary.is_valid:
                print ('geomA boundary is  valid ')
                geomA = Polygon(geomA.boundary)
                if not geomA.is_valid:
                    print ('geomA is still not valid')
                    print (explain_validity(geomA))
                    


if __name__ == "__main__":
    
        
    srcDS = '/Volumes/GRASS2020/GRASSproject/SRTM/region/basin/amazonia_drain_filldem1cell/0/stage2/basin_mouth2_000125c.shp'
    PolyInspect(srcDS)