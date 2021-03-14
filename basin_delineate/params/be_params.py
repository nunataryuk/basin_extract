'''
Created on 31 Oct 2020

@author: thomasgumbricht
'''

from pprint import pprint
import params.be_xml as be_xml
import os
import sys
from ds_manage.datasource import DS

class Params(DS):
    ''' Convert xml tree dict co class parameters
    '''
    def __init__(self, xmlFPN):

        DS.__init__(self)

        # Read params from XML
        d = be_xml.ReadXML(xmlFPN)

        # Loop over the dict with params
        for root in d:

            self.system = d[root]['userproj']['@system']

            self.region = d[root]['userproj']['@tractid']
            
            self.alias = d[root]['userproj']['@alias']

            delete = d[root]['process']['delete']

            if delete[0] in ['T','t','Y','y']:
                self.delete = True
            else:
                self.delete = False

            overwrite = d[root]['process']['overwrite']

            if overwrite[0] in ['T','t','Y','y']:
                self.overwrite = True
            else:
                self.overwrite = False
                
            filecheck = d[root]['process']['filecheck']

            if filecheck[0] in ['T','t','Y','y']:
                self.filecheck = True
                self.overwrite = False
                self.delete = False
            else:
                self.filecheck = False

            paramsD = d[root]['process']['parameters']

            self.stage = int(paramsD['@stage'])
            self.invertedFill = int(paramsD['@invertedFill'])
            self.adjacentdist = float(paramsD['@adjacentdist'])
            self.flowaccReclassPath = paramsD['@flowaccReclassPath']
            self.outlet = paramsD['@outlet']
            self.distill = paramsD['@distill']
            self.clusteroutlet = paramsD['@clusteroutlet']
            self.verbose = int(paramsD['@verbose'])
            self.proj4CRS = paramsD['@proj4CRS']
            self.watershed = paramsD['@watershed']
            self.basinCellThreshold = int(paramsD['@basinCellThreshold'])
            self.minBasinAreaKm2 = float(paramsD['@minBasinAreaKm2'])
            self.thresholdTerminalClumps = float(paramsD['@thresholdTerminalClumps'])
            self.filldirsingle = float(paramsD['@filldirsingle'])
            self.pitfillmaxarea = float(paramsD['@pitfillmaxarea'])
            self.tilesize = int(paramsD['@tilesize'])
            self.tilecelloverlap = int(paramsD['@tilecelloverlap'])
            self.grassDEM = paramsD['@grassDEM']

            self.stage = int(d[root]['process']['parameters']["@stage"])

            self.srcpath = '/Volumes/%s' %(d[root]['process']['srcpath']['@volume'])

            self.dstpath = '/Volumes/%s' %(d[root]['process']['dstpath']['@volume'])

            self.srcCompD = d[root]['process']['srccomp']

            self.dstCompD = d[root]['process']['dstcomp']
            

            for key in self.srcCompD:
                self.src = self.srcCompD[key]['@source']

                self.folder = self.srcCompD[key]['@folder']
                
                self.product = self.srcCompD[key]['@product']
                
                self.suffix = self.srcCompD[key]['@suffix']
                
            self.srcfp = os.path.join(self.srcpath, 'GRASSproject', self.src, 'region', self.folder, self.alias, '0')
            

            if self.filecheck:
                
                self.verbose = 2

            if self.verbose:
                infostr = '    Reading xml file: %s' %(xmlFPN)
                print (infostr)
                print ('        system:',self.system)
                print ('        region:',self.region)
                
                print ('        alias:',self.alias)
                print ('        delete flag:', self.delete)
                print ('        overwrite flag:',self.overwrite)
                print ('        paramsD:',paramsD)
                print ('        srcpath:',self.srcpath)
                print ('        dstpath:',self.dstpath)
                print ('        srcCompD:',self.srcCompD)
                print ('        dstCompD:',self.dstCompD)

        if self.verbose > 1:

            pprint(d)

        # Set and control source Datasets
        #self.SourceDS()

        # Arrange processing setup dependent on stage

        #if self.stage == 0: #Stage = is the preparation of the initial GRASS processing, only generates a bare-bone script file
            # Not yet implemented
            #pass

        #elif self.stage == 2: # Stage 1 is the distillation of outlet points and the scripting of the extraction of basins
            # Not yet implemented
            #pass

        # Set the source dataset(s) to use
        #self.SourceDS()

        # Check for source data sets given the defined method
        #self.CheckSrcMethod()

        # Set the dst files, and check if it exists
        self.DestinationDS()
        
        print ('self.dstCompD')
        print (self.dstCompD)
        
        print ('self.srcCompD')
        print (self.srcCompD)
        
        FISK

        '''
        if stage == 1:
            return (stage, overwrite, paramsD, region, dstfpn, dstremfpn, srcfp, SFDsrcfpn, MFDsrcfpn, dsttxtfn, dstPolygonfpn_s1, basinRoot, dsttxtfpn)

        if stage == 2: # Stage 2 removes overlapping basins and optionally removes incomplete basins
            return (stage, overwrite, paramsD, region, dstfpn, dstremfpn, srcfp, SFDsrcfpn, MFDsrcfpn, dsttxtfn, dstPolygonfpn_s1, basinRoot, dsttxtfpn)
        '''

    def SourceDS(self):
        ''' Set and check the paths to source files
        '''

        # The first part of the data source name (prefix) is given

        if 'basin-mouth-outlet-pt' in self.srcCompD:

            self.outletSrcFpn_s2 = self.Mouthsrcfpn = self.SetSrcDSFPN('basin-mouth-outlet-pt')

        if 'basin-SFD-outlet-pt' in self.srcCompD:

            self.outletSrcFpn_s2 = self.SFDsrcfpn = self.SetSrcDSFPN('basin-SFD-outlet-pt')

        if 'basin-MFD-outlet-pt' in self.srcCompD:

            self.outletSrcFpn_s2 = self.MFDsrcfpn = self.SetSrcDSFPN('basin-MFD-outlet-pt')

        if 'shorewall-pt'  in self.srcCompD:

            self.shorewallsrcfpn = self.SetSrcDSFPN('shorewall-pt')

    def SetSrcDSFPN(self, name):
        ''' Construct full path to source data set
        '''

        self.src = self.srcCompD[name]['@source']

        self.folder = self.srcCompD[name]['@folder']

        self.srcfp = os.path.join(self.srcpath, 'GRASSproject', self.src, 'region', self.folder, self.alias, '0')

        self.prefix = self.srcCompD[name]['@prefix']

        self.product = self.srcCompD[name]['@product']

        self.suffix = self.srcCompD[name]['@suffix']

        srcfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':self.prefix, 'prod':self.product,
                    'reg':self.region, 'sf':self.suffix}

        srcfpn = os.path.join(self.srcfp,srcfn)

        if not os.path.exists(srcfpn):

            exitstr = 'Source file %s does not exist' %(srcfpn)

            sys.exit(exitstr)

        return srcfpn

    def CheckSrcMethod(self):

        if self.outlet.upper() == 'MOUTH' and not self.Mouthsrcfpn:
            exitstr = 'to run a full width mouth outlet analysis you must include a Mouth point vector file'
            sys.exit(exitstr)

        if self.outlet.upper() == 'SFD' and not self.SFDsrcfpn:
            exitstr = 'to run an SFD outlet analysis you must include an SFD point vector file'
            sys.exit(exitstr)

        if self.outlet.upper() == 'MFD' and not self.MFDsrcfpn:
            exitstr = 'to run an MFD outlet analysis you must include an MFD point vector file'
            sys.exit(exitstr)

        if self.outlet.upper() == 'MOUTH' and not self.Mouthsrcfpn:
            exitstr = 'to run an SFD outlet analysis you must include an SFD point vector file'
            sys.exit(exitstr)

        if self.outlet.upper() == 'SFD' and not self.SFDsrcfpn:
            exitstr = 'to run a full width mouth outlet analysis you must include a Mouth point vector file'
            sys.exit(exitstr)

        if self.outlet.upper() == 'MFD' and not self.MFDsrcfpn:
            exitstr = 'to run an MFD outlet analysis you must include an MFD point vector file'
            sys.exit(exitstr)

    def DestinationDS(self):
        '''
        '''
        # The destination DS name is by default completely defined from the source DS
        # In case there are no attributes in the tag <basin-outlet> a dummy attribute must be added
        if self.dstCompD['basin-outlet'] == '':
            self.dstCompD['basin-outlet'] = {'dummy':0}

        if '@source' in self.dstCompD['basin-outlet'].keys():
            src = self.dstCompD['basin-outlet']['@source']
        else:
            src = self.src

        if '@product' in self.dstCompD['basin-outlet'].keys():
            product = self.dstCompD['basin-outlet']['@product']
        else:
            product = self.product

        if '@folder' in self.dstCompD['basin-outlet'].keys():
            folder = self.dstCompD['basin-outlet']['@folder']
        else:
            folder = self.folder

        # Note that for all intermediate files, the path is GRASSProject, i.e. the 'srcpath'
        self.dstfp = os.path.join(self.srcpath, 'GRASSproject', src, 'region', folder, self.alias, '0')
        
        '''
        if '@prefix' in self.dstCompD['basin-outlet'].keys():
            prefix = self.dstCompD['basin-outlet']['@prefix']

        else:
            
            #Force prefix
            prefix = 'basin-%s-outlet' %(self.outlet.lower())
            if self.outlet.lower() == 'sfd' and self.distill.lower() == 'mfd':
                prefix = 'basin-sfdxmfd-outlet'
            if self.outlet.lower() == 'mfd' and self.distill.lower() == 'mfd':
                prefix = 'basin-mfd2-outlet'
            if self.outlet.lower() == 'mouth' and self.distill.lower() == 'mouth':
                prefix = 'basin-mouth2-outlet'
                prefix_s4 = 'basin-mouth2-poly'
                prefix_s4_omitted = 'omitted-mouth2-poly'
                prefix_s4_remain_outlets = 'remain-mouth2-outlet'
                prefix_s4_duplicate_outlets = 'duplicate-mouth2-outlet'
        
        
        pfL = prefix.split('-')
        pfpoly = pfL[0]
        for p in range(1,len(pfL)-1):
            pfpoly = '%s_%s' %(pfpoly, pfL[p])
        self.basinRoot = pfpoly

        pfpoly = '%s_%s' %(pfpoly, 'area')
        '''
        if '@suffix' in self.dstCompD['basin-outlet'].keys():
            suffix = self.dstCompD['basin-outlet']['@suffix']
        else:
            suffix = self.suffix
            
        # Set detination file for stage 0
        
        self.Stage0_dstDS(product, suffix)
        
        self.Stage1_dstDS(prefix, product, suffix)

        self.Stage2_dstDS(prefix, product, suffix)
        
        #basin-combined-outlet-pt_drainage_amazonia_0_cgiar-250.shp
        #basin-mouth2-outlet_drainage_amazonia_0_cgiar-250.shp

        dstfn_s3 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s-s3.shp' %{'pf':prefix, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.dstfpn_s3 = os.path.join(self.dstfp,dstfn_s3)
        
        dstfpn_s4_basin_polygons = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s-s4.shp' %{'pf':prefix_s4, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.dstfpn_s4_basin_polygons = os.path.join(self.dstfp,dstfpn_s4_basin_polygons)
        
        dstfn_s4_omitted_polygons = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s-s4.shp' %{'pf':prefix_s4_omitted, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.dstfpn_s4_omitted_polygons = os.path.join(self.dstfp,dstfn_s4_omitted_polygons)
        
        dstfn_s4_remain_outlets = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s-s4.shp' %{'pf':prefix_s4_remain_outlets, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.dstfpn_s4_remain_outlets = os.path.join(self.dstfp,dstfn_s4_remain_outlets)
        
        dstfn_s4_duplicate_outlets = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s-s4.shp' %{'pf':prefix_s4_duplicate_outlets, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.dstfpn_s4_duplicate_outlets = os.path.join(self.dstfp,dstfn_s4_duplicate_outlets)

        self.grassbasefn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s' %{'pf':'replaceme', 'prod':product,
                    'reg':self.region, 'sf':suffix}

        self.grassdstfn = prefix.replace('-','_')
        
        if self.filecheck:
            return
            
        if self.stage == 2 and (self.overwrite or self.delete) and os.path.exists(self.dstPolygonfpn_s2):
            #self.DelDs(self.dstPolygonfpn_s2)
            pass
            
        if self.stage == 4 and (self.overwrite or self.delete) and os.path.exists(self.dstfpn_s4_basin_polygons):
            self.DelDs(self.dstfpn_s4_basin_polygons)
            self.DelDs(self.dstfpn_s4_omitted_polygons)
            self.DelDs(self.dstfpn_s4_remain_outlets)
            self.DelDs(self.dstfpn_s4_duplicate_outlets)

    '''
    FIX THE REST
        dsttxtfn = '%(pf)s-pt_%(prod)s_%(reg)s_0_%(sf)s.txt' %{'pf':prefix, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.dsttxtfpn = os.path.join(self.dstfp,dsttxtfn)

        self.dstremfpn = False
        if self.outlet.lower() == 'sfd' and self.distill.lower() == 'mfd':
            dstremfn = 'rem-%s' %(dstfn)
            self.dstremfpn = os.path.join(self.dstfp,dstremfn)

        if self.outlet.lower() == 'mfd' and self.distill.lower() == 'mfd':
            dstremfn = 'rem-%s' %(dstfn)
            self.dstremfpn = os.path.join(self.dstfp,dstremfn)

        if self.stage == 1 and os.path.exists(self.dstfpn):
            if self.overwrite or self.delete:
                if self.verbose:
                    infostr = '    Deleting existing shape file %s' %(self.dstfpn)
                    print (infostr)
                self.DelDs(self.dstfpn)
                if self.delete:
                    sys.exit()

            elif not os.path.exists(self.dstfp):
                os.makedirs(self.dstfp)
    '''
            
    def Stage0_dstDS(self, product, suffix):
        '''
        '''
        
        terminalClumps_s0 = 'terminal-clumps_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.terminalClumpsFPN_s0 = os.path.join(self.dstfp,terminalClumps_s0)
                
        terminalClumpsLarge_s0 = 'terminal-clumps-large_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.terminalClumpsLargeFPN_s0 = os.path.join(self.dstfp,terminalClumpsLarge_s0)
        
        terminalClumpsSmall_s0 = 'terminal-clumps-small_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.terminalClumpsSmallFPN_s0 = os.path.join(self.dstfp,terminalClumpsSmall_s0)
        
        inlandCompDEM_s0 = 'inland-comp-DEM_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.inlandCompDEMFPN_s0 = os.path.join(self.dstfp,inlandCompDEM_s0)
        
        hydoFillPt_s0 = 'hydro-fill-pt_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydoFillPtFPN_s0 = os.path.join(self.dstfp,hydoFillPt_s0)
        
        invertedFillPt_s0 = 'inverted-fill-pt_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.invertedFillPtFPN_s0 = os.path.join(self.dstfp,invertedFillPt_s0)
        '''
        hydoFillPtDEM_s0 = 'hydro-fill-pt-dem_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydoFillPtDEMFPN_s0 = os.path.join(self.dstfp,hydoFillPtDEM_s0)
        '''
        hydoFillArea_s0 = 'hydro-fill-area_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydoFillAreaFPN_s0 = os.path.join(self.dstfp,hydoFillArea_s0)
        '''
        hydoFillAreaDEM_s0 = 'hydro-fill-area-dem_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydoFillAreaDEMFPN_s0 = os.path.join(self.dstfp,hydoFillAreaDEM_s0)
        '''
        hydoFillDEM_s0 = 'hydro-fill-dem_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydoFillDEMFPN_s0 = os.path.join(self.dstfp,hydoFillDEM_s0)
        
        invertedFillDEM_s0 = 'inverted-fill-dem_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.invertedFillDEMFPN_s0 = os.path.join(self.dstfp,invertedFillDEM_s0)
        
        fillDEMtiles_s0 = 'fill-dem-tiles_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.fillDEMtiles_s0 = os.path.join(self.dstfp,"stage0",fillDEMtiles_s0)
        
    def Stage1_dstDS(self, prefix, product, suffix):
        
        upstreamLnSFD_s1 = 'upstream-ln-SFD_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.upstreamLnSFDfpn_s1 = os.path.join(self.dstfp,'stage1',upstreamLnSFD_s1)
        
        upstreamLnMFD_s1 = 'upstream-ln-MFD_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.upstreamLnMFDfpn_s1 = os.path.join(self.dstfp,'stage1',upstreamLnMFD_s1)
        
        lowlevelOutletClumps_s1 = 'lowlevel-outlet-clumps_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.lowlevelOutletClumpsfpn_s1 = os.path.join(self.dstfp,'stage1',lowlevelOutletClumps_s1)
        
        shorelineOutletClumps_s1 = 'shoreline-outlet-clumps_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.shorelineOutletClumpsfpn_s1 = os.path.join(self.dstfp,'stage1',shorelineOutletClumps_s1)
        
        allOutlets_s1 = 'all-outlets-pt_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.allOutletsfpn_s1 = os.path.join(self.dstfp,'stage1',allOutlets_s1)
        
        SFDOutlets_s1 = 'SFD-outlets-pt_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.SFDOutletsfpn_s1 = os.path.join(self.dstfp,'stage1',SFDOutlets_s1)
        
        MFDOutlets_s1 = 'MFD-outlets-pt_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.MFDOutletsfpn_s1 = os.path.join(self.dstfp,'stage1',MFDOutlets_s1)
        
        thickwall_s1 = 'thickwall_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.thickwallfpn_s1 = os.path.join(self.dstfp,'stage1',thickwall_s1)
        
        mouthshoreline_s1 = 'mouthshoreline_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.mouthshorelinefpn_s1 = os.path.join(self.dstfp,'stage1',mouthshoreline_s1)
        
        shorewall_s1 = 'shorewall_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.shorewallfpn_s1 = os.path.join(self.dstfp,'stage1',shorewall_s1)
        
        shorewallpt_s1 = 'shorewall-pt_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.shorewallptfpn_s1 = os.path.join(self.dstfp,'stage1',shorewallpt_s1)
        
        shorefillDEM_s1 = 'shorefill-DEM_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.shorefillDEMfpn_s1 = os.path.join(self.dstfp,'stage1',shorefillDEM_s1)
        
    def Stage2_dstDS(self, prefix, product, suffix): 
    
        dstfn_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.dstfpn_s2 = os.path.join(self.dstfp,dstfn_s2)
        
        dstPolygonfn_s2 = '%(pf)s-area_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.dstPolygonfpn_s2 = os.path.join(self.dstfp,dstPolygonfn_s2)
 
        print ('self.dstfpn_s2',self.dstfpn_s2)
        print ('self.dstPolygonfpn_s2',self.dstPolygonfpn_s2)

