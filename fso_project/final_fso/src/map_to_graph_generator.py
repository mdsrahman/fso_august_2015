
from xml.etree.ElementTree import iterparse
import geopy.distance
from shapely import geometry as shgm
import numpy as np
import logging
from my_util import MyGridBin
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


class MapToGraphGenerator:
  '''
  takes a file path to the .osm file as input and 
  saves the generated adjacency graph in the output file path
  major steps:
  1. laod and parse the osm file
  2. bin the buildings
  3.  calculate the LOS and generate the graph
  4. save the graph in the output file
  fields:
  --------------
  building_lat: dictionary of latitudes of all buildings, 
                keyed by building id as found in the osm file
  self.building_lon: dictionary of longitudes of all buildings, 
                keyed by building id as found in the osm file
  max_lon, max_lat: maximum longitude and latitude, used to map to Euclidian coordinate
  min_lon, self.min_lat: minimum longitude and latitude, used to map to Euclidian coordinate
  node_x, node_y: list of node x,y coordinates, list index being the node index
  buidling_bin, node_bin: bins for the building and nodes respectively,
                         a node belongs to only one grid bin, but a building (bid) may belong to multiple bins
  '''
  def __init__(self, mapFilePath, 
               outputFilePath, 
               short_edge_length, 
               long_edge_length, 
               small_map_debug = False):
    
    self.mapFilePath = mapFilePath
    self.outputFilePath = outputFilePath
    self.long_edge_length = long_edge_length
    self.short_edge_length = short_edge_length
    self.small_map_debug = small_map_debug
    self.small_graph_edge_list = []
    
    self.building_lat = None 
    self.building_lon = None 
    self.min_lat = None 
    self.min_lon = None 
    self.max_lat = None 
    self.max_lon =  None

    self.building_bin = None
    self.node_bin = None
    self.node_x = None
    self.node_y = None
    
    self.building_bin_x_interval =  100 #meter
    self.building_bin_y_interval = 100 #meter
    self.builidng_bin_pair_intersections =  None
    self.edge_counter = 0
    
    self.logger = logging.getLogger('logger')
    self.logger.addHandler(logging.StreamHandler())
    self.logger.setLevel(logging.DEBUG)
    
    self.last_time = time.clock()
  
  def getElapsedTime(self):
    '''
    return the time difference since the last of call of this method
    used for performance evaluation of method calls
    '''
    e_time =  time.clock() - self.last_time  
    self.last_time =  time.clock()
    return e_time 
  
  def parseMapFile(self):
    '''
    parse the osm file and save the corner points (lat,lon) of buildings, also min and max lat+lon
    '''
    ways = {}
    tags = {}
    refs = []
    way_id = None
    lat= {}
    lon = {}
    self.building_lat = {}
    self.building_lon = {}
    self.max_lat = self.max_lon = -float('inf')
    self.min_lat = self.min_lon = float('inf')
    building_counter = 0
    
     
    context = iterparse(self.mapFilePath, events=("start", "end"))
  
    context = iter(context)
  
    # get the root element
    event, root = context.next()
    
    for event, elem in context:
        if event == 'start': continue
        if elem.tag == 'tag':
            tags[elem.attrib['k']] = elem.attrib['v']
        #-----------------------------------------------------#
        #              node processing
        #-----------------------------------------------------#
        elif elem.tag == 'node':
            osmid = int(elem.attrib['id'])
            lat[osmid] = float(elem.attrib['lat'])
            lon[osmid] = float(elem.attrib['lon'])
            tags = {}
        #-----------------------------------------------------#
        #              node ref i.e nd processing
        #-----------------------------------------------------#          
        elif elem.tag == 'nd':
            refs.append(int(elem.attrib['ref']))
        #-----------------------------------------------------#
        #              way_id  processing
        #-----------------------------------------------------# 
        elif elem.tag == 'way_id':
          if elem.attrib['role'] == 'outer':
            way_id = int(elem.attrib['ref'])
            #members.append((int(elem.attrib['ref']), elem.attrib['type'], elem.attrib['role']))
        #-----------------------------------------------------#
        #              way processing
        #-----------------------------------------------------# 
        elif elem.tag == 'way':
            osm_id = int(elem.attrib['id'])
            ways[osm_id] = refs
            if 'building' in tags.keys():
              blat_list = [lat[nid] for nid in refs]
              del blat_list[-1]
              blon_list = [lon[nid] for nid in refs]
              del blon_list[-1]
              self.building_lat[osm_id] = blat_list
              self.building_lon[osm_id] = blon_list
              self.max_lat = max(blat_list+[self.max_lat])
              self.max_lon = max(blon_list+[self.max_lon])
              self.min_lat = min(blat_list+[self.min_lat])
              self.min_lon = min(blon_list+[self.min_lon])
              building_counter +=1
              #print "DEBUG:building# ",building_counter
              #ways.append((osm_id, tags, refs)) 
            refs = []
            tags = {}
        #-----------------------------------------------------#
        #              relation processing
        #-----------------------------------------------------# 
        elif elem.tag == 'relation':
            osm_id = int(elem.attrib['id'])
            if 'building' in tags.keys() and way_id:
              #<-----process the ways right here
              blat_list = [lat[nid] for nid in ways[way_id]]
              del blat_list[-1]
              blon_list = [lon[nid] for nid in ways[way_id]]
              del blon_list[-1]
              self.building_lat[osm_id] = blat_list
              self.building_lon[osm_id] = blon_list
              self.max_lat = max(blat_list+[self.max_lat])
              self.max_lon = max(blon_list+[self.max_lon])
              self.min_lat = min(blat_list+[self.min_lat])
              self.min_lon = min(blon_list+[self.min_lon])
              building_counter +=1
              #print "DEBUG:building# ",building_counter
              #relations.append((osm_id, tags, members))
            way_id = None
            tags = {}
        root.clear()
  
  def transformToCartesianCoord(self):
    '''
    transform the lat, lon of buildings into Cartesian Coord (x,y) 
    and save it in building_x, building_y dicts indexed by building id as set in osm file
    '''
    self.building_x = {}
    self.building_y = {}
    self.max_x = self.max_y = -float('inf')
    
    for bindx, blats in self.building_lat.iteritems():
      self.building_x[bindx] = []
      self.building_y[bindx] = []
      blons = self.building_lon[bindx]
      for blat, blon in zip(blats, blons):
        xy = geopy.Point(blat,blon)
        x_ref = geopy.Point(self.min_lat, blon)
        y_ref = geopy.Point(blat, self.min_lon)
        x = geopy.distance.distance(xy, x_ref).m
        y = geopy.distance.distance(xy, y_ref).m
        self.building_x[bindx].append(x)
        self.building_y[bindx].append(y)
        self.max_x = max(self.max_x, x)
        self.max_y = max(self.max_y, y)
    del self.building_lat
    del self.building_lon
  
      
  def binBuilding(self):
    '''
    bin the building for easier calculation of LOS between nodes
    bin the nodes in the node_bin
    '''
    self.logger.info("Binning the buildings...")
    self.building_bin = MyGridBin(x_interval = self.building_bin_x_interval,
                                  y_interval = self.building_bin_y_interval,
                                  max_x = self.max_x,
                                  max_y = self.max_y)
    
    self.node_bin = MyGridBin(x_interval = self.building_bin_x_interval,
                                  y_interval = self.building_bin_y_interval,
                                  max_x = self.max_x,
                                  max_y = self.max_y)
    self.node_x = []
    self.node_y = []
    node_count = 0
    for bid in self.building_x:
      bxs = list(self.building_x[bid])
      bys = list(self.building_y[bid])
      for x,y in zip(bxs,bys):
        self.building_bin.put(bid, x, y)
        self.node_bin.put(node_count,x,y)
        node_count += 1
        self.node_x.append(x)
        self.node_y.append(y)
        
  def gridRayTrace(self, x0, y0, x1, y1):
    '''
    returns grids intersected by line (x0,y0,x1,y1)
    '''
    visited_grids = []
    dx = abs(x0 - x1)
    dy = abs(y0 - y1)
    x = x0
    y = y0
    n = 1 + dx + dy
    x_inc = 1 if x1 > x0 else -1
    y_inc = 1 if y1 > y0 else -1
    error = dx - dy
    dx *= 2
    dy *= 2
    
    while n>0:
      visited_grids.append((x,y))
      if error >=0:
        x += x_inc
        error -= dy
      else:
        y += y_inc
        error += dx
      n -= 1 
    return visited_grids

  def isIntersecting(self, x1,y1,x2,y2,xs,ys):
    line = shgm.LineString([(x1, y1), (x2, y2)])
    polygon = shgm.Polygon(zip(xs,ys))
    if line.crosses(polygon): # and not line.touches(polygon):
      return True
    else: 
      return False
    
  
  def isEdge(self,x1,y1,x2,y2,building_list):
    '''
    given two cooridinates (x1,y1) and (x2,y2), checks wheter the line between these two points
    intersect with any of the buildings in building_list
    '''
    for bid in building_list:
      bxs = self.building_x[bid]
      bys = self.building_y[bid]
      if self.isIntersecting(x1, y1, x2, y2, bxs, bys):
        return False
    return True
  
  def getEdgeType(self,x1,y1,x2,y2):
    dist_sqrt = (x1-x2)**2 + (y1-y2)**2
    if dist_sqrt>self.long_edge_length*self.long_edge_length:
      return 'non-edge'
    if dist_sqrt <= self.short_edge_length * self.short_edge_length:
      return 'short' 
    else:
      return 'long'
    
  def calculateLOS(self):
    '''
    for each bin:
    1. find two sets of nodes A and B as follows:
      A: Nodes of the current bins
      B: Nodes of the bins extending upto 1000 meter left and down
    2. run the following steps for
       I) every pair of nodes in set A
       II) every pair of nodes between set A and B  
       ---------------------------------------------------
        i) find the node distance, if its greater than long edge distance, skip
        ii) if its smaller than short edge distance, mark potenial edge-type either long or short
        iii) find the bin index of n1 and n2
        iv) using indices of iii, find all ray traced grids and hence buildings for intersections
        v) if does not intersect any buildings, save the edge with type
    '''
    self.edge_counter = 0
    stat_node_pairs = 0
    f = open(self.outputFilePath,"w")
    #save the node coords first
    f.write("Node Positions")
    for i,x in enumerate(self.node_x):
      y = self.node_y[i]
      f.write("\n"+str(i)+"\t"+str(x)+"\t"+str(y))
    
    f.write("\nEdges")
    
    
    max_node_gx, max_node_gy = self.node_bin.getMaxGridCoords()
    
    for gx in range(max_node_gx+1):
      for gy in range(max_node_gy+1):
        
        #step 1: find edges between pairs of nodes of this bin
        
        node_A = self.node_bin.getbyGridCoords(gx, gy)
        grid_node_count_A = len(node_A)
        
        building_list = self.building_bin.getbyGridCoords(gx, gy) #the buildings of this grid
        
        for i in xrange(grid_node_count_A - 1):
          x1 = self.node_x[ node_A[i] ]
          y1 = self.node_y[ node_A[i] ]
          
          for j in xrange(i+1, grid_node_count_A):
            x2 = self.node_x[ node_A[j] ]
            y2 = self.node_y[ node_A[j] ]
            stat_node_pairs +=1
            if stat_node_pairs%100000 == 0:
              self.logger.info("processed node pairs:"+str(stat_node_pairs)+
                           " total edges so far:"+str(self.edge_counter)+
                           " time:"+str(self.getElapsedTime()))
            
            if not self.isEdge(x1, y1, x2, y2, building_list):
              continue
            #else check edge type----
            edge_type = self.getEdgeType(x1, y1, x2, y2)
            if  edge_type != 'non-edge':
              f.write("\n"+str(node_A[i])+"\t"+str(node_A[j])+"\t"+edge_type)
              self.edge_counter += 1
              if self.small_map_debug:
                self.small_graph_edge_list.append((node_A[i],node_A[j],edge_type))
                
        #step 2: find edges between pairs of nodes of this bin and bins left and down of it
        #-------------------------------------------------------------------------------------------------------------
        grid_offset_x = int(self.long_edge_length/self.building_bin_x_interval)
        grid_offset_y = int(self.long_edge_length/self.building_bin_y_interval)
         
        max_gx = min(max_node_gx, gx+grid_offset_x)
        max_gy = min(max_node_gy, gy+grid_offset_y)
         
        for ngx in xrange(gx, max_gx+1):
          for ngy in xrange(max_gy+1):
            if ngx== gx and ngy<=gy:
              continue
             
            node_B = self.node_bin.getbyGridCoords(ngx, ngy)
            grid_node_count_B = len(node_B)
             
            intersecting_grids = list(self.gridRayTrace(gx, gy, ngx, ngy))
            bids = []
            for bx,by in intersecting_grids:
              bids.extend(self.building_bin.getbyGridCoords(bx, by))
            bids = list( set(bids) ) #remove duplicates
             
            for i in xrange(grid_node_count_A):
              x1 = self.node_x[ node_A[i] ]
              y1 = self.node_y[ node_A[i] ]
              for j in xrange(grid_node_count_B):
                if node_B[j]>=len(self.node_x):
                  print j
                  print node_B
                  print len(self.node_x)
                  print "!!!!!!"
                x2 = self.node_x[ node_B[j] ]
                y2 = self.node_y[ node_B[j] ]
                stat_node_pairs +=1
                if stat_node_pairs%100000 == 0:
                  self.logger.info("processed node pairs:"+str(stat_node_pairs)+
                               " total edges so far:"+str(self.edge_counter)+
                               " time:"+str(self.getElapsedTime()))
                 
                if not self.isEdge(x1, y1, x2, y2, bids):
                  continue
                #else check edge type----
                edge_type = self.getEdgeType(x1, y1, x2, y2)
                if  edge_type != 'non-edge':
                  f.write("\n"+str(node_A[i])+"\t"+str(node_B[j])+"\t"+edge_type)
                  self.edge_counter += 1
                  if self.small_map_debug:
                    self.small_graph_edge_list.append((node_A[i],node_B[j],edge_type))
    f.close()    
  
  def debugGenerateVisualGraph(self):
    patches = [] 
    for bid, bxs in self.building_x.iteritems():
      #print "bid:",i,"------------------------------" 
      bys = self.building_y[bid]
      pcoord = np.asarray(zip(bxs, bys), dtype = float)
      polygon = Polygon(pcoord, fc='grey')
      #polygon.set_facecolor('none')
      patches.append(polygon) 
      
    
    fig, ax = plt.subplots()
    
    p = PatchCollection(patches, match_original=True)
    ax.add_collection(p)     
    
    #if self.small_map_debug:
    for u,v,edge_type in self.small_graph_edge_list:
        edge_color = 'g'
        if edge_type == 'short':
          edge_color = 'r'
        plt.plot([ self.node_x[u], self.node_x[v] ],\
                 [ self.node_y[u], self.node_y[v] ], color = edge_color ,  ls ='dotted')
    
    ax.autoscale(enable=True, axis = 'both', tight= True)
    ax.set_aspect('equal', 'box')
    plt.show()
    return 
  
  def debugPrintSummary(self): 
    self.logger.info("Map File Path:"+str(self.mapFilePath))
    self.logger.info("Total Area:"+str(self.max_x)+" X "+str(self.max_y)+" meters")
    self.logger.info("Total Buildings:"+str(len(self.building_x)))
    self.logger.info("Total Corner Points:"+str(len(self.node_x)))
    self.logger.info("Total Edges:"+str(self.edge_counter))
    
  def generateMapToGraph(self):
    '''
    complete all the steps of generating graph form osm map in this method
    all other methods are called in here
    Major Steps:
    1. parse the xml osm file map and store the building corner points lat, lon
    2. transform the building corner points(lat, lon) to Euclidian coord (x,y)
    3. select nodes, from building corner points (x,y), either all of those or some subset of those
      for example discarding all points within 20m of a chosen points etc.
    4. bin the buildings according to their corner points and interval values
    5. calculate the LOS among selected points and save the graph
    '''
    self.parseMapFile()
    self.transformToCartesianCoord()
    self.binBuilding()
    self.debugPrintSummary()
    self.calculateLOS()
    if self.small_map_debug:
      self.debugGenerateVisualGraph()
    self.debugPrintSummary()
    #self.generateBinPairIntersections()

    
    
    
    