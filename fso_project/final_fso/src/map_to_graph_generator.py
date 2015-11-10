
from xml.etree.ElementTree import iterparse
import geopy.distance
from shapely import geometry as shgm
import numpy as np
import logging
#from my_util import MyGridBin
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
               buildingOutputFilePath ,
               graphOutputFilePath, 
               short_edge_length, 
               long_edge_length, 
               edge_append_mode = False,
               small_map_debug = False):
    
    self.mapFilePath = mapFilePath
    self.buildingOutputFilePath = buildingOutputFilePath 
    self.graphOutputFilePath = graphOutputFilePath
    
    self.long_edge_length = long_edge_length
    self.short_edge_length = short_edge_length
    self.edge_append_mode = edge_append_mode
    self.small_map_debug = small_map_debug

    self.isEdgeCallCounter = 0
    self.isIntersectingCallCounter = 0
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
    self.minimum_proximity = 5 #meter
    
    self.building_bin_x_interval =  50 #meter
    self.building_bin_y_interval =  50 #meter
    self.max_bin_x  = -1
    self.max_bin_y  = -1
    self.builidng_bin_pair_intersections =  None
    self.edge_counter = 0
    
    self.logger = logging.getLogger('logger')
    self.logger.addHandler(logging.StreamHandler())
    self.logger.setLevel(logging.INFO)
    
    self.last_time = time.clock()
    self.intersectionCallDuration = 0
    
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
    #self.building_lat = {}
    #self.building_lon = {}
    self.building_lat = []
    self.building_lon = []
    
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
              blat_list = [float(i) for i in blat_list]
              blon_list = [float(i) for i in blon_list]
              self.building_lat.append(list(blat_list))
              self.building_lon.append(list(blon_list))
              
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
              
              blat_list = [float(i) for i in blat_list]
              blon_list = [float(i) for i in blon_list]
              self.building_lat.append(list(blat_list))
              self.building_lon.append(list(blon_list))
              
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
    total_buildings = len(self.building_lat)
    self.building_x = []
    self.building_y = []
    for i in range(total_buildings):
      self.building_x.append([])
      self.building_y.append([])
    
    self.max_x = self.max_y = -float('inf')
    f = open(self.buildingOutputFilePath,"w")
    f.write("BuildingID: lat,lon;....")
    for i in range(total_buildings):#self.building_lat.iteritems():
      f.write("\n"+str(i)+" : ")
      blats = self.building_lat[i]
      blons = self.building_lon[i]
      for blat, blon in zip(blats, blons):
        f.write(str(blat)+" , "+str(blon)+" ; ")
        xy = geopy.Point(blat,blon)
        x_ref = geopy.Point(self.min_lat, blon)
        y_ref = geopy.Point(blat, self.min_lon)
        x = geopy.distance.distance(xy, x_ref).m
        y = geopy.distance.distance(xy, y_ref).m
        self.building_x[i].append(x)
        self.building_y[i].append(y)
        self.max_x = max(self.max_x, x)
        self.max_y = max(self.max_y, y)
    f.close()
    #task 1: save the coordinates in a file, don't delete the following until nodes are saved too
    #del self.building_lat
    #del self.building_lon
  
  
  def loadNodeCoord(self):
    '''
    loads the nodes from the file and hence populate the following lists
    
    '''
   
    #--**************load nodes-***********************
    with open(self.graphOutputFilePath,"r") as f:
      self.node_lat = []
      self.node_lon = []
      next(f)
      file_read_mode = 'load_node'
      for line in f:
        #21 : 40.7149461 , -74.0066863 ; 
        if file_read_mode == 'load_node':
          if line =='Edges\n' or line == 'Edges':
            break
          node_id, coords = line.split(':')
          lat_lon =  coords.split(';')
          lat, lon = lat_lon[0].split(',')
          lat = float(lat)
          lon = float(lon)
          self.node_lat.append(lat)
          self.node_lon.append(lon)
          self.min_lat = min(self.min_lat, lat)
          self.min_lon = min(self.min_lon, lon)
    #now calculate the cartesian node-coords
    self.node_x = []
    self.node_y = []
    for nlat,nlon in zip(self.node_lat,self.node_lon):
      xy = geopy.Point(nlat,nlon)
      x_ref = geopy.Point(self.min_lat, nlon)
      y_ref = geopy.Point(nlat, self.min_lon)
      x = geopy.distance.distance(xy, x_ref).m
      y = geopy.distance.distance(xy, y_ref).m
      self.node_x.append(x)
      self.node_y.append(y)
      self.max_x = max(self.max_x, x)
      self.max_y = max(self.max_y, y)
      
  def binBuildingAndNodes(self):
    '''
    bin the building for easier calculation of LOS between nodes
    bin the nodes in the node_bin
    '''
    self.logger.info("Binning the buildings...")
#     self.building_bin = MyGridBin(x_interval = self.building_bin_x_interval,
#                                   y_interval = self.building_bin_y_interval,
#                                   max_x = self.max_x,
#                                   max_y = self.max_y)
#     
#     self.node_bin = MyGridBin(x_interval = self.building_bin_x_interval,
#                                   y_interval = self.building_bin_y_interval,
#                                   max_x = self.max_x,
#                                   max_y = self.max_y)
    
    self.max_bin_x  = int(self.max_x // self.building_bin_x_interval)
    self.max_bin_y  = int(self.max_y // self.building_bin_y_interval)
    
    self.building_bin = []
    for i in range(self.max_bin_x+1):
      self.building_bin.append([])
      for j in range(self.max_bin_y+1):
        self.building_bin[i].append([])
        
    self.node_bin = []
    for i in range(self.max_bin_x+1):
      self.node_bin.append([])
      for j in range(self.max_bin_y+1):
        self.node_bin[i].append([])
    
    
    node_count = 0
    if not self.edge_append_mode:
      self.node_x = []
      self.node_y = []
      f = open(self.graphOutputFilePath,"w")
      f.write("Node Positions")
      
    else:
      self.loadNodeCoord()
      for x,y in zip(self.node_x, self.node_y):
        #self.node_bin.put(node_count,x,y)
        bin_x = int(x // self.building_bin_x_interval)
        bin_y = int(y // self.building_bin_y_interval)
        self.node_bin[bin_x][bin_y].append(node_count)
        node_count += 1
      
    for bid,bxs in enumerate(self.building_x):
      #bxs = list(self.building_x[bid])
      bys = list(self.building_y[bid])
      last_x = -100
      last_y = -100 
      corner_point_index = -1
      for x,y in zip(bxs,bys):
        #self.building_bin.put(bid, x, y)
        bin_x = int(x // self.building_bin_x_interval)
        bin_y = int(y // self.building_bin_y_interval)
        if bid not in self.building_bin[bin_x][bin_y]:
          self.building_bin[bin_x][bin_y].append(bid)
        corner_point_index += 1
        if not self.edge_append_mode:
          dist_sq = (last_x - x)**2 + (last_y - y)**2
          if dist_sq < self.minimum_proximity * self.minimum_proximity:
            continue
          
          self.node_x.append(x)
          self.node_y.append(y)
          #self.node_bin.put(node_count,x,y)
          if node_count not in self.node_bin[bin_x][bin_y]:
            self.node_bin[bin_x][bin_y].append(node_count)
          node_lat = self.building_lat[bid][corner_point_index]
          node_lon = self.building_lon[bid][corner_point_index]
          f.write("\n"+str(node_count)+" : "+str(node_lat)+" , "+str(node_lon)+" ; ")
          
          node_count += 1
          last_x = x
          last_y = y
        #task 2: save the node latitude longitude in a file according to the index in node_x
    for bid,bxs in enumerate(self.building_x):
      bys = list(self.building_y[bid])
      #p1 = shgm.Polygon(zip(bxs,bys))
      bx_min =  min(bxs)
      by_min =  min(bys)
      bx_max =  max(bxs) 
      by_max  = max(bys)
      x_bin_min = int(bx_min // self.building_bin_x_interval) 
      y_bin_min = int(by_min // self.building_bin_y_interval)
      x_bin_max = int(bx_max // self.building_bin_x_interval)
      y_bin_max = int(by_max // self.building_bin_y_interval)
      for i in range(x_bin_min,x_bin_max+1):
        for j in range(y_bin_min, y_bin_max+1):
          if bid not in self.building_bin[i][j]:
            self.building_bin[i][j].append(bid)
      
    if not self.edge_append_mode:
      f.close()
    del self.building_lat
    del self.building_lon
        
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

  def isOverlapping(self, bbox_1, bbox_2):
    '''
    check if two rectangles are overlapping:
    bbox formats must be [xmin, ymin, xmax, ymax]
    '''
    l1_x = bbox_1[0]
    l2_x = bbox_2[0]
    r1_x = bbox_1[2]
    r2_x = bbox_2[2]
    
    l1_y = bbox_1[3]
    l2_y = bbox_2[3]
    r1_y = bbox_1[1]
    r2_y = bbox_2[1]
    
    
    if l1_x > r2_x or l2_x > r1_x:
      return False
    if l1_y < r2_y or l2_y < r1_y:
      return False
    return True
     

  def isIntersecting(self, x1,y1,x2,y2,xs,ys): 
  
    line = shgm.LineString([(x1, y1), (x2, y2)])
    polygon = shgm.Polygon(zip(xs,ys))
    
#     bbox_1 = line.bounds
#     bbox_2 = polygon.bounds
#     l1_x = bbox_1[0]
#     l2_x = bbox_2[0]
#     r1_x = bbox_1[2]
#     r2_x = bbox_2[2]
#     
#     l1_y = bbox_1[3]
#     l2_y = bbox_2[3]
#     r1_y = bbox_1[1]
#     r2_y = bbox_2[1]
#     
# #     hasOverlap = self.isOverlapping(line.bounds, polygon.bounds)
# #      
# #      
# #     if not hasOverlap:
# #       return False
#     
#     if l1_x > r2_x or l2_x > r1_x:
#       return False
#     if l1_y < r2_y or l2_y < r1_y:
#       return False
      
    if line.crosses(polygon): # and not line.touches(polygon):
      return True
    else: 
      return False
    
  
#   def isEdge(self,x1,y1,x2,y2,building_list):
#     '''
#     given two cooridinates (x1,y1) and (x2,y2), checks wheter the line between these two points
#     intersect with any of the buildings in building_list
#     '''
#     #-----******************************-chage it back
#     self.isEdgeCallCounter += 1
#     line = shgm.LineString([(x1, y1), (x2, y2)])
#     #return False
#     #-----******************************-----------------
#     for bid in building_list:
#       bxs = self.building_x[bid]
#       bys = self.building_y[bid]
#       self.isIntersectingCallCounter += 1
#       if self.isIntersectingCallCounter%10000000==0:
#         #self.logger.setLevel(logging.DEBUG)
#         self.logger.info("Total intersection call made so far:"+str(1.0*self.isIntersectingCallCounter/1000000)+" million")
#         #self.logger.setLevel(logging.INFO)
#       l_time =  time.clock()
#       polygon = shgm.Polygon(zip(bxs,bys))
#       intersection_status = line.crosses(polygon)
#       e_time = time.clock() - l_time
#        
#       self.intersectionCallDuration += e_time
#       if self.isIntersectingCallCounter%100000==0:
#         avg_intersection_duration =  1.0*self.intersectionCallDuration/self.isIntersectingCallCounter
#         self.logger.info("avg intersection call duration: "+str(avg_intersection_duration))
#       if intersection_status:
#         return False
#     return True #<--
#     #return False
  

    
  def isEdge_V2(self,p0_x, p0_y, p1_x, p1_y,building_list):
    '''
    ...
    '''
    s1_x = p1_x - p0_x     
    s1_y = p1_y - p0_y
    #hasIntersected = False
    
    for bid in building_list:
      self.isIntersectingCallCounter += 1
      s_time =  time.clock()
      
#       if hasIntersected:
#         self.intersectionCallDuration = time.clock() - s_time
#         if self.isIntersectingCallCounter%10000==0:
#           avg_intersection_duration =  1.0*self.intersectionCallDuration/self.isIntersectingCallCounter
#           self.logger.info("avg intersection call duration: "+str(avg_intersection_duration))
#        return False
      bxs = list(self.building_x[bid])+[self.building_x[bid][0]]
      bys = list(self.building_y[bid])+[self.building_y[bid][0]]
      total_corners = len(bxs)-1
      for i in range(total_corners):
#         if hasIntersected:
#           self.intersectionCallDuration = time.clock() - s_time
#           if self.isIntersectingCallCounter%10000==0:
#             avg_intersection_duration =  1.0*self.intersectionCallDuration/self.isIntersectingCallCounter
#             self.logger.info("avg intersection call duration: "+str(avg_intersection_duration))
#           return False
        p2_x = bxs[i]
        p2_y = bys[i]
        p3_x = bxs[i+1]
        p3_y = bys[i+1]
        
        s2_x = p3_x - p2_x     
        s2_y = p3_y - p2_y
        l = -s2_x * s1_y + s1_x * s2_y

        if l==0:
          continue
        s = 1.0*(-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / l
        t = 1.0*( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / l
        
        if s >= 0 and s <= 1 and t >= 0 and t <= 1:
          i_x = p0_x + 1.0*t * s1_x
          i_y = p0_y + 1.0*t * s1_y
          if (i_x==p0_x and i_y==p0_y) or (i_x==p1_x and i_y==p1_y):
            continue
          else:
  #           hasIntersected = True
            self.intersectionCallDuration += time.clock() - s_time
            if self.isIntersectingCallCounter%1000000==0:
              avg_intersection_duration =  1.0*self.intersectionCallDuration/self.isIntersectingCallCounter
              self.logger.info("avg intersection call duration: "+str(avg_intersection_duration))
            #self.logger.info("intersection detected!!! "+"s:"+str(s)+" t:"+str(t))
            return False
          #break #crossing detected
    self.intersectionCallDuration += time.clock() - s_time
    if self.isIntersectingCallCounter%1000000==0:
      avg_intersection_duration =  1.0*self.intersectionCallDuration/self.isIntersectingCallCounter
      self.logger.info("avg intersection call duration: "+str(avg_intersection_duration))
    return True
#     char get_line_intersection(float p0_x, float p0_y, float p1_x, float p1_y, 
#     float p2_x, float p2_y, float p3_x, float p3_y, float *i_x, float *i_y)
# {
#     float s1_x, s1_y, s2_x, s2_y;
#     s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
#     s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;
# 
#     float s, t;
#     s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
#     t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);
# 
#     if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
#     {
#         // Collision detected
#         if (i_x != NULL)
#             *i_x = p0_x + (t * s1_x);
#         if (i_y != NULL)
#             *i_y = p0_y + (t * s1_y);
#         return 1;
#     }
# 
#     return 0; // No collision
# }
#     
  
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
    self.logger.info("@calculateLOS....")
    self.edge_counter = 0
    stat_node_pairs = 0
    f = open(self.graphOutputFilePath,"a")
    
    if not self.edge_append_mode:
      f.write("\nEdges")
    
    
    #max_node_gx, max_node_gy = self.node_bin.getMaxGridCoords() #task
    total_node_bin_count = (self.max_bin_x+1)* (self.max_bin_y+1)
    processed_node_bin_count = 0
    gx_start = gy_start = 0
    gx_end =  self.max_bin_x+1
    gy_end =  self.max_bin_y+1
    #if edge_append_mode, then retrieve the last node under processing
    if self.edge_append_mode:
      with open(self.mapFilePath+'.grid_log.txt', 'r') as f1:
        for line in f1:
          if line =='\n':
            break
          gx_start, gy_start = line.split()
          gx_start = int(gx_start)
          gy_start = int(gy_start)
      
      
    f2 = open(self.mapFilePath+'.grid_log.txt', 'a')
    
    self.logger.info("Grid-bin size:"+str(self.max_bin_x)+"X"+str(self.max_bin_y))
    for gx in range(gx_start, gx_end): 
      for gy in range(gy_start, gy_end):
        #-----the following for stat purpose---***_------
        self.logger.info("processing node grid "+str(gx)+" "+str(gy))
        
        self.logger.debug("=================================================")
        
        total_building_considered = 0
        total_intersection_call_made = 0
        
        #step 1: find edges between pairs of nodes of this bin
        
        f2.write(str(gx)+"\t"+str(gy)+"\n")
        #f2.flush()
        processed_node_bin_count += 1
        #node_A = self.node_bin.getbyGridCoords(gx, gy)
        node_A = list(self.node_bin[gx][gy])
        self.logger.debug("\t#nodes in this bin:"+str(len(node_A)))
        grid_node_count_A = len(node_A)
        if grid_node_count_A == 0:
          continue
        self.logger.info("Current node bin index: "+str(gx)+" "+str(gy)+ " progress:"
                            +str(processed_node_bin_count)
                            +"/"
                            +str(total_node_bin_count)
                            + " #nodes in this bin: "+str(grid_node_count_A)
                            + " elapsed-time(s): "+str(self.getElapsedTime())
                            +" #edges so far:"+str(self.edge_counter))
                          
        
        #building_list = self.building_bin.getbyGridCoords(gx, gy) #the buildings of this grid
        building_list = list(self.building_bin[gx][gy])
        self.logger.debug("\t#of buildings in this bin:"+str(len(building_list)))
        debug_this_bin_pair_count = 0
        for i in xrange(grid_node_count_A - 1):
          x1 = self.node_x[ node_A[i] ]
          y1 = self.node_y[ node_A[i] ]
          
          for j in xrange(i+1, grid_node_count_A):
            x2 = self.node_x[ node_A[j] ]
            y2 = self.node_y[ node_A[j] ]
            stat_node_pairs +=1
            debug_this_bin_pair_count += 1
            if debug_this_bin_pair_count%1000==0:
              self.logger.debug("\t#of node-pair from this bin processed so far:"+str(debug_this_bin_pair_count))
            if stat_node_pairs%100000 == 0:
              self.logger.info("processed node pairs:"+str(stat_node_pairs)+
                           " total edges so far:"+str(self.edge_counter)+
                           " time:"+str(self.getElapsedTime()))
            total_intersection_call_made += 1 
            if not self.isEdge_V2(x1, y1, x2, y2, building_list):
              continue
            #else check edge type----
            edge_type = self.getEdgeType(x1, y1, x2, y2)
            if  edge_type != 'non-edge':
              f.write("\n"+str(node_A[i])+"\t"+str(node_A[j])+"\t"+edge_type)
              self.edge_counter += 1
              if self.small_map_debug:
                self.small_graph_edge_list.append((node_A[i],node_A[j],edge_type))
                
        #step 2: find edges between pairs of nodes of this bin and bins right and down of it
        #-------------------------------------------------------------------------------------------------------------
        #grid_offset_x = int(self.long_edge_length/self.building_bin_x_interval)
        #grid_offset_y = int(self.long_edge_length/self.building_bin_y_interval)
         
        neighboring_grids_x = []
        neighboring_grids_y = []
        
        grid_offset = int(self.long_edge_length/self.building_bin_x_interval)
        for i in range(gy+1, gy+grid_offset):
          if i> self.max_bin_y  :#max_node_gy:
            break
          neighboring_grids_x.append(gx)
          neighboring_grids_y.append(i)
          
        for i in range(gx+1, gx+grid_offset):
          if i > self.max_bin_x:#max_node_gx:
            break
          for j in range(gy-grid_offset,gy+grid_offset):
            if j<0:
              continue
            if j> self.max_bin_y:
              break
            neighboring_grids_x.append(i)
            neighboring_grids_y.append(j)
        self.logger.debug("\t#of neighboring grids considered:"+str(len(neighboring_grids_x)))
        #self.logger.setLevel(logging.INFO)
        for ngx, ngy in zip(neighboring_grids_x, neighboring_grids_y):
          #check if the distance to the middle is less than long-edge  
          self.logger.debug("\t*******************currently considering grid:"+str(ngx)+" "+str(ngy))
          dx1 = (gx+0.5)*self.building_bin_x_interval
          dy1 = (gy+0.5)*self.building_bin_x_interval
          
          dx2 = (ngx+0.5)*self.building_bin_x_interval
          dy2 = (ngy+0.5)*self.building_bin_x_interval
          
          dist_sq = (dx1 - dx2)**2 + (dy1 - dy2)**2
          if dist_sq>self.long_edge_length*self.long_edge_length:
            continue
          
          
          node_B = list(self.node_bin[ngx][ngy])
          grid_node_count_B = len(node_B)
          self.logger.debug("\t\t#nodes in this neighboring bin:"+str(len(node_B)))
          if grid_node_count_B == 0:
            continue
          intersecting_grids = list(self.gridRayTrace(gx, gy, ngx, ngy))
          self.logger.debug("\t\t#no of grids intersected by grid-pair:"+str(len(intersecting_grids)))
          bids = []
          for bx,by in intersecting_grids:
            #bids.extend(self.building_bin.getbyGridCoords(bx, by))
            bids.extend(self.building_bin[bx][by])
          bids = list( set(bids) ) #remove duplicates
          total_building_considered = len(bids) #<--for stat purpose
          self.logger.debug("\t\tTotal building contained in the LOS-grids:"+str(total_building_considered))
          for i in xrange(grid_node_count_A):
            x1 = self.node_x[ node_A[i] ]
            y1 = self.node_y[ node_A[i] ]
            for j in xrange(grid_node_count_B):
              x2 = self.node_x[ node_B[j] ]
              y2 = self.node_y[ node_B[j] ]
              stat_node_pairs +=1
#               if stat_node_pairs%100000 == 0:
#                 self.logger.info("processed node pairs:"+str(stat_node_pairs)+
#                              " total edges so far:"+str(self.edge_counter)+
#                              " time:"+str(self.getElapsedTime()))
                 
              total_intersection_call_made += 1
              if not self.isEdge_V2(x1, y1, x2, y2, bids):
                continue
              #else check edge type----
              edge_type = self.getEdgeType(x1, y1, x2, y2)
              if  edge_type != 'non-edge':
                f.write("\n"+str(node_A[i])+"\t"+str(node_B[j])+"\t"+edge_type)
                self.edge_counter += 1
                if self.small_map_debug:
                  self.small_graph_edge_list.append((node_A[i],node_B[j],edge_type))
    f.close() 
    f2.close()   
    self.logger.info("Current node bin index: "+str(gx)+" "+str(gy)+ " progress:"
                            +str(processed_node_bin_count)
                            +"/"
                            +str(total_node_bin_count)
                            + " #nodes in this bin: "+str(grid_node_count_A)
                            + " elapsed-time(s): "+str(self.getElapsedTime())
                            +" #edges so far:"+str(self.edge_counter)
                            +" #buildings considered: "+str(total_building_considered)
                            +" #intersection call made: "+str(total_intersection_call_made)
                            +" #total intersection call so far:"+str(self.isIntersectingCallCounter)
                            )                       
  def debugGenerateVisualGraph(self):
    self.logger.info("@debugGeneateVisualGRaph....")
    patches = [] 
    for bid, bxs in enumerate(self.building_x):
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
    self.logger.info("Total Nodes:"+str(len(self.node_x)))
    self.logger.info("Total Edges:"+str(self.edge_counter))
    self.logger.info("Total #isintersecting calls made:"+str(self.isIntersectingCallCounter))
    for i in range(self.max_bin_x):
      for j in range(self.max_bin_y):
        node_count = len(self.node_bin[i][j])
        building_count = len(self.building_bin[i][j])
        if node_count>200:
          self.logger.info("Bin index:"+str(i)+" "+str(j)+" #node:"+str(node_count)+" #building:"+str(building_count))
    
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
    self.binBuildingAndNodes()
    self.debugPrintSummary()
    
    self.calculateLOS()
    self.logger.info("isEdgeCallCounter:"+str(self.isEdgeCallCounter))
    self.logger.info("isIntersectingCounter:"+str(self.isIntersectingCallCounter))
    self.debugPrintSummary()
    if self.small_map_debug:
      self.debugGenerateVisualGraph()
    
    #self.generateBinPairIntersections()

    
    
    
    