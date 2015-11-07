
import geopy.distance
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

class VisualizeGraph:
  
  def __init__(self,
               buildlingFilePath,
               adjFilePath
               ):
    print "initializing VisualizeGraph(..)"
    self.buildlingFilePath = buildlingFilePath
    self.adjFilePath = adjFilePath
    self.building_lat = None
    self.building_lon = None
    self.node_lat = None
    self.node_lon = None
    self.node_x =  None
    self.node_y =  None
    
    self.min_lat = self.min_lon = float('inf')
    self.adj_graph_edge_list = []
    

  def loadData(self):
    '''
    loads the following from the file:
    1. buildig lat-lon
    2. node coordinates and edges
    '''
    #--**************load building coordinates first-***********************
    self.min_lat = self.min_lon = float('inf')
    with open(self.buildlingFilePath,"r") as f:
      #278112888 : 40.7131171 , -74.0022105 ; 40.713108 , -74.0022527 ; 40.7133529 , -74.0020722 ; 40.713338 , -74.0020477 ;
      next(f) #skip the first  line
      self.building_lat = {}
      self.building_lon = {}
      for line in f:
        bid, lat_lon  = line.split(':')
        coords = lat_lon.split(';')
        lat = []
        lon = []
        for coord in coords[:-1]:
          x,y = coord.split(',')
          x = float(x)
          y = float(y)
          lat.append(x)
          lon.append(y)
          self.min_lat = min(self.min_lat, x)
          self.min_lon = min(self.min_lon, y)
        self.building_lat[bid] = list(lat)
        self.building_lon[bid] = list(lon)
    #--**************load nodes and  edges-***********************
    self.adj_graph_edge_list = []
    with open(self.adjFilePath,"r") as f:
      self.node_lat = []
      self.node_lon = []
      next(f)
      file_read_mode = 'load_node'
      for line in f:
        #21 : 40.7149461 , -74.0066863 ; 
        if file_read_mode == 'load_node':
          if line =='Edges\n':
            file_read_mode = 'load_edge'
            continue
          node_id, coords = line.split(':')
          lat_lon =  coords.split(';')
          lat, lon = lat_lon[0].split(',')
          lat = float(lat)
          lon = float(lon)
          self.node_lat.append(lat)
          self.node_lon.append(lon)
          self.min_lat = min(self.min_lat, lat)
          self.min_lon = min(self.min_lon, lon)
          
        elif file_read_mode == 'load_edge':
          n1, n2, edge_type = line.split()
          self.adj_graph_edge_list.append((
                                          int(n1),
                                          int(n2),
                                          str(edge_type[0])
                                          ))
        
   
  def transformToCartesianCoord(self):
    '''
    transform the lat, lon of buildings into Cartesian Coord (x,y) 
    and save it in building_x, building_y dicts indexed by building id as set in osm file
    '''
    self.building_x = {}
    self.building_y = {}
    

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
        
    self.node_x = []
    self.node_y = []
    for node_lat, node_lon in zip(self.node_lat, self.node_lon):
        xy = geopy.Point(node_lat, node_lon)
        x_ref = geopy.Point(self.min_lat, node_lon)
        y_ref = geopy.Point(node_lat, self.min_lon)
        x = geopy.distance.distance(xy, x_ref).m
        y = geopy.distance.distance(xy, y_ref).m    
        self.node_x.append(x)
        self.node_y.append(y)  


          
  def plotFigure(self):
    #self.logger.info("@debugGeneateVisualGRaph....")
      
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
    

    for u,v,edge_type in self.adj_graph_edge_list:
        edge_color = 'g'
        if edge_type == 's':
          edge_color = 'r'
        plt.plot([ self.node_x[u], self.node_x[v] ],\
                 [ self.node_y[u], self.node_y[v] ], color = edge_color ,  ls ='dotted')
    
    ax.autoscale(enable=True, axis = 'both', tight= True)
    ax.set_aspect('equal', 'box')
    plt.show()
    return 
   
  def draw(self):
    '''
    i) load the building from files
    ii) load the nodes  from file
    iii) also load the edges into networkx graphs
    iv) calculate the cartesian coords
    v) draw the graph and save it in file
    '''
    self.loadData()
    self.transformToCartesianCoord()
    self.plotFigure()
    
    
    
    