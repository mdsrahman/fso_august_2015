# This Module contains the class GenerateInput
# that handles 
# i) parsting the configuration file and reading the params
# ii) generates all input graph whether synthetic or from the map file i.e osm file

import ConfigParser
#import numpy as np
import networkx as nx
#import matplotlib.pyplot as plt
import random 
import logging
import time
 


class InputGenerator:
  '''
  description of class-level public fields:

  **-----------------------------------------------------------------------------**
        the following fields are read from file titled configFile(*class-field)
  **-----------------------------------------------------------------------------**
  configFile(str):stores the name of the input config-file supplied to the constructor
  
  graphType(str):either of the values 'map' or 'synthetic', 
                used to determine whether to read osm file or 
                generate synthetic random graph for input
                
  seed(int):used as seed for the random class of python class
  
  experiment_name(str):used to prefix filenames in all graph output files i.e 
                      if the value is 'san_francisco_north', then 
                      adjacency_graph: san_francisco_north.adj.txt
                      background_graph: san_francisco_north.bg.txt
                      static_graph: san_francisco_north.static.txt
                      dynamic_graph: san_francisco_north.dynamic.txt
                      dynamic_graph_for_java_input: san_francisco_north.dyn_java.txt 
                      
 
  gateway_to_node_ratio(float): >0 and <=1, ratio of nodes to be randomly 
                              selected as gateways from the total nodes
  link_capacity(float): >0, capacity of the fso-links in Mbps
  
  fso_per_node(int): >0, maximum number of fso-links (pairs of up and down-links) per node, 
                      is the max. node degree in the undirected static graph,
                    also, the number of directed fso-link (up or down link) is twice this number
  fso_per_gateway(int): >0, same as fso_per_node above, except applies to gateways only               
  max_x_coord(float): >0, maximum value of x-coordinate in the rectangular map area,
                    if graphType(*class-field) is 'map', it is retrieved from the
                    latitude/longitude of the buildings in the input osm-file,
                    otherwise, this value is read from input configFile(*class-field)
                    
  max_y_coord(float): >0, maximum value of y-coordinate in the rectangular map area,
                    if graphType(*class-field) is 'map', it is retrieved from the
                    latitude/longitude of the buildings in the input osm-file,
                    otherwise, this value is read from input configFile(*class-field)
                    
  target_spacing(float): >0 and <= max_x_coord (*class-field), 
                      the targets to be covered are placed on the corners
                      of square grids, this is the length of the side of 
                      the square in meter
  
  max_short_edge_length(float): >0, maximum length of short fso-links, 
                              only these short links are used in backbone graph building
                              
  max_long_edge_length(float): >max_short_edge_length(*class-field), 
                              maximum length of the long fso-links, 
                              these links are added to both dynamic and static graph 
                              after step 3 of the heuristic placement algorithm
                                                   
  output_statistics_file(str): the full-path to the file that the output statistics are
                      appended to, see output_stat_fields, for full description & order
                      of the these values, important info for graph plots
                      
  graph_output_folder: the folder to store all the graphs 
                    like .adj.txt, .bg.txt, .static.txt, .dynamic.txt, .dyn_java.txt etc. 
  
  coverage_radius(float): coverage radius of target-covering nodes in meter,
                           used in subclass step_1 for target-coverage
  
  ratio_of_max_added_nodes_in_step_4(float): >=0, ratio of the maximum no. of nodes to be added in
                                  step 4 (both dynamic and static) to the existing nodes in 
                                  backbone_graph built after step 3
  percent_of_pattern_nodes_in_avg_flow_calculation(int): <=0 and <=100, as the name implied
                          used in the subclass related to stat-collection
  number_of_pattern_in_avg_flow_calculation(int): <=0 and <=number_of_nodes, as the name
                          implies, this number of iterations are used for avg. flow calculations
                          involving the pattern nodes in the subclass related to stat-collection
  -----the following fields are meaningful for only synthetic graph---------
                  
  no_of_samples(int): >0, number of synthetic random graphs 
                    to be generated for the experiment, 
                    each sample has 
                    the same number of nodes i.e number_of_nodes (*class-field)
                    and the same number of max. edges i.e. max_no_of_edges(*class-field)
                    and the same ratio of short-to-long edges 
                    i.e. short_to_long_edge_ratio (*class-field)
                    
  number_of_nodes(int): >0, total number of nodes in each synthetic graph
  
  max_no_of_edges(int): >0, maximum number of edges in each synthetic graph
  
  short_to_long_edge_ratio(float): >=0 and <=1, ratio between total short and long
                                  edges in each synthetic graph 
  
  **-----------------------------------------------------------------------------**
        the following class-fields get values internally
  **-----------------------------------------------------------------------------**  
  config(ConfigParser obj): used to read configuration values 
                          from the configFile(*class-field)  
                          
  node_x(list of float): stores the x-coordinate of all nodes        
  
  node_y(list of float): stores the y-coordinate of all nodes     
  
  max_short_edge_per_node(int): max. allowed short edges 
                        i.e. edges having length<=max_short_edge_length,
                        calculated in method generateEdges(...) as follows: 
                        expected_no_of_edge: max_no_of_edges/number_of_nodes
                        then applying short_to_long_edge_ratio to expected_no_of_edge accordingly
                        
  max_long_edge_per_node(int): max. allowed long edges 
                        i.e. edges having length l such that,
                        max_short_edge_length< l <=max_long_edge_length,
                        calculated in method generateEdges(...) as follows: 
                        expected_no_of_edge: max_no_of_edges/number_of_nodes
                        then applying short_to_long_edge_ratio to expected_no_of_edge accordingly
  
  short_edge_counter(nparray of type int): initialized only in method generateEdges(...),
                      keeps track of number of short edges added to node i where i is the index
  
  long_edge_counter(nparray of type int): initialized only in method generateEdges(...),
                      keeps track of number of long edges added to node i where i is the index 
  
  adj (networkx.Graph obj): the adjacency matrix generated either synthetically or from the map
                        initialized in generateSyntheticGraph(..) method or 
                        #TODO: document corresponding function for map-calls
  short_edge_adj(networkx.Graph obj): a subgraph of adj (*class-field) with only the 'short'
                        fso links the corresponding nodes, initialized and built up in method 
                        generateEdges(..)
  gateways(list of int): the list of gateway/sink nodes, set at the method selectGateways(..)
  logger(logging refernce, used for debugging in all the sub-classes and itself)
  self.no_of_connected_components(int): total connected components in self.adj 
                                        value set in selectGateways(), used for debugging and stat
  
  last_time (time_t obj): used in getting elapsed_time() method
  ''' 
  
  def __init__(self, configFile):
    #----------------------------------------------------------
    #           Class Fields
    #----------------------------------------------------------
    #--------the following class-fields get initial values externally-----
    
    self.configFile = configFile
    self.graphType =  None
    self.seed = None
    self.experiment_name = None

    self.gateway_to_node_ratio = None
    self.link_capacity = None
    self.fso_per_node = None
    self.fso_per_gateway =  None
    self.max_short_edge_length = None
    self.max_long_edge_length = None
    self.target_spacing = None
    self.coverage_radius = None
    self.output_statistics_file = None
    self.graph_output_folder = None
    
    self.ratio_of_max_added_nodes_in_step_4 = None
      #-----the following are initialized only for synthetic graphs------
    self.no_of_samples = None
    self.number_of_nodes = None
    self.max_short_edge_per_node = None
    self.max_long_edge_per_node = None
    self.max_no_of_edges = None
    self.max_x_coord = None
    self.max_y_coord = None
    
      #-----the following are initialized only for map files------
    #TODO: Declare the param variables for map input
    
    #------the following class-fields get values internally--------
    self.config = None
    self.node_x = None
    self.node_y = None
    #self.max_short_edge_per_node = None
    #self.max_long_edge_per_node = None
    self.short_edge_counter =  None
    self.long_edge_counter = None
    self.adj = None
    self.short_edge_adj = None
    self.gateways = None
    self.no_of_connected_components = None
    self.percent_of_pattern_nodes_in_avg_flow_calculation =  None
    self.number_of_pattern_in_avg_flow_calculation =  None
    self.last_time = time.clock()
    #----------------------------------------------------------
    #           End of Class Fields
    #----------------------------------------------------------
    
    self.logger = logging.getLogger('logger')
    self.logger.addHandler(logging.StreamHandler())
    
    self.parseConfigFile()
    random.seed(self.seed)
  
  def getElapsedTime(self):
    '''
    return the time difference since the last of call of this method
    used for performance evaluation of method calls
    '''
    e_time =  time.clock() - self.last_time  
    self.last_time =  time.clock()
    return e_time
  
  def parseConfigFile(self):
    '''
    populates class fields with values from configFile(*class-field)
    all parameter values under 'global' section are stored in the class-fields,
    if graphType(*class-field) is 'map',
      then all parameter values under 'map' section are stored in the class-fields,
    else if graphType(*class-field) is 'synthetic'
      then all parameter values under 'synthetic' section are stored in the class-fields
    if any expected section or option not found, then uncaught exception are raised.. 
    no default value is assumed for any option
    
    Args:
      None
      
    Returns:
      None
    '''
    self.config = ConfigParser.ConfigParser()
    self.config.read(self.configFile)
    
    self.graphType = self.config.get('global','graphType')
    self.seed = self.config.getint('global','seed')
    self.experiment_name = self.config.get('global','experiment_name')

    self.gateway_to_node_ratio = self.config.getfloat('global','gateway_to_node_ratio')
    self.link_capacity = self.config.getfloat('global','link_capacity')
    self.fso_per_node = self.config.getint('global','fso_per_node')
    self.backup_fso_per_node = self.fso_per_node
    self.fso_per_gateway = self.config.getint('global','fso_per_gateway')
    self.backup_fso_per_gateway = self.fso_per_gateway
    self.target_spacing =  self.config.getfloat('global','target_spacing')
    self.coverage_radius =  self.config.getfloat('global','coverage_radius')
    self.max_short_edge_length =  self.config.getfloat('global','max_short_edge_length')
    self.max_long_edge_length = self.config.getfloat('global','max_long_edge_length')
    self.output_statistics_file =  self.config.get('global','output_statistics_file')
    self.graph_output_folder =  self.config.get('global','graph_output_folder')
    self.ratio_of_max_added_nodes_in_step_4 =  self.config.getfloat('global','ratio_of_max_added_nodes_in_step_4')
    self.percent_of_pattern_nodes_in_avg_flow_calculation = \
        self.config.getint('global','percent_of_pattern_nodes_in_avg_flow_calculation')
    self.number_of_pattern_in_avg_flow_calculation = \
        self.config.getint('global','number_of_pattern_in_avg_flow_calculation')
    if self.graphType=='synthetic':
      self.no_of_samples = self.config.getint('synthetic','no_of_samples')
      self.number_of_nodes = self.config.getint('synthetic','number_of_nodes')
      self.max_short_edge_per_node = self.config.getint('synthetic','max_short_edge_per_node')
      self.max_long_edge_per_node = self.config.getint('synthetic','max_long_edge_per_node')
      self.max_no_of_edges = self.config.getint('synthetic','max_no_of_edges')
      self.max_x_coord = self.config.getfloat('synthetic','max_x_coord')
      self.max_y_coord = self.config.getfloat('synthetic','max_y_coord')
    
    #TODO: populate the values for map input
    #elif self.graphType == 'map':
    #-----------------??
  
  def runInputGenerator(self):
    '''effectively executes all the steps related to tasks of this class
      processing:
         
        i) set the adjacency matrix adj by either of the following:
           if graphtType is 'synthetic', calls generateSyntheticGraph() method
           else if graphType is 'map', calls generateMapGraph() method 
                                       #TODO:complete method generateMapGraph()
        ii) selects the gateways by calling selectGateways() method
    '''
    if self.graphType == 'synthetic':
      self.generateSyntheticGraph()
    self.logger.debug("completed generating sythetic graph:time:"+str(self.getElapsedTime()))
    self.selectGateways()
    self.logger.debug("completed selecting sinks"+str(self.getElapsedTime()))
    
  def generateSyntheticGraph(self):
    '''
    Pre-condition: None
    Post-condition:
      i) Generates the random node position by calling method self.generateNodePositions(..)
      ii) Initializes the self.adj networkx object by creating self.number_of_nodes in it
      iii) Add edges to self.adj object by calling method self.generateEdges(..)
      Therefore, the synthetic graph self.adj is ready with nodes & edges after this method-call
    Args: None
    Returns: None
    '''
    #task (i)
    self.generateNodePositions()
    self.logger.debug("Completed generating node positions")
    #task (ii)
    self.adj = nx.Graph()
    self.adj.graph['name'] = 'Adjacency Graph'
    for i in xrange(self.number_of_nodes):
      self.adj.add_node(i) 
    self.logger.debug("Completed adding all nodes to Adjacency Graph")
    #task (iii)
    self.generateEdges()
    return
  
  def generateNodePositions(self):
    '''
    generates the random positions for all number_of_nodes, 
    stores the position values in list node_x, node_y,
    here, the list index is the node number starting from 0 through (number_of_nodes - 1)
    Args: None
    Returns: None
    '''
    self.node_x = []
    self.node_y = []
    for i in xrange(self.number_of_nodes):
      x = random.uniform(0, self.max_x_coord)
      y = random.uniform(0, self.max_y_coord)
      self.node_x.append(x)
      self.node_y.append(y)
      
    return
  
  def getDistanceSquare(self,n1,n2):
    '''
    Args: 
      n1: node number or index (0<= and <=self.number_of_nodes - 1)
      n2: node number or index (0<= and <=self.number_of_nodes - 1)
    Returns:
      d: the square (for performance reaseon) of Euclidian distance
         between node n1 and node n2,
         after reading the coordinate values from self.node_x and self.node_y
    '''
    x1 = self.node_x[n1]
    y1 = self.node_y[n1]
    x2 = self.node_x[n2]
    y2 = self.node_y[n2]
    d = (x1- x2)**2 + (y1-y2)**2
    return d 
  
  def hasRoomForShortEdge(self,n):
    '''
    Args: n: node number, 0<=n<self.number_of_nodes
    Returns: True if n has total short edges < self.max_short_edge_per_node, False otherwise
    '''
    if self.short_edge_counter[n] < self.max_short_edge_per_node:
      return True
    else: 
      return False
    
  def hasRoomForLongEdge(self,n):
    '''
    Args: n: node number, 0<=n<self.number_of_nodes
    Returns: True if n has total long edges < self.max_long_edge_per_node, False otherwise
    '''
    if self.long_edge_counter[n] < self.max_long_edge_per_node:
      return True
    else: 
      return False
    
  def addEdge(self,n1,n2,edge_type):
    '''
    Args:
      n1: node index,  0<=n1<self.number_of_nodes
      n2: node index,  0<=n2<self.number_of_nodes
      edge_type: either 'long' or 'short'
    Returns: None
    
    pre-condition: 
      i) the adjacency networkx graph self.adj must be instantiated
      ii) the nparrays self.short_edge_counter and self.long_edge_counter must be allocated
          and initialized to zeros
    post-condition: 
      i) adds undirected edge(n1,n2) to networkx graph self.adj
      ii) increases the edge counter (either long or short based on edge_type) of nodes
           n1 and n2
    '''
    self.adj.add_edge(n1,n2,con_type=edge_type)
    if edge_type == 'short':
      self.short_edge_counter[n1] += 1
      self.short_edge_counter[n2] += 1
    elif edge_type == 'long':
      self.long_edge_counter[n1] += 1
      self.long_edge_counter[n2] += 1
    return
  
  def generateEdges(self): 
    '''
      pre-condition: 
              the adjacency networkx graph i.e. self.adj must be initialize with all the nodes
              before any call to this method
      post-condition:
              edges are added to self.adj such that
              i) total number of edges does <= self.max_no_of_edges
              ii) for any node:
                  total_short_edge <= self.max_short_edge_per_node
                  total_long_edge <= self.max_long_edge_per_node
              iii) self.short_edge_adj is also initialized and 'short' edges are added 
                  simultaneously as with adj
      Args: None
      Returns: None
    '''
    #initialize the max short and long edge values
    self.short_edge_adj = nx.Graph()
    self.short_edge_adj.graph['name'] =  'Short Edge Adjacency Graph'
    for n in xrange(self.number_of_nodes):
      self.short_edge_adj.add_node(n)
      
#     expected_edge_per_node =  int(round( 1.0*self.max_no_of_edges/ self.number_of_nodes, 0))
#     self.max_long_edge_per_node = int(round(1.0* expected_edge_per_node / (self.short_to_long_edge_ratio+1),0))
#     self.max_short_edge_per_node = int(round(1.0*self.max_long_edge_per_node*self.short_to_long_edge_ratio,0))
#     
    #initialize the long and short edge counter per node
    #self.short_edge_counter = np.zeros(self.number_of_nodes, dtype = int)
    #self.long_edge_counter = np.zeros(self.number_of_nodes, dtype = int)
    self.short_edge_counter = []
    self.long_edge_counter = []
    
    for i in xrange(self.number_of_nodes):
      self.short_edge_counter.append(0)
      self.long_edge_counter.append(0)
    
    #now consider pairs of nodes and assign edges both long and short between them
    total_edges_added = 0

    for n1 in xrange(self.number_of_nodes - 1):
      if n1%100==0:
        self.logger.debug("n1:elapsed_time:"+str(n1)+":"+str(self.getElapsedTime()))
      if self.short_edge_counter[n1]>= self.max_short_edge_per_node and\
         self.long_edge_counter[n1] >= self.max_long_edge_per_node:
        continue
      
      for n2 in xrange(n1+1, self.number_of_nodes):
        if self.short_edge_counter[n2]>= self.max_short_edge_per_node and\
           self.long_edge_counter[n2] >= self.max_long_edge_per_node:
          continue
        
        distance_sq =  self.getDistanceSquare(n1, n2) 
        if distance_sq  <= self.max_short_edge_length**2:
          if self.short_edge_counter[n1] < self.max_short_edge_per_node and\
             self.short_edge_counter[n2] < self.max_short_edge_per_node:
            self.addEdge(n1, n2, 'short')
            self.short_edge_adj.add_edge(n1,n2)
            total_edges_added += 1
        elif distance_sq  <= self.max_long_edge_length**2:
          if self.long_edge_counter[n1] < self.max_long_edge_per_node and\
              self.long_edge_counter[n2] < self.max_long_edge_per_node:
            self.addEdge(n1, n2, 'long')
            total_edges_added += 1
        if total_edges_added >= self.max_no_of_edges:
          return
    return
  
  def selectGateways(self):
    '''
    pre-condition: self.adj and self.short_edge_adj must be initialized
    post-condition:
      i) selects one gateway from each connected component from self.short_edge_adj (IMP!!) graph
      ii) selects the rest number of gateways from the rest of the nodes which has degree>0
    Args: None
    Returns: None
    '''
    
    self.gateways = []

    connected_components = nx.connected_components(self.short_edge_adj) 
        #IMP! finding connected components for short-edge-graph
    self.no_of_connected_components = 0
    
    for l in connected_components:
      self.no_of_connected_components += 1
      n = random.choice(l)
      self.gateways.append(n)
      
    num_gateways = int(round(1.0 * self.gateway_to_node_ratio*self.number_of_nodes,0))
     
    self.logger.debug("No. of connected components:"+str(self.no_of_connected_components))
    self.logger.debug("no. of gateways as permitted by the config ratio:"+str(num_gateways))
    
    if num_gateways<1:
      num_gateways =  1
      
    selection_attempts = 0
    while (num_gateways > len(self.gateways)):
      selection_attempts += 1
      if selection_attempts > self.short_edge_adj.number_of_nodes():
        break
      candidate_gateways = random.sample(self.short_edge_adj.nodes(), num_gateways - len(self.gateways))
      for i in candidate_gateways :     
        if i not in self.gateways and self.short_edge_adj.degree(i)>0:  
          self.gateways.append(i)
    self.logger.debug("total gateways:"+str(len(self.gateways)))
    return
  
  def isShortEdge(self,n1,n2):
    '''
    Args: n1,n2: node no, must be 0 <= and <self.number_of_nodes
    Returns: True if (n1,n2) is a short edge in self.adj graph, otherwise False
    '''
    if self.adj.has_edge(n1, n2):
      con_type =  self.adj[n1][n2]['con_type']
      if con_type == 'short':
        return True
    return False

  def getEdgeType(self, n1, n2):
    con_type = 'non_edge'
    if self.adj.has_edge(n1, n2):
      con_type =  self.adj[n1][n2]['con_type']
    return con_type
  
  def visualizeGraph(self, g, show_edges = True):
    '''
    pre-condition: must be called after self.adj and self.gateways are set or initialized
    visualize any graph that is subgraph of self.adj using pyplot
    Args: g, must be a networkx.Graph object and must be a subgraph of self.adj
    Returns: None
    '''
    #build node positions and node colors:
    node_positions = {}
    node_colors = []
    for n in g.nodes():
      node_positions[n] = (self.node_x[n], self.node_y[n])
      if n in self.gateways:
        node_colors.append('r')
      else:
        node_colors.append('w')
        
    
    #build edge_colors:
    edge_list = []
    if show_edges:
      edge_list = g.edges()
      
    edge_colors=[]
    for u,v in g.edges():
      if self.isShortEdge(u, v):
        edge_colors.append('r')
      else:
        edge_colors.append('k')
        
    nx.draw_networkx(G = g, 
                     pos = node_positions , 
                     with_labels = True, 
                     edgelist =  edge_list,
                     node_color = node_colors,
                     edge_color = edge_colors)

  
  
  
  