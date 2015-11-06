
import networkx as nx
from step_4_dynamic import Step_4_dynamic

class Step_4_static(Step_4_dynamic):
  '''
  This class handles all processing related to the static graph
  #-----class-field description------
  self.static_graph : the static graph, output of this class, built in runStep_4_static method
  ''' 
  def __init__(self,configFile):
    Step_4_dynamic.__init__(self, configFile)
    #---class fields------
    self.static_graph = None
    #---end of class fields
  
  def reset(self):
    Step_4_dynamic.reset(self)
    #---class fields------
    self.static_graph = None
    #---end of class fields
  def checkFSOPerNode(self):
    '''
    It can so happen that the existing degree of a non-gateway node in backbone graph
    already exceeds that max value specified in the config file and later stored in 
    self.fso_per_node, in that case, the field self.fso_per_node is updated and 
    and info-log message is displayed
    '''
    gateway_node_degrees =  self.backbone_graph.degree(self.gateways)
    max_gateway_node_degree = max(gateway_node_degrees.values())
    if max_gateway_node_degree > self.fso_per_gateway:
      self.logger.info("Max. Degree of gateway nodes exceed fso_per_gateway specified in config")
      self.logger.info("resetting fso_per_gateway from "
                       +str(self.fso_per_gateway)
                       +" to "
                       +str(max_gateway_node_degree))
      self.fso_per_node = max_gateway_node_degree
    
    if not self.non_gateway_backbone_nodes:
      return #there is no non-gateway nodes in the graph, so no need to update fso_per_node

    non_gateway_degrees =  self.backbone_graph.degree(self.non_gateway_backbone_nodes)
    max_non_gateway_degree = max(non_gateway_degrees.values())
    if max_non_gateway_degree > self.fso_per_node:
      non_gateway_node_with_max_degree = []
      for n,d in non_gateway_degrees.iteritems():
        if d==max_non_gateway_degree:
          non_gateway_node_with_max_degree.append(n)
    
      self.logger.info("Max. Degree of non-gateway nodes exceed fso_per_node specified in config")
      self.logger.info("Non gateway nodes with max. degree:"+str(non_gateway_node_with_max_degree))
      self.logger.info("resetting fso_per_node from "
                       +str(self.fso_per_node)
                       +" to "
                       +str(max_non_gateway_degree))
      self.fso_per_node = max_non_gateway_degree
  
  def hasRoomForNewEdgeInStaticGraph(self,n, no_of_new_edges = 1):
    '''
    Args:n: must be a node of static graph
         no_of_new_edges: no. of new edges to connect to this node n
    if n is a gateway and degree of n in backbone_graph<fso_per_gateway
      or n is a non-gatway and degree of n in backbone_graph<fso_per_node:
        Returns True
    else:
        Returns False
    '''
    #self.logger.debug("@hasRoomForNewEdgeInStaticGraph: n:"+str(n))
    if n in self.gateways:
      if self.static_graph.degree(n)<=self.fso_per_gateway - no_of_new_edges:
        return True
    else:
      if self.static_graph.degree(n)<=self.fso_per_node - no_of_new_edges:
        return True
    return False
  
  def initStaticGraph(self):
    '''
    check whether the max-degree of non-gateway nodes is <= self.fso_per_node as read from configfile
    if not, update it
    
    initializes the static graph from self.backbone_graph
    adds all nodes and edges from backbone_graph from 
    also, adds all possible edges between every node pair u,v as in the self.adj graph
      without violating(either the node is a gateway so no bound on its degree or the 
      degree of the node is < fso_per_node
  
    '''
    self.checkFSOPerNode()
    self.static_graph = nx.Graph(self.backbone_graph)
    self.static_graph.graph['name'] = 'Static graph'
    
    for u,v in self.adj.edges():
      if not self.backbone_graph.has_edge(u, v) and \
         self.backbone_graph.has_node(u) and \
         self.backbone_graph.has_node(v): #new edge
        if self.hasRoomForNewEdgeInStaticGraph(u) and\
           self.hasRoomForNewEdgeInStaticGraph(v):
          self.static_graph.add_edge(u,v)
  
  def isPathValidForStaticGraph(self,p): 
      ''' 
      for path p:
      the first node and last node must have room for one new edge
      all intermediate node that are in static graph must have degree permitting two new edges 
      '''
      if len(p)<3: #has to be with at least one new node 
        return  False
      
      if not self.hasRoomForNewEdgeInStaticGraph(p[0]):
        return False
      if not self.hasRoomForNewEdgeInStaticGraph(p[-1]):
        return False
      
      for i in p[1:-1]: #all intermediate nodes
        if self.static_graph.has_node(i) and \
        not self.hasRoomForNewEdgeInStaticGraph(i,2):
          return False
        
      return True

  def getVaidShortestPathsForStaticGraph(self, source, existing_nodes, new_nodes):
    '''
    Args: 
      source: the source of singel_source_shortest_path, must be a valid node in self.adj
      existing_nodes: a list of nodes, should exist in self.adj
      new_nodes: a list of nodes, should exist in self.adj and there souldn't be any intersection
                between existing_nodes and new_nodes
    solves single_source_shortest_path from 'source' to all nodes in self.adj
    for each path: check the followig:
      i) src and dest node are either gateways or each has degree<= fso_per_node - 1
      ii) all intermediate node n on the path, such that n is in static_graph:
            degree of n <= fso_per_node  - 2
      only then add the count of new nodes on the path to dict, add the path itself to 
      another dict
      otherwise the both the dict are 0 and empty list respectively
    Returns: 
        i) dict of of count of new nodes on the path
        ii) dict of the path itself
    '''

      
    length_of_shortest_path_from_source, shortest_path_from_source =\
             nx.single_source_dijkstra(self.adj, source)
    reachable_nodes = length_of_shortest_path_from_source.keys()
    
    count_of_new_nodes_on_paths = {}
    shortest_paths = {}
    
    for n in existing_nodes:
      count_of_new_nodes_on_paths[n] = 0
      shortest_paths[n] = []
      if n in reachable_nodes:
        shortest_path_to_n = shortest_path_from_source[n]
        #check path validity
        if self.isPathValidForStaticGraph(shortest_path_to_n ):
          for u in shortest_path_to_n:
            if u in new_nodes:
              count_of_new_nodes_on_paths[n]+=1
          if count_of_new_nodes_on_paths[n]>0:
            shortest_paths[n] =  list(shortest_path_to_n)
    return count_of_new_nodes_on_paths, shortest_paths
   
  def runStep_4_static(self):
    '''
    entry point for processing all tasks related to this class
    almost modeled on the steps of runStep_4_dynamic() of Step_4_dynamic
    except few exceptions

    processing:
      i) make an initial static graph self.static_graph as follows:
        1) initialize self.static_graph with all the nodes in self.backbone_graph
        2) add more edges between the existing nodes without violating the degree constraint 
          i.e. self.fso_per_node
      
      ii) create pool of nodes i.e. available nodes
           from new_nodes_for_step_4, that exist in self.adj but not in self.backbone_graph 
           IMP: available nodes will be reduced as the processing goes on
           
      iii) while (available nodes list is non-empty) and 
          (number of nodes in static graph < max_no_of_node for static graph)
          1. create a new flow-graph, adding super-source and super-sink and unit capacity to each
            link from the existing static graph
          2. solve max-flow on this flow-graph and compute residual graph
          3. for every node n in the residual graph:
            a) compute source-potential of n by treating n as the supersource
            b) compute sink-potential of n by treating n as the supersink
            
          4. for every pair of nodes u,v in the current static graph
            such that- degree of u < (fso_per_node-1) and degree of v <= (fso_per_node - 1)
            a) compute shortest path P between u and v 
               such that for each node n other than u,v that exists on path u-->v and
               self.static_graph, the degree of n <= (fso_per_node - 2)
            b)if
              there is new node on this path then compute path-benefit(u-->v) and path-benefit(v-->u)
              and save the path which has max-benefit and also the path
          5. if there is one such-max-benefit path, add the PATH (IMP: not all possible edges) 
              to the static graph and remove the new nodes on the path from available node list
          6. if no path was added in the last iteration, break the loop
    *most steps uses methods of class Step_4_dynamic except getNewNodesOnShortestPaths()
    rather this method uses getValidShortesPathForStaticGraph() which returns the full path (not just
    the list of new nodes) along with the count of new nodes on this path (default zero)
    '''
    #task (i)
    self.initStaticGraph()
    
    #task (ii)
    available_nodes = list(self.new_nodes_for_step_4)
    no_of_nodes_added_in_step_4 = 0
    #task (iii)
    while available_nodes and no_of_nodes_added_in_step_4<self.max_extra_nodes_for_step_4:
      #task (1) and (2)
      residual_graph = self.generateResidualGraph(self.static_graph)
      #task (3)
      source_potentials = self.getSourcePotentialOfAllNodes(residual_graph)
      sink_potentials = self.getSinkPotentialOfAllNodes(residual_graph)
      
      max_path_benefit = 0
      max_beneficial_path_new_node_count = 0
      max_beneficial_path = []
      
      for u in self.static_graph.nodes():
        if u not in self.gateways and self.static_graph.degree(u) > self.fso_per_node-1:
          continue
        #task (4)
        count_of_new_nodes_on_paths, shortest_paths = \
            self.getVaidShortestPathsForStaticGraph(u, self.static_graph.nodes(), available_nodes)
        for v in self.static_graph.nodes():
          if u==v:
            continue
          if v not in self.gateways and self.static_graph.degree(v) > self.fso_per_node-1:
            continue
        
          count_of_new_nodes_on_path_u_v =  count_of_new_nodes_on_paths[v]
          self.logger.debug("u,v,shortest_path:"+str(u)+" "+str(v)+" "+str(shortest_paths[v]))
          self.logger.debug("\tcount_of_new_nodes_on_path_u_v:"+str(count_of_new_nodes_on_path_u_v))
          if count_of_new_nodes_on_path_u_v >0:
            path_benefit_u_v = \
              min( source_potentials[u], sink_potentials[v])/(1.0*count_of_new_nodes_on_path_u_v)
            self.logger.debug("\tsource_potentials[u]:"+str(source_potentials[u]))
            self.logger.debug("\tsink_potentials[v]:"+str(sink_potentials[v]))
            self.logger.debug("\tpath_benefit_u_v:"+str(path_benefit_u_v))
            if path_benefit_u_v>max_path_benefit:
              max_path_benefit = path_benefit_u_v
              max_beneficial_path = list(shortest_paths[v])
              max_beneficial_path_new_node_count =  count_of_new_nodes_on_path_u_v
              self.logger.debug("max_path_benefit changed to:"+str(max_path_benefit))
              self.logger.debug("max_beneficial_path changed to:"+str(max_beneficial_path))
      #task (5)
      if max_beneficial_path:
        self.logger.info("Adding path to static graph:"+str(max_beneficial_path))
        self.static_graph.add_path(max_beneficial_path)
        no_of_nodes_added_in_step_4 += max_beneficial_path_new_node_count
        available_nodes = list( set(available_nodes) - set(max_beneficial_path) ) 
      else:
        break #couldn't add any new node in the last iteration, so processing is over
    
    
    
    
    
    
    
    
    