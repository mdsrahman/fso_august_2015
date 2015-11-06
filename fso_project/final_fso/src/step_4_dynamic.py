
from step_3 import Step_3
import networkx as nx
import networkx.algorithms.flow as flow
import logging

class Step_4_dynamic(Step_3):
  '''
  class that adds extra edges/nodes/paths to the backbone graph to build dynamic graph
  self.max_extra_nodes_for_step_4 (int): maximum number of nodes to be added in step 4 for
                    both dynamic and static graphs using self.ratio_of_max_added_nodes_in_step_4
                    which is set by config spec
  new_nodes_for_step_4(list of int): list of nodes that exist 
                        in self.adj but not in self.backbone_graph 
  dynamic_graph (networkx.Graph obj): the dynamic graph, output of this class,
                                    initialized in initDynamicGraph() method and changed/updated/built
                                    in run_Step_4_dynamic() method
  '''
  def __init__(self,configFile):
    Step_3.__init__(self,configFile)
    #---class fields----
    self.max_extra_nodes_for_step_4 = None
    self.new_nodes_for_step_4 = None
    self.dynamic_graph =  None
    #----end of class fields----
    self.logger.setLevel(logging.INFO) #just for debugging
    self.logger.debug("Step_4 constructor")
  
  def reset(self):
    Step_3.reset(self)
    #---class fields----
    self.max_extra_nodes_for_step_4 = None
    self.new_nodes_for_step_4 = None
    self.dynamic_graph =  None
    #----end of class fields----
    self.logger.setLevel(logging.INFO) #just for debugging

    
  def initDynamicGraph(self):
    '''
    creates an initial dynamic graph from backbone_graph
    '''
    self.dynamic_graph = nx.Graph()
    self.dynamic_graph.graph['name'] = 'Dynamic Graph'
    
    self.dynamic_graph.add_nodes_from(self.backbone_graph.nodes())
    
    for u,v in self.adj.edges():
      if self.dynamic_graph.has_node(u) and self.dynamic_graph.has_node(v):
        self.dynamic_graph.add_edge(u,v)
    self.logger.debug("dynamic graph initialized")
    
  def updateDynamicGraph(self, node_list):
    '''Args: node_list: list of nodes to be added to the dynamic graph
      adds all the nodes n in the node_list and then all edges between n and all other nodes of
      the dynamic graph as those exist in self.adj graph
    '''
    #task: add new nodes
    self.dynamic_graph.add_nodes_from(node_list)
      
    #task: add new edges, may already exist in the graph but just overwritten, no side-effect
    for n in node_list:
      nbrs =  self.adj.neighbors(n)
      for u in nbrs:
        if self.dynamic_graph.has_node(u):
          self.dynamic_graph.add_edge(n,u)
          
  def initMaxAllowedNodesForStep_4(self):
    '''
    initializes the extra number of nodes and hence, the pool of new nodes to be added in step 4
    these fields are used in both dynamic and static graph building
    ''' 
    self.max_extra_nodes_for_step_4 = \
        int(round(self.ratio_of_max_added_nodes_in_step_4*self.backbone_graph.number_of_nodes(),0))
    
    self.new_nodes_for_step_4 = \
      list(set(self.adj.nodes()) -  set(self.backbone_graph.nodes())) 
      #finds nodes that are in adj graph but not in backbone graph
  
  def generateResidualGraph(self, g):
    '''
    processing:
      given graph g(networkx.Graph instance):
      i) create a new Digraph obj  instance max_flow_g from g
      ii) add 'capacity' of 1 in each edge of g
      ii) add super source 'src' and add a directed edge to each non-gateway node n in g with 
          capacity infinite i.e. src->n with capacity infinity
      iii) add super sink 'snk' and add a directed edge from each gateway node n in g with to 'snk'
          with capacity infinite i.e.n->snk with capacity infinity 
      iv) run max-flow on the max_flow_g
      
      v) compute the residual graph in max_flow_g by setting the capacity of each edge 
        as (capacity - flow) and the flow value to zero
        and return max_flow_g (note, 'src' and 'snk' exists in the residual graph)
    '''   
    #task (i)
    non_gateway_nodes_in_g =  list(set(g.nodes()) -  set(self.gateways))
    
    max_flow_g =  nx.DiGraph()
    max_flow_g.add_nodes_from(g.nodes())
    for u,v in g.edges():
      max_flow_g.add_edge(u, v, capacity = 1.0)
      max_flow_g.add_edge(v, u, capacity = 1.0)
      
    #task (ii)
    max_flow_g.add_node('src')
    for n in non_gateway_nodes_in_g:
      max_flow_g.add_edge('src', n, capacity = float('inf'))
    
    #task (iii)
    max_flow_g.add_node('snk')
    for n in self.gateways:
      max_flow_g.add_edge(n, 'snk', capacity = float('inf'))
    
    #task (iv)
    residual_max_flow_g = flow.shortest_augmenting_path(G=max_flow_g, s='src', t='snk')
    
    #task(v)
    for u,v in residual_max_flow_g.edges():
        residual_max_flow_g[u][v]['capacity'] -= residual_max_flow_g[u][v]['flow']
        residual_max_flow_g[u][v]['flow'] = 0.0
        
    return residual_max_flow_g
  
  def getSourcePotentialOfAllNodes(self,g):
    '''
    Args: g is a residual di-graph with 'src' and 'snk' node
    processing:
      i) create a new graph max_flow_g from g and remove node 'src' from max_flow_g
      ii) for each node n in max_flow_g:
          if n is gateway: set source potential of n as infinity in the dict
          else:
            solve max-flow treating n as source and 'snk' as sink 
            (there shouldn't be any edge n->'snk')
            and store the max_flow value in the dict
      return the dict
    '''
    source_potentials = {}
    max_flow_g = nx.DiGraph(g)
    max_flow_g.remove_node('src')
    for n in max_flow_g.nodes():
      if n in self.gateways:
        source_potentials[n]=float('inf')
      elif n!='snk':
        residual_max_flow_g = flow.shortest_augmenting_path(G=max_flow_g, s=n, t='snk')
        source_potentials[n] = residual_max_flow_g.graph['flow_value']
    return source_potentials
   
  def getSinkPotentialOfAllNodes(self,g):
    '''
    Args: g is a residual di-graph with 'src' and 'snk' node
    processing:
      i) create a new graph max_flow_g from g and remove node 'snk' from max_flow_g
      ii) for each node n in max_flow_g:
          if n is gateway: set sink potential of n as infinity in the dict
          else:
            temporarily remove the edge 'src'->n (there should be  edge 'src'->n)
            solve max-flow treating n as sink and,'src' as source 
            and store the max_flow value in the dict
            restore the edge 'src'->n with capacity as before and flow_val as zero
      return the dict
    '''
    sink_potentials = {}
    max_flow_g = nx.DiGraph(g)
    max_flow_g.remove_node('snk')
    
    for n in max_flow_g.nodes():
      if n in self.gateways:
        sink_potentials[n]=float('inf')
      elif n!='src':
        backup_capacity_src_to_n = max_flow_g['src'][n]['capacity']
        backup_capacity_n_to_src = max_flow_g[n]['src']['capacity']
        max_flow_g.remove_edge('src', n)
        max_flow_g.remove_edge(n, 'src')
        
        residual_max_flow_g = flow.shortest_augmenting_path(G=max_flow_g, s='src', t=n)
        sink_potentials[n] = residual_max_flow_g.graph['flow_value']
        
        max_flow_g.add_edge('src', n, capacity =  backup_capacity_src_to_n, flow=0.0)
        max_flow_g.add_edge(n, 'src', capacity =  backup_capacity_n_to_src, flow=0.0)
    return sink_potentials
  
  def getNewNodesOnShortestPaths(self, source, existing_nodes, new_nodes):
    '''
    Args: 
      source: the source of singel_source_shortest_path, must be a valid node in self.adj
      existing_nodes: a list of nodes, should exist in self.adj
      new_nodes: a list of nodes, should exist in self.adj and there souldn't be any intersection
                between existing_nodes and new_nodes
    solves single_source_shortest_path from 'source' to all nodes in self.adj
    then gets the  nodes from 'new_nodes' on each path (if any) from 'source'
    to each of the nodes in 'existing_nodes', if no path exist,then defaults to empty list
    Returns: the dict of list of nodes from existing_nodes on shortest path 
    to each node in existing_nodes(dict key) 
    '''

      
    length_of_shortest_path_from_source, shortest_path_from_source =\
             nx.single_source_dijkstra(self.adj, source)
    reachable_nodes = length_of_shortest_path_from_source.keys()
    
    new_nodes_on_shortest_paths = {}
    
    for n in existing_nodes:
      new_nodes_on_shortest_paths[n] = []
      
    reachable_existing_nodes = list(set(existing_nodes) & set(reachable_nodes))
    
    for n in reachable_existing_nodes :
      shortest_path_to_n = shortest_path_from_source[n]
      new_nodes_on_path_to_n =  set(new_nodes) & set(shortest_path_to_n)
      new_nodes_on_shortest_paths[n] =  list(new_nodes_on_path_to_n)
      
    return new_nodes_on_shortest_paths
  
  def runStep_4_dynamic(self):
    '''
    method that completes all the task related to step 4 dynamic
    processing:
      i) make an initial dynamic graph self.dynamic_graph as follows:
        1) initialize self.dynamic_graph with all the nodes in self.backbone_graph
        2) add all edges between pairs of nodes in self.dynamic_graph according to self.adj  
      
      ii) compute the max. number of nodes to be added in step 4, i.e. self.max_extra_nodes_for_step_4
          using ratio_of_max_added_nodes_in_step_4 and hence max_no_of_node for dynamic graph and
          static graph
      iii) create pool of nodes i.e. available_nodes from new_nodes_for_step_4, 
          that exist in self.adj but not in self.backbone_graph 
          IMP: available nodes will be reduced in size as the processing goes on
      iv) while (available_nodes is non-empty) and 
          (number of nodes in dynamic graph < max_no_of_node for dynamic graph)
          1. create a new flow-graph, adding super-source and super-sink and unit capacity to each
            link from the existing dynamic graph
          2. solve max-flow on this flow-graph and compute residual graph
          3. for every node n in the residual graph:
            a) compute source-potential of n by treating n as the supersource
            b) compute sink-potential of n by treating n as the supersink
          4. for every pair of nodes u,v in the current dynamic graph:
            a) compute shortest path between u and v and no. of nodes from new_nodes_for_step_4 
              in the path
            b)if
              there is new node on this path then compute path-benefit(u-->v) and path-benefit(v-->u)
              and save the path which has max-benefit and also the path
          5. if there is one such-max-benefit path, add the nodes to the dynamic graph along with
             all possible edges and remove the new nodes from available node list
          6. if no path was added in the last iteration, break the loop
      v) prepare dynamic graph spec for TM's java code
      vi) call TM's java code and retrieve upper and lower-bound of dynamic graph, 
          specifying the number of patterns as percent of total nodes in dynamic graph
          and the path to the dynamic graph file

    '''
    #task (i)
    self.initDynamicGraph()
    self.initMaxAllowedNodesForStep_4()
    #task (ii)+(iii)
    available_nodes = list(self.new_nodes_for_step_4)
    no_of_nodes_added_in_step_4 = 0
    #task (iv)
    while available_nodes and no_of_nodes_added_in_step_4<self.max_extra_nodes_for_step_4:
      self.logger.debug("no_of_nodes_added_in_step_4:"+str(no_of_nodes_added_in_step_4)
                       +":elapsed time:"+str(self.getElapsedTime()))
      #task (1) and (2)  
      residual_graph = self.generateResidualGraph(self.dynamic_graph)
      #task (3)
      source_potentials = self.getSourcePotentialOfAllNodes(residual_graph)
      sink_potentials = self.getSinkPotentialOfAllNodes(residual_graph)
      
      max_path_benefit = 0
      max_beneficial_new_nodes = []
      
      #debug_processed_node = 0
      debug_total_node = self.dynamic_graph.number_of_nodes()
      
      for u in self.dynamic_graph.nodes():
        #task (4).a
#         debug_processed_node += 1
#         if debug_processed_node%100==0:
#           self.logger.info("processing node/total_node:"+str(debug_processed_node)+"/"+str(debug_total_node))
        if source_potentials[u]<=0:
          continue
        new_nodes_on_shortest_path_from_u = \
            self.getNewNodesOnShortestPaths(u, self.dynamic_graph.nodes(), available_nodes)
        for v in self.dynamic_graph.nodes():
          if u==v:
            continue
          if sink_potentials[v]<=0:
            continue
          #task (4).b
          new_nodes_on_path_u_v = new_nodes_on_shortest_path_from_u[v]
          count_of_new_nodes_on_path_u_v = len(new_nodes_on_path_u_v)

          if count_of_new_nodes_on_path_u_v>0:
            path_benefit_u_v = \
             min( source_potentials[u], sink_potentials[v])/(1.0*count_of_new_nodes_on_path_u_v)
            if path_benefit_u_v>max_path_benefit:
              max_path_benefit = path_benefit_u_v
              max_beneficial_new_nodes = list(new_nodes_on_path_u_v)
              self.logger.debug("new max_path_benefit:"+str(max_path_benefit))
              self.logger.debug("new max_beneficial_new_nodes:"+str(max_beneficial_new_nodes))
      if max_beneficial_new_nodes: #at least this list contains one node, so has been set in the loops above
        #task (5)
        self.updateDynamicGraph(max_beneficial_new_nodes)
        no_of_nodes_added_in_step_4 += len(max_beneficial_new_nodes)
        self.logger.info("New Nodes added to dynamic graph:"+str(max_beneficial_new_nodes))
        available_nodes = list( set(available_nodes) - set(max_beneficial_new_nodes) ) 
      else:
        break #couldn't add any new node in the last iteration, so processing is over
        
      
      
      
      
      
      
    
    
    
    