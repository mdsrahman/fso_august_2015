from final_fso.src.step_2 import Step_2
import networkx as nx

import logging
from matplotlib import pyplot as plt

class Step_3(Step_2):
  '''
  class to reduce node degrees in the backbone graph
  ---------------class-fields---------------------------------------------
  backbone_graph_after_step_2(networkx.Graph() obj): a backup version of backbone
                            graph before degree reduction, set in runStep_3() method
  self.non_gateway_backbone_nodes(list of int): non-gateway nodes in self.backbone_graph so far,
                            initialized in runStep_3() method and used by getMaxDegreeNodes() method
  '''
  def __init__(self,configFile):
    Step_2.__init__(self, configFile) #constructor for super-class must top everything
    #--class-fields------
    self.backbone_graph_after_step_2 = None
    self.non_gateway_backbone_nodes = [] 
    self.logger.debug("Step_3 initialized")
    #---end-of-class fields----
    self.logger.setLevel(logging.INFO)
  
  def reset(self):
    Step_2.reset(self) 
    #--class-fields------
    self.backbone_graph_after_step_2 = None
    self.non_gateway_backbone_nodes = [] 
    self.logger.debug("Step_3 initialized")
    #---end-of-class fields----
    self.logger.setLevel(logging.INFO) 
    
  def getMaxDegreeNodes(self):
    '''
    Args: None
    Processing: Returns the nodes in self.backbone_graph with maximum degrees (>2)
      i) find the dict  node_degrees of self.backbone_graph
      ii) find the max_degree from values of of the node degrees
      iii) if max_degree <=2: return max_degree, empty list
      iv) append all nodes (node_degree keys) having degree=max_degree to maxDegreeNodes list
      v) return max_degree, maxDegreeNodes list
    '''
    max_degree_nodes = []
    
    node_degrees =  self.backbone_graph.degree()
    max_degree = max(node_degrees.values())
    if max_degree<=2:
      return max_degree, max_degree_nodes #max_degree_nodes is empty now
    
    for n,d in node_degrees.iteritems():
      if d==max_degree:
        max_degree_nodes.append(n)
      
    return max_degree, max_degree_nodes
  
  def getCandidateNodes(self, max_degree):
    '''
    Args: max_degree: int such that <2
    Processing: Returns all in self.backbone_graph with degree d such that 1<=d<=max_degree-2
      i) set candidateNodes as all the gateways (must be present in self.backbone_graph)
      i) if max_degree<2, return candidateNodes i.e. []
      ii) find the dict  node_degrees of self.backbone_graph
      iii) scan node_degrees for every  node having degree d such that
                1<= d <=max_degree-2  *for non-gateway nodes, 1<=d ensures that it is connected
                and 0<=d<=max_degree_2 for gateway nodes,
          and append it to candidateNodes list
      iv) return candidateNodes list
    '''
    candidateNodes = []
    if max_degree<2:
      return candidateNodes
    node_degrees =  self.backbone_graph.degree()
    for n,d in node_degrees.iteritems():
      if d<=max_degree-2:
        if n in self.gateways or 1<=d:
          candidateNodes.append(n)
    return candidateNodes
  
  def isValidForConnection(self,n,u,v):
    '''
    called to check whether it is valid to remove edge n-u and add edge u-v, 
    Args: n, u, v: nodes such that each 0<= and <self.number_of_nodes
    Returns:
      True if 
            there is a 'short' edge u-v in self.adj graph
            v is either a gateway or still connected to one of the gateways in self.backbone_graph
            after the edge n-u is removed
    '''

    if not self.isShortEdge(u,v): #if self.adj does not have a short edge u-v
      return False
    
    if v in self.gateways: 
      return True
    
    #temporarily remove edge n-u from self.backbone_graph
    self.backbone_graph.remove_edge(n, u)
    
    v_is_gateway_connected = False
    for g in self.gateways:
      if nx.has_path(self.backbone_graph, g, v):
        v_is_gateway_connected = True
        break
    
    #add back the edge n-u to self.backbone_graph
    self.backbone_graph.add_edge(n, u)
    return v_is_gateway_connected
  
  def getSuccessorInBFSTree(self,n):
    '''
    Args: n: node
    Returns: the list of successors nodes bfs_successor_n of node n in a tree rooted at
           a super-gateway (artificial gateway connected to only gateways)
           in the graph self.backbone_graph
    '''
    bfs_successors_n = []
    if not self.backbone_graph.has_node(n):
      return bfs_successors_n
    #add super gateway node sg temporarily
    for g in self.gateways:
      self.backbone_graph.add_edge('sg', g)
      
    all_bfs_successors = nx.bfs_successors(self.backbone_graph,'sg')
    #self.logger.debug("all_bfs_successors.keys()"+str(sorted(all_bfs_successors.keys())))
    if n in all_bfs_successors.keys():
      bfs_successors_n = list(all_bfs_successors[n])
    #remove super gateway node sg, this will effectively remove all its edges too 
    self.backbone_graph.remove_node('sg')
    return bfs_successors_n
  
  def runStep_3(self):
    '''
    processing:
    only in self.backbone_graph:
      ia) intialize non_gate_way_backbone_nodes
      i) check every  non-gateway node n with max_degree d> 2,
      ii)  find a list of nodes cnodes with degree at most d-2:
      iii)    run bfs starting at every gateway g
      iv)    find n in one of the bfs-paths
      v)    find n's successor in bfs-path i.e. n_bfs_successors
      vi)    for each u in n_bfs_successors:
      vii)      if there is an edge between u-c in self.adj graph for any node c in cnodes
                and (c is either a gateway or still has a path to one of the gateways 
                      even if n-u is disconnected):
      viii)        disconnect n-u and connect u-c
      ix)        continue to step (i) if step (viii) worked this iteration
      x)   if no such u-c found for any node n, then stop the process
    '''
    self.backbone_graph_after_step_2 = nx.Graph(self.backbone_graph)
    self.backbone_graph_after_step_2.graph['name']='Backbone Graph After Step 2'
    self.non_gateway_backbone_nodes =  list(set(self.backbone_graph.nodes()) - set(self.gateways))
    
    
    degree_reduced_in_last_iteration = True
    
    while degree_reduced_in_last_iteration:
      degree_reduced_in_last_iteration = False
      #task (i)
      max_degree, max_degree_nodes =  self.getMaxDegreeNodes()
      self.logger.debug("max_degree:"+str(max_degree))
      self.logger.debug("max_degree_nodes:"+str(max_degree_nodes))
      if not max_degree_nodes:
        break # no nodes with degree>2 or no non-gateway nodes found, so terminate processing now
      #task (ii) 
      candidateNodes = self.getCandidateNodes(max_degree)
      self.logger.debug("candidateNodes:"+str(candidateNodes))
      for n in max_degree_nodes:
        #task (iii)+task (iv)+task(v)
        debug_is_n_gateway = n in self.gateways
        self.logger.debug("gateways:"+str(self.gateways))
        self.logger.debug("self.backbone_graph.nodes():"+str(sorted(self.backbone_graph.nodes())))
        self.logger.debug("max_degree_nodes:"+str(sorted(max_degree_nodes)))
        self.logger.debug("n,is_n_gateway,nbrs,:"+str(n)+","+str(debug_is_n_gateway)+\
                           ", "+str(sorted(self.backbone_graph.neighbors(n))) ) 
        bfs_successors_n  =  self.getSuccessorInBFSTree(n)
        self.logger.debug("bfs_successors_n:"+str(bfs_successors_n))
        #task (vi)
        for u in bfs_successors_n:
        #task (vii)
          for v in candidateNodes:
            n_u_v_connectionValidity = self.isValidForConnection(n, u, v)
            self.logger.debug("n,u,v,n_u_v_connectionValidity"
                          +str(n)+" "
                          +str(u)+" "
                          +str(v)+" "
                          +str(n_u_v_connectionValidity))
            if n_u_v_connectionValidity:
              #task (viii)
              self.backbone_graph.remove_edge(n, u)
              self.backbone_graph.add_edge(u, v)
              
              degree_reduced_in_last_iteration = True
              break # no more searching in the candidateNodes list
          if degree_reduced_in_last_iteration:
            break #no more searching in the bfs_successor list of node n
        if degree_reduced_in_last_iteration:
            break #no more searching for nodes in max_degree_nodes list
        
      
    
    
    
    
    