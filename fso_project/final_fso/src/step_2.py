import networkx as nx


from step_1 import Step_1
from my_util import MyHeap

class Step_2(Step_1):
  '''
  description of the class fields
  backbone_graph(networkx.Graph obj): background graph generated after running the BFS,
                          initialized and set in runStep_2() method
                          
  bfs_nodes(list of int): each element is node index n, 0 <= n <self.number_of_nodes,
                        these list includes all the gateway nodes plus the target-covering 
                        nodes generated in the previous step in list self.node_cover
  
  bfs_node_index(list of int): of length self.number_of_nodes, keeps index of every node n
                        in self.bfs_nodes, if n is not in self.bfs_nodes, then value is -1 
                        i.e. invalid index, otherwise values is 0<= and < len(self.bfs_nodes)
                        initialized or set in self,buildBFSPaths() method
                        also used by self.getBFSPath() and self.getBFSPathDistance() methods
  
  bfs_paths(2D list of path P(another list)): bfs_paths[u][v] is the list of nodes on the bfs_path
                        from u to v
  
  bfs_path_distance(2D list of int): bfs_path_distance[u][v] is the no. of hops in the bfs-path
                                    between node u and v, default 0 means no path exists
  
  bfs_heap(MyHeap obj): implements a priority queue, queue-object is tuple (u,v) where
                         u and v are nodes from self.bfs_nodes and queue-priority is the bfs-distance
                        between u and v, i.e. self.bfs_path_distance[u][v]
                        IMPORTANT: never push (u,v) with ZERO BFS-PATH length, it is interpreted 
                                  as non-existent path, 
                                  so check priority>0 before every bfs_heap.push() method call
  remaining_bfs_node_cover(list of int): initialized, popped and emptied  in self.runStep_2() method
                          initially all the nodes from self.node_cover, but gradually popped for
                          each node visited in the BFS
  '''
  def __init__(self, configFile):
    Step_1.__init__(self, configFile) 
        #super class initializer must be called at the top
        
    #---class fields------------   
    self.backbone_graph = None
    self.bfs_nodes = None
    self.bfs_node_index = None
    self.bfs_paths =  None
    self.bfs_path_distance = None
    self.bfs_heap = None
    self.remaining_bfs_node_cover = None
    #---end of class fields----- 
    
  def reset(self):
    Step_1.reset(self)
        #super class initializer must be called at the top
        
    #---class fields------------   
    self.backbone_graph = None
    self.bfs_nodes = None
    self.bfs_node_index = None
    self.bfs_paths =  None
    self.bfs_path_distance = None
    self.bfs_heap = None
    self.remaining_bfs_node_cover = None
    
  def getBFSPath(self, u, v):
    '''
    Args:
      u, v: node no, must be 0<= and < self.number_of_nodes
    Returns:
      returns the corresponding bfs-path between u and v from self.bfs_paths
    '''
    u_index = self.bfs_node_index[u]
    v_index = self.bfs_node_index[v]
    bfs_path_uv =  list(self.bfs_paths[u_index][v_index]) 
    return bfs_path_uv
  
  def getBFSPathDistance(self, u, v):
    '''
    Args:
      u, v: node no, must be 0<= and < self.number_of_nodes
    Returns:
      returns the corresponding bfs-path distance between u and v from self.bfs_path_distance
    '''
    u_index = self.bfs_node_index[u]
    v_index = self.bfs_node_index[v]
    bfs_path_distance_uv =  self.bfs_path_distance[u_index][v_index] 
    return bfs_path_distance_uv
  
  def buildBFSPaths(self):
    '''
    processing:
      i) initializes self.bfs_nodes 
          and self.bfs_node_index (for fast access of node-index in self.bfs_nodes)
          
      ii) push nodes from self.gateways and self.node_cover into self.bfs_nodes and
          at the same time
            generate index-value for each n in self.bfs_nodes for entry in self.bfs_node_index
      
      iii) initializes self.bfs_paths (2D array of list of int as nodes on the bfs-path)
          and self.bfs_path_distance
          
      iv) for each pair of nodes (u,v) in self.bfs_nodes:
            a) self.bfs_paths[u][v]<- bfs_path (essentially list of nodes) between u and v
            b) self.bfs_path_distance[u][v]<- bfs_path length between u and v 
              in graph IMPORTANT!!: self.short_edge_adj, graph with only edges=short-link-fso
    '''
    #task (i)
    self.bfs_nodes = []
    self.bfs_node_index = []
    for i in xrange(self.number_of_nodes):
      self.bfs_node_index.append(-1)
    
    #task (ii)
    index_counter = 0
    for n in self.gateways:
      self.bfs_nodes.append(n)
      self.bfs_node_index[n] = index_counter
      index_counter += 1
    
    node_cover_without_gateways = list(set(self.node_cover) - set(self.gateways)) 
      #exclude the gateways to prevent duplication
    for n in node_cover_without_gateways:
      self.bfs_nodes.append(n)
      self.bfs_node_index[n] = index_counter
      index_counter += 1
    
    #task (iii)
    total_bfs_nodes = len(self.bfs_nodes)
    
    self.bfs_path_distance=[]
    for i in xrange(total_bfs_nodes):
      self.bfs_path_distance.append([])
      for j in xrange(total_bfs_nodes):
        self.bfs_path_distance[i].append(0) 
    
    
    self.bfs_paths = []
    for i in xrange(total_bfs_nodes):
      self.bfs_paths.append([])
      for j in xrange(total_bfs_nodes):
        self.bfs_paths[i].append([])
    
    #task (iv)
    for i_index,i in enumerate(self.bfs_nodes):
      shortest_paths= nx.single_source_shortest_path(self.short_edge_adj, i)
      reachable_nodes_on_shortest_paths = shortest_paths.keys() 
        #reachable_nodes_on_shortest_paths: all nodes reachable from node i
      reachable_bfs_nodes =  list(set(reachable_nodes_on_shortest_paths) & set(self.bfs_nodes))
        #reachable_bfs_nodes: all nodes in self.bfs_nodes that are reachable from node i
      for j in reachable_bfs_nodes:
        j_index = self.bfs_node_index[j]
        self.bfs_path_distance[i_index][j_index] = len(shortest_paths[j])
        self.bfs_paths[i_index][j_index] = list(shortest_paths[j])
        

  def runStep_2(self):
    '''
    processing:
      i) call self.buildBFSPaths() to set self.bfs_paths and self.bfs_distance
      ii) initialize the backbone_graph self.backbone_graph with nodes from self.bfs_nodes,
          so self.backbone_graph will contain all the gateways as well
      iii) create min-heap self.bfs_heap and push all the bfs-paths 
            between every gateway node g to all other node n [discard same pairs (a,a)]
      iv) initialize bfs_cover_node <- self.node_cover
      
      v) while heap is not empty and bfs_cover_node is non-empty
        a)  pop (i_index,j_index) from heap and hence i and j
        b) if j in bfs_cover_node:
          1. add bfs_path from i to j to self.backbone_graph
          2. remove j from bfs_cover_node
          3. for all node k in self.bfs_nodes:
              if j!=k
              push (j,k) into heap with priority <-self.bfs_distance[j][k]
              
    '''
    #task (i)
    self.buildBFSPaths()
    
    #task (ii)
    self.backbone_graph=nx.Graph()
    self.backbone_graph.graph['name'] = 'Backbone Graph'
    
    for n in self.bfs_nodes: #already pushing all gateways and all node_covers found in step i
      self.backbone_graph.add_node(n)
    
    #task (iii)
    self.bfs_heap = MyHeap()
    
    for g in self.gateways:
      for n in self.bfs_nodes:
        if n in self.gateways:
          continue #no point in finding shortest path between set of gateways
        bfs_distance_gn = self.getBFSPathDistance(g, n)
        if g!=n and bfs_distance_gn>0: 
            #don't push (g,g) pairs i.e loop-edge and don't push non-path/zero-len-path pairs
          self.bfs_heap.push((g,n), bfs_distance_gn )
    
    #task (iv)
    self.remaining_bfs_node_cover = list(set(self.node_cover) - set(self.gateways)) 
          #IMP: only non-gateway node-covers should be covered
    
    #task (v)
    while not self.bfs_heap.isHeapEmpty() and self.remaining_bfs_node_cover:
      i, j = self.bfs_heap.pop()
      if j in self.remaining_bfs_node_cover:
        bfs_path_ij =  self.getBFSPath(i, j)
        #bfs_path_distance_ij =  self.getBFSPathDistance(i, j)
        self.backbone_graph.add_path(bfs_path_ij)
        self.remaining_bfs_node_cover.remove(j)
        for k in self.bfs_nodes:
          bfs_distance_jk = self.getBFSPathDistance(j, k)
          if j!=k and bfs_distance_jk>0: 
              #don't push (g,g) pairs i.e loop-edge and don't push non-path/zero-len-path pairs
            self.bfs_heap.push((j,k), bfs_distance_jk)
      
          
      
               
            
            
            
