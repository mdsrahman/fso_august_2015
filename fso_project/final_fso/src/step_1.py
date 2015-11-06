
from my_util import MyHeap
from my_util import MyGridBin
from input_generator import InputGenerator
import networkx as nx
 
class Step_1(InputGenerator): 
  '''
  #TODO: complete field description
  node_bin(MyGridBin obj): bins the nodes in 2D grid, 
    used for quick retrieval of nodes in heuristic set cover method
  node_hp(MyHeap obj): node priority queue, is built in buildNodeHeap(..) method,
                      here priority value is the negative of total number of nodes covered by that node
  max_tx, max_ty: just like max_x, max_y except these are the maximum grid index for targets
                  i.e relative coordinate derived by dividing the whole area by target_spacing
  total_targets(int): total number of targets to be covered
  no_of_targets_covered(int): targets covered by nodes popped from the node heap node_hp
  node_target_count(list of int): keeps track of total targets covered by node n(the index), 
                                initialized in buildNodeHeap() method
  covered_target(2D array of bool): keeps track of which target has already been covered by an existing node
                      popped from node_hp, initialized in buildNodeHeap() method
  node_cover(list of int): list of target-covering nodes, set in the method heuristicTargetCover(),
                          in essence, the output of this whole class
                          might include gateways as well
  '''


  def __init__(self, configFile):
    '''post-condition:
      i) bins the nodes
      ii) build the heap
      iii) call the heuristic set cover
    '''
    InputGenerator.__init__(self, configFile) #super class initializer must be called at the top
    #------------class fields--------------------
    self.node_bin = None
    self.node_hp = None
    
    self.max_tx = int(self.max_x_coord // self.target_spacing)
    self.max_ty = int(self.max_y_coord // self.target_spacing)
    self.total_targets = (self.max_tx+1)*(self.max_ty+1)
    
    self.no_of_targets_covered = 0
    self.node_target_count = None
    self.covered_target = None
    self.node_cover = []
    #---end of class fields
    
  def reset(self): 
    self.node_bin = None
    self.node_hp = None
    
    self.max_tx = int(self.max_x_coord // self.target_spacing)
    self.max_ty = int(self.max_y_coord // self.target_spacing)
    self.total_targets = (self.max_tx+1)*(self.max_ty+1)
    
    self.no_of_targets_covered = 0
    self.node_target_count = None
    self.covered_target = None
    self.node_cover = []
    
  def runStep_1(self):
    '''runs all the tasks for step 1
    '''
    
    self.binNodes()
    self.logger.debug("Completed binning nodes")
    self.buildNodeHeap()
    self.logger.debug("Completed building node heap")
    self.heuristicTargetCover()
    self.logger.debug("Completed heuristic set cover and hence step 1")
    
  def binNodes(self):
    '''
    bins the nodes in to a 2D bin according to their coordinate, max_short_edge_length or short-fso-link-distance
    is used for bin_interval along both directions
    Args: None
    Returns: None
    '''
    self.node_bin = MyGridBin(x_interval = self.coverage_radius,
                              y_interval = self.coverage_radius,
                              max_x = self.max_x_coord,
                              max_y = self.max_y_coord)
    
    for n in xrange(self.number_of_nodes):
      self.node_bin.put(n, self.node_x[n], self.node_y[n])
  
  def isWithinDistance(self,x1,y1,x2,y2, d):
    '''
    Args:
      x1,x2,y1,y2: x and y cooridnates of point 1 and 2
      d: the distance to be tested against
    Returns:
      True: if (x1,y1) and (x2,y2) are within (Euclidian) distance d of each other
      False: otherwise
    '''
    current_distance = (x1 - x2)**2 + (y1 - y2)**2
    if current_distance <= d*d:
      return True
    else:
      return False
    
  def getTargetsCoveredByNode(self,n):
    '''
    Args:n, node index, 0<= and <self.number_of_nodes
    Returns, covered_target_grid_index:
        list of tuples (tx,ty) such that tx and ty are the target-index in the rectangular grid
    
    processing description:
        
    i) find 4 tx,ty coordinate sets for n_min_tx,n_max_tx, n_min_ty, n_max_ty
        make sure they are within the range [0,tx_max] and [0,ty_max]
        
    ii) append every target in the rectangular grid [n_tx_min,n_tx_max][n_ty_min,n_ty_max] to the return_target_list
      but do the following test only for grids such that it's tx,ty are the extremities
        test: find the absolute coordinate of the middle point of target-grid, find the sq-distance
        of this point to the center of the node-coordinate, it this less than the sq-radius, only then add it 
        to the return list
    '''  
    n_x = self.node_x[n]
    n_y = self.node_y[n]
    
    n_min_tx =  int(max(0, n_x - self.coverage_radius) // self.target_spacing)
    n_max_tx =  int(min(self.max_x_coord, n_x + self.coverage_radius) // self.target_spacing)
    
    n_min_ty =  int(max(0, n_y - self.coverage_radius) // self.target_spacing)
    n_max_ty =  int(min(self.max_y_coord, n_y + self.coverage_radius) // self.target_spacing)
    
    covered_target_grid_index = []
    for i in xrange(n_min_tx, n_max_tx+1):
      for j in xrange(n_min_ty, n_max_ty+1):
        '''
        check for rigorous coverage for targets lying at the extremities
        '''
        if i in [n_min_tx,n_max_tx] or j in [n_min_ty,n_max_ty]:
          tx = (i + 0.5)*self.target_spacing
          ty = (j + 0.5)*self.target_spacing
          if self.isWithinDistance(n_x, n_y, tx, ty, self.coverage_radius):
            covered_target_grid_index.append((i,j))
        else: #otherwise add it to the list anyway
          covered_target_grid_index.append((i,j))
    return covered_target_grid_index
  
  def isCoveredByNode(self,n, ti,tj):
    '''
    Args:
      n: node index, 0<= n < self.number_of_nodes
      tx: target grid_x_Coordinate
      ty: target grid_y_Coordinate
    Returns:
      if the node n covers the target positioned at grid (ti,tj) then True
      otherwise, False
      
    Processing:
      i) find the x,y coordinate of node n
      ii) find the x,y cooridinate of target at grid (ti,tj)
          here the x<- (ti+0.5)*self.target_spacing, similarly for y
      iii) if the above two points are within range of self.coverage_radius, then True
    '''
    n_x = self.node_x[n]
    n_y = self.node_y[n]
    t_x =  (ti+0.5)*self.target_spacing
    t_y =  (tj+0.5)*self.target_spacing
    return self.isWithinDistance(n_x, n_y, t_x, t_y, self.coverage_radius) 
  
  def buildNodeHeap(self):
    '''
      i) initialize covered_target , node_target_count and node_hp
      
      ii) for each node n:
        finds the list of targets covered by n and hence the total number of targets n_t
        push it into the node heap node_hp with priority -n_t
        
    '''
    self.node_target_count = []
    for n in xrange(self.number_of_nodes):
      self.node_target_count.append(0)
    
    #following is a 2D array for boolean values
    self.covered_target = []
    for i in xrange(self.max_tx+1):
      self.covered_target.append([])
      for j in xrange(self.max_ty+1):
        self.covered_target[i].append(False)
      
    self.node_hp = MyHeap()
    for n in xrange(self.number_of_nodes):
      number_of_targets_covered = 0
      target_grid_coord_list = self.getTargetsCoveredByNode(n) 
        #target_grid_coord is a list of tuples  i.e. [(t1_x,t1_y),(t2_x,t2_y),...]
      for target_grid_coord in target_grid_coord_list:
        ti,tj = target_grid_coord #target_grid_coord is a tuple i.e. (tx,ty)
        if self.isCoveredByNode(n, ti, tj):
          number_of_targets_covered += 1
      self.node_target_count[n] = number_of_targets_covered
      if self.node_target_count[n]>0:
        self.node_hp.push(n, -self.node_target_count[n])
    return
  
  def heuristicTargetCover(self):
    '''
    pre-condition: must be called after calling buildNodeHeap() method 
    processing:
      First add all the gateway nodes to node_cover and update the affected node keys accordingly
      while either node heap node_hp is not empty AND all the targets are not covered:
          i) pop node n form heap, if gateway node, then continue 
          (gateways have been handled previously)
          ii) find the targets covered by this node n calling getTargetsCoveredByNode() method
          iii) find all the nearby nodes binned by node_bin
          iv) for each target, 
            a) update its covering status i.e self.covered_target
            b) update the self.node_target_count for neighboring nodes
            c) push back those nodes into the heap essentially overwriting the old priority
            d) update self.no_of_targets_covered, 
                  based on the number of new targets covered by this node n
            e) and finally, if at least one new target is covered, add that node to the node_cover
    '''
    #STEP: add only gateway nodes-------------------------------------------------------
    '''
    processing:
      for all gateway node n:
      add n to covered_node
       find targets covered by n
         for each target t:
            if t is yet uncovered:
               change it to covered
               find all nodes covering t, decrement their hash-key
    '''
    for n in self.gateways:
      self.node_cover.append(n)
      n_x = self.node_x[n]
      n_y = self.node_y[n]
      #get the list of targets
      target_grid_coord_list = self.getTargetsCoveredByNode(n)
      
      #get potential nodes covering targets in the above list
      nbr_nodes = self.node_bin.get(n_x, n_y, self.coverage_radius)
      if n in nbr_nodes:
        nbr_nodes.remove(n)
        
      affected_nbr_nodes = set()
      
      for target_grid_coord in target_grid_coord_list:
        ti,tj = target_grid_coord #target_grid_coord is a tuple i.e. (tx,ty)
        if not self.covered_target[ti][tj]:
          #if gateway is covering this target, 
            #update the target's status and count of total covered targets
          self.covered_target[ti][tj] = True
          self.no_of_targets_covered += 1
          for nbr in nbr_nodes:
            if self.isCoveredByNode(nbr, ti, tj):
              self.node_target_count[nbr] -= 1
              affected_nbr_nodes.add(nbr)
              
      while affected_nbr_nodes:
        pushback_node = affected_nbr_nodes.pop()
        if self.node_target_count[pushback_node]>0: 
                #no need to push those node that does not cover any new targets
          self.node_hp.push(pushback_node, -self.node_target_count[pushback_node])  
        
    self.logger.debug("completed processing targets covered by gateway nodes") 
    #STEP: add only non-gateway nodes-------------------------------------------------------
    while not self.node_hp.isHeapEmpty() and self.no_of_targets_covered < self.total_targets:
      n = self.node_hp.pop()
      if n is None: #safety
        break
      if n in self.gateways: #IMP: gateways have been already handled, no need to consider them
        continue
      total_new_targets_covered = 0
      
      n_x = self.node_x[n]
      n_y = self.node_y[n]
      
      nbr_nodes = self.node_bin.get(n_x, n_y, self.coverage_radius)
      if n in nbr_nodes:
        nbr_nodes.remove(n)
        
      affected_nbr_nodes = set()
      
      target_grid_coord_list = self.getTargetsCoveredByNode(n) 
        #target_grid_coord is a list of tuples  i.e. [(t1_x,t1_y),(t2_x,t2_y),...]
        
      for target_grid_coord in target_grid_coord_list:
        ti,tj = target_grid_coord #target_grid_coord is a tuple i.e. (tx,ty)
        if not self.covered_target[ti][tj]:
          #task (a)
          self.covered_target[ti][tj] = True
          total_new_targets_covered += 1
          #task (b)
          for nbr in nbr_nodes:
            if self.isCoveredByNode(nbr, ti, tj):
              self.node_target_count[nbr] -= 1
              affected_nbr_nodes.add(nbr)
              
      #task (c)
      while affected_nbr_nodes:
        pushback_node = affected_nbr_nodes.pop()
        if self.node_target_count[pushback_node]>0: 
                #no need to push those node that does not cover any new targets
          self.node_hp.push(pushback_node, -self.node_target_count[pushback_node])  
        
      #task (d)+(e)  
      if total_new_targets_covered>0:
        self.no_of_targets_covered += total_new_targets_covered
        self.node_cover.append(n)
        
  def visualizeGraph(self, g, show_edges = True):
    '''
    over-ridden to show nodes in self.node_cover in green
    pre-condition: must be called after self.adj and self.gateways and self.node_cover
               are set or initialized
    visualize any graph that is subgraph of self.adj using pyplot
    Args: g, must be a networkx.Graph object and must be a subgraph of self.adj
    Returns: None
    '''
    #build node positions and node colors:
    node_positions = {}
    node_colors = []
    for n in g.nodes():
      node_positions[n] = (self.node_x[n], self.node_y[n])
      if n in self.gateways and n in self.node_cover:
        node_colors.append('b')
      elif n in self.node_cover:
        node_colors.append('g')
      elif n in self.gateways:
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
    
  
  
  
  
  
  