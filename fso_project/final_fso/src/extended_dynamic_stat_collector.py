from stat_collector import StatCollector
import networkx as nx
import random
import networkx.algorithms.flow as flow
import csv
from matplotlib import pyplot as plt
import logging

class ExtendedDynamicStatCollector(StatCollector):
  def __init__(self, configFile):
    StatCollector.__init__(self, configFile)
    self.weigh_increment = 5
    self.current_run_no = 0
    self.static_graph_spec_file_path = None
    self.onlyGeneratePlot = False
    
    self.static_fso_per_node = None
    self.static_fso_per_gateway = None
    self.static_eq_dynamic = None
    self.static_eq_avg_max_flow = None
    self.static_eq_avg_upperbound_flow = None
    
    self.stat_header.append('Static_eq_fso_per_node')
    self.stat_header.append('Static_eq_fso_per_gateway')
    self.stat_header.append('Static_eq_avg_max_flow')
    self.stat_header.append('Static_eq_avg_upperbound_flow')
    
  def reset(self):
    StatCollector.reset(self)
    self.weigh_increment = 5
    self.current_run_no = 0
    self.static_graph_spec_file_path = None
    
    self.static_fso_per_node = None
    self.static_fso_per_gateway = None
    self.static_eq_dynamic = None
    self.static_eq_avg_max_flow = None
    self.static_eq_avg_upperbound_flow = None
    
    
  def createStaticGraphSpecOutputFile(self):
    '''
    save the static graph according to TM's spec 
    must set self.static_graph_spec_file_path correctly before calling
    '''
    
  
    with open(self.static_graph_spec_file_path, 'w') as f:
      f.write("trafficSources:\n")
      for n in self.static_graph.nodes():
        f.write(str(n)+"\n")
      f.write('\n')
      
      
      #---FSONodes:
      fso_counter ={}
      gateway_degree_in_static_graph = self.static_graph.degree()
      static_nodes = self.static_graph.nodes() 
      
      f.write("FSONodes:\n")
      for n in static_nodes:
        for fso in range(1, gateway_degree_in_static_graph[n]+1):
          f.write(str(n)+"_fso"+str(fso)+"\n")
          fso_counter[n]=0
      f.write('\n')   
      
      #---FSOLinks:---
      f.write("FSOLinks:\n")   
      edges = self.static_graph.edges()
      for u,v in edges:
        fso_counter[u] += 1
        fso_counter[v] += 1
        fso_u =  fso_counter[u]
        fso_v = fso_counter[v]
        
        f_text = str(u)+"_fso"+str(fso_u)+"To"\
                    +str(v)+"_fso"+str(fso_v)+" "+str(self.link_capacity)+"Mbps\n" 
        f.write(f_text)
        f_text = str(v)+"_fso"+str(fso_v)+"To"\
                    +str(u)+"_fso"+str(fso_u)+" "+str(self.link_capacity)+"Mbps\n" 
        f.write(f_text)
      f.write('\n') 
  
      #------gateways--:
      static_gateways =  list(set(self.gateways) & set(static_nodes))
      f.write('gateways:\n')
      for n in static_gateways:
        f.write(str(n)+"\n")
  
  def getEquivalendDynamicGraph(self, in_patternNodes):
    '''
    Args: patternNodes: list of non-gateway source nodes
    Returns: an equivalent dynamic graph
    Processing:
    intialize equiDynamic graph as empty
    initalize curDynamic graph as directed dynamic graph
    phase 1 (connectivity phase):
    -------------------------------
    set-up:
      i. add super-sink and connect it to all the gateway nodes in curDynamic graph
      ii. add 1 as edge weight to all edges
      iii. initialize remaining_in_degree and remaning_out_degree 
          counter for each node to max. permitted degree
        (as determined by fso_per_node or fso_per_gateway params)
      
    1. while there is remaining nodes in patternNodes:
    2. find Dijsktra paths to all nodes taking super-sink as the source
    
    3. add the min. cost path (if any in curDynamic graph) to one of the pattern nodes n
       (both edges and incident nodes on the path except super-sink) to equiDynamic graph, 
       then increment weight of each edge on this path by given param c in curDynamic graph,
       decrement the remaining_degree counter of the nodes on this path, 
        and freeze** on curDynamic graph the edges of all nodes on this path having 0 value
         in remaining_degree counter, 
        also remove from patternNode any node existing on this path 
        (this node thus got a path to gatways already) and continue to step 1
        
    4. if no such path found, find dijsktra paths taking current pattern node n as the source
        and find the min-cost path (if any in curDynamic graph) to one of 
        the pattern nodes already added to 
        equiDynamic graph (thus having a path to the gateways), then add the path to equiDynamic
        graph, increment the weights of the edges on the path in curDynamic graph by given param c,
        decrement the remaining_degree counter of the nodes on this path, 
        and freeze** on curDynamic graph the edges of all nodes on this path having 0 value
         in remaining_degree counter, 
        also remove from patternNode any node existing on this path 
        (this node thus got path already)
        
    5. if a full scan of the patternNodes list was made but no node could be connected to gateways,
    then stop the iteration, the remaining nodes are un-connectable nodes
    
    phase 2 (flow enhancement phase):
    ------------------------------------
    1. add a super-source connecting all the nodes in patterNodes in curDynamic graph (as it is
      with the weights of the edges accumulated in the previous phase)
      
    2. compute Dijsktra path from super-source to super-sink in curDynamic graph, 
      then add the min-cost path (if any) to equiDynamic graph,
      then increment weight of each edge on this path by given param c in curDynamic graph
       and freeze** on curDynamic graph the edges of all nodes on this path having 0 value
         in remaining_degree counter,
       
    3. if no such path found in step 2, then break else continue to step 2
    
    Note: freeze**: if remaining_in_degree is zero for a node n, then remove all other 
        in-coming edges of n in curDynamic graph except those existing in equiDynamic graph
        similarly for remaining_out_degree and out going edges...
        
    Note: edges to suepr-sink and super-source do not contribute to degrees of the nodes
    return the equiDynamic graph
    
    '''
    #initialize phase
    patternNodes = list(in_patternNodes)
    cur_dynamic_graph = nx.DiGraph()
    for u,v in self.dynamic_graph.edges():
      cur_dynamic_graph.add_edge(u,v, weight = 1)
      cur_dynamic_graph.add_edge(v,u, weight = 1)
      
    cur_dynamic_graph.add_node('snk')
    for g in self.gateways:
      cur_dynamic_graph.add_edge(g,'snk',weight = 1)
    
    cur_dynamic_graph.add_node('src')
    for n in patternNodes:
      cur_dynamic_graph.add_edge('src',n,weight = 1)
      
    eq_dynamic_graph = nx.DiGraph()
    
    remaining_in_degree = []
    remaining_out_degree = []
    
    for i in xrange(self.number_of_nodes):
      remaining_in_degree.append(0)
      remaining_out_degree.append(0)
      
    for n in self.dynamic_graph.nodes():
      if n in self.gateways:
        remaining_in_degree[n] = self.fso_per_gateway
        remaining_out_degree[n] = self.fso_per_gateway
      else:
        remaining_in_degree[n] = self.fso_per_node
        remaining_out_degree[n] = self.fso_per_node
    
    #Phase 1   
    while True:
      path = []
      try:
        path = nx.dijkstra_path(cur_dynamic_graph, 'src', 'snk', weight='weight')
      except nx.NetworkXNoPath:
        break
      
      # add path to equiDynamic graph
      eq_dynamic_graph.add_path(path[1:-1])
      self.logger.debug("adding path to eq_dynamic:"+str(path[1:-1]))
      # update degree of each node on this path and freeze degree saturated node, also
      # update the weigh on edges of this path
      u = path[1]
      for v in path[2:-1]:
        cur_dynamic_graph[u][v]['weight'] += self.weigh_increment
        if remaining_out_degree[u]>=1:
          remaining_out_degree[u] -= 1
          self.logger.debug("\t u,remaining_out_degree[u]:"+str(u)+","+str(remaining_out_degree[u]))
          if remaining_out_degree[u] == 0:
            self.logger.debug("freezing(out) node:"+str(u))
            for x,y in cur_dynamic_graph.out_edges(u):
              if not eq_dynamic_graph.has_edge(x, y):
                cur_dynamic_graph.remove_edge(x, y)
                self.logger.debug("\t removing edge:"+str(x)+"-"+str(y))
        
        if remaining_in_degree[v]>=1: 
          remaining_in_degree[v] -= 1
          self.logger.debug("\t v,remaining_in_degree[v]:"+str(v)+","+str(remaining_in_degree[v]))
          if remaining_in_degree[v] == 0:
            self.logger.debug("freezing(in) node:"+str(u))
            for x,y in cur_dynamic_graph.in_edges(v):
              if not eq_dynamic_graph.has_edge(x,y):
                cur_dynamic_graph.remove_edge(x,y)
                self.logger.debug("\t removing edge:"+str(x)+"-"+str(y))
        u=v
        
      # remove nodefrom patternNode
      patternNodes.remove(path[1])
      
      #disconenct the edge from src to this pattern node in cur_dynamic_graph
      cur_dynamic_graph.remove_edge('src', path[1])
      
    #phase 2
    #add back the src-to-pattern node edges
    for n in in_patternNodes:
      cur_dynamic_graph.add_edge('src',n,weight = 1)
      
    empty_iteration_count = 0
    
    while empty_iteration_count<100:
      self.logger.debug("empty_iteration_count:"+str(empty_iteration_count))
      added_new_edge = False
      path = []
      try:
        path = nx.dijkstra_path(cur_dynamic_graph, 'src', 'snk', weight='weight')
      except nx.NetworkXNoPath:
        break
      
      # update degree of each node on this path and freeze degree saturated node, also
      # update the weigh on edges of this path
      u = path[1]
      for v in path[2:-1]:
        cur_dynamic_graph[u][v]['weight'] += self.weigh_increment
        if not eq_dynamic_graph.has_edge(u,v):
          added_new_edge = True
        
        if remaining_out_degree[u]>=1:
          remaining_out_degree[u] -= 1
          self.logger.debug("\t u,remaining_out_degree[u]:"+str(u)+","+str(remaining_out_degree[u]))
          if remaining_out_degree[u] == 0:
            self.logger.debug("freezing(out) node:"+str(u))
            for x,y in cur_dynamic_graph.out_edges(u):
              if not eq_dynamic_graph.has_edge(x, y):
                cur_dynamic_graph.remove_edge(x, y)
                self.logger.debug("\t removing edge:"+str(x)+"-"+str(y))
        
        if remaining_in_degree[v]>=1: 
          remaining_in_degree[v] -= 1
          self.logger.debug("\t v,remaining_in_degree[v]:"+str(v)+","+str(remaining_in_degree[v]))
          if remaining_in_degree[v] == 0:
            self.logger.debug("freezing(in) node:"+str(u))
            for x,y in cur_dynamic_graph.in_edges(v):
              if not eq_dynamic_graph.has_edge(x,y):
                cur_dynamic_graph.remove_edge(x,y)
                self.logger.debug("\t removing edge:"+str(x)+"-"+str(y))
        u=v
            # add path to equiDynamic graph
      if added_new_edge:
        eq_dynamic_graph.add_path(path[1:-1]) 
        self.logger.debug("phase 2: adding path to eq_dynamic:"+str(path[1:-1]))
      #instead of disconnecting, increment the src-node weight
      cur_dynamic_graph['src'][path[1]]['weight'] += self.weigh_increment
      if added_new_edge:
        empty_iteration_count = 0
      else:
        empty_iteration_count += 1
      
    cur_dynamic_graph.remove_node('src')
    cur_dynamic_graph.remove_node('snk')
    return eq_dynamic_graph, cur_dynamic_graph, patternNodes
  
     
      
      
  def computeExtDynamicAvgFlow(self):
    '''
    runs this stat colletor
    processing:
      for each iteration 
      i) select pattern nodes
      ii) find an equivalent dynamic graph, 
          if its the first iteration, save it as equivalen static graph 
      iii) run max flow and save it
      for all iterations, find the avg max flow using the saved values
    '''
    
    max_flow_list_dyn = []
    upperbound_flow_list_dyn = []
    
    candidate_pattern_nodes = list(set(self.dynamic_graph.nodes()) - set(self.gateways))
    if not candidate_pattern_nodes: #no non-gateway nodes in the static graph, so this avg finding is mute
      self.ext_dynamic_avg_flow = 0.0
      self.ext_dynamic_upperbound_flow = 0.0
      self.logger.info("No non-gateway node found in dynamic graph, avg and upperbound max flows set to 0!!")
      return
    number_of_candidate_pattern_nodes = len(candidate_pattern_nodes)
    number_of_pattern_nodes = \
      int(round(1.0*number_of_candidate_pattern_nodes * self.percent_of_pattern_nodes_in_avg_flow_calculation/100.0, 0))
    if number_of_pattern_nodes<1: #check if the ratio was too low
      number_of_pattern_nodes = 1
    
    number_of_iterations = self.number_of_pattern_in_avg_flow_calculation
    
    if number_of_pattern_nodes>= number_of_candidate_pattern_nodes:
      self.logger.debug("number_of_pattern_nodes>= number_of_candidate_pattern_nodes!!")
      number_of_iterations =  1
      number_of_pattern_nodes = number_of_candidate_pattern_nodes
    
    for i in xrange(number_of_iterations):
      pattern_nodes = random.sample(candidate_pattern_nodes, number_of_pattern_nodes)
      #self.logger.debug("pattern nodes:"+str(pattern_nodes))
      eq_dyamic_graph, self.cur_dynamic_graph, remaining_pattern_nodes = \
            self.getEquivalendDynamicGraph(pattern_nodes)
      #self.logger.debug("un-connectable pattern nodes:"+str(remaining_pattern_nodes))
      
        
      for u,v in eq_dyamic_graph.edges():
        eq_dyamic_graph[u][v]['capacity'] =  1.0

      eq_dyamic_graph.add_node('src')
      for n in pattern_nodes:
        eq_dyamic_graph.add_edge('src',n, capacity = float('inf'))
      
      eq_dyamic_graph.add_node('snk')
      for g in self.gateways:
        eq_dyamic_graph.add_edge(g, 'snk', capacity = float('inf'))
        
      residual_max_flow_dyn = flow.shortest_augmenting_path(G=eq_dyamic_graph, s='src', t='snk')
      max_flow_list_dyn.append(residual_max_flow_dyn.graph['flow_value'])

      pattern_nodes_included_in_eq_dynamic_graph =\
              list(set(pattern_nodes) - set(remaining_pattern_nodes))
      pattern_node_degrees = self.dynamic_graph.degree(pattern_nodes_included_in_eq_dynamic_graph)
      
      upper_bound_flow_dyn = 0
      for n,d in pattern_node_degrees.iteritems():
        upper_bound_flow_dyn += min(self.fso_per_node, d)

      upperbound_flow_list_dyn.append(upper_bound_flow_dyn)
      
    self.ext_dynamic_avg_flow =\
        1.0*self.link_capacity*sum(max_flow_list_dyn)/len(max_flow_list_dyn)
    self.ext_dynamic_upperbound_flow =\
        1.0*self.link_capacity*sum(upperbound_flow_list_dyn)/len(upperbound_flow_list_dyn)
  
  def saveStatInFile(self): #overriding parent method
    '''
    overriding the method in stat_collector,
    everything same except now saving shaifur's avg. dyn flow and upper. dyn flow insteadm of TM's
    save the statistics in a file
    the order of the fields are very important
    ''' 
    stat_row={}
    stat_row['number_of_nodes_in_input_graph'] = self.adj.number_of_nodes()
    stat_row['number_of_edges_in_input_graph'] = self.adj.number_of_edges()
    stat_row['number_of_gateways'] = len(self.gateways)
    stat_row['number_of_nodes_in_static_graph'] = self.static_graph.number_of_nodes()
    stat_row['number_of_nodes_in_dynamic_graph'] = self.dynamic_graph.number_of_nodes()
    stat_row['number_of_fso_per_node'] = self.fso_per_node
    stat_row['number_of_fso_per_gateway'] = self.fso_per_gateway
    stat_row['percent_of_pattern_nodes'] = self.percent_of_pattern_nodes_in_avg_flow_calculation
    stat_row['number_of_patterns'] = self.number_of_pattern_in_avg_flow_calculation
    
    stat_row['statc_upperbound_max_flow'] = self.static_upperbound_flow
    stat_row['total_gateway_capacity_static'] = self.total_gateway_capacity_static
    stat_row['static_avg_max_flow'] = self.static_avg_flow
    
    #the following stats are changed now
    stat_row['dynamic_upperbound_max_flow'] = self.ext_dynamic_upperbound_flow
    stat_row['total_gateway_capacity_dynamic'] = self.total_gateway_capacity_dynamic
    stat_row['dynamic_avg_max_flow'] = self.ext_dynamic_avg_flow
    
    #the following were added for computation of eq static graph
    stat_row['Static_eq_fso_per_node'] = self.static_fso_per_node
    #self.logger.info("value-check inside SaveStatInFile:self.static_fso_per_gateway"+str(self.static_fso_per_gateway))
    stat_row['Static_eq_fso_per_gateway'] = self.static_fso_per_gateway
    stat_row['Static_eq_avg_max_flow'] = self.static_eq_avg_max_flow
    stat_row['Static_eq_avg_upperbound_flow'] = self.static_eq_avg_upperbound_flow
    
    f= open(self.output_statistics_file, 'ab')
    writer = csv.DictWriter(f, self.stat_header)
    writer.writerow(stat_row)
    f.close()
    if self.printSummaryAfterSavingStat:
      print"--------------Stat-Summary--------------"
      print "Experiment Name:",self.experiment_name," run#:",self.current_run_no
      print "Number of Nodes in Input Graph:",stat_row['number_of_nodes_in_input_graph']
      print "Number of Edges in Input Graph:",stat_row['number_of_edges_in_input_graph']
      print "Number of Gateways in Input Graph:",stat_row['number_of_gateways']
      print "Number of Nodes in Static Graph:",stat_row['number_of_nodes_in_static_graph']
      print "Number of Nodes in Dynamic Graph:",stat_row['number_of_nodes_in_dynamic_graph']
      print "Number of FSO-per-node:",stat_row['number_of_fso_per_node']
      print "Number of FSO-per-gateway:",stat_row['number_of_fso_per_gateway']
      print "Percent of total nodes used in patterns for avg. max. flow calculation:",stat_row['percent_of_pattern_nodes']
      print "Number of patterns used in avg. max. flow calculation:",stat_row['number_of_patterns']
      print "Static:.................."
      print "Static upperbound max. flow:",stat_row['statc_upperbound_max_flow']
      print "Static total gateway capacity:",stat_row['total_gateway_capacity_static']
      print "Static avg. max. flow:",stat_row['static_avg_max_flow'], "!!"
      
      print "Ratio of avg. max. flow to upperbound max. flow for Static graph:",
      if stat_row['statc_upperbound_max_flow']>0:
        print 1.0*stat_row['static_avg_max_flow']/stat_row['statc_upperbound_max_flow']
      else:
        print "N/A"
        
      print "Ratio of avg. max. flow to total gateway capacity for Static graph:",
      if stat_row['total_gateway_capacity_static']>0:
        print 1.0*stat_row['static_avg_max_flow']/stat_row['total_gateway_capacity_static']
      else: 
        print "N/A"
        
      print "Dynamic:.................."
      print "Dynamic upperbound max. flow:",stat_row['dynamic_upperbound_max_flow']
      print "Dynamic total gateway capacity:",stat_row['total_gateway_capacity_dynamic']
      print "Dynamic avg. max. flow:",stat_row['dynamic_avg_max_flow'], "!!"
      
      print "Ratio of avg. max. flow to upperbound max. flow for Dynamic graph:",
      if stat_row['dynamic_upperbound_max_flow']>0:
        print 1.0*stat_row['dynamic_avg_max_flow']/stat_row['dynamic_upperbound_max_flow']
      else:
        print "N/A"
      
      print "Ratio of avg. max. flow to total gateway capacity for Dynamic graph:",
      if stat_row['total_gateway_capacity_dynamic']>0:
        print  1.0*stat_row['dynamic_avg_max_flow']/stat_row['total_gateway_capacity_dynamic']
      else:
        print 'N/A'
      
      print "Static equivalent fso-per-node:",stat_row['Static_eq_fso_per_node'] 
      print "Static equivalent fso-per-gateway:", stat_row['Static_eq_fso_per_gateway']
      print "Static equivalent avg max. flow:",stat_row['Static_eq_avg_max_flow']
      print "Static equivalent avg upperboud flow:",stat_row['Static_eq_avg_upperbound_flow'] 
      print "==============End of stat-summary==========="
  
  def runExtendedDynamicStatCollector(self):
    self.logger.info("in runExtendedDynamicStatCollector method...")
    
    self.getStaticAvgFlow()
    self.setTotalGatewayCapacities()
    #instead of java folder dump this file into graph output folder
    self.dynamic_graph_spec_file_path =\
     self.graph_output_folder+"/"+str(self.experiment_name)+"-"+str(self.current_run_no)+"_dyn.txt"
    
    self.static_graph_spec_file_path =\
     self.graph_output_folder+"/"+str(self.experiment_name)+"-"+str(self.current_run_no)+"_stat.txt" 
     
    self.createDynamicGraphSpecOutputFile()
    self.createStaticGraphSpecOutputFile()
    
    #------save adj graph--------------------------
    filePath = \
      self.graph_output_folder+"/"+str(self.experiment_name)+"-"+str(self.current_run_no)+"_adj.txt"
    self.saveGraphInFile(self.adj, filePath)
    
    #------save backbone graph--------------------------
    filePath = \
      self.graph_output_folder+"/"+str(self.experiment_name)+"-"+str(self.current_run_no)+"_backbone.txt"
    self.saveGraphInFile(self.backbone_graph, filePath)
    
    #------save static graph--------------------------
    filePath = \
      self.graph_output_folder+"/"+str(self.experiment_name)+"-"+str(self.current_run_no)+"_static.txt"
    self.saveGraphInFile(self.static_graph, filePath)
    #------save dynamic graph--------------------------
    filePath = \
      self.graph_output_folder+"/"+str(self.experiment_name)+"-"+str(self.current_run_no)+"_dynamic.txt"
    self.saveGraphInFile(self.dynamic_graph, filePath)
    
    self.computeExtDynamicAvgFlow()
    self.findFlowEquivalentStaticGraph()
    self.saveStatInFile()
    
  
  def getAverageAndUpperboundFlow(self, g, no_of_pattern_nodes, no_of_iterations):
    '''
    returns the avg and upperbound flow of arbitrary graph g
    follows the spec of number of patterns and iteratios specified in th
          
    '''  
    avg_max_flow_val = 0.0
    avg_upperbound_flow_val = 0.0

    candidate_pattern_nodes = list( set(g.nodes()) - set(self.gateways))
    if not candidate_pattern_nodes: 
              #no non-gateway nodes in the static graph, so this avg finding is mute
      return avg_max_flow_val, avg_upperbound_flow_val
    
    number_of_candidate_pattern_nodes = len(candidate_pattern_nodes)
    number_of_pattern_nodes = \
      int(round(1.0*number_of_candidate_pattern_nodes * no_of_pattern_nodes/100.0, 0))
      
    if number_of_pattern_nodes<1: #check if the ratio was too low
      number_of_pattern_nodes = 1
    
    number_of_iterations = no_of_iterations
    
    if number_of_pattern_nodes>= number_of_candidate_pattern_nodes:
      number_of_iterations =  1
      number_of_pattern_nodes = number_of_candidate_pattern_nodes
      
    max_flow_vals = []
    upperbound_flow_vals = []
    
    for i in xrange(number_of_iterations):
      pattern_nodes = random.sample(candidate_pattern_nodes, number_of_pattern_nodes)
      max_flow_val = self.getMaxFlowForSourceNodes(g, pattern_nodes)
      pattern_node_degrees =  g.degree(pattern_nodes)
      upper_bound_flow = sum(pattern_node_degrees.values()) 
      
      max_flow_vals.append(max_flow_val)
      upperbound_flow_vals.append(upper_bound_flow)
      
    avg_max_flow_val = 1.0*self.link_capacity*sum(max_flow_vals)/len(max_flow_vals)
    avg_upperbound_flow_val = 1.0*self.link_capacity*sum(upperbound_flow_vals)/len(upperbound_flow_vals)
    if avg_upperbound_flow_val is None:
      avg_upperbound_flow_val = 0.0
    return avg_max_flow_val, avg_upperbound_flow_val
  
  
  
  
  
  def findFlowEquivalentStaticGraph(self):
    '''
     increase the fso-per-node and fso-per-gateway by an increment step
      upto some max, until the avg flow value for the static graph increases
    Processing:
      i) increase the fso-per-node and fso-per-gateway
      ii) add more edges to the static graph between the existing nodes
      iii) find the avg. max flow for the static graph
      iv) if it just exceeds that of dynamic max flow stop and save the current
      fso-per-node fso-per-gateway
      v) else continue to step i if no. of iterations<=max_allowed or new edges were added
      in the last step
    '''
    #initialize the variables
    self.static_eq_dynamic = nx.Graph(self.static_graph)
    self.static_eq_avg_max_flow = self.static_avg_flow
    self.static_eq_avg_upperbound_flow = self.static_upperbound_flow
    self.static_fso_per_node = self.fso_per_node
    #self.logger.info("value check: fso_per_gatway"+str(self.fso_per_gateway))
    self.static_fso_per_gateway = self.fso_per_gateway
    #self.logger.info("value check: static_fso_per_gatway"+str(self.static_fso_per_gateway))
    
    edge_added_in_prev_iterations = True
    debug_iter_counter = 0
    while edge_added_in_prev_iterations and (self.static_eq_avg_max_flow < self.ext_dynamic_avg_flow):
      edge_added_in_prev_iterations = False
      debug_iter_counter += 1
      #task (i)
      self.static_fso_per_gateway += 1
      #self.logger.info("value check: static_fso_per_gatway"+str(self.static_fso_per_gateway))
      self.static_fso_per_node += 1
      #task (ii)
      max_allowed_degree = {}
      current_degree = self.static_eq_dynamic.degree()
      for n in self.static_eq_dynamic.nodes():
        max_allowed_degree[n] = self.static_fso_per_node
        if n in self.gateways:
          max_allowed_degree[n] = self.static_fso_per_gateway      
      for u in self.static_eq_dynamic.nodes():
        for v in self.static_eq_dynamic.nodes():
          if current_degree[u]>=max_allowed_degree[u]:
            break
          if current_degree[v]>=max_allowed_degree[v]:
            continue
          #else, explore the possibility to add edge u-v
          if not self.static_eq_dynamic.has_edge(u, v) and self.adj.has_edge(u, v):
            self.static_eq_dynamic.add_edge(u,v)
            self.logger.debug("Added edge (u,v) to static_eq_dynamic graph:"+str(u)+","+str(v))
            current_degree[u] += 1
            current_degree[v] += 1
            edge_added_in_prev_iterations = True
      
      #task (iii)
      #self.logger.info("value check inside findFlowEquivalentStaticGraph():self.static_fso_per_gateway:"+str(self.static_fso_per_gateway))
      if edge_added_in_prev_iterations:
        self.static_eq_avg_max_flow, self.static_eq_avg_upperbound_flow = \
          self.getAverageAndUpperboundFlow(self.static_eq_dynamic, 
                                           self.percent_of_pattern_nodes_in_avg_flow_calculation, 
                                           self.number_of_pattern_in_avg_flow_calculation)
  def plotFlowStat(self, in_plot_keys, ylabel): 
    '''
    opens the current self.output_statistics_file and plots the required plotType
    processing:
      i) opens the file and stores the columns
      ii) make the plots for the required files
      iii) show the plot as output or save it in a file
    '''
    x=[]
    y={}
    plot_keys = list( set(self.stat_header) & set(in_plot_keys) )
    
    for k in plot_keys:
        y[k]= []
        
    max_flow_val = 0
    min_flow_val = float('inf')
    
    y['number_of_fso_per_node']=[]   
    y['Static_eq_fso_per_node']=[]
    counter = 0
    with open(self.output_statistics_file) as csvfile:
      reader = csv.DictReader(csvfile, fieldnames=self.stat_header)
      for row in reader:
        counter += 1
        x.append(counter)
        for k in plot_keys:
          y[k].append(float(row[k]))
          max_flow_val = max(max_flow_val, float(row[k]))
          min_flow_val = min(min_flow_val, float(row[k]))
        y['number_of_fso_per_node'].append(int(row['number_of_fso_per_node']))  
        y['Static_eq_fso_per_node'].append(int(row['Static_eq_fso_per_node']))
    #fig, ax = plt.subplots()
    for k in plot_keys:
      #print "x:", x
      #print "y[k]:", y[k]
      ls = '-'
      if k.find('upper') >-1:
        ls = '--'
      plt.plot(x,y[k], label= k, ls= ls, lw=1.5, marker = 's')
    
    legend = plt.legend(loc='upper right',bbox_to_anchor=[0.5, 1.2])
    plt.grid(True)
    plt.xlabel('Sample No')
    plt.ylabel(ylabel)
    plt.xticks(xrange(0,len(x)+2))
    plt.yticks(xrange(int(round(min_flow_val,-3))-2000,int(round(max_flow_val,-3))+ 2000,10000))
    plt.show()
    
  
  def plotFSOStat(self): 
    '''
    opens the current self.output_statistics_file and plots the required plotType
    processing:
      i) opens the file and stores the columns
      ii) make the plots for the required files
      iii) show the plot as output or save it in a file
    '''
    x=[]
    y={}
    
    y['number_of_fso_per_node']=[]   
    y['Static_eq_fso_per_node']=[]
    
    counter = 0
    with open(self.output_statistics_file) as csvfile:
      reader = csv.DictReader(csvfile, fieldnames=self.stat_header)
      for row in reader:
        counter += 1
        x.append(counter)
        y['number_of_fso_per_node'].append(int(row['number_of_fso_per_node']))  
        y['Static_eq_fso_per_node'].append(int(row['Static_eq_fso_per_node']))

    plt.plot(x, y['number_of_fso_per_node'], label= 'Dynamic FSO-per-node', ls='-', lw=1.5, marker = 's')
    plt.plot(x, y['Static_eq_fso_per_node'], label = 'Static FSO-per-node', ls='-', lw=1.5, marker = 's')
    
    legend = plt.legend(loc='upper left')
    plt.grid(True)
    plt.xlabel('Sample No')
    plt.ylabel('No. of FSO per node')
    plt.xticks(xrange(0,len(x)+2))
    plt.yticks(xrange(0, max(y['Static_eq_fso_per_node'])+2))
    plt.show()
   
  def runAll(self):
    static_avg_max_flows = []
    static_avg_upperbound_flows = []
    dynamic_avg_max_flows = []
    dynamic_avg_upperbound_flows = []
    static_eq_fso_per_node = []
    
    if not self.onlyGeneratePlot:
      for i in xrange(1, self.no_of_samples+1):
        self.current_run_no = i
        self.fso_per_node = self.backup_fso_per_node 
        self.fso_per_gateway = self.backup_fso_per_gateway 
        self.runInputGenerator()
        self.runStep_1()
        self.runStep_2()
        self.runStep_3()
        self.runStep_4_dynamic()
        self.runStep_4_static()
        self.runILPSolver()
        self.printSummaryAfterSavingStat = True #will now print on screen summary
        self.runExtendedDynamicStatCollector()
        #backup the following stat in master_stat
        #number_of_fso_per_node, 
#         number_of_edges_per_node, 
#         percent_of_pattern_nodes,
#         static_avg_max_flow, 
#         static_up_flow, 
#         dyn_avg_max_flow, 
#         dyn_up_flow,
#         static_eq_fso_per_node
        static_avg_up_flow = min(self.static_upperbound_flow, self.total_gateway_capacity_static)
        dynamic_avg_up_flow = min(self.ext_dynamic_upperbound_flow, self.total_gateway_capacity_dynamic)
        
        static_avg_max_flows.append(self.static_avg_flow)
        dynamic_avg_max_flows.append(self.ext_dynamic_avg_flow)
        static_avg_upperbound_flows.append(static_avg_up_flow)
        dynamic_avg_upperbound_flows.append(dynamic_avg_up_flow)
        static_eq_fso_per_node.append(self.static_fso_per_node) #TODO:
        self.reset()
        
    aggr_static_avg_max_flow = 1.0*sum(static_avg_max_flows)/len(static_avg_max_flows)
    aggr_dynamic_avg_max_flow = 1.0*sum(dynamic_avg_max_flows)/len(dynamic_avg_max_flows)
    aggr_static_avg_up_flow =  1.0*sum(static_avg_upperbound_flows)/len(static_avg_upperbound_flows)
    aggr_dynamic_avg_up_flow =  1.0*sum(dynamic_avg_upperbound_flows)/len(dynamic_avg_upperbound_flows)
    aggr_static_fso_per_node = 1.0*sum(static_eq_fso_per_node)/len(static_eq_fso_per_node)
    f_text = str(self.fso_per_node)+","\
                +str(self.max_short_edge_per_node)+","\
                +str(self.max_long_edge_per_node)+","\
                +str(self.percent_of_pattern_nodes_in_avg_flow_calculation)+","\
                +str(aggr_static_avg_max_flow)+","\
                +str(aggr_static_avg_up_flow)+","\
                +str(aggr_dynamic_avg_max_flow)+","\
                +str(aggr_dynamic_avg_up_flow)+","\
                +str(aggr_static_fso_per_node)+"\n"
    with open("./output/master_stat.txt","a") as f:
      f.write(f_text)
      
    '''
    plot_keys = [
                 'static_avg_max_flow',
                 'dynamic_avg_max_flow',
                 'statc_upperbound_max_flow',
                 'dynamic_upperbound_max_flow'
                 ]
    self.plotFlowStat( plot_keys, 'Flow Value (Mbps)' )
    self.plotFSOStat()'''
    
    
    
    
    
    
    
    