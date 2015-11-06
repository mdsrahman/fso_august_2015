import subprocess
from ilp_solver import ILPSolver
import random
import networkx as nx
import networkx.algorithms.flow as flow
import csv

class StatCollector(ILPSolver):
  '''
  Responsible for saving all the graphs, generating and saving all the experimental data
  #---descriptions of class-fields----
  
  #-----end of descriptions of class fields-----
  '''
  def __init__(self,configFile):
    ILPSolver.__init__(self, configFile)
    #---class fields-------
    self.dynamic_upperbound_flow = None
    self.dynamic_avg_flow = None
    self.static_upperbound_flow = None
    self.static_avg_flow = None
    self.total_gateway_capacity_static = None
    self.total_gateway_capacity_dynamic = None
    self.path_to_java_code_for_avg_calc = './java/tm.jar'
    self.dynamic_graph_spec_file_path = './java/temp_dynamic_spec.txt'
    self.java_code_stat_file_path = './java/temp_java_stat.txt'
    self.java_code_stdout_file_path = './java/temp_java_stdout.txt'
    self.skipJavaCodeCall = False
    self.stat_header=[
                      'number_of_nodes_in_input_graph',
                      'number_of_edges_in_input_graph',
                      'number_of_gateways',
                      'number_of_nodes_in_static_graph',
                      'number_of_nodes_in_dynamic_graph',
                      'number_of_fso_per_node',
                      'number_of_fso_per_gateway',
                      'percent_of_pattern_nodes',
                      'number_of_patterns',
                      'statc_upperbound_max_flow',
                      'total_gateway_capacity_static',
                      'static_avg_max_flow',
                      'dynamic_upperbound_max_flow',
                      'total_gateway_capacity_dynamic',
                      'dynamic_avg_max_flow'
                      ]
    self.printSummaryAfterSavingStat = False
    #---end of class fields------
  
  
  def reset(self):
    ILPSolver.reset(self)
    #---class fields-------
    self.dynamic_upperbound_flow = None
    self.dynamic_avg_flow = None
    self.static_upperbound_flow = None
    self.static_avg_flow = None
    self.total_gateway_capacity_static = None
    self.total_gateway_capacity_dynamic = None
    self.path_to_java_code_for_avg_calc = './java/tm.jar'
    self.dynamic_graph_spec_file_path = './java/temp_dynamic_spec.txt'
    self.java_code_stat_file_path = './java/temp_java_stat.txt'
    self.java_code_stdout_file_path = './java/temp_java_stdout.txt'
    self.skipJavaCodeCall = False

  def saveGraphInFile(self, g, filePath):
    '''
    save the graph g in 'fileName' file as follows:
    Nodes:
      n x y
      ...
    Edges:
      n1 n2 eType
      ...
    
    Gateways:
      g1
      ...
    
    Here, x,y are the coordinates of node n
    eType is the edge type (either long or short) of edge n1-n2
    '''
    with open(filePath, 'w') as f:
      gateway_list = []
      f.write('Nodes:\n')
      for n in g.nodes():
        if n in self.gateways:
          gateway_list.append(n)
        x = self.node_x[n]
        y = self.node_y[n]
        f_text = str(n)+" "+str(x)+" "+str(y)+"\n"
        f.write(f_text)
        
      f.write('Edges:\n')
      for u,v in g.edges():
        edge_type = self.getEdgeType(u, v)
        if edge_type == 'long' or edge_type == 'short':
          f_text =  str(u)+" "+str(v)+" "+edge_type+"\n"
          f.write(f_text)
      
      f.write('Gateways:\n')
      for g in gateway_list:
        f.write(str(g)+"\n")
      
      
  def createDynamicGraphSpecOutputFile(self):
    '''
    creates the temporary dynamic spec file and save the self.dynamic_graph according
    to TM's spec
    '''
    with open(self.dynamic_graph_spec_file_path, 'w') as f:
      f.write("trafficSources:\n")
      for n in self.dynamic_graph.nodes():
        f.write(str(n)+"\n")
      f.write('\n')
      
      
      #---FSONodes:
      f.write("FSONodes:\n")
      gateway_degree_in_static_graph = self.static_graph.degree(self.gateways)
      nodes = self.dynamic_graph.nodes() 
      for n in nodes:
        no_of_fso_for_n = self.fso_per_node
        if n in self.gateways:
          no_of_fso_for_n =  self.fso_per_gateway #!!!TODO: change back to following if needed
          #no_of_fso_for_n = min(self.fso_per_gateway, gateway_degree_in_static_graph[n])
        for fso in range(1, no_of_fso_for_n+1):
          f.write(str(n)+"_fso"+str(fso)+"\n")
      f.write('\n')   
      
      #---FSOLinks:---
      f.write("FSOLinks:\n")   
      edges = self.dynamic_graph.edges()
      for u,v in edges:
        no_of_fso_for_u = self.fso_per_node
        if u in self.gateways:
          no_of_fso_for_u = self.fso_per_gateway #!!!TODO: change back to following if needed
          #no_of_fso_for_u = min(self.fso_per_gateway, gateway_degree_in_static_graph[u])
        for fso1 in range(1, no_of_fso_for_u+1):
          no_of_fso_for_v = self.fso_per_node
          if v in self.gateways: 
            no_of_fso_for_v = self.fso_per_gateway #!!!TODO: change back to following if needed
            #no_of_fso_for_v = min(self.fso_per_gateway, gateway_degree_in_static_graph[v])
          for fso2 in range(1, no_of_fso_for_v+1):  
            #"0_0_fso2To2_0_fso1 10Gbps"
            f_text = str(u)+"_fso"+str(fso1)+"To"\
                    +str(v)+"_fso"+str(fso2)+" "+str(self.link_capacity)+"Mbps\n" 
            f.write(f_text)
            f_text = str(v)+"_fso"+str(fso2)+"To"\
                    +str(u)+"_fso"+str(fso1)+" "+str(self.link_capacity)+"Mbps\n" 
            f.write(f_text)
      f.write('\n')   
      #------gateways--:
      dynamic_gateways =  list(set(self.gateways) & set(nodes))
      f.write('gateways:\n')
      for n in dynamic_gateways:
        f.write(str(n)+"\n")
        
  def callJavaCodeToGetDynamicAvgFlow(self):
    '''
    i) prepare the dynamic graph spec and save it in ./java/temp_dynamic_spec.txt
    ii) call the program ./java/tm.jar
    iii) after the program returns, read the file ./java/temp_java_output.txt and collect
       the upper-bound and avg-flow and save those as floats
    '''
    self.createDynamicGraphSpecOutputFile()
    #open the file for streaing java stdout
    if not self.skipJavaCodeCall:
      with  open(self.java_code_stdout_file_path, 'w') as f_java_stdout:
        subprocess.call(['java','-jar',
                         self.path_to_java_code_for_avg_calc,
                         self.dynamic_graph_spec_file_path,
                         self.java_code_stat_file_path,
                         str(self.percent_of_pattern_nodes_in_avg_flow_calculation),
                         str(self.number_of_pattern_in_avg_flow_calculation) 
                         ],
                        stdout = None) #TODO: change back to stdout=None if any slower for file
                        #stdout = f_java_stdout)
  
      with open(self.java_code_stat_file_path,'r') as f:
        for line in f:
          vals = line.split(',')
          self.dynamic_upperbound_flow = float(vals[0])
          self.dynamic_avg_flow = float(vals[1])
  
  def getMaxFlowForSourceNodes(self, g, sourceList):
    '''
    processing:
      given graph g(networkx.Graph instance) and list of nodes sourceList in g:
      i) create a new Digraph obj  instance max_flow_g from g
      ii) add 'capacity' of 1 in each edge of g
      iii) add super source 'src' and add a directed edge to each non-gateway node n in sourceList with 
          capacity infinite i.e. src->n with capacity infinity
      iv) add super sink 'snk' and add a directed edge from each gateway node n in g with to 'snk'
          with capacity infinite i.e.n->snk with capacity infinity 
      v) run max-flow on the max_flow_g, retrieve the max_flow and return it
    ''' 
      
    #task (i)+(ii)
    max_flow_g =  nx.DiGraph()
    max_flow_g.add_nodes_from(g.nodes())
    for u,v in g.edges():
      max_flow_g.add_edge(u, v, capacity = 1.0)
      max_flow_g.add_edge(v, u, capacity = 1.0)
      
    #task (iii)
    max_flow_g.add_node('src')
    for n in sourceList:
      if n not in self.gateways:
        max_flow_g.add_edge('src', n, capacity = float('inf'))
    
    #task (iv)
    max_flow_g.add_node('snk')
    for n in self.gateways: #assuming all gateways are already there in the graph g, otherwise error!!
      max_flow_g.add_edge(n, 'snk', capacity = float('inf'))
    
    #task (v)
    residual_max_flow_g = flow.shortest_augmenting_path(G=max_flow_g, s='src', t='snk')
    
    return residual_max_flow_g.graph['flow_value']
  
  def getStaticAvgFlow(self):
    '''
    This method runs max-flow for pattern of non-gateway nodes number of times as specified
    in the config file
    processing:
      i) find the number of nodes for each pattern
      ii) find the list of non-gateway nodes in static graph
      iii) for number of iteration (it's 1 if pattern_node>=non-gateway nodes)
          1. find max_flow connecting only the pattern node to super-src and 
            the gateways to super-sink
          2. store the sum of degrees and max-flow in separate lists
      iv) after all these iterations, find the mean of the above list values 
          and store those as static upper bound and static avg flow
          
    '''  
    #task (i)
    candidate_pattern_nodes = list( set(self.static_graph.nodes()) - set(self.gateways))
    self.logger.debug("candidate_pattern_nodes:"+str(candidate_pattern_nodes))
    if not candidate_pattern_nodes: #no non-gateway nodes in the static graph, so this avg finding is mute
      self.static_avg_flow = 0.0
      self.static_upperbound_flow = 0.0
      self.logger.info("No non-gateway node found in static graph, avg and upperbound max flows set to 0!!")
      return
    number_of_candidate_pattern_nodes = len(candidate_pattern_nodes)
    number_of_pattern_nodes = \
      int(round(1.0*number_of_candidate_pattern_nodes * self.percent_of_pattern_nodes_in_avg_flow_calculation/100.0, 0))
    self.logger.debug("number_of_pattern_nodes:"+str(number_of_pattern_nodes))
    self.logger.debug("number_of_candidate_pattern_nodes:"+str(number_of_candidate_pattern_nodes))
    if number_of_pattern_nodes<1: #check if the ratio was too low
      number_of_pattern_nodes = 1
    
    number_of_iterations = self.number_of_pattern_in_avg_flow_calculation
    
    if number_of_pattern_nodes>= number_of_candidate_pattern_nodes:
      self.logger.debug("number_of_pattern_nodes>= number_of_candidate_pattern_nodes!!")
      number_of_iterations =  1
      number_of_pattern_nodes = number_of_candidate_pattern_nodes
      
    static_max_flow_vals = []
    static_upper_bound_vals = []
    
    for i in xrange(number_of_iterations):
      pattern_nodes = random.sample(candidate_pattern_nodes, number_of_pattern_nodes)
      max_flow_val = self.getMaxFlowForSourceNodes(self.static_graph, pattern_nodes)
      pattern_node_degrees =  self.static_graph.degree(pattern_nodes)
      upper_bound_flow = sum(pattern_node_degrees.values()) 
      
      static_max_flow_vals.append(max_flow_val)
      static_upper_bound_vals.append(upper_bound_flow)
      self.logger.debug("\tmax_flow_val:"+str(max_flow_val))
      self.logger.debug("\tupper_bound_flow:"+str(upper_bound_flow))
      self.logger.debug("\tpattern_node_degrees:"+str(pattern_node_degrees))
      
    self.static_avg_flow = 1.0*self.link_capacity*sum(static_max_flow_vals)/len(static_max_flow_vals)
    self.static_upperbound_flow = 1.0*self.link_capacity*sum(static_upper_bound_vals)/len(static_upper_bound_vals)
    
  def saveStatInFile(self):
    '''
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
    
    stat_row['dynamic_upperbound_max_flow'] = self.dynamic_upperbound_flow
    stat_row['total_gateway_capacity_dynamic'] = self.total_gateway_capacity_dynamic
    stat_row['dynamic_avg_max_flow'] = self.dynamic_avg_flow
    
    
    f= open(self.output_statistics_file, 'a')
    writer = csv.DictWriter(f, self.stat_header)
    writer.writerow(stat_row)
    f.close()
    if self.printSummaryAfterSavingStat:
      print"--------------Stat-Summary--------------"
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
      print "Ratio of avg. max. flow to upperbound max. flow for Static graph:",\
            1.0*stat_row['static_avg_max_flow']/stat_row['statc_upperbound_max_flow']
      print "Ratio of avg. max. flow to total gateway capacity for Static graph:",\
            1.0*stat_row['static_avg_max_flow']/stat_row['total_gateway_capacity_static']
      print "Dynamic:.................."
      print "Dynamic upperbound max. flow:",stat_row['dynamic_upperbound_max_flow']
      print "Dynamic total gateway capacity:",stat_row['total_gateway_capacity_dynamic']
      print "Dynamic avg. max. flow:",stat_row['dynamic_avg_max_flow'], "!!"
      print "Ratio of avg. max. flow to upperbound max. flow for Dynamic graph:",
      if stat_row['dynamic_upperbound_max_flow']>0:
        print 1.0*stat_row['dynamic_avg_max_flow']/stat_row['dynamic_upperbound_max_flow']
      else:
        print "N/A"
      print "Ratio of avg. max. flow to total gateway capacity for Dynamic graph:",\
            1.0*stat_row['dynamic_avg_max_flow']/stat_row['total_gateway_capacity_dynamic']
      print "==============End of stat-summary==========="
  
  def setTotalGatewayCapacities(self):
    '''
    computes the sum of degrees of gateway nodes in both static and dynamic graph
    and store the total capacities
    '''
    static_gateway_degrees = sum(self.static_graph.degree(self.gateways).values())
    self.total_gateway_capacity_static =  1.0*self.link_capacity*static_gateway_degrees
    #dynamic_gateway_degrees = sum(self.dynamic_graph.degree(self.gateways).values())
    self.total_gateway_capacity_dynamic = 1.0*self.link_capacity*\
            self.fso_per_gateway*len(self.gateways)
    
  def runStatCollector(self):
    '''
    all the processing related to this class are done here
    tasks:
      i) generate the dynamic graph and save it in the appropriate format
        in the current directory as temp_dynamic_graph.txt (overwrite mode)
      ii) run method to get upper-bound flow and avg. max flow for static graph
      iii) call TM's code (jar file in the same directory) and wait for its termination
      iv) retrieve the avg. max flow and upper bound flow from the file java_output.txt 
        (overwritten) file and save it
      v) save all the experiment stats in the self.output_statistics_file file (append mode)
      vi) save all the graphs with proper title in the self.graph_output_folder folder
    '''
    self.logger.info("in runStatCollector method...")
    
    self.getStaticAvgFlow()
    
    self.setTotalGatewayCapacities()
    
    self.logger.info("Running java code for dynamic avg flow calculation...")
    self.callJavaCodeToGetDynamicAvgFlow()
    self.saveStatInFile()
    
    
    