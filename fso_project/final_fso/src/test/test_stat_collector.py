import unittest
from final_fso.src.stat_collector import StatCollector
import matplotlib.pyplot as plt
import logging
import networkx as nx

class TestStatCollector(unittest.TestCase):
  def setUp(self):
    unittest.TestCase.setUp(self)
    self.stat_collector = StatCollector('./config.txt')
    #self.gi.generateNodePositions()
    print "Test for StatCollector class.." 
  
  def visualizeAllGraphs(self):
    covering_node_graph = nx.Graph()
    covering_node_graph.graph['name'] = 'Target-covering nodes only graph'
    covering_node_graph.add_nodes_from(self.stat_collector.node_cover)
    
    graphs = [self.stat_collector.adj, 
              self.stat_collector.short_edge_adj, 
              covering_node_graph,
              self.stat_collector.backbone_graph_after_step_2, 
              self.stat_collector.backbone_graph,
              self.stat_collector.dynamic_graph,
              self.stat_collector.static_graph]
    
    #for i in xrange(len(graphs)):
    for i in xrange(len(graphs)):
      plt_graph = graphs[i]
      plt_title = graphs[i].graph['name']
      plt.figure(i+1)
      
      plt.xticks( xrange(0,int(self.stat_collector.max_x_coord)+1,int(self.stat_collector.target_spacing)) )
      plt.yticks( xrange(0,int(self.stat_collector.max_y_coord)+1,int(self.stat_collector.target_spacing)) )
      plt.axis([-10, self.stat_collector.max_x_coord+10, -10, self.stat_collector.max_y_coord+10])
      plt.grid(True)
      plt.title(plt_title)
      
      self.stat_collector.visualizeGraph(plt_graph)
    
    plt.show()

  def getNodesNotConnectedToGatewaysInDynamicGraph(self):
    '''
    returns the list of nodes not directly connected to the gateways
    must be called after self.adj has been set/built
    '''
    list_of_nodes = []
    for n in self.stat_collector.dynamic_graph.nodes():
      if n not in self.stat_collector.gateways:
        isConnectedToGateways = False
        for g in self.stat_collector.gateways:
          if n==g or self.stat_collector.dynamic_graph.has_edge(n, g):
            isConnectedToGateways = True
            break
        if not isConnectedToGateways:
          list_of_nodes.append(n)
    return list_of_nodes
      
    
  def test_StatCollectorMethods(self):
    '''
    selective calls to methods...
    '''
    #self.stat_collector.logger.setLevel(logging.DEBUG)
    self.stat_collector.runInputGenerator()
    #self.stat_collector.logger.setLevel(logging.INFO)
    print "Input Generator Complete"
    
    #self.stat_collector.logger.setLevel(logging.DEBUG)
    self.stat_collector.runStep_1()
    #self.stat_collector.logger.setLevel(logging.INFO)
    print "step 1 complete"
    
    #self.stat_collector.logger.setLevel(logging.DEBUG)
    self.stat_collector.runStep_2()
    #self.stat_collector.logger.setLevel(logging.INFO)
    print "step 2 complete"
    
    #self.stat_collector.logger.setLevel(logging.DEBUG)
    self.stat_collector.runStep_3()
    #self.stat_collector.logger.setLevel(logging.INFO)
    print "step 3 complete"
    
    #self.stat_collector.logger.setLevel(logging.DEBUG)
    self.stat_collector.runStep_4_dynamic()
    print "step 4-dynamic complete"
    #self.stat_collector.logger.setLevel(logging.INFO)
    
    #self.stat_collector.logger.setLevel(logging.DEBUG)
    self.stat_collector.runStep_4_static()
    #self.stat_collector.logger.setLevel(logging.INFO)
    print "step 4-static complete"
    
    self.stat_collector.runILPSolver()
     
    
    
    #need to reset the relative path as python code is called at './src/test' rather than at './src'
    self.stat_collector.path_to_java_code_for_avg_calc = '../java/tm.jar'
    self.stat_collector.dynamic_graph_spec_file_path = '../java/temp_dynamic_spec.txt'
    self.stat_collector.java_code_stat_file_path = '../java/temp_java_stat.txt'
    self.stat_collector.java_code_stdout_file_path = '../java/temp_java_stdout.txt' 
    
    self.stat_collector.printSummaryAfterSavingStat = True #will now print on screen summary
    
    #self.stat_collector.logger.setLevel(logging.DEBUG)
    self.stat_collector.runStatCollector()
    #self.stat_collector.logger.setLevel(logging.INFO)
    
    print "Non-gateway nodes not connected to the gateway in the dynamic graph:",
    print self.getNodesNotConnectedToGatewaysInDynamicGraph()
    #----------------------------------
    self.visualizeAllGraphs()
    #----------------------------------
    
    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  