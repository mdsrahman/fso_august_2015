import unittest
from final_fso.src.extended_dynamic_stat_collector import ExtendedDynamicStatCollector
import matplotlib.pyplot as plt
import logging
import networkx as nx

class TestStatCollector(unittest.TestCase):
  def setUp(self):
    unittest.TestCase.setUp(self)
    configFileName = './config.txt'
    self.ext_stat_collector = ExtendedDynamicStatCollector(configFileName)
    #self.gi.generateNodePositions()
    print "Test for ExtDynamicStatCollector class.." 
  
  def visualizeAllGraphs(self):
    covering_node_graph = nx.Graph()
    covering_node_graph.graph['name'] = 'Target-covering nodes only graph'
    covering_node_graph.add_nodes_from(self.ext_stat_collector.node_cover)
    
    graphs = [self.ext_stat_collector.adj, 
              self.ext_stat_collector.short_edge_adj, 
              covering_node_graph,
              self.ext_stat_collector.backbone_graph_after_step_2, 
              self.ext_stat_collector.backbone_graph,
              self.ext_stat_collector.dynamic_graph,
              self.ext_stat_collector.static_graph]
    
    #for i in xrange(len(graphs)):
    for i in xrange(len(graphs)):
      plt_graph = graphs[i]
      plt_title = graphs[i].graph['name']
      plt.figure(i+1)
      
      plt.xticks( xrange(0,int(self.ext_stat_collector.max_x_coord)+1,int(self.ext_stat_collector.target_spacing)) )
      plt.yticks( xrange(0,int(self.ext_stat_collector.max_y_coord)+1,int(self.ext_stat_collector.target_spacing)) )
      plt.axis([-10, self.ext_stat_collector.max_x_coord+10, -10, self.ext_stat_collector.max_y_coord+10])
      plt.grid(True)
      plt.title(plt_title)
      
      self.ext_stat_collector.visualizeGraph(plt_graph)
    
    plt.show()

  def getNodesNotConnectedToGatewaysInDynamicGraph(self):
    '''
    returns the list of nodes not directly connected to the gateways
    must be called after self.adj has been set/built
    '''
    list_of_nodes = []
    for n in self.ext_stat_collector.dynamic_graph.nodes():
      if n not in self.ext_stat_collector.gateways:
        isConnectedToGateways = False
        for g in self.ext_stat_collector.gateways:
          if n==g or self.ext_stat_collector.dynamic_graph.has_edge(n, g):
            isConnectedToGateways = True
            break
        if not isConnectedToGateways:
          list_of_nodes.append(n)
    return list_of_nodes
      
    
  def test_StatCollectorMethods(self):
    self.ext_stat_collector.onlyGeneratePlot = False
    self.ext_stat_collector.runAll()
    #----------------------------------
    #self.visualizeAllGraphs()
    #----------------------------------
    
    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  