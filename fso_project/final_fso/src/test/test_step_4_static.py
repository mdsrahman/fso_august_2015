import unittest
from final_fso.src.step_4_static import Step_4_static
import matplotlib.pyplot as plt
import logging
import networkx as nx

class TestStep_4_static(unittest.TestCase):
  def setUp(self):
    unittest.TestCase.setUp(self)
    self.s4_static = Step_4_static('./config.txt')
    #self.gi.generateNodePositions()
    print "Test for Step_4_static class.."
  
  

  def test_Step_4__static_graph_plot(self):
    '''
    selective calls to methods...
    '''
    #self.s4_static.logger.setLevel(logging.DEBUG)
    self.s4_static.runInputGenerator()

    #self.s4_static.logger.setLevel(logging.INFO)
    self.s4_static.runStep_1()
    self.s4_static.runStep_2()
    self.s4_static.runStep_3()
    #self.s4_static.logger.setLevel(logging.DEBUG)
    self.s4_static.runStep_4_dynamic()
    self.s4_static.fso_per_node = 2
    #self.s4_static.logger.setLevel(logging.DEBUG)
    self.s4_static.runStep_4_static()
    print "No. of gateways:" ,len(self.s4_static.gateways)
    
    graphs = [self.s4_static.adj, 
              self.s4_static.short_edge_adj, 
              self.s4_static.backbone_graph_after_step_2, 
              self.s4_static.backbone_graph,
              self.s4_static.dynamic_graph,
              self.s4_static.static_graph]
    
    for i in xrange(len(graphs)):
      plt_graph = graphs[i]
      plt_title = graphs[i].graph['name']
      plt.figure(i+1)
      plt.grid(True)
      plt.xticks( xrange(0,int(self.s4_static.max_x_coord)+1,int(self.s4_static.target_spacing)) )
      plt.yticks( xrange(0,int(self.s4_static.max_y_coord)+1,int(self.s4_static.target_spacing)) )
      plt.title(plt_title)
      self.s4_static.visualizeGraph(plt_graph)
    
    plt.show()
  
    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  