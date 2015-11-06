import unittest
from final_fso.src.step_4_dynamic import Step_4_dynamic
import matplotlib.pyplot as plt
import networkx as nx
import logging

class TestStep_4_dynamic(unittest.TestCase):
  def setUp(self):
    unittest.TestCase.setUp(self)
    self.s4_dynamic = Step_4_dynamic('./config.txt')
    #self.gi.generateNodePositions()
    print "Test for Step_4_dynamic class.."
  
  

  def test_Step_4__dynamic_graph_plot(self):
    '''
    selective calls to methods...
    '''
    self.s4_dynamic.runInputGenerator()
    self.s4_dynamic.runStep_1()
    self.s4_dynamic.runStep_2()
    self.s4_dynamic.runStep_3()
    #self.s4_dynamic.logger.setLevel(logging.DEBUG)
    self.s4_dynamic.runStep_4_dynamic()
    
    graphs = [self.s4_dynamic.adj, 
              self.s4_dynamic.short_edge_adj, 
              self.s4_dynamic.backbone_graph_after_step_2, 
              self.s4_dynamic.backbone_graph,
              self.s4_dynamic.dynamic_graph]
    
    for i in xrange(len(graphs)):
      print "DEBUG:i",i
      plt_graph = graphs[i]
      plt_title = graphs[i].graph['name']
      plt.figure(i+1)
      plt.grid(True)
      plt.xticks( xrange(0,int(self.s4_dynamic.max_x_coord)+1,int(self.s4_dynamic.target_spacing)) )
      plt.yticks( xrange(0,int(self.s4_dynamic.max_y_coord)+1,int(self.s4_dynamic.target_spacing)) )
      plt.title(plt_title)
      self.s4_dynamic.visualizeGraph(plt_graph)
    
    plt.show()
  
    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  