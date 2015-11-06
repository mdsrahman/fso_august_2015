import unittest
from final_fso.src.step_2 import Step_2
import matplotlib.pyplot as plt

class TestTargetHandler(unittest.TestCase):
  def setUp(self):
    unittest.TestCase.setUp(self)
    self.s2 = Step_2('./config.txt')
    #self.gi.generateNodePositions()
    print "Test for Step_2 class.."
    

  def test_Step2_graph_plot(self):
    '''
    selective calls to methods...
    '''
    self.s2.runInputGenerator()
    self.s2.runStep_1()
    self.s2.runStep_2()
    
    graphs = [self.s2.adj, self.s2.short_edge_adj, self.s2.backbone_graph]
    graph_title = ['Adjacency graph', 'Short Edge Adjacency graph', 'Backbone graph']
    
    for i in xrange(len(graphs)):
      plt_graph = graphs[i]
      plt_title = graph_title[i]
      plt.figure(i+1)
      plt.grid(True)
      plt.xticks( xrange(0,int(self.s2.max_x_coord)+1,int(self.s2.target_spacing)) )
      plt.yticks( xrange(0,int(self.s2.max_y_coord)+1,int(self.s2.target_spacing)) )
      plt.title(plt_title)
      self.s2.visualizeGraph(plt_graph)
    
    plt.show()
  
    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  