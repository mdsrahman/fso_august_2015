import unittest
from final_fso.src.step_3 import Step_3
import matplotlib.pyplot as plt
import networkx as nx

class TestStep_3(unittest.TestCase):
  def setUp(self):
    unittest.TestCase.setUp(self)
    self.s3 = Step_3('./config.txt')
    #self.gi.generateNodePositions()
    print "Test for Step_3 class.."
  
  def test_runStep_3(self):
    '''
    processing:
    before calling runStep_3:
      i) initialize self.gateways
      ii) initialize self.backbone_graph
      iii) initialize self.adj
      iv) then call self.runStep_3
      v) verify by visualization
    '''
      
    backbone_edges=[
                    ( 0, 1),
                    ( 0, 2),
                    ( 0, 3),
                    ( 3, 4),
                    ( 3, 5),
                    ( 3, 6),
                    ( 6, 7),
                    ( 6, 8),
                    ( 9, 10),
                    ( 9, 11),
                    ( 10, 12),
                    ( 10, 13),
                    ( 10, 14)
                    ]
    
      
    
    extra_edges_in_adj = [
                          ( 8, 12),
                          ( 6, 11),
                          ( 0, 12)
                         ]
    
    self.s3.gateways = [0, 9]
    
    self.s3.backbone_graph =  nx.Graph()
    self.s3.backbone_graph.graph['name'] = 'Backbone Graph'
    
    self.s3.adj = nx.Graph()
    self.s3.adj.graph['name'] = 'Adjacency Graph'
    
    for i in xrange(15):
      self.s3.backbone_graph.add_node(i)
      self.s3.adj.add_node(i)
      
    for (u,v) in backbone_edges:
      self.s3.backbone_graph.add_edge(u,v)
      self.s3.adj.add_edge(u,v,con_type='short')
      
    for (u,v) in extra_edges_in_adj:
      self.s3.adj.add_edge(u,v,con_type = 'short')
    
    self.s3.runStep_3()

    plt.figure(1)
    nx.draw_networkx(self.s3.backbone_graph_after_step_2)
    plt.figure(2)
    nx.draw_networkx(self.s3.backbone_graph)
    plt.show()

  def mtest_Step3_graph_plot(self):
    '''
    selective calls to methods...
    '''
    self.s3.runInputGenerator()
    self.s3.runStep_1()
    self.s3.runStep_2()
    self.s3.runStep_3()
    
    graphs = [self.s3.adj, 
              self.s3.short_edge_adj, 
              self.s3.backbone_graph_after_step_2, 
              self.s3.backbone_graph]
    graph_title = ['Adjacency graph', 
                   'Short Edge Adjacency graph', 
                   'Backbone graph after step 2', 
                   'Backbone graph step 3']
    
    for i in xrange(len(graphs)):
      plt_graph = graphs[i]
      plt_title = graph_title[i]
      plt.figure(i+1)
      plt.grid(True)
      plt.xticks( xrange(0,int(self.s3.max_x_coord)+1,int(self.s3.target_spacing)) )
      plt.yticks( xrange(0,int(self.s3.max_y_coord)+1,int(self.s3.target_spacing)) )
      plt.title(plt_title)
      self.s3.visualizeGraph(plt_graph)
    
    plt.show()
  
    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  