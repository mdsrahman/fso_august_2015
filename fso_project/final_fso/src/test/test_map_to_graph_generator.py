import unittest
from final_fso.src.map_to_graph_generator import MapToGraphGenerator

class TestMapToGraphGenerator(unittest.TestCase):
  def setUp(self):
    filename = 'manhattan'
    mapFilePath = './maps/'+filename+'.osm' 
    building_outputFilePath = './maps/graphs/'+filename+'.bldg.txt'
    graph_outputFilePath = './maps/graphs/'+filename+'.adj.txt'
    short_edge_length = 100
    long_edge_length = 1000  
    edge_append_mode = False
    visualize_graph = True #change this one if needed
    unittest.TestCase.setUp(self)
    self.m2gg = MapToGraphGenerator(mapFilePath, 
                                    building_outputFilePath,
                                    graph_outputFilePath,
                                    short_edge_length, 
                                    long_edge_length, 
                                    edge_append_mode,
                                    visualize_graph)
    print "Test for MapToGraphGenerator....."
    

  def test_Step2_graph_plot(self):
    '''
    call the method of MapToGraphGenerator...
    '''
    self.m2gg.generateMapToGraph()
    #self.m2gg.debugGenerateVisualGraph()
    self.m2gg.debugPrintSummary()
    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  