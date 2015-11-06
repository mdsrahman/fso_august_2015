import unittest
from final_fso.src.map_to_graph_generator import MapToGraphGenerator

class TestMapToGraphGeneratorVisualizer(unittest.TestCase):
  def setUp(self):
    mapFilePath = './maps/tiny_nyc.osm'
    outputFilePath = './maps/graphs/dummy.txt'
    short_edge_length = 100
    long_edge_length = 1000 
    visualize_graph = True
    unittest.TestCase.setUp(self)
    self.losFilePath = './maps/graphs/tiny_nyc_graph.txt'
    self.m2ggv = MapToGraphGenerator(mapFilePath, 
                                    outputFilePath, 
                                    short_edge_length, 
                                    long_edge_length, 
                                    visualize_graph)
    print "Test for MapToGraphGenerator....."
    

  def test_Step2_graph_plot(self):
    '''
    call the method of MapToGraphGenerator...
    '''
    self.m2ggv.parseMapFile()
    self.m2ggv.transformToCartesianCoord()
    self.m2ggv.loadLOS(self.losFilePath)
    self.m2ggv.debugPrintSummary()
    self.m2ggv.debugGenerateVisualGraph()
    
    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  