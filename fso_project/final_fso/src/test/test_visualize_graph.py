import unittest
from final_fso.src.visualize_graph import VisualizeGraph

class TestMapToGraphGenerator(unittest.TestCase):
  def setUp(self):
    filename = 'manhattan_o'
    buildlingFilePath = './maps/graphs/'+filename+'.bldg.txt'
    adjFilePath = './maps/graphs/'+filename+'.adj.txt'

    unittest.TestCase.setUp(self)
    self.vz = VisualizeGraph(buildlingFilePath,
                             adjFilePath)
    print "Test for VisualizeGraph....."
    

  def test_Step2_graph_plot(self):
    '''
    call the method draw.
    '''
    self.vz.draw()

    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  