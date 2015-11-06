import unittest
from final_fso.src.step_1 import Step_1


class TestStep_1(unittest.TestCase):
  def setUp(self):
    unittest.TestCase.setUp(self)
    self.s1 = Step_1('./config.txt')
    #self.gi.generateNodePositions()
    print "Test for class Step_1 .."
  

  def test_heuristicTargetCover(self):
    self.s1.runInputGenerator()
    self.s1.runStep_1()
    print "Set Covers:", sorted(self.s1.node_cover)
    print "Total Set Cover nodes:",len(self.s1.node_cover)
    print "Total Nodes:",self.s1.number_of_nodes
    print "Targets Covered:",self.s1.no_of_targets_covered
    print "Total Targets:",self.s1.total_targets
    print "Target max_x:",self.s1.max_tx
    print "Target max_y:",self.s1.max_ty
    print "Is Node-heap empty:",self.s1.node_hp.isHeapEmpty()
    
    
  
    
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  