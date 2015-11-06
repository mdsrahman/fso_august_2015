import unittest
from final_fso.src.my_util import MyGridBin , MyHeap
import random
 
class TestMyUtil(unittest.TestCase):
  
  def test_MyGridBin(self):
    max_x =  100
    max_y =  100
    x_interval = 10
    y_interval = 10
    gb = MyGridBin(
                   x_interval,
                   y_interval,
                   max_x,
                   max_y
                   )
    
    gb.put(0,12,12)
    gb.put(1,22,22)
    gb.put(2,25,10)
    gb.put(3,0,5)
    
    print"DEBUG:bin content:"
    for i in xrange(gb.max_bx+1):
      for j in xrange(gb.max_by+1):
        print gb.bin[i][j],",",
      print ""
    
    ret_list =  sorted(gb.get(9.9,9.9,10))
    exp_list = sorted([0,3])
    
    self.assertEqual(ret_list, exp_list, 'wrong object list'+
                        str(ret_list)+' vs '+str(exp_list)+
                        ' bin-content:'+ str(gb.bin) )
  def test_MyHeap(self):
    hp = MyHeap()
    self.assertEqual(hp.isHeapEmpty(),True,'error in the method isHeapEmpty')
    hp.push(0,-9)
    hp.push(2,-90)
    hp.push(100,-100)
    hp.push(2, 90)
    hp.push(150,150)
    
    self.assertEqual(hp.pop(),100,'heap pop error')
    self.assertEqual(hp.pop(),0,'heap pop error')
    self.assertEqual(hp.pop(),2,'heap pop error')
    self.assertEqual(hp.isHeapEmpty(),False,'error in the method isHeapEmpty')
    
    self.assertEqual(hp.pop(),150,'heap pop error')
    self.assertEqual(hp.pop(),None,'heap pop error')
    self.assertEqual(hp.isHeapEmpty(),True,'error in the method isHeapEmpty')
    hp.push(150,-150)
    self.assertEqual(hp.pop(),150,'heap pop error')
    self.assertEqual(hp.pop(),None,'heap pop error')
    self.assertEqual(hp.isHeapEmpty(),True,'error in the method isHeapEmpty')
if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  
  
  
  
  
  
  
  
  
  
  