import unittest
from final_fso.src.input_generator import InputGenerator
import networkx as nx
from matplotlib import pyplot as plt

class TestInputGenerator(unittest.TestCase):
  def setUp(self):
    unittest.TestCase.setUp(self)
    self.gi = InputGenerator('./config.txt')
    #self.gi.generateNodePositions()
    print "Test for GenerateInput.."
    self.gi.runInputGenerator()

  def mtest_params_synthetic_graph(self):
    self.assertIsInstance(self.gi.graphType, str, 'type(str) error for:graphType')
    self.assertEqual(self.gi.graphType, 'synthetic', 'value error for:graphType')
     
    self.assertIsInstance(self.gi.seed, int, 'type(int) error for:seed')
    self.assertEqual(self.gi.seed, 101349, 'value error for:seed')
     
    self.assertIsInstance(self.gi.experiment_name, str, 'type(str) error for:experiment_name')
    self.assertEqual(self.gi.experiment_name, 'Test', 'value error for:experiment_name')
     
    self.assertIsInstance(self.gi.gateway_to_node_ratio, float, 'type(float) error for:gateway_to_node_ratio')
    self.assertEqual(self.gi.gateway_to_node_ratio, 0.01, 'value error for:gateway_to_node_ratio')
     
    self.assertIsInstance(self.gi.link_capacity, float, 'type(float) error for:link_capacity')
    self.assertEqual(self.gi.link_capacity, 1000.0, 'value error for:link_capacity')
     
    self.assertIsInstance(self.gi.fso_per_node, int, 'type(int) error for:fso_per_node')
    self.assertEqual(self.gi.fso_per_node, 4, 'value error for:fso_per_node')
     
    self.assertIsInstance(self.gi.target_spacing, float, 'type(float) error for:target_spacing')
    self.assertEqual(self.gi.target_spacing, 10, 'value error for:target_spacing')
     
    self.assertIsInstance(self.gi.max_short_edge_length, float, 'type(float) error for:max_short_edge_length')
    self.assertEqual(self.gi.max_short_edge_length, 100.0, 'value error for:max_short_edge_length')
    
    self.assertIsInstance(self.gi.max_long_edge_length, float, 'type(float) error for:max_long_edge_length')
    self.assertEqual(self.gi.max_long_edge_length, 1000.0, 'value error for:max_long_edge_length')
    
    self.assertIsInstance(self.gi.output_statistics_file, str, 'type(str) error for:output_statistics_file')
    self.assertEqual(self.gi.output_statistics_file, './stat/synthetic_graph_stat.txt', 'value error for:output_statistics_file')
      
    self.assertIsInstance(self.gi.graph_output_folder, str, 'type(str) error for:graph_output_folder')
    self.assertEqual(self.gi.graph_output_folder, './output/', 'value error for:graph_output_folder')
     
    if self.gi.graphType == 'synthetic':
      
      self.assertIsInstance(self.gi.no_of_samples, int, 'type(int) error for:no_of_samples')
      self.assertEqual(self.gi.no_of_samples, 1, 'value error for:no_of_samples')
        
      self.assertIsInstance(self.gi.number_of_nodes, int, 'type(int) error for:number_of_nodes')
      self.assertEqual(self.gi.number_of_nodes, 500, 'value error for:number_of_nodes')
        
      self.assertIsInstance(self.gi.max_no_of_edges, int, 'type(int) error for:max_no_of_edges')
      self.assertEqual(self.gi.max_no_of_edges, 1000, 'value error for:max_no_of_edges')
        
      self.assertIsInstance(self.gi.short_to_long_edge_ratio, float, 'type(float) error for:short_to_long_edge_ratio')
      self.assertEqual(self.gi.short_to_long_edge_ratio, 1.5, 'value error for:short_to_long_edge_ratio')
        
      self.assertIsInstance(self.gi.max_x_coord, float, 'type(float) error for:max_x_coord')
      self.assertEqual(self.gi.max_x_coord, 1000, 'value error for:max_x_coord')
        
      self.assertIsInstance(self.gi.max_y_coord, float, 'type(float) error for:max_y_coord')
      self.assertEqual(self.gi.max_y_coord, 500, 'value error for:max_y_coord')
        
  #TODO: test case for params read for map input
  #elif self.gi.graphType == 'map':
    #todo...

  def mtest_generateNodePositions(self):
    '''
    ensures that:
      i) both node_x and node_y has exactly  nodes = number_of_nodes
      ii node_x and node_y has maximum value of max_x_coord and max_y_coord respectivley
    '''
    self.assertEqual(len(self.gi.node_x), self.gi.number_of_nodes, "Error in the size of node_x")
    self.assertEqual(len(self.gi.node_y), self.gi.number_of_nodes, "Error in the size of node_y")
    self.assertLessEqual(max(self.gi.node_x), self.gi.max_x_coord, "node_x value exceeds max_x_coord")
    self.assertLessEqual(max(self.gi.node_y), self.gi.max_y_coord, "node_y value exceeds max_y_coord")
  
 
  def mtest_generateSyntheticGraph(self):
    '''ensure that:
        i) total nodes in adj equals self.number_of_nodes
        ii) total_edges in adj <= self.max_no_of_edges
        iii) for each node, total short edges <= self.max_short_edges_per_node
                            total long edges <= self.max_long_edges_per_node
    '''
    self.assertEqual(self.gi.adj.number_of_nodes(), self.gi.number_of_nodes, 'Incorrect no. of nodes in adj')
    self.assertLessEqual(self.gi.adj.number_of_edges(), self.gi.max_no_of_edges, 'max_no_of_edges exceeded for adj')
    
    for u in self.gi.adj.nodes():
      u_short_edge_counter = 0
      u_long_edge_counter = 0
      nbr_u = self.gi.adj.neighbors(u)
      for v in nbr_u:
        edge_data = self.gi.adj.get_edge_data(u, v)
        if edge_data['con_type'] == 'short':
          u_short_edge_counter += 1
        elif edge_data['con_type'] == 'long':
          u_long_edge_counter += 1
      self.assertLessEqual(u_short_edge_counter, self.gi.max_short_edge_per_node,'max_short_edge_per_node exceeded')
      self.assertLessEqual(u_long_edge_counter, self.gi.max_long_edge_per_node,'max_long_edge_per_node exceeded')
  
  
  
  def mtest_selectGateways(self):
    '''
    ensures that
      i) if number of gateways is higher than specified as ratio in the config file,
          then each comes from separate connected components, i.e. there is no path 
          between any pairs of gateways
      ii) else, the number of gateways must be almost equal to the ratio specified
    '''
    
    current_gateway_to_node_ratio  = 1.0*len(self.gi.gateways)/self.gi.adj.number_of_nodes()
    if current_gateway_to_node_ratio > self.gi.gateway_to_node_ratio:
      for u in self.gi.gateways:
        for v in self.gi.gateways:
          if u != v:
            self.assertFalse(nx.has_path(self.gi.adj,u,v),
              'Gateways have path between them even though the gateway ratio is higher i.e they are from the same connected component')
    else:
      self.assertAlmostEqual(current_gateway_to_node_ratio, 
                           self.gi.gateway_to_node_ratio, places = 4, 
                           msg = 'gateway_to_node ratio is erroneous by 4 decimal places:'
                           +str(current_gateway_to_node_ratio)+' vs '
                           +str(self.gi.gateway_to_node_ratio))
      
      for g in self.gi.gateways:
        self.assertGreaterEqual( self.gi.adj.degree(g), 0, 'warning! Disconnected node selected as gateway')
    

 
  def mtest_visualizeGraph(self):
    '''
    for now, just visualize the two graphs self.adj and self.short_edge_adj
    '''
    plt.figure(1)
    plt.title("Adjacency Graph")
    self.gi.visualizeGraph(self.gi.adj)
    plt.figure(2)
    plt.title("Short-edge Adjacency Graph Graph")
    self.gi.visualizeGraph(self.gi.short_edge_adj)
    plt.show()
  
  def test_summmary(self):
    '''
    calls the mtest_* methods that is only the methods of interest
    '''
    print "running inputGenerator Processes"
    
    self.mtest_generateSyntheticGraph()
    self.mtest_selectGateways()
    self.mtest_visualizeGraph()

if __name__=='__main__':
  unittest.main(verbosity = 2)
  
  