'''
This module contains custom built utilities classes

'''
import heapq as hq

class MyHeap: 
  def __init__(self, removed_obj_marker='!'):
    self.pq = []                         # list of entries arranged in a heap
    self.object_finder = {}               # mapping of objects to entries
    self.REMOVED = removed_obj_marker      # placeholder for a removed obj
    
  def push(self, obj, priority):
      if obj in self.object_finder: #this might be costly for a large heap
        entry = self.object_finder.pop(obj)
        entry[1] = self.REMOVED
          
      entry = [priority, obj]
      self.object_finder[obj] = entry
      
      hq.heappush(self.pq, entry)
  
  def pop(self):
      while self.pq:
          priority, obj = hq.heappop(self.pq)
          if obj != self.REMOVED:
              del self.object_finder[obj]
              return obj
      return None #return None if trying to pop from empty heap
    
  def debugPrintHeapContent(self):
    
    print "Heap Content:-------"
    for obj in self.object_finder:
      entry = self.object_finder[obj]
      if entry[1] != self.REMOVED:
        print "Obj:",entry[1]," priority:",entry[0]
        
  def isHeapEmpty(self):
    if not self.object_finder:
      return True
    else:
      return False
  
class MyGridBin():
  '''
  bin class-sketch:
  creates a 2D bin, constructors take bin x-interval, bin y-interval, max_x, max_y
  method: set(id, x,y), x,y absolute coordinate, bin the object id to the approp.bin
  method: get(x,y,range=0): x,y absolute coordinate, range optionally 0,
          returns the ids of objects belonging to the bin that the object (x,y) would belong
          also if range>0, returns the ids of all objects belonging to the bins of rectangular
          area that are range distance away from x,y positions, may return some bin's object that might
          not intersect with the circle of radius centered at x,y
  '''
  def __init__(self, x_interval, y_interval, max_x, max_y):
    self.x_interval = x_interval
    self.y_interval = y_interval
    self.max_x = max_x
    self.max_y = max_y
    
    
    self.max_bx  = int(self.max_x // self.x_interval)
    self.max_by  = int(self.max_y // self.y_interval)
    
    #initialize the bin-array
    self.bin = []
    for i in xrange(self.max_bx+1):
      self.bin.append([])
      for j in xrange(self.max_by+1):
        self.bin[i].append([])
  
  def put(self, id, x, y):
    bx = int(x//self.x_interval)
    by = int(y//self.y_interval)
    
    self.bin[bx][by].append(id)
    
  def getMaxGridCoords(self):
    return self.max_bx, self.max_by
  
  def getGridCoords(self,x,y):
    bx = int(x//self.x_interval)
    by = int(y//self.y_interval)
    return bx,by
  
  def getbyGridCoords(self, bx, by):
    if bx<=self.max_bx and by<=self.max_by:
      return list(self.bin[bx][by])
    else:
      return []
  
  def get(self, x, y, radius = 0.0):
    return_list = []
    bx = int(x//self.x_interval)
    by = int(y//self.y_interval)
    
    return_list.extend(self.bin[bx][by])
    
    if radius>0.0:
      r_min_bx = int(max(0, x -  radius)//self.x_interval)
      r_max_bx = int(min(self.max_x, x + radius)//self.x_interval)
      
      r_min_by = int(max(0, y -  radius)//self.y_interval)
      r_max_by = int(min(self.max_y, y + radius)//self.y_interval)
      
      for i in xrange(r_min_bx, r_max_bx+1):
        for j in xrange(r_min_by, r_max_by+1):
          if i==bx and j==by:
            continue
          return_list.extend(self.bin[i][j])
      
    return return_list
  
  
  
  
  
    
