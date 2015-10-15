import math
from platypus.core import POSITIVE_INFINITY

def euclidean_dist(solution1, solution2):
    o1 = solution1.objectives
    o2 = solution2.objectives
    return math.sqrt(sum([math.pow(o1[i]-o2[i], 2.0) for i in range(len(o1))]))

class DistanceMatrix(object):
    """Maintains pairwise distances between solutions.
    
    The distance matrix, used by SPEA2, maintains the pairwise distances
    between solutions.  It also provides convenient routines to lookup the
    distance between any two solutions, find the most crowded solution, and
    remove a solution.
    """
    
    def __init__(self, solutions, distance_fun=euclidean_dist):
        super(DistanceMatrix, self).__init__()
        self.distances = []
        
        for i in range(len(solutions)):
            distances_i = []
            for j in range(len(solutions)):
                if i != j:
                    distances_i.append((j, distance_fun(solutions[i], solutions[j])))
                      
            self.distances.append(sorted(distances_i, cmp=lambda x, y : cmp(x[1], y[1])))                       
    
    def find_most_crowded(self):
        """Finds the most crowded solution.
        
        Returns the index of the most crowded solution, which is the solution
        with the smallest distance to the nearest neighbor.  Any ties are
        broken by looking at the next distant neighbor.
        """
        minimum_distance = POSITIVE_INFINITY
        minimum_index = -1
        
        for i in range(len(self.distances)):
            distances_i = self.distances[i]
            
            if distances_i[0][1] < minimum_distance:
                minimum_distance = distances_i[0][1]
                minimum_index = i
            elif distances_i[0][1] == minimum_distance:
                for j in range(len(distances_i)):
                    dist1 = distances_i[j][1]
                    dist2 = self.distances[minimum_index][j][1]
                    
                    if dist1 < dist2:
                        minimum_index = i
                        break
                    if dist2 < dist1:
                        break
        
        return minimum_index
    
    def remove_point(self, index):
        """Removes the distance entries for the given solution.
        
        Parameters
        ----------
        index : int
            The index of the solution
        """
        del self.distances[index]
        
        for i in range(len(self.distances)):
            self.distances[i] = [(x if x < index else x-1, y) for (x, y) in self.distances[i] if x != index]
    
    def kth_distance(self, i, k):
        """Returns the distance to the k-th nearest neighbor.
        
        Parameters
        ----------
        i : int
            The index of the solution
        k : int
            Finds the k-th nearest neightbor distance
        """
        return self.distances[i][k][1]
    
    def __getitem__(self, key):
        if isinstance(key, tuple) and len(key) == 2:
            if key[0] == key[1]:
                return 0.0
            else:
                for i, d in self.distances[key[0]]:
                    if i == key[1]:
                        return d
                
                raise ValueError("key not found")
        else:
            raise ValueError("key must be a tuple")