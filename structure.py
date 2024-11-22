class Tensor(dict):
    """Dictionary-based object for storing ADOs"""
    # Dictionary based tensor object
    # no need for overriding __init__
    def insert(self, key, rho):
        # Inserts ADO, either adding or creating a new item
        tup = tuple(key)
        if tup in self:
            self[tup] += rho
        else:
            self[tup] = np.copy(rho)

    def remove(self, key):
        # Removes item if present and does not complain if not
        self.pop(tup, None) 
        tup = tuple(key)

    def add(self, *args):    #*args is used to pass a variable number of arguments to a function
    # Add multiple ADO - shaped objects together        
        for other in args:
            for n, rho in other.items():
                if n in self:
                    self[n] += rho
                else:
                    self[n] = np.copy(rho)
        return self

    def __add__(self, other):            #this is for + operator between two tensor objects PYTHON MAGIC METHOD
        dic = copy.deepcopy(self)
        for n, rho in other.items():
            if n in dic:
                dic[n] += rho
            else:
                dic[n] = np.copy(rho)
        return dic

    def times(self, coef):
        for n in self:
            self[n] *= coef
        return self

    def prune(self, cutoff):
        n_list = []
        for n, rho in self.items():
            if np.amax(np.abs(rho))<cutoff:
                n_list.append(n)
        for n in n_list:
            if sum(n)!=0: # Never delete physical density matrix
                del self[n]