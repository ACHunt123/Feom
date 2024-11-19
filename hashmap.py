
import numpy as np
import matplotlib.pyplot as plt
#
#   Hashmap.py generates a hashmap that maps the index of the BCF to the index of the ADO
#

# INPUTS
# K       #number of exponential terms in the BCF - either the truncation of the number of beads + 1
# max_N   #maximum depth of the ADO expansion

# OUTPUTS
# I_to_index #hash map from the index of the ADO to the index of the BCF
# index_to_I #hash map from the index of the BCF to the index of the ADO


def generateHashmap(K,max_N,write_to_file = False):

    # number of ADOs in the n'th tier with k elements
    def length(n,k):
        return np.math.factorial(n+k-1) // np.math.factorial(k-1) // np.math.factorial(n) # = [ (n+k-1) C (k-1)]

    #generates the set of indicies corresponding to the n'th tier with K elements
    def generatenumbers(n, k): 

        # output is a list of the strings of the indicies separated by commas
        numbers = np.empty(length(n,k),dtype=object) 
        # numbers = ['xxx' for i in range(length(n,k))]

        def generate_tuples_of_numbers(n, k):
            max_element = n
            allowed = range(max_element, -1, -1) #the allowed elements in the set

            #recursive function to generate the set of indicies. tuple t of length k whose elements sum to n
            def helper(n, k, t): 

                # k: # of elements left to choose
                # n: sum of these elements left to choose
                # t: tuple of the elements so far

                # if there are no elements left to choose and the sum of the elements is n (ie the sum of elements left to choose = 0), return the tuple
                if k == 0:
                    if n == 0: 
                        yield t 

                # if there is one element left to choose and the sum of the elements left to choose is n, return the tuple with the last element added
                elif k == 1:
                    if n in allowed:
                        yield t + (n,) 

                elif k>1 : #if there are more than one elements left to choose
                    for v in allowed: # select first digit v from the allowed digits
                        yield from helper(n - v, k - 1, t + (v,)) #recursively call the function with the sum reduced by v, the number of elements reduced by 1 and the tuple with the new element added

            return helper(n, k, ()) #call with the empty tuple ()
        tups = generate_tuples_of_numbers(n, k) #generate the tuples of indicies

        # Reformat the tuples to strings
        for i, tup in enumerate(tups):
            numbers[i] = str(tup).replace('(','').replace(')','').replace(' ','')

        return numbers

    # OUTPUTS
    # I is a natural number that gives the position of the ADO in the  density
    # index are the indicies of the ADO

    allnums = np.array([])
    for i in range(0,max_N+1):
        allnums= np.concatenate((allnums,generatenumbers(i,K))) #concatenate the set of indicies to the list of all indicies

    I_to_index = {index:I for index,I in enumerate(allnums) } #hash map from the index of the ADO to the index of the BCF
    index_to_I = {I:index for index,I in enumerate(allnums) } #hash map from the index of the BCF to the index of the ADO

    if write_to_file:
        with open('hashmap.txt','w') as f:
            f.write(str(I_to_index))
            f.write('\n')
            f.write(str(index_to_I))

    return I_to_index, index_to_I

if __name__ == '__main__':
    # index_to_I = dict() #hash map from the index of the BCF to the index of the ADO

    # INPUTS
    K = 2      #number of exponential terms in the BCF - either the truncation of the number of beads + 1
    max_N = 5   #maximum depth of the ADO expansion

    I_to_index, index_to_I = generateHashmap(K,max_N)

    print(I_to_index)
    print(index_to_I)

