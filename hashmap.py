
import numpy as np
import matplotlib.pyplot as plt
#
#   Hashmap.py generates a hashmap that maps the index of the BCF to the index of the ADO
#

# INPUTS
# K       #number of exponential terms in the BCF - either the truncation of the number of beads + 1
# L   #maximum depth of the ADO expansion

# OUTPUTS
# I_to_index #hash map from the index of the ADO to the index of the BCF
# index_to_I #hash map from the index of the BCF to the index of the ADO

# The formatting of the Hashmaps are as follows:
# I2ind[I] : int I -> list of ints corresponding to the BCF indecies
# ind2I[ind] :tuple of ints -> int I
# This is done because lists are not hashable, and tuples are

# number of ADOs in the n'th tier with k elements
def length(n,k):
    return np.math.factorial(n+k-1) // np.math.factorial(k-1) // np.math.factorial(n) # = [ (n+k-1) C (k-1)]
    # return np.math.factorial(n) // (np.math.factorial(k_)* np.math.factorial(n-k_))

# total number of ADOs for a given K and L
def total_length(K,L,N_nonmats):
    Ktot=K+N_nonmats
    return sum([length(n,Ktot) for n in range(0,L+1)])
    


def generateHashmap(K,L,N_nonmats,write_to_file = False):
    Ktot=K+N_nonmats
    # K here is the number of matsubara exponentials in the BCF
    # N_nonmats is the number of non-matsubara terms in the BCF
    # Ktot is the total number of exponential terms in the BCF (Ktot = K + N_nonmats)
    # L is the maximum depth of the ADO expansion

    # generates the set of indicies corresponding to the n'th tier with K elements
    def generatenumbers(n, k): 

        # output is a list of the strings of the indicies separated by commas
        numbers = np.empty(length(n,k),dtype=object) 
        # numbers = ['xxx' for i in range(length(n,k))]

        if k==1: #if there is only one element in the set
            return [f'{n}']

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

    # Collect all of the indices
    allnums = np.array([])
    for i in range(0,L+1):
        allnums= np.concatenate((allnums,generatenumbers(i,Ktot))) #concatenate the set of indicies to the list of all indicies

    def tup2list(tup):
        return [int(i) for i in tup.split(',')]

    # Create the hashmaps and format them as described above
    I_to_index = {I:tup2list(index) for I,index in enumerate(allnums) } #hash map from the index of the ADO to the index of the BCF
    index_to_I = {index:I for I,index in enumerate(allnums) } #hash map from the index of the BCF to the index of the ADO

    if write_to_file:
        with open('hashmap.txt','w') as f:
            f.write(str(I_to_index))
            f.write('\n')
            f.write(str(index_to_I))

    return I_to_index, index_to_I

### Convert the hashmaps to a list - FOR FORTRAN
# also find the first indices of each tier for faster indexing
def Convert_to_list(I_to_index):
    I_to_index_list = np.array([I_to_index[i] for i in range(len(I_to_index))])         # List of each ADO index
    tiers = np.array([np.sum(I_to_index_list[i]) for i in range(len(I_to_index_list))]) # Tier of each ADO index
    I0s =[0]
    tiernow = 0
    for index, tier in enumerate(tiers):
        if tier != tiernow:
            tiernow = tier
            I0s.append(index)
    I0s.append(len(I_to_index_list))
    I0s = np.array(I0s)
    return I_to_index_list, I0s

if __name__ == '__main__':

    # INPUTS
    K = 3   #number of exponential terms in the BCF - either the truncation of the number of beads + 1
    L = 3   #maximum depth of the ADO expansion
    N_nonmats = 1

    # Get hashmaps
    I_to_index, index_to_I = generateHashmap(K,L,N_nonmats)   
    print('\n')
    ADO_index=Convert_to_list(I_to_index) # this list has dimensions [I, K] (where K is the number of exponential terms in the BCF)
    print(ADO_index)
    print('\n')

    tier=2
    n=2
    pm=-1
    I0s = np.array([np.sum([length(n,K+1) for n in range(0,tier)]) for tier in range(0,L+1)],dtype=int)
    I0 = I0s[tier]
    I1 = I0s[tier+1]
    print(I0,I1)
    print(ADO_index[I0:I1,:])
    sys.exit()

    print(lst[1,:])  

