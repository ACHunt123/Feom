
import numpy as np
import sys
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
    def generatenumbersOLD(n, k): 

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

    def generatenumbers(n, k): #faster version
        numbers = np.empty(length(n,k),dtype=object) 
        index=np.zeros(k,dtype=int) # temporary array to store the index
        index[0]=n      # set the first element to n (largest number)
        numbers[0]= ','.join(map(str, index))
        if n==0 or k==1:  return numbers # exit if n=0 as [x] is the answer
        for I in range(1,length(n,k)):
            nend = index[-1]
            locs=np.where(index!=0)[0] # find the locations of the non-zero elements
            index[np.max(locs)]-=1  # decrement the right most non-zero element
            if np.max(locs)==Ktot-1: #if the right most element is the last element

                # remove the index from the list of indicies
                locs=np.delete(locs,np.argmax(locs))
                if locs.size==0: print('ERROR: locs is empty')

                index[Ktot-1]=0
                index[np.max(locs)]-=1
                index[np.max(locs)+1]=nend+1

            else: 
                index[np.max(locs)+1]+=1
            numbers[I]= ','.join(map(str, index))
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
    Ktot=N_nonmats+K

    # Get hashmaps
    I_to_index, index_to_I = generateHashmap(K,L,N_nonmats) 
    sys.exit()  
    print('\n')
    ADO_index,I0s=Convert_to_list(I_to_index) # this list has dimensions [I, K] (where K is the number of exponential terms in the BCF)
    print(ADO_index)
    print(I0s)

    # ADO_tier=ADO_index[I0s[tier]:I0s[tier+1]]
    # print(ADO_tier)

    ## now to make the algo
    index = [1,0,0,2,0,0]
    # starting min and max indices
    tier= np.sum(index)
    I0=I0s[tier]
    sn=0 # Running total of indices that have been found so far (left to right)
    Lengths = np.zeros((tier+1,Ktot),dtype=int)
    for n in range(0,tier+1):
        for k in range(1,Ktot):
            Lengths[n,k]=length(n,k)

    for p in range(1,Ktot): # Loop over the digit to focus on [x,.,.,.,.,.] then [.,x,.,.,.,.] etc.
        ni = index[p-1]     # The number of the digit
        sn+=ni              # Add this to the running total
        for n in range(0,tier-sn):
            I0+=Lengths[n,Ktot-p] # Move to the first instance where the leading digit is x = ni

        if sn==tier:        # If the running total of all the digits is equal to the tier, then we have found the correct index
            break
    print('\n')
    print(ADO_index[I0:I0+1],'final')





    sys.exit()



