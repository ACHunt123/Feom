
import numpy as np
import sys
import matplotlib.pyplot as plt
'''    +---------------------------------------+
       |   Hashmap: Generates the hashmaps     |  
       |   for the ADOs in the FEOM            |
       |           By A. C. Hunt 2025          |
       +---------------------------------------+
'''


def length(n,k): # number of ADOs in the n'th tier with k elements
    return np.math.factorial(n+k-1) // np.math.factorial(k-1) // np.math.factorial(n) # = [ (n+k-1) C (k-1)]

def total_length(K,L,N_nonmats): # total number of ADOs for a given K and L
    Ktot=K+N_nonmats
    return sum([length(n,Ktot) for n in range(0,L+1)])
    
def generateHashmap(K,L,N_nonmats,write_to_file = False):
    ''' Generate the hashmaps for the ADOs in the FEOM
    INPUTS
    K           Number of Matsubara exponentials in the BCF
    N_nonmats   Number of non-matsubara terms in the BCF
    Ktot        Total number of exponential terms in the BCF (Ktot = K + N_nonmats)
    L           Maximum depth of the ADO expansion

    OUTPUTS
    I2ind       List of ADO indices, where each index is a list of integers corresponding to the BCF indices
    I0s         List of indices of the first ADO in each tier
    '''
    Ktot=K+N_nonmats

    def generatenumbers(n, k):  # Generate all of the numbers in the n-th tier with k elements
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

    # Collect all of the indices, going through each tier
    allnums = np.array([])      # this will hold all of the indices of the ADOs
    I0s =[0]                    # this will hold the index of the first ADO in the newest tier
    for i in range(0,L+1):
        allnums= np.concatenate((allnums,generatenumbers(i,Ktot))) # concatenate the set of indicies to the list of all indicies
        I0s.append(len(allnums))
    I0s = np.array(I0s) # convert to numpy array for easier indexing

    def tup2list(tup): # Re-format the tuple of ints into a list of ints
        return [int(i) for i in tup.split(',')]

    # Create the hashmaps and format them as described above
    I2ind = np.array([tup2list(index) for index in allnums])         # List of each ADO index

    if write_to_file:
        with open('hashmap.txt','w') as f:
            f.write(str(I_to_index))
            f.write('\n')
            f.write(str(index_to_I))

    return I2ind, I0s

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



