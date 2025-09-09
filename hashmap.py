
import numpy as np
import sys,math
import matplotlib.pyplot as plt
'''    +---------------------------------------+
       |   Hashmap: Generates the hashmaps     |  
       |   for the ADOs in the FEOM            |
       |           By A. C. Hunt 2025          |
       +---------------------------------------+
'''


def length(n,k): # number of ADOs in the n'th tier with k elements
    # return np.math.factorial(n+k-1) // np.math.factorial(k-1) // np.math.factorial(n) # = [ (n+k-1) C (k-1)]
    return math.factorial(n+k-1) // math.factorial(k-1) // math.factorial(n) # = [ (n+k-1) C (k-1)]

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


    return I2ind, I0s

if __name__ == '__main__':

    ### Generate figure showing the number of ADOs in each tier
    Kmax=10
    Lmax=10
    data= np.zeros((Lmax+1,Kmax+1),dtype=int)
    for K in range(0,Kmax+1):
        for L in range(0,Lmax+1):
            data[L,K]=total_length(K,L,1)   

        
    # Bold numeric labels only
    col_labels = [fr'$\mathbf{{{k}}}$' for k in range(Kmax + 1)]
    row_labels = [fr'    $\mathbf{{{l}}}\!\!$' for l in range(Lmax + 1)]

    # Create plot
    scale=2
    xx=6 * scale
    yy=4 *scale
    fsz= (xx, yy)
    fig, ax = plt.subplots(figsize=fsz)
    # ax.axis('tight')
    ax.axis('off')
    # ax.set_aspect('equal')

    # Create table
    table = ax.table(cellText=data,
                    colLabels=col_labels,
                    rowLabels=row_labels,
                    cellLoc='center',
                    loc='center')

    # Adjust table appearance

    table.scale(1, 2.1)
    for key, cell in table.get_celld().items():
        cell.set_fontsize(20)

    # Add axis labels manually using text
    # Positioning may need adjustment depending on scale
    ax.text(0.5, 0.875, 'K', transform=ax.transAxes,
            ha='center', va='center', fontsize=16, fontweight='bold')

    ax.text(-0.12, 0.5, 'L', transform=ax.transAxes,
            ha='center', va='center', fontsize=16, fontweight='bold')

    plt.subplots_adjust(top=1.025, bottom=0.05)


    fig.suptitle("Number of ADOs for Debye bath HEOM", fontsize=18)
    fig.savefig("/home/ach221/data/Feom/telluride2025/figures/ADO_count_debye.pdf")

    plt.show()



    sys.exit()



