'''
Auxillary library with functions for
a) projecting objects down to their respective local (linear) subspace (top l-dimensions) : ENCODE
b) projecting the objects back to full dimensionality : DECODE
'''

from sklearn.decomposition import PCA
import numpy as np
import itertools

np.random.seed(42)

'''
auxillary routine to compute subspace E for a given cluster
the subspace is constructed by choosing those eigenvectors
which have the l-biggest corresponding eigenvalues 
output sorted by eigval descending order
'''
def compute_subspace(C, l):
    # the len(C[0]) is a 'hack' to perform PCA obtaining all components and then
    # choosing only those with the l-lowest singular values
    #12.01.20: set the line below from n_components=len(C[0]) to min(len(C[0]),len(C))
    ncomp = min(l,len(C))
    pca = PCA(n_components=ncomp, svd_solver='full')
    pcs = pca.fit_transform(C) #we encoded shit
    #now decoding shit
    pcd = pca.inverse_transform(pcs)

    return [pcd, pca.components_, pca]



'''
routine for encode for outlier:
'''
def compute_subspace_out(C, l):
    # the len(C[0]) is a 'hack' to perform PCA obtaining all components and then
    # choosing only those with the l-lowest singular values
    #12.01.20: set the line below from n_components=len(C[0]) to min(len(C[0]),len(C))
    ncomp = min(len(C[0]),len(C))
    pca = PCA(n_components=ncomp, svd_solver='full')
    pca.fit(C)
    pcs = pca.components_[:l]
    sval = pca.explained_variance_[:l]
    #print('SVAL: ',sval)
    #print('PCS: ',pcs)
    return [pcs, sval]


'''
function takes a cluster and a dimensionality of the subspace to be
projected to, and returns the objects projected to an
l-dimensional embedding, the l-eigenvectors and their corr. eigenvals
'''
def encode(C, l):
    E, evals = compute_subspace_out(C, l)
    Cproj = []
    for e in C:
        cp = np.dot(e, E.T)
        Cproj.append(cp)
    #print('Transformed:', Cproj)

    return [Cproj, E, evals]


'''
function takes downprojected objects from a cluster c and an eigenvector matrix
which has been previously used to downproject and returns the 
reconstructed data based on the top l-eigenvectors E
'''
def decode(Cpj, E, means):
    Crec = []
    for e in Cpj:
        cr = np.dot(e, E) + means
        Crec.append(cr)
    #print('Reconstructed: ', Crec)

    return Crec





'''
function that takes entire clustering and computes
for each cluster the reconstruction by further taking:
- penalty factor delta on the dimensionality --> the lower dim, the lower the penalty
- penalty factor gamma on the number of clusters --> the more clusters the higher the penalty

loss function is:
MSE(X,X^)+delta*dim_enc+gamma*clust

X and Xrec are list of lists (clusters), for X the clusters are full dimensional for Xrec the reconstructed full dim.
'''
from sklearn.metrics import mean_squared_error


def compute_loss_old(X,Xrec,delta,gamma,l):
    clusnum = len(X)
    #lines below to make from [[],[],[]] --> []
    Xflat = [] #list(itertools.chain.from_iterable(X))
    Xrecflat = [] #list(itertools.chain.from_iterable(Xrec))
    for clust in X:
        for e in clust:
            Xflat.append(e)

    for clust in Xrec:
        for e in clust:
            Xrecflat.append(e)
    MSE = mean_squared_error(Xflat,Xrecflat)
    #print('MSE : : ',MSE,' DELTA*DIM: ',delta*l)
    regloss = MSE + delta*l + gamma*clusnum
    return regloss


'''
MSE Loss invariant of dim and clustsize
'''

'''
auxillary subroutine to compute dimdist:
1/d sum_i=0 ^d (xî - x^i)^2 
'''
def compute_dimdist(x,xrec):
    sumoverdim = 0
    for e in zip(x,xrec):
        xi, xirec = e
        sumoverdim = sumoverdim + (xi-xirec)**2
        #print('xi, ',xi,' xirec ',xirec)
    return sumoverdim/len(x)

'''
auxillary subroutine to compute dimdist over all objects of a SINGLE CLUSTER
normalized by the CLUSTER SIZE!

1/|ci| sum_j=1^|ci| dimdist(xij,xirecj)
'''
def compute_clust_dimdist(ci, cirec):
    sumoverclust = 0
    for e in zip(ci,cirec):
        x, xrec = e
        sumoverclust = sumoverclust + compute_dimdist(x,xrec)
    return sumoverclust/len(ci)


#print(compute_clust_dimdist([(1,2,2),(3,1,3)],[(2,2,2),(4,1,3)]))

def compute_loss(X, Xrec, delta, gamma, l):
    clusnum = len(X)
    MSE = 0
    for clus in zip(X,Xrec):
        ci, cirec = clus
        MSE = MSE + compute_clust_dimdist(ci, cirec)
    #print('MSE : : ',MSE,' DELTA*DIM: ',delta*l)
    regloss = MSE + delta*l + gamma*clusnum
    return regloss



'''
data = [[1,1,1],[2,2,2],[3,3,3],[4,4,4]]
enc, E, evals, means = encode(data,1)
datarec = decode(enc, E, means)
print('------')
print(compute_loss_old(data, datarec,0,0,0))


c1 = [[1,1],[2,2],[3,3],[4,4]]
c2 = [[0.5,1],[1,2],[2,4]]

c1rec = [[1,0],[2,2],[3,1],[4,4]]
c2rec = [[1,1],[2,2],[2,4]]
print(compute_loss_old([c1,c2],[c1rec,c2rec],0,0,1))
'''