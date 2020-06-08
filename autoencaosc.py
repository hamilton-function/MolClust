import hyperpipeline as hy
import encodecode as ende
import numpy as np

'''
Setup for ORCLUS in autoencoder setting as submitted in PKDD:
CAUTION: k,l as hyperparameters of ORCLUS, l represents the dimensionality of the subspace
The subspace definition in ORCLUS is, that we compute the l-WEAKEST eigenvectors based on their corr. eigenvalues
which means for the autoencoder-fashioned approach: l --> d-l dimensions to which the objects shall be projected to!
'''


def generate_orclus_hypergrid(res,maxk,maxdim):
    #modified for explicit k=16!
    k_range = [n for n in range(4, maxk+1) if n % res == 0]
    #for exp. 2 only, ...set back to original:
    dimrange = [n for n in range(275, maxdim + 1) if n % res == 0]
    #dimrange = [n for n in range(290, maxdim+1) if n % res == 0]
    #dimrange = [n for n in range(maxdim, maxdim+1) if n % res == 0]
    orcgrid = [(x,y) for x in k_range for y in dimrange]
    return orcgrid




'''
routine for...
1. Applying ORCLUS on 
2. a hyperparameter grid (k,dim)
3. evaluating the loss
(4. storing the best setting and clustering (according to loss) via pickle)

-repeat for each hyperparameter setting 10 (n-)times since issue with centroid initialization
'''

def scan_orclus_loss_on_grid(res, dim, maxk, maxdim, path2elki, path2data, path2out, delta, gamma):
    orcgrid = generate_orclus_hypergrid(res, maxk, maxdim)
    besthyperset = []
    bestloss = np.inf
    bestpwlinearsubspace = []
    bestdim = np.inf
    bestclusnum = np.inf

    itercount = 0
    for e in orcgrid:
        print('----------- ITER ',itercount,'/',len(orcgrid),'----------')
        itercount += 1

        k, l = e
        hy.elki_run_orclus(path2elki, path2data, path2out, k, l)
        clustering = hy.get_all_results(path2out)

        if(clustering != 'NORESULT'):
            #now encode
            original = []
            encoded = []
            decoded = []
            pwlinearsubspace = []
            for key in clustering.keys():
                if len(clustering[key]) != 0:
                    clusterdata = clustering[key]

                    #TODO: should be in case of ORCLUS d-l as argument, since routine chooses d-l vectors with largest eigenvalue
                    projdim = dim-l
                    Crec, E, pcaobj = ende.compute_subspace(clusterdata, projdim)


                    #computed subspace
                    pwlinearsubspace.append(E)
                    #reconstructed data
                    decoded.append(Crec)
                    #original data
                    original.append(clusterdata)

            projdim = dim - l
            #compute loss
            loss = ende.compute_loss(original, decoded, delta, gamma, projdim)

            if (loss == bestloss) and (projdim <= bestdim) and (k <= bestclusnum):
                bestloss = loss
                besthyperset = [k,l]
                bestpwlinearsubspace = pwlinearsubspace
                bestdim = projdim
                bestclusnum = k

            #check if current clustering surpasses previous results
            if(loss < bestloss):
                bestloss = loss
                besthyperset = [k,l]
                bestpwlinearsubspace = pwlinearsubspace
                bestdim = projdim
                bestclusnum = k

            print('k: ',k,' l: ',l,' current loss: ',loss)

    print('Best loss: ',bestloss,' best hyperparameter: ', besthyperset,' best piecewise-linear-subspace: ',bestpwlinearsubspace)
    return (bestloss, bestdim, bestclusnum, besthyperset, bestpwlinearsubspace)



elki = "E:\\DesktopBackup\\elki-0.7.1\\elki-bundle-0.7.1.jar"
datain = "E:\\MolClust\\cov2vect.csv"
resoutORC = "E:\\DesktopBackup\\elki-0.7.1\\ORC_resrec"

#scan_orclus_loss_on_grid(res, dim, maxk, maxdim, path2elki, path2data, path2out, delta, gamma)
#scan_orclus_loss_on_grid(2, 300, 10, 299, elki, datain, resoutORC, 0.5, 0.5)

#15.04.2020, 15:15 till 15:31 15:46 16:02 16:17 16:32
#we have dim >> n := num. of samples --> can not exploit full range of dimensions! , 100/10 ~ 10 objects per cluster --> maxdim = 10


#scan_orclus_loss_on_grid(2, 300, 6, 292, elki, datain, resoutORC, 0.5, 0.5)
#scan_orclus_loss_on_grid(2, 300, 2, 296, elki, datain, resoutORC, 0.5, 0.5)

#ATTENTION, modified generate_orclus_hypergrid(res,maxk,maxdim) k beginning from 16!
scan_orclus_loss_on_grid(1, 300, 4, 275, elki, datain, resoutORC, 0.5, 0.5)

#scan_orclus_loss_on_grid(1, 300, 2, 296, elki, datain, resoutORC, 0.5, 0.5)

'''
27.05.20
Interesting observations: 
(l,k)
(4,2) DONE
(4,4) DONE
(8,8) DONE
(8,10) DONE
(10,10) DONE
(4,16) DONE

compute spectral clustering -->  DONE

next experiment: k=4, l ranges from 5 to 25, in steps of 5
do we observe differences with same k but different number of subspaces? if yes: 
can we observe certain special characteristics learned with different l? --> DONE


compute median per cluster in k-means and slink --> DONE


maybe compute medoid in subspaces? with that: compute k-closest objects to medoid as orbiting objects
--> augments in visually recognizing the subspace? For that proceed as follows:
1. project objects of a cluster to their respective subspace
2. compute medoid of projected objects and return the CID
3. compute k-closest objects to medoid and return their CIDs
4. profit!
'''