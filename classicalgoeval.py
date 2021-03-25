'''
in this setting we check for the best clustering using k-means from sklearn
the best result is determined by computing the silhouette coefficient for each k
the emerging 'best' clustering has to be checked w.r.t. the semantics (what has been clustered? structure, target etc.)
'''

from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import SpectralClustering
from sklearn import metrics
from sklearn.decomposition import PCA

import numpy as np
import csv
import pickle
import pubchempy as pcp
import glob
import os


# read from csv
def read_csv(csvfile):
    X = []
    with open(csvfile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
        for row in csv_reader:
            X.append(row)
    return X


'''
routine that compute a k-means clustering + silh coefficient, default: k-means++ initialization
each run is by default repeated 10 times and the best is taken in terms of inertia
'''


def get_best_clustering_by_silhouette(X, krange):
    best_silhouette = -np.inf
    bestk = -np.inf
    bestcentroids = []
    best_labels = []

    for k in krange:
        kmeans_model = KMeans(n_clusters=k).fit(X)
        labels = kmeans_model.labels_
        silscore = metrics.silhouette_score(X, labels, metric='euclidean')
        if silscore >= best_silhouette:
            best_silhouette = silscore
            bestk = k
            bestcentroids = kmeans_model.cluster_centers_
            best_labels = kmeans_model.labels_

    return bestk, best_silhouette, bestcentroids, best_labels


def get_best_clustering_by_silhouette_slink(X, krange):
    best_silhouette = -np.inf
    bestk = -np.inf
    best_labels = []

    for k in krange:
        slink_model = AgglomerativeClustering(n_clusters=k).fit(X)
        labels = slink_model.labels_
        silscore = metrics.silhouette_score(X, labels, metric='euclidean')
        if silscore >= best_silhouette:
            best_silhouette = silscore
            bestk = k
            best_labels = slink_model.labels_

    return bestk, best_silhouette, 'nocentroids', best_labels


def get_best_clustering_by_silhouette_spectral(X, krange):
    best_silhouette = -np.inf
    bestk = -np.inf
    best_labels = []

    for k in krange:
        spectral_model = SpectralClustering(n_clusters=k, assign_labels="discretize", random_state=0).fit(X)
        labels = spectral_model.labels_
        silscore = metrics.silhouette_score(X, labels, metric='euclidean')
        if silscore >= best_silhouette:
            best_silhouette = silscore
            bestk = k
            best_labels = spectral_model.labels_

    return bestk, best_silhouette, 'nocentroids', best_labels


'''
routine to fetch names for a given vector
'''


def vect2namedesc(vect):
    # load mapping
    infile = open('covmapping.pkl', 'rb')
    mapping = pickle.load(infile)
    infile.close()

    queryres = []

    queryvec = np.asarray(vect)

    for key in mapping:
        keyvec = np.asarray(((mapping[key])[1])[0])
        if np.array_equal(keyvec, queryvec):
            queryres = (mapping[key])[0]
            break

    return queryres


'''
routine to fetch cid for a given vector
'''


def vect2cid(vect):
    # load mapping
    infile = open('covmapping.pkl', 'rb')
    mapping = pickle.load(infile)
    infile.close()

    queryres = []

    queryvec = np.asarray(vect)

    for key in mapping:
        keyvec = np.asarray(((mapping[key])[1])[0])
        if np.array_equal(keyvec, queryvec):
            queryres = key
            break

    return queryres


'''
routine to fetch vector for a given cid
'''


def cid2vect(cid):
    # load mapping
    infile = open('covmapping.pkl', 'rb')
    mapping = pickle.load(infile)
    infile.close()

    ''' 
    #interesting case where two compounds have the same vectorization but different stereochemistry
    #this is however not possible to be captured by 2D SDF information and would require 3D SDFs which are not
    #generally available due to multiple degrees of freedom

    v1 = ((mapping[60749 ])[1])[0]
    v2 = ((mapping[57336515])[1])[0]
    eq = np.array_equal(v1, v2)
    print('EQ: ',eq)
    '''

    print(cid)

    return ((mapping[cid])[1])[0]


'''
routine to fetch mol structure as png files and store them in a given folder
'''


def get_structure_image(compname, targetdir, clusterid):
    molres = pcp.get_cids(compname, 'name', 'substance', list_return='flat')
    if len(molres) == 0:
        print('NO COMPOUND FOUND FOR ', compname)
    else:
        toprelcid = molres[0]
        pcp.download('PNG', targetdir + str(clusterid) + '_' + compname + '.png', str(toprelcid), 'cid', overwrite=True)


'''
routine to fetch mol structure  BY CID as png files and store them in a given folder
'''


def get_structure_image_by_cid(cid, compname, targetdir, clusterid):
    # hack to avoid filenames above 255 chars...FUCK YOU WINDOWS, for limiting in the year 2020 to 255...
    truncatedcomp = compname[0:245]
    pcp.download('PNG', targetdir + str(clusterid) + '_' + str(cid) + '.png', str(cid), 'cid', overwrite=True)


'''
get structure images of clustering
'''


def get_structure_image_clustering(clusterdict, targetdir):
    for key in clusterdict:
        print(clusterdict[key])
        for e in clusterdict[key]:
            cname = e[0]
            cid = e[1]
            get_structure_image_by_cid(cid, cname, targetdir, key)
    print('DONE - STRUCTURAL IMAGES RETRIEVED.')


'''
routine to compute the medoid of a given cluster
input: set of vectors
output: medoid vector
'''
from scipy.spatial import distance_matrix


def get_medoid_compound(cluster):
    # compute distance matrix
    dmx = distance_matrix(cluster, cluster)

    # obtain index of medoid vector
    medidx = np.argmin(dmx.sum(axis=0))
    # obtain vector with index medidx
    medvect = cluster[medidx]

    # obtain compound name and corresponding
    idnamedesc = vect2namedesc(medvect)[0]
    cid = vect2cid(medvect)

    # obtain compound image
    return [cid, idnamedesc, medvect]


def extract_clusters_from_clusdirectory(clusdir):
    # retrieve all files and therefore compounds from directory
    allcomp = glob.glob(clusdir + '*.png')

    cidarray = []

    # processing to get CID
    # X_12345.png will be first split to X, 12345.png, then 12345.png will be split to 12345 and png.
    for e in allcomp:
        cid = int((e.split('_')[1]).split('.')[0])
        cidarray.append(cid)

    vecarr = []
    # retrieve 300d embedding vectors for each cid
    for e in cidarray:
        vec = cid2vect(e)
        vecarr.append(vec)

    # now compute the medoid with existing routine
    medoidres = get_medoid_compound(vecarr)

    return medoidres


'''
Routine to get for a clustering directory the medoids for all clusters
returned as a dictionary with key: cluster id, and value, the medoid cid and name
'''


def get_all_cluster_medoids(clusteringdir):
    medoiddict = {}
    direct = [name for name in os.listdir(clusteringdir) if os.path.isdir(os.path.join(clusteringdir, name))]
    for e in direct:
        r = extract_clusters_from_clusdirectory(clusteringdir + e + '\\')[0:2]
        medoiddict[e] = r
    return medoiddict


'''
Routine to obtain medoid of subspace + k orbits
'''


def get_all_subspace_cluster_medoids(clusteringdir, l, korbit):
    medoiddict = {}
    direct = [name for name in os.listdir(clusteringdir) if os.path.isdir(os.path.join(clusteringdir, name))]
    for e in direct:
        r = get_subspace_medoid_and_orbit(clusteringdir + e + '\\', l, korbit)
        medoiddict[e] = r
    return medoiddict


'''
Auxillary routine to obtain medoid of subsapce
'''


def get_subspace_medoid_compound(cluster, mapping):
    # compute distance matrix
    dmx = distance_matrix(cluster, cluster)

    # obtain index of medoid vector
    medidx = np.argmin(dmx.sum(axis=0))
    # obtain vector with index medidx
    medvect = cluster[medidx]

    # obtain compound name and corresponding
    idnamedesc = ''
    cidtarget = np.inf
    for e in mapping:
        cid, fullvec, subspacevec = e
        if np.array_equal(subspacevec, medvect):
            idnamedesc = vect2namedesc(fullvec)[0]
            cidtarget = cid

    # obtain compound image
    return [cidtarget, idnamedesc, medvect]


'''
Subroutine to obtain kNN of a given medoid compound via cid
'''
from scipy.spatial import distance


def get_medoid_kNN(korbit, cluster, medoid, mapping):
    # compute distances of all objects in cluster against the medoid
    medoiddistarr = (distance.cdist(cluster, [medoid]))

    flatmedoiddistarr = [item for sublist in medoiddistarr for item in sublist]

    korbitwithmedoid = korbit  # +1

    # obtain index of k-smallest distances to medoid; korbit+1 is due to the fact that d(medoid,medoid)=0
    lowestkidx = np.argpartition(flatmedoiddistarr, korbitwithmedoid)

    korbitarr = []

    # retrieve subspace vectors --> vectors --> cid --> compound name
    for e in lowestkidx[:korbitwithmedoid + 1]:
        # get subspace vector
        subsvec = cluster[e]
        cid = (mapping[e])[0]
        fullvec = (mapping[e])[1]
        cname = vect2namedesc(fullvec)[0]
        meddist = medoiddistarr[e]

        # uncommented line below contains the lower dimensional subspace cluster vect
        # korbitarr.append([cid,cname,meddist,subsvec])
        korbitarr.append([cid, cname, meddist])

    return korbitarr


'''
Routine to get medoid + orbiting of a subspace cluster
input: directory with png of compounds, with their cids in their names
output: medoid + orbiting objects (= k closest objects to medoid) within its subspace
Approach:
1. get cid names and retrieve their respective vectors
2. compute the subspace of dimensionality l and project objects down to this subspace
3. compute for lower-dimensional objects their medoid [and their orbiting objects]
'''


def get_subspace_medoid_and_orbit(clusdir, l, korbit):
    '''
    Part 1
    '''
    # retrieve all files and therefore compounds from directory
    allcomp = glob.glob(clusdir + '*.png')

    cidarray = []

    # processing to get CID
    # X_12345.png will be first split to X, 12345.png, then 12345.png will be split to 12345 and png.
    for e in allcomp:
        cid = int((e.split('_')[1]).split('.')[0])
        cidarray.append(cid)

    vecarr = []
    # retrieve 300d embedding vectors for each cid
    for e in cidarray:
        vec = cid2vect(e)
        vecarr.append(vec)

    '''
    Part 2
    '''
    ncomp = min(l, len(vecarr))
    pca = PCA(n_components=ncomp, svd_solver='full')
    vecarr_proj = pca.fit_transform(vecarr)
    # print(vecarr_proj)

    subspacemapping = zip(cidarray, vecarr, vecarr_proj)
    materializedsubspacemapping = list(subspacemapping)

    '''
    Part 3
    '''
    medoidres = get_subspace_medoid_compound(vecarr_proj, materializedsubspacemapping)

    medsubvec = medoidres[2]
    korbitres = get_medoid_kNN(korbit, vecarr_proj, medsubvec, materializedsubspacemapping)
    # for e in korbitres:
    #    print(e)
    # print('-----')

    return korbitres


'''
auxillary function:
Project and compute distance
1. compute l-dim. subspace of a given cluster
2. project compound vectors to l-dim. subspace
3. project query vector to l-dim. subspace
4. compute distance between two objects in l-dim. subspace
5. profit
'''


def proj_comp_dist_to_query(clusdir, l, p, q):
    '''
    Part 1
    '''
    # retrieve all files and therefore compounds from directory
    allcomp = glob.glob(clusdir + '*.png')

    cidarray = []

    # processing to get CID
    # X_12345.png will be first split to X, 12345.png, then 12345.png will be split to 12345 and png.
    for e in allcomp:
        cid = int((e.split('_')[2]).split('.')[0])
        cidarray.append(cid)

    vecarr = []
    # retrieve 300d embedding vectors for each cid
    for e in cidarray:
        vec = cid2vect(e)
        vecarr.append(vec)

    '''
    Part 2
    '''
    ncomp = min(l, len(vecarr))
    pca = PCA(n_components=ncomp, svd_solver='full')
    vecarr_proj = pca.fit_transform(vecarr)
    # print(vecarr_proj)

    ptrans = pca.transform([p])
    qtrans = pca.transform([q])
    subspacedist = np.linalg.norm(ptrans - qtrans)

    return subspacedist


'''
Visualization routine:
- computes for a given cluster c_i the medoid and ALL distances to the medoid, yielding triplets in form of: [cid, cname, meddist]
- the list of the triplets is sorted in ascending order by distance to medoid (including the medoid itself)
- either on the horizontal axis or at each 'dot' in the connected scatter plot we see the CID
'''
import matplotlib.pyplot as plt


def visualize_allNN_dist(clusdir, l):
    allcomp = os.listdir(clusdir)
    dirsize = len(allcomp) - 1
    # obtain distances of all objects of the cluster to the cluster medoid
    medoiddist = get_subspace_medoid_and_orbit(clusdir, l, dirsize)

    # sort compounds by distance in ascending order
    medoiddistsorted = sorted(medoiddist, key=lambda x: x[2])

    xlis = [i for i in range(0, len(medoiddistsorted))]
    ylis = []
    llis = []

    for e in medoiddistsorted:
        dist = (e[2])[0]
        cid = e[0]
        ylis.append(dist)
        llis.append(cid)

    fig, ax = plt.subplots()

    ax.plot(xlis, ylis, '-o')
    ax.set_ylim(bottom=-40)
    ax.yaxis.set_ticks(np.arange(min(ylis), max(ylis), 10))

    for i, txt in enumerate(llis):
        ax.annotate(txt, (xlis[i], ylis[i] - 10), rotation=90)

    ax.set_ylabel('Distance in subspace to medoid')
    ax.set_xlabel('kNN')

    # set title, composed of kXlY and cI
    figurelabel = clusdir.split('\\')[-3] + ' ' + clusdir.split('\\')[-2]
    ax.set_title(figurelabel)

    plt.savefig(clusdir + figurelabel + '.png', bbox_inches='tight', dpi=300)
    # plt.show()


'''
Routine to obtain kNN distance plots
'''


def get_all_subspace_cluster_kNNdist_plots(clusteringdir, l):
    medoiddict = {}
    direct = [name for name in os.listdir(clusteringdir) if os.path.isdir(os.path.join(clusteringdir, name))]
    for e in direct:
        r = visualize_allNN_dist(clusteringdir + e + '\\', l)
    return None


'''
takes a directory (=cluster)
returns to a user defined target directory the clusters sorted
'''


def dir_sorted_allKNN(clusdir, l):
    allcomp = os.listdir(clusdir)
    dirsize = len(allcomp) - 1
    # obtain distances of all objects of the cluster to the cluster medoid
    medoiddist = get_subspace_medoid_and_orbit(clusdir, l, dirsize)

    # sort compounds by distance in ascending order
    medoiddistsorted = sorted(medoiddist, key=lambda x: x[2])

    xlis = [i for i in range(0, len(medoiddistsorted))]
    ylis = []
    llis = []

    for e in medoiddistsorted:
        dist = (e[2])[0]
        cid = e[0]
        ylis.append(dist)
        llis.append(cid)

    cluscomp = glob.glob(clusdir + '*.png')

    # processing to get CID
    # X_12345.png will be first split to X, 12345.png, then 12345.png will be split to 12345 and png.
    for e in cluscomp:
        cid = int((e.split('_')[1]).split('.')[0])
        # get corresponding rank (=index) in llis
        idx = llis.index(cid)

        # print(e)

        # old filename only
        oldfilename = e

        # new filename only, by appending 'RANK_oldfilename from idx variable
        newfilename = e[0:e.rfind('\\')] + '\\' + str(idx) + '_' + e[e.rfind('\\') + 1:]
        print(newfilename)

        os.rename(oldfilename, newfilename)

    print('ALL RENAMED')


'''
takes a directory of directories (=clusterING)
returns to a user defined target directory ALL clusters sorted
'''


def all_dir_sorted_allKNN(clusteringdir, l):
    direct = [name for name in os.listdir(clusteringdir) if os.path.isdir(os.path.join(clusteringdir, name))]
    for e in direct:
        r = dir_sorted_allKNN(clusteringdir + e + '\\', l)
    return None


'''
routine to compute clustermap of compounds -within- a l-dimensional projection (subspace)
Steps:
1. obtain compounds from directory
2. project compounds to l-dimensional representations
3. compute clustermap (seaborn?)
'''
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd
import pandas as pd
import seaborn as sns


def compute_pairwise_subspacedist_clustermap(clusdir, l):
    '''
    Part 1
    '''
    # retrieve all files and therefore compounds from directory
    allcomp = glob.glob(clusdir + '*.png')

    cidarray = []

    # processing to get CID
    # X_12345.png will be first split to X, 12345.png, then 12345.png will be split to 12345 and png.
    for e in allcomp:
        # set ot [1] instead of [0] since we have now ordering per cluster folder

        cid = int((e.split('_')[2]).split('.')[0])
        cidarray.append(cid)

    vecarr = []
    # retrieve 300d embedding vectors for each cid
    for e in cidarray:
        vec = cid2vect(e)
        vecarr.append(vec)

    '''
    Part 2
    '''
    ncomp = min(l, len(vecarr))
    pca = PCA(n_components=ncomp, svd_solver='full')
    vecarr_proj = pca.fit_transform(vecarr)
    # print(vecarr_proj)

    subspacemapping = zip(cidarray, vecarr_proj)
    materializedsubspacemapping = list(subspacemapping)

    '''
    Part 3
    '''
    dmx = distance_matrix(vecarr_proj, vecarr_proj)
    dmx_df = pd.DataFrame(dmx)

    idxdict = dict(zip([i for i in range(len(cidarray))], cidarray))

    dmx_df.rename(index=idxdict, inplace=True)
    dmx_df.rename(columns=idxdict, inplace=True)

    cmp = sns.clustermap(dmx_df, xticklabels=True, yticklabels=True)
    # plt.show()

    figurelabel = clusdir.split('\\')[-3] + ' ' + clusdir.split('\\')[-2]
    cmp.savefig(clusdir + figurelabel + '.png')

    return None


'''
takes a directory of directories (=clusterING)
returns to a user defined target directory ALL clusters their clustermap
'''


def all_compute_pairwise_subspacedist_clustermap(clusteringdir, l):
    medoiddict = {}
    direct = [name for name in os.listdir(clusteringdir) if os.path.isdir(os.path.join(clusteringdir, name))]
    for e in direct:
        r = compute_pairwise_subspacedist_clustermap(clusteringdir + e + '\\', l)
    return None



'''
Code fragment to plot hierarchy of compounds
'''
# input array:
covdata = read_csv('E:\MolClust\cov2vect.csv')
krangelis = [n for n in range(24, 25) if n % 2 == 0]

# get arary of labels/compound names
compnamelabels = []
for comp in covdata:
    compname = vect2cid(comp)
    compnamelabels.append(compname)

print(compnamelabels)

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering


def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram ; 120
    dendrogram(linkage_matrix, color_threshold=110, labels=compnamelabels, **kwargs)


# setting distance_threshold=0 ensures we compute the full tree.
model = AgglomerativeClustering(distance_threshold=0, n_clusters=None)

model = model.fit(covdata)
plt.title('Hierarchical Clustering of Mol2Vec Embedded Compounds')
# plot the top three levels of the dendrogram
plot_dendrogram(model)
plt.xlabel("CID")
plt.ylabel("Linkage Distance")
plt.show()



'''
Code fragment to compute k-means clustering, project compounds down to 2d representation using t-SNE,
and plotting the 2d-clustering

color schema from xkcd: https://xkcd.com/color/rgb/
'''

from sklearn import cluster, datasets
from sklearn.manifold import TSNE
from matplotlib import pyplot as plt

k=16

k_means = cluster.KMeans(n_clusters=k)
k_means.fit(covdata)

#rs=42
tsne = TSNE(n_components=2, random_state=42)
X_2d = tsne.fit_transform(covdata)


target_ids = range(k)
plt.figure(figsize=(6, 5))
colors = 'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'xkcd:sky blue', 'xkcd:spring green'
for i, c in zip(target_ids, colors):
    plt.scatter(X_2d[k_means.labels_ == i, 0], X_2d[k_means.labels_ == i, 1], c=c)
plt.show()


'''
#switch here for slink vs. kmeans
bestres = get_best_clustering_by_silhouette(covdata, krangelis)
#print('Best k: ',bestres[0],'\n best silh.score: ',bestres[1],'\n best centroids ',bestres[2],'\n best labels: ',bestres[3])
#collect vectors by cluster label
clusresdict = {}
for e in zip(covdata,bestres[3]):
    vec, label = e
    if label not in clusresdict.keys():
        clusresdict[label] = [vec]
    else:
        clusresdict[label].append(vec)
'''
# testrun for prototype
# c1 = clusresdict[1]
# mdc = get_medoid_compound(c1)
# print('medoid: ',mdc)

# now fetch the names + descriptions:
# invec = covdata[1]
# print(vect2namedesc(invec))

'''
get compounds + description for all elements of each cluster:
'''

'''
clusrescomp = {}
for key in clusresdict:
    for vec in clusresdict[key]:
        if key not in clusrescomp.keys():
            namedesc = vect2namedesc(vec)[0]
            cid = vect2cid(vec)
            clusrescomp[key] = [[namedesc, cid]]
        else:
            namedesc = vect2namedesc(vec)[0]
            cid = vect2cid(vec)
            clusrescomp[key].append([namedesc,cid])
for key in clusrescomp:
    for e in clusrescomp[key]:
        print('KEY ',key,'  :  ',e)
get_structure_image_clustering(clusrescomp, 'E:\\MolClust\\expres\\kmeans\\k24\\')
'''

# meddic = get_all_cluster_medoids('E:\\MolClust\\expres\\kmeans\\k24\\')
# for key in meddic:
#    print(key,' medoid: ',meddic[key])


# test for single orclus subspace cluster
# r2 = get_all_subspace_cluster_medoids('E:\\MolClust\\expres\\k2l4\\', 4, 3)
# for key in r2:
#    print(key, ' : ', r2[key])


# t = visualize_allNN_dist('E:\\MolClust\\expres\\k4l10\\c1\\', 10)
# get_all_subspace_cluster_kNNdist_plots('E:\\MolClust\\expres\\k16l4\\',4)


# p = cid2vect(3005572)
# q = cid2vect(2733525)
# pqdist = proj_comp_dist_to_query('E:\\MolClust\\expres\\k10l8\\c9\\', 8, p, q)
# print('P-Q-DIST:',pqdist)


# dir_sorted_allKNN('E:\\MolClust\\expres\\k16l4\\c2\\', 4)
# all_dir_sorted_allKNN('E:\\MolClust\\expres\\k16l4\\',4)


# compute_pairwise_subspacedist_clustermap('E:\\MolClust\\expres\\k4l4\\c2\\',4)
# all_compute_pairwise_subspacedist_clustermap('E:\\MolClust\\expres\\k4l4\\',4)
