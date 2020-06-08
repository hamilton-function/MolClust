import subprocess
import numpy as np
import os
import shutil

'''
Pipeline STAGE 1:
applying ORCLUS
correlation clustering algorithm 
on a given dataset via python-script
'''

elki = "E:\\DesktopBackup\\elki-0.7.1\\elki-bundle-0.7.1.jar"
datain = "E:\\DesktopBackup\\elki-0.7.1\\ORCCAS_test.csv"
resout = "E:\\DesktopBackup\\elki-0.7.1\\ORCCAS_res"


'''
Function takes the following parameters and executes the CLIQUE algorithm in ELKI:
path2elki: path to the single all-in-one ELKI jar file, get it from: 
https://elki-project.github.io/releases/release0.7.1/elki-bundle-0.7.1.jar

path2data: full path to the csv (input) file
as the output it generates a folder with the cluster files

CAUTION: use subprocess.call instead of subprocess.popen, since call is blocking while popen is not, leading to
exceptions
'''
def elki_run_orclus(path2elki, path2data, path2out, k, l):
    if os.path.exists(path2out):
        shutil.rmtree(path2out)
    process = subprocess.call('java -jar '+path2elki+' KDDCLIApplication     -dbc.in '+path2data+'     -time     -algorithm clustering.correlation.ORCLUS     -projectedclustering.k '+str(k)+'     -projectedclustering.l '+str(l)+'   -resulthandler ResultWriter     -out '+path2out+'')
    print("Elki job done")


def elki_run_cash(path2elki, path2data, path2out, m, s, j):
    if os.path.exists(path2out):
        shutil.rmtree(path2out)
    #process = subprocess.call('java -jar '+path2elki+' KDDCLIApplication     -dbc.in '+path2data+'     -time     -algorithm clustering.correlation.CASH     -cash.minpts '+str(m)+'     -cash.maxlevel '+str(s)+'   -cash.mindim '+str(1)+'   -cash.jitter '+str(j)+'   -resulthandler ResultWriter     -out '+path2out+'')
    #ELKI-CASH seems to be choking on certain hyperparameter configurations - not reproducable so far. thus a temporal barrier is introduced to terminate the CASH run after 10 seconds
    try:
        process = subprocess.run('java -jar '+path2elki+' KDDCLIApplication     -dbc.in '+path2data+'     -time     -algorithm clustering.correlation.CASH     -cash.minpts '+str(m)+'     -cash.maxlevel '+str(s)+'   -cash.mindim '+str(1)+'   -cash.jitter '+str(j)+'   -resulthandler ResultWriter     -out '+path2out+'', timeout=60)
    except subprocess.TimeoutExpired:
        print('CASH run took too long - TERMINATED.')
    print("Elki job done")


def elki_run_copac(path2elki, path2data, path2out, knn, eps, minpts):
    if os.path.exists(path2out):
        shutil.rmtree(path2out)
    process = subprocess.call('java -jar '+path2elki+' KDDCLIApplication     -dbc.in '+path2data+'     -time     -algorithm clustering.correlation.COPAC     -copac.knn '+str(knn)+'     -dbscan.epsilon '+str(eps)+'   -dbscan.minpts '+str(minpts)+'   -resulthandler ResultWriter     -out '+path2out+'')
    print("Elki job done")


def elki_run_4c(path2elki, path2data, path2out, eps, minpts, delta, lambd):
    if os.path.exists(path2out):
        shutil.rmtree(path2out)
    process = subprocess.call('java -jar '+path2elki+' KDDCLIApplication     -dbc.in '+path2data+'     -time     -algorithm clustering.correlation.FourC     -dbscan.epsilon '+str(eps)+'     -dbscan.minpts '+str(minpts)+'   -pca.filter.absolute -pca.filter.delta '+str(delta)+'     -predecon.lambda '+str(lambd)+'   -resulthandler ResultWriter     -out '+path2out+'')
    print("Elki job done")


def elki_run_lmclus(path2elki, path2data, path2out, maxdim, minpts, thresh):
    if os.path.exists(path2out):
        shutil.rmtree(path2out)
    process = subprocess.call('java -jar '+path2elki+' KDDCLIApplication     -dbc.in '+path2data+'     -time     -algorithm clustering.correlation.LMCLUS     -lmclus.maxdim '+str(maxdim)+'     -lmclus.minsize '+str(minpts)+'   -lmclus.threshold '+str(thresh)+'   -resulthandler ResultWriter     -out '+path2out+'')
    print("Elki job done")


'''
Pipeline STAGE 2:
extract objects from each of the cluster_X.txt files
'''


def extract_res(cluster_i, clid):
    with open(cluster_i) as f:
        content = f.readlines()
        content = [x.strip() for x in content]

    clustermembers = []
    for e in content:
        # check if it is a line containing a coordinate, indicated with "ID="
        if ("ID=" in e):
            # split line by whitespace, discard first entry, this is the ID= stuff...
            coorvec = (np.array((e.split(" "))[1:]).astype(np.float)).tolist()
            clustermembers.append(coorvec)
    # print('coordinates: ',clustermembers)

    return clustermembers



'''
Auxillary function to extract the quantitative models from CASH clustering results
'''

def extract_model(cluster_i, clid):
    with open(cluster_i) as f:
        content = f.readlines()
        content = [x.strip() for x in content]


        #get dimensionality of original data:
        oridim = np.inf
        for e in content:
            if ("ID=" in e):
                # split line by whitespace, discard first entry, this is the ID= stuff...
                coorvec = (np.array((e.split(" "))[1:]).astype(np.float)).tolist()
                oridim = len(coorvec)
                break

        clustermodels = []
        modeldim = 0
        for e in content:
            #check if hashtag and x_0 is in line, if yes: it is a model line
            if("#" in e and "x_" in e):
                clustermodels.append(e)
                modeldim += 1
        # modeldim = 1 corresponds to d-1 dimensional subspace
        subspacedim = oridim - modeldim

    #print([modeldim,clustermodels])
    return [subspacedim,clustermodels]

#extract_model("E:\\DesktopBackup\\elki-0.7.1\\CASHTEST_REC\\cluster.txt",0)



'''
Auxillary function to extract the dimensions from COPAC clustering results
'''

def extract_dimension(cluster_i, clid):
    with open(cluster_i) as f:
        content = f.readlines()
        content = [x.strip() for x in content]

        clusterdim = 0
        for e in content:
            #check if hashtag and Dimension is in line, if yes: it is a dimensionality line
            if("# Dimension: " in e):

                #split by whitespace and get last item which is the dimensionality, append it as int
                dimlis = e.split()
                dim = int(dimlis[-1])
                clusterdim = dim
    return clusterdim



'''
Auxillary function to extract the dimensions from LMCLUS clustering results
which are strangely embedded in the filename like e.g. cluster_4d_1
'''

def extract_lmclus_dimension(cluster_i, clid):
    #slice filename from where 'cluster' starts until excl. .txt
    #first remove file ending .txt by index (endidx)
    #then split by underline and obtain the 1st position which includes the dimensionality Xd
    #then split by d and get the 0th position and cast it to integer to get dimensionality
    endidx = len(cluster_i)-4
    dim = int((((cluster_i[0:endidx]).split('_'))[1]).split('d')[0])
    return dim




'''
Function iterates over all files in the specified path2out (output folder of ELKI)
and returns a dictionary which assigns for each cluster file a cluster id as key
and as value the list with the equation of the linear function and data points belonging to the cluster
'''


def get_all_results(path2out):
    # all clusters are stored in a dictionary with the structure:
    # key: cluster_id, mean vector, subspaces, clustermembers
    clusdict = {}

    directory = os.fsencode(path2out)
    #print(directory)

    if os.path.exists(path2out):
        clid = 1
        noid = -1
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if "cluster" in filename:
                # print(filename)
                clusdict[clid] = extract_res(path2out + '\\' + filename, clid)
                clid += 1
            elif "noise" in filename:
                clusdict[noid] = extract_res(path2out + '\\' + filename, noid)
                noid -= 1

        return clusdict
    else:
        return "NORESULT"



'''
Auxillary function to extract all quantitative models from CASH run
'''

def get_all_models(path2out):
    # all clusters are stored in a dictionary with the structure:
    # key: cluster_id, mean vector, subspaces, clustermembers
    clusdict = {}

    directory = os.fsencode(path2out)

    if os.path.exists(path2out):
        clid = 1
        noid = -1
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if "cluster" in filename:
                clusdict[clid] = extract_model(path2out + '\\' + filename, clid)
                clid += 1
            elif "noise" in filename:
                clusdict[noid] = extract_model(path2out + '\\' + filename, noid)
                noid -= 1

        return clusdict
    else:
        return "NORESULT"


#testrun for cluster model extraction
#modelres = get_all_models("E:\\DesktopBackup\\elki-0.7.1\\CASHTEST_REC\\")
#for key in modelres.keys():
#    print(key,modelres[key])

#elki_run_orclus(elki, datain, resout, 3, 1)
#extract_res("E:\\DesktopBackup\\elki-0.7.1\\ORCCAS_res\\cluster_0.txt",0)


'''
Auxillary function to extract all quantitative models from COPAC run
'''

def get_all_dimensions(path2out):
    clusdict = {}

    directory = os.fsencode(path2out)

    if os.path.exists(path2out):
        clid = 1
        noid = -1
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if "cluster" in filename:
                clusdict[clid] = extract_dimension(path2out + '\\' + filename, clid)
                clid += 1
            elif "noise" in filename:
                clusdict[noid] = extract_dimension(path2out + '\\' + filename, noid)
                noid -= 1

        return clusdict
    else:
        return "NORESULT"


'''
Auxillary function to extract all quantitative models from LMCLUS run
'''

def get_all_lmclus_dimensions(path2out):
    clusdict = {}

    directory = os.fsencode(path2out)

    if os.path.exists(path2out):
        clid = 1
        noid = -1
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if "cluster" in filename:
                clusdict[clid] = extract_lmclus_dimension(filename, clid)
                clid += 1
            elif "noise" in filename:
                clusdict[noid] = 0
                noid -= 1

        return clusdict
    else:
        return "NORESULT"