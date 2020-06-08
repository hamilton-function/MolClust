'''
file to extract results and get names + descriptions of detected clusters
'''


import hyperpipeline as hyp
import classicalgoeval as cge
import pubchempy as pcp
import pickle



'''
pickle clustering dictionary
'''
def pickle_clustering(fname,cdict):
    filename = fname+'.pkl'
    outfile = open(filename, 'wb')
    pickle.dump(cdict, outfile)
    outfile.close


'''
load pickled clustering dictionary
'''
def load_pickeled_clustering(filename):
    infile = open(filename+'.pkl','rb')
    mapping = pickle.load(infile)
    infile.close()
    return mapping


'''
routine to fetch mol structure  BY CID as png files and store them in a given folder
'''
def get_structure_image_by_cid(cid, compname, targetdir, clusterid):
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






#activate block comment when using classicalgoeval.py!
'''
resdict = hyp.get_all_results('E:\\DesktopBackup\\elki-0.7.1\\ORC_resrec')


clusrescomp = {}
for key in resdict:
    for vec in resdict[key]:
        if key not in clusrescomp.keys():
            namedesc = cge.vect2namedesc(vec)[0]
            cid = cge.vect2cid(vec)
            clusrescomp[key] = [[namedesc, cid]]
        else:
            namedesc = cge.vect2namedesc(vec)[0]
            cid = cge.vect2cid(vec)
            clusrescomp[key].append([namedesc,cid])


for key in clusrescomp:
    for e in clusrescomp[key]:
        print('KEY ',key,'  :  ',e)


pickle_clustering('new_orcdictk4l25',clusrescomp)

clusterdict = load_pickeled_clustering('new_orcdictk4l25')

get_structure_image_clustering(clusterdict, 'E:\\MolClust\\expres\\k4l25\\')
'''
