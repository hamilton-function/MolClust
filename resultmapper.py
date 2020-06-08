'''
Two dictionaries to retrieve compound information

dict 1:
key: compound id, value: vector

dict 2:
key: compound id, value: compound name + description
'''

import glob
import mol2vecpipe as m2vp
import pubchempy as pcp
import csv
import pickle



'''
get compound name + description
'''
def get_compoundnamearray(csvfile):
    namearray = []
    with open(csvfile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            namearray.append([row[0],row[1]])
    return namearray



'''
def generate_mappings(csvvectdir, rawcsv):
    cid2vect_dict = {}
    cid2comp_dict = {}

    csvarray = glob.glob(csvvectdir+'*.csv')

    print('POPULATING CID2VECT LOOKUP-TABLE:')
    for csvfile in csvarray:
        #get the compound id from the csv filename:
        compid = int((csvfile.split('\\')[-1])[:-4])

        molvec = m2vp.molvec_extractor(csvfile)
        cid2vect_dict[compid] = list(molvec)

    namearray = get_compoundnamearray(rawcsv)

    print('POPULATING CID2COMP LOOKUP-TABLE:')

    for entry in namearray:
        e,descr = entry
        molres = pcp.get_cids(e, 'name', 'substance', list_return='flat')
        if len(molres) == 0:
            print('NO COMPOUND FOUND FOR ',e)
        else:
            print('FETCHING CID FOR COMPOUND: ',e)
            toprelcid = molres[0]
            cid2comp_dict[toprelcid] = [e,descr]

    ultradict = {}
    for key in cid2vect_dict.keys():
        ultradict[key] = (cid2comp_dict[key],cid2vect_dict[key])

    return ultradict
'''


def generate_mappings(csvvectdir, rawcsv):
    cid2vect_dict = {}
    cid2comp_dict = {}

    csvarray = glob.glob(csvvectdir+'*.csv')

    print('POPULATING CID2VECT LOOKUP-TABLE:')
    for csvfile in csvarray:
        #get the compound id from the csv filename:
        compid = int((csvfile.split('\\')[-1])[:-4])

        molvec = m2vp.molvec_extractor(csvfile)
        cid2vect_dict[compid] = list(molvec)

    namearray = get_compoundnamearray(rawcsv)

    print('POPULATING CID2COMP LOOKUP-TABLE:')

    for entry in namearray:
        cid,name = entry
        cid2comp_dict[int(cid)] = [name,'']


    ultradict = {}
    for key in cid2vect_dict.keys():
        print('key: ',key)
        ultradict[key] = (cid2comp_dict[key],cid2vect_dict[key])

    return ultradict


#create mapping and pickle it!
'''
covmapping = generate_mappings('E:\\MolClust\\csvvect\\', 'E:\\MolClust\\covidmancandidates.csv')

filename = 'covmapping.pkl'
outfile = open(filename,'wb')
pickle.dump(covmapping,outfile)
outfile.close


for key in covmapping:
    print(covmapping[key])
'''