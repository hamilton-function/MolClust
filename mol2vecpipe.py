import subprocess
import glob
import pandas as pd
import csv


'''
Routine takes a sdf file as input
executes a featurization according to https://mol2vec.readthedocs.io/en/latest/
and returns a csv file including smiles information and a 300d mol2vec representation of the compound

* for the mapping the model_300dim.pkl is taken which relies on 19.9 million molecules from ChEMBL (ver. 23)
and ZINC 15 databases according to https://www.blopig.com/blog/2018/08/mol2vec-finding-chemical-meaning-in-300-dimensions/

* requires RDKit for the mol2vec library: https://github.com/samoturk/mol2vec/

* executes CLI routine: mol2vec featurize -i infile.sdf -o outfile.csv -m fullpath2model -r 1 --uncommon UNK
'''


def sdf2vec(sdffile):
    outfile = (sdffile.split('\\')[-1])[:-4]
    csvdir = 'E:\\MolClust\\csvvect\\'
    print('FEATURIZING SDF')
    process = subprocess.call('mol2vec featurize -i '+sdffile+' -o '+csvdir+outfile+'.csv -m E:\\MolClust\\model_300dim.pkl -r 1 --uncommon UNK')
    print('FEATURIZATION DONE')



'''
Auxillary routine to do a mapping from sdf to vec representation for all objects within a directory
input: directory consisting only of .sdf files
output: directory with corresponding csv files
'''

def mass_featurize_sdf2vec(inputdir):
    sdfarray = glob.glob(inputdir+'*.sdf')
    for e in sdfarray:
        print('PROCESSING FILE: ',e)
        sdf2vec(e)




'''
Auxillary routine which extracts 300d-vector for each compound
in: csv file with smiles information and vector representation
out: 300d-array / vector of compound
'''

def molvec_extractor(incsv):
    data = pd.read_csv(incsv)
    vec_cols = data.filter(regex='mol2vec')
    molvec = vec_cols.values

    return molvec



'''
Routine reads in all csv files with smiles and mol2vec information from a given directory
returns 
* dictionary with compound id as key and vector as values
* csv file with only the vectors which can then be used for autoenc-clustering
'''

def extract_vectorized_molecules(inputdir,csvout):
    molvec_dict = {}

    csvarray = glob.glob(inputdir+'*.csv')
    for csvfile in csvarray:
        #get the compound id from the csv filename:
        compid = int((csvfile.split('\\')[-1])[:-4])

        molvec = molvec_extractor(csvfile)
        molvec_dict[compid] = list(molvec)

    #routine to write all vectors to a csv file
    vecarr = list(molvec_dict.values())

    with open(csvout+'.csv','w', newline='') as file:
        writer = csv.writer(file)
        for vec in vecarr:
            for scalar in vec:
                writer.writerow(scalar)

    return molvec_dict


#mass_featurize_sdf2vec('E:\\MolClust\\sdfs\\')
#extract_vectorized_molecules('E:\\MolClust\\csvvect\\','cov2vect')