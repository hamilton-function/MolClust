import pubchempy as pcp
import csv


'''
Routine to fetch array of compounds from csv file
The compound names are located in the first column 
'''
def get_compoundnamearray(csvfile):
    namearray = []
    with open(csvfile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            namearray.append(row[0])
    return namearray



'''
reads in a csv of structure:
col0: cid 
col1: compound name

returns: dictionary with cid as key and compound name as value
'''

def get_cid_compound_dict(csvfile):
    ciddict = {}
    with open(csvfile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            ciddict[row[0]] = row[1]
    return ciddict




'''
Pipeline for PubChem interface (pubchem.ncbi.nlm.nih.gov)

Takes array of substance names
returns SDF or SMILES files

In the case that we obtain a list of compounds like for Emetine dihydrochloride hydrate
[Compound(201899), Compound(3068143), Compound(10219), Compound(21127119)]
we can simply take the first element of the array, since the results are sorted "by relevance" 
according to PubChem
'''

def get_sdf_from(arrayofnames,subdirpath):
    #dictionary with key: compound name, value: cid
    compounddict = {}

    for e in arrayofnames:
        molres = pcp.get_cids(e, 'name', 'substance', list_return='flat')
        if len(molres) == 0:
            print('NO COMPOUND FOUND FOR ',e)
        else:
            print('FETCHING SDF FOR COMPOUND: ',e)
            toprelcid = molres[0]
            pcp.download('SDF', subdirpath+str(toprelcid)+'.sdf', str(toprelcid), 'cid', overwrite=True)
            compounddict[e] = toprelcid

    print('DONE')
    return compounddict


#compnames = get_compoundnamearray('cov19candidates.csv')
#compound_dict = get_sdf_from(compnames,'E:\\MolClust\\sdfs\\')


'''
special function which downloads for a keyword: covid-19
all available compounds in sdf format

pcp.Compound.from_cid(5090)
'''
def get_all_pubchem_covid_compounds_sdf(ciddict,subdirpath):
    compounddict = {}

    for e in ciddict.keys():
        compound = ciddict[e]
        print('CID: ',e,' Compound: ',compound)
        pcp.download('SDF', subdirpath+str(e)+'.sdf', str(e), 'cid', overwrite=True)

        #line to download 3d structures
        '''
        res = pcp.get_compounds(e, 'cid', record_type='3d')
        if len(res) > 0:
            pcp.download('SDF', subdirpath + str(e) + '.sdf', str(e), 'cid', record_type = '3d', overwrite=True)
        '''

        compounddict[compound] = e

    print('DONE')
    return compounddict


#cidict = get_cid_compound_dict('E:\\MolClust\\PubChem_compound_text_covid-19.csv')
#compound_dict = get_all_pubchem_covid_compounds_sdf(cidict,'E:\\MolClust\\sdfs_covid19pubchem\\')


#cidict = get_cid_compound_dict('E:\\MolClust\\covidmancandidates.csv')
#compound_dict = get_all_pubchem_covid_compounds_sdf(cidict,'E:\\MolClust\\sdfs\\')