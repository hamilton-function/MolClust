from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
import glob
import subprocess


'''
def desalt(sdfname):
    # load sdf
    suppl = Chem.ForwardSDMolSupplier(sdfname, sanitize=False)

    # create salt eliminator
    remover = SaltRemover()

    wr = Chem.SDWriter(sdfname[0:-4] + '_ds.sdf')

    print('DE-SALTING...')

    for m in suppl:
        m.UpdatePropertyCache()
        print(m.GetNumAtoms())
        clean_mol = remover.StripMol(m)
        print(clean_mol.GetNumAtoms())

        #since rdkit eliminates for whatever reasons ALL hydrogens, we try this out...
        #hydmol = Chem.AddHs(clean_mol)
        wr.write(clean_mol)

    wr.flush()
    wr.close()

    print('NO MORE SALT - ENJOY!')
'''
# the rdkit de-salting is pure crap and a massive pain-in-the-ass ...if you are going to remove the salts, it also eliminates the damn hydrogens ...
# even if there is a solution - I go now for openbabel - seems to be by far more reliable by doing the -ONE- thing it should do, and not other side-bullshit.


allsdfs = glob.glob("E:\\MolClust\\saltycandidates\\*.sdf")

#obabel xxx.sdf -O xxx_ds.sdf -r
for sdf in allsdfs:
    outdir = sdf[0:-4]+'_ds.sdf'
    print(outdir)
    process = subprocess.run('obabel '+sdf+' -O '+outdir+' -r', timeout=60)
