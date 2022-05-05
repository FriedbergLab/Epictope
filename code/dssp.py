import os
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB import PDBList

arr = os.listdir("./data/pdb/")
for element in arr:
    file = element.split(".")[0]
    p = PDBParser(PERMISSIVE = True, QUIET = True)
    structure = p.get_structure(file, "data/pdb/"+file+".pdb")
    # print(structure)
    model = structure[0]
    # print(model)
    dssp = DSSP(model, "data/pdb/"+file+".pdb", dssp='mkdssp')
    dssp = DSSP(model, "data/pdb/"+file+".pdb", dssp='mkdssp', acc_array='Sander', file_type='PDB')
    a_key = list(dssp.keys())[2]
    # print(list(dssp.keys()))
    n = 0
    for n in range(0, len(list(dssp.keys()))):
        print(dssp.keys()[n], dssp[list(dssp.keys())[n]])
        n += 1
    # print(dssp[a_key])
    dssp_tuple = dssp_dict_from_pdb_file("data/pdb/"+file+".pdb")
    dssp_dict = dssp_tuple[0]
    # print(dssp_dict.keys())
    fhw = open("output/AlphaFoldPredictions_Features/"+file+".txt", "w")
    for j in dssp.keys():
        fhw.write(str(j[0])+"\t"+str(j[1][1])+"\t"+str(dssp[j][1])+"\t"+str(dssp[j][2])+"\t"+str(dssp[j][3])+"\n")
    fhw.close()
