#import Bio.PDB
#import Bio.SVDSuperimposer

import Bio.PDB.PDBParser
import Bio.SVDSuperimposer
#import Bio.PDB.Polypeptide
import numpy as np

pdb_code = "1JOY"
pdb_filename = "%s.pdb" % pdb_code
pdb_out_filename = "%s_aligned.pdb" % pdb_code

seq_str = 'MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEECNAIIEQFIDYLR'
use_str = 'MA---------RTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDI---------------'
use = [(letter != "-") for letter in use_str]
assert len(use) == len(seq_str)


print ("Loading PDB file %s" % pdb_filename)
structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)

print(structure[0])

print ("Everything aligned to first model...")
ref_model = structure[0]

for alt_model in structure :
    #Build paired lists of c-alpha atoms, ref_atoms and alt_atoms
    ref_atoms = []
    alt_atoms = []
    for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
        for ref_res, alt_res, amino, allow in zip(ref_chain, alt_chain, seq_str, use) :
            assert ref_res.resname == alt_res.resname
            assert ref_res.id      == alt_res.id
            assert amino == Bio.PDB.Polypeptide.three_to_one(ref_res.resname)
            if allow :
                #CA = alpha carbon
                #ref_atoms.append(ref_res['CA'])                
                #alt_atoms.append(alt_res['CA'])
                ref_atoms.append(ref_res['CA'].get_coord())
                alt_atoms.append(alt_res['CA'].get_coord())
                
    #Align these paired atom lists:
    #super_imposer = Bio.PDB.Superimposer()
    #super_imposer.set_atoms(ref_atoms, alt_atoms)
    
    sup = Bio.SVDSuperimposer.SVDSuperimposer()
    ref_atoms = np.asarray(ref_atoms)
    alt_atoms = np.asarray(alt_atoms)
    #print(ref_atoms)
    sup.set(ref_atoms, alt_atoms)
    
    print ("model 1/model %i" % alt_model.id+1)
    #Print RMSD for pair of structures before alignment
    print ("RMS(bef) = %0.2f" % (sup._rms(ref_atoms, alt_atoms)))
    
    sup.run()

    '''
    if ref_model.id == alt_model.id :
        #Check for self/self get zero RMS, zero translation
        #and identity matrix for the rotation.
        #assert numpy.abs(super_imposer.rms) < 0.0000001
        #assert numpy.max(numpy.abs(super_imposer.rotran[1])) < 0.000001
        #assert numpy.max(numpy.abs(super_imposer.rotran[0]) - numpy.identity(3)) < 0.000001
        assert numpy.abs(sup._rms(ref_atoms, alt_atoms)) < 0.0000001
        assert numpy.max(numpy.abs(sup.get_rotran()[1])) < 0.000001
        assert numpy.max(numpy.abs(sup.get_rotran()[0]) - numpy.identity(3)) < 0.000001
    else :
        #Update the structure by moving all the atoms in
        #this model (not just the ones used for the alignment)
        #super_imposer.apply(alt_model.get_atoms())
    '''
    #Print RMSD for pair of structures after alignment
    print ("RMS(aft) = %0.2f" % (sup.get_rms()))
    print (sup.get_rotran()[0], sup.get_rotran()[1])
    
    for res_atom in alt_model.get_atoms():
        res_atom.transform(sup.get_rotran()[0],sup.get_rotran()[1])
    
    #structure[alt_model] = alt_model
    
io=Bio.PDB.PDBIO()
io.set_structure(structure)
io.save(pdb_out_filename)
    