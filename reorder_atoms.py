# -*- coding: utf-8 -*-
"""
script for reordering atoms in molecules
according to the names of a reference molecule

Originally meant to reorder POPC atoms from CHARMM-GUI to Slipids ordering
"""

import pmx


class MoleculesCompared:
    """
    Class for comparing molecules (pmx.Model instances)
    (meant mostly for two molecules, but written for any number)

    Initiated with a list of molecules to be compared.
    Molecules have to be of the same size (no. atoms),
    and they must be the instance of pmx.Molecule.

    Residue name is also checked, but only warning is given if it differs.

    1st molecule in the list is the reference molecule

    especially: meant for sorting CHARMM-GUI POPC lipids
    to match the order of Slipids

    See the code for further details
    """
    def __init__(self, mols):
        self.natoms = None
        self.resname = None

        if not isinstance(mols, (tuple, list)):
            raise RuntimeError, "given object is neither a molecule nor a list"

        for item in mols:
            if not isinstance(item, pmx.molecule.Molecule):
                raise RuntimeError, "at least one item in the list is not the instance of pmx.Molecule"
            else:
                if self.resname == None:
                    # should be done only for the first item
                    self.resname = item.resname
                elif self.resname != item.resname:
                    print "WARNING: molecules don't have the same resname! ", self.resname, "  -vs-  ", item.resname
                if self.natoms == None:
                    # should be done only for the first item
                    self.natoms = len(item.atoms)
                elif self.natoms != len(item.atoms):
                    raise RuntimeError, "the number of atoms in the given molecules don't agree -- probably different types of molecules"

        self.mols = mols

    def get_eq_at_index(self, atom):
        """
        Returns the index of the equivalent atom in the ref-molecule
        (ref-mol: the first one in the list-of-molecules self.mols)
        """
        refmol=self.mols[0]
        # fetch atom instance(s) with the same name
        eq_at = refmol.fetch(atom.name)
        if len(eq_at) != 1:
            raise RuntimeError, "Atom not found or multiple atoms with the same name exist!"
        else:
            return eq_at[0].id

    def get_ref_sorted_atoms(self):
        """
        Returns self.mols (list of mols) with name-sorted atoms
        after the first molecule in the list
        """
        sorted_list_mols = [self.mols[0]]
        for mol in self.mols[1:]:
            # sort atoms in the current molecule mol
            sorted_list_atoms = sorted(mol.atoms, key=self.get_eq_at_index)
            # turn the list of atoms into a pmx.Molecule instance
            sorted_mol = pmx.molecule.Molecule(resname=mol.resname,id=mol.id)
            for atom in sorted_list_atoms:
                sorted_mol.append(atom)
            # append the sorted_mol (pmx.Molecule instance)
            # at the end of the list of molecules
            sorted_list_mols.append(sorted_mol)
        return sorted_list_mols

#%%

if __name__ == '__main__':

    # load-in the files
    charmm = pmx.Model(filename="step5_charmm2gmx.pdb")
    slipids = pmx.Model(filename="confout_memb_slipids.gro")

    #%%

    # set the reference POPC residue as the 1st item in the list-of-mols
    list_POPC_mols = [slipids.residue(1)]

    # add all other POPC molecules to-be-sorted to the list
    for mol in charmm.residues:
        if mol.resname == 'POPC':
            list_POPC_mols.append(mol)

    #%%

    # sort the molecules! -- get the list of reference-sorted molecules
    refsorted = MoleculesCompared(list_POPC_mols).get_ref_sorted_atoms()

    # renumber atoms
    for mol in refsorted:
        mol.renumber_atoms()

    #%%

    # ref-sorted residues/molecules keep the original id (resid)
    # --> can replace the original residues easily
    for mol in refsorted[1:]:
        charmm.replace_residue(charmm.residue(mol.id), mol)

    # write it out in PDB (or GRO) format
    charmm.write("reordered_Charmm-Slipids_POPC_memb.pdb")
