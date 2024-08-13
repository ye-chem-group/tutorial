import numpy as np

from pyscf import gto, df


if __name__ == '__main__':
    atom = '''
    O          0.00000        0.00000        0.11779
    H          0.00000        0.75545       -0.47116
    H          0.00000       -0.75545       -0.47116
    '''
    basis = 'cc-pvdz'
    auxbasis = 'cc-pvdz-ri'
    charge = 0
    spin = 0

    mol = gto.M(atom=atom, basis=basis, charge=charge, spin=spin)
    auxmol = df.addons.make_auxmol(mol, auxbasis)

    j3c = mol.intor('int3c2e')
    j2c = mol.intor_symmetric('int2c2e')
    print(j3c.shape)
    print(j2c.shape)
