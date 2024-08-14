import os
import numpy as np

from pyscf import gto, scf, mp
from pyscf import ao2mo


def dump_ao_integrals(mol, out_path):
    enuc = mol.energy_nuc()
    T = mol.intor_symmetric('int1e_kin')
    vnuc = mol.intor_symmetric('int1e_nuc')
    h = T + vnuc
    S = mol.intor_symmetric('int1e_ovlp')
    V = ao2mo.restore('s1', mol.intor('int2e', aosym='s8'), mol.nao_nr())

    # save integrals to file
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    np.save(f'{out_path}/h.npy', h)
    np.save(f'{out_path}/V.npy', V)
    np.save(f'{out_path}/S.npy', S)
    np.save(f'{out_path}/enuc.npy', enuc)

if __name__ == '__main__':
    fxyz = 'geom/h2o_eq.xyz'
    basis = 'sto-3g'
    charge = 0
    spin = 0

    atom = open(fxyz, 'r').read().rstrip('\n').lstrip('\n').split('\n')[2:]

    mol = gto.M(atom=atom, basis=basis, charge=charge, spin=spin)
    mol.verbose = 4

    print('Number of electrons: %d' % (mol.nelectron))
    print('Number of occupied MOs: %d' % (mol.nelectron//2))

    out_path = 'ints'
    dump_ao_integrals(mol, out_path)

    mf = scf.RHF(mol)
    mf.kernel()

    mmp = mp.MP2(mf)
    mmp.kernel()
