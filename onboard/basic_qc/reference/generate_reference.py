import numpy as np

from pyscf import gto, scf, mp
from pyscf import ao2mo


def dump_ao_integrals(mf):
    h = mf.get_hcore()
    S = mf.get_ovlp()
    V = ao2mo.restore('s1', mf._eri, mf.mol.nao_nr())

    np.save('h', h)
    np.save('V', V)
    np.save('S', S)

if __name__ == '__main__':
    atom = '''
    O          0.00000        0.00000        0.11779
    H          0.00000        0.75545       -0.47116
    H          0.00000       -0.75545       -0.47116
    '''
    basis = 'sto-3g'
    charge = 0
    spin = 0

    mol = gto.M(atom=atom, basis=basis, charge=charge, spin=spin)

    mf = scf.RHF(mol)
    mf.kernel()

    print('Nuclear repulsion: %.15f' % (mf.energy_nuc()))
    print('Electronic part of the HF energy: %.15f' % (mf.energy_elec()[0]))

    s1e = mf.get_ovlp()
    dm = mf.make_rdm1()
    PS = np.dot(dm,s1e)
    print('Mulliken populations:', np.diag(PS))
    print('Sum of Mulliken populations:', np.trace(PS))


    dump_ao_integrals(mf)

    mmp = mp.MP2(mf)
    mmp.kernel()
