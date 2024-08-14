from pyscf import gto
from pyscf.data.nist import BOHR

from zflow.pyscf_helper import LAT


if __name__ == '__main__':
    tasks = [
        (0.9, 'h2o_eq'),
        (1.8, 'h2o_2eq'),
        (2.7, 'h2o_3eq'),
    ]

    for R,prefix in tasks:
        atom = f'''
        O
        H 1 {R}
        H 1 {R} 2 104.5
        '''
        mol = gto.M(atom=atom, basis='sto-3g')
        atms = [mol.atom_symbol(i) for i in range(mol.natm)]
        rs = mol.atom_coords() * BOHR
        atom = [(atm,r) for atm,r in zip(atms,rs)]
        lat = LAT().init_from_pyscf_atom(atom)
        lat.dump_xyz(fout=f'{prefix}.xyz')
