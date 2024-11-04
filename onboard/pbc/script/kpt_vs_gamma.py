import numpy as np

from pyscf.pbc import gto, scf
from pyscf.pbc.tools import super_cell
from pyscf import lib


if __name__ == '__main__':
    alat = 3.567    # cubic lattice constant for diamond at equilibrium
    basis = 'cc-pvdz'
    pseudo = None   # all-electron
    kmesh = (2,2,2) # 2x2x2 supercell/kpoint mesh
    Nk = int(np.prod(kmesh))

    atom = f'''
    C    0.0000000000    0.0000000000    0.0000000000
    C    {alat*0.25}    {alat*0.25}    {alat*0.25}
    '''
    a = np.asarray([
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5]
    ]) * alat

    cell = gto.Cell(atom=atom, a=a, basis=basis, pseudo=pseudo) # fill in other settings
    cell.verbose = 5
    cell.max_memory = 4000  # unit: MB; adjust for your need
    cell.build()

    log = lib.logger.new_logger(cell)

    # unit cell + kpts
    kpts = cell.make_kpts(kmesh)
    kmf = scf.KRHF(cell, kpts=kpts).density_fit()
    kmf.kernel()
    ehf_kpt = kmf.e_tot

    # supercell + Gamma point (i.e., 1x1x1 kpoint mesh)
    scell = super_cell(cell, kmesh)
    smf = scf.RHF(scell).density_fit()
    smf.kernel()
    ehf_sup = smf.e_tot / Nk    # normalize to 1 unit cell

    log.info('E(kpoint)   = %.10f', ehf_kpt)
    log.info('E(sup cell) = %.10f', ehf_sup)
