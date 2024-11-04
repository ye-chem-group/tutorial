import sys
import numpy as np

from pyscf.pbc import gto, scf
from pyscf import lib

log = lib.logger.Logger(sys.stdout, 6)


if __name__ == '__main__':
    try:
        kmesh = [int(x) for x in sys.argv[1]]
        max_memory = float(sys.argv[2]) * 950   # GB to MB with a conservative factor 0.95
    except:
        log.info('Usage: kmesh, max_memory (in GB)')
        sys.exit(1)

    log.info('Input args:')
    log.info('kmesh = %s', kmesh)
    log.info('max_memory = %s', max_memory)
    log.info('\n')

    alat = 3.567    # cubic lattice constant for diamond at equilibrium
    basis = 'cc-pvdz'
    pseudo = None   # all-electron

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
    cell.max_memory = max_memory
    cell.build()

    # unit cell + kpts
    kpts = cell.make_kpts(kmesh)
    kmf = scf.KRHF(cell, kpts=kpts).density_fit()
    kmf.kernel()
    ehf_kpt = kmf.e_tot

    log.info('Final SCF energy is %.10f', ehf_kpt)
