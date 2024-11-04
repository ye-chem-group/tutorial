import sys
import numpy as np

from pyscf.pbc import gto, scf
from pyscf import lib

log = lib.logger.Logger(sys.stdout, 6)


if __name__ == '__main__':
    try:
        fchk = sys.argv[1]
        basis = sys.argv[2]
        alat = float(sys.argv[3])   # cubic lattice constant for diamond at equilibrium
        pseudo = sys.argv[4]
        if pseudo.lower() == 'none': pseudo = None  # all-electron
        kmesh = [int(x) for x in sys.argv[5]]
        max_memory = float(sys.argv[6]) * 950   # GB to MB with a conservative factor 0.95
    except:
        log.info('Usage: fchk, basis, alat, pseudo, kmesh, max_memory (in GB)')
        sys.exit(1)

    log.info('Input args:')
    log.info('fchk = %s', fchk)
    log.info('basis = %s', basis)
    log.info('alat = %s', alat)
    log.info('kmesh = %s', kmesh)
    log.info('max_memory = %s', max_memory)
    log.info('\n')

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
    kmf.chkfile = fchk
    kmf.kernel()
    ehf_kpt = kmf.e_tot

    log.info('Final SCF energy is %.10f', ehf_kpt)
