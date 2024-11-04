import sys
import numpy as np

from pyscf.pbc import gto, scf, mp
from pyscf.data.elements import chemcore
from pyscf import lib

log = lib.logger.Logger(sys.stdout, 6)


if __name__ == '__main__':
    try:
        fchk = sys.argv[1]
        max_memory = float(sys.argv[2]) * 950   # GB to MB with a conservative factor 0.95
    except:
        log.info('Usage: fchk, max_memory (in GB)')
        sys.exit(1)

    log.info('Input args:')
    log.info('fchk = %s', fchk)
    log.info('max_memory = %s', max_memory)
    log.info('\n')

    cell, scf_res = scf.chkfile.load_scf(fchk)
    cell.max_memory = max_memory    # The new input `max_memory` may be
                                    # different from previously used

    mf = scf.KRHF(cell, kpts=scf_res['kpts']).density_fit()
    for k,v in scf_res.items():
        setattr(mf, k, v)
    mf.converged = True

    frozen = chemcore(cell)
    log.info('Freezing %d orbitals per cell', frozen)

    mymp = mp.KMP2(mf, frozen=frozen)
    mymp.kernel(with_t2=False)

    log.info('Final MP2 correlation energy is %.10f', mymp.e_corr)
    log.info('Final MP2 total energy is %.10f', mymp.e_tot)
