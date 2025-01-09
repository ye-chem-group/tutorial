import sys
import numpy as np

from pyscf.pbc import gto, scf
from pyscf import lib

log = lib.logger.Logger(sys.stdout, 6)


def handle_cderi(mf, fcderi, log=None):
    import os
    if log is None:
        from pyscf.lib import logger
        log = logger.new_logger(mf)
    if fcderi is None:
        return mf
    if os.path.isfile(fcderi):
        mf.with_df._cderi = fcderi
    else:
        mf.with_df._cderi_to_save = fcderi
    return mf


if __name__ == '__main__':
    try:
        fxyz = sys.argv[1]
        falat = sys.argv[2]
        basis = sys.argv[3]
        fchk = sys.argv[4]
        kmesh = list(map(int, sys.argv[5]))
        mem = float(sys.argv[6]) * 950
        fcderi = sys.argv[7]
    except:
        log.warn('Usage: fxyz, falat, basis, fchk, kmesh, mem, fcderi')
        sys.exit(1)

    log.info('Command line arguments')
    log.info('fxyz = %s', fxyz)
    log.info('falat = %s', falat)
    log.info('basis = %s', basis)
    log.info('fchk = %s', fchk)
    log.info('kmesh = %s', kmesh)
    log.info('mem = %s', mem)
    log.info('fcderi = %s', fcderi)
    log.info('\n')

    atom = '\n'.join(open(fxyz, 'r').read().lstrip('\n').rstrip('\n').split('\n')[2:])
    a = np.loadtxt(falat)

    cell = gto.M(atom=atom, basis=basis, a=a).set(verbose=5)
    kpts = cell.make_kpts(kmesh)

    kmf = scf.KRHF(cell, kpts=kpts).density_fit()
    kmf.chkfile = fchk
    handle_cderi(kmf, fcderi)
    kmf.kernel()
