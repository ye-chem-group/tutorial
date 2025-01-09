import sys
import numpy as np

from pyscf.pbc import gto, scf
from pyscf import lo
from pyscf.data.elements import chemcore
from pyscf import lib

from lno.tools.kpts2supcell import k2s_scf, k2s_aoint
from lno.tools.tools import sort_orb_by_cell

log = lib.logger.Logger(sys.stdout, 6)


if __name__ == '__main__':
    try:
        fchk = sys.argv[1]
        flo = sys.argv[2]
        loname = sys.argv[3]
        mem = float(sys.argv[4]) * 950
    except:
        log.warn('Usage: fchk, flo, loname, mem')
        sys.exit(1)

    log.info('Command line arguments')
    log.info('fchk = %s', fchk)
    log.info('flo = %s', flo)
    log.info('loname = %s', loname)
    log.info('mem = %s', mem)
    log.info('\n')

    cell, scf_res = scf.chkfile.load_scf(fchk)

    kpts = scf_res['kpts']
    nkpts = Ncell = len(kpts)
    kmf = scf.KRHF(cell, kpts=kpts).density_fit()
    for k,v in scf_res.items():
        setattr(kmf, k, v)
    kmf.converged = True
    nocc_per_kpt = [np.count_nonzero(kmf.mo_occ[k] > 1e-6) for k in range(nkpts)]
    Nocc = sum(nocc_per_kpt)

    frozen_per_cell = chemcore(cell)
    frozen = frozen_per_cell * nkpts

    mf = k2s_scf(kmf)
    scell = mf.cell

    loname0 = f'{loname}_unsorted'

    occ_coeff = mf.mo_coeff[:,frozen:Nocc]
    mlo = lo.PipekMezey(scell, occ_coeff).set(verbose=4)
    lo_coeff = mlo.kernel()

    while True:
        lo_coeff1 = mlo.stability_jacobi()[1]
        if lo_coeff1 is lo_coeff:
            break
        mlo = lo.PipekMezey(scell, lo_coeff1).set(verbose=4)
        mlo.init_guess = None
        lo_coeff = mlo.kernel()
    lib.chkfile.dump(flo, loname0, lo_coeff)

    s1e = k2s_aoint(cell, kpts, kmf.get_ovlp())
    lo_coeff_sorted = sort_orb_by_cell(scell, lo_coeff, nkpts, s=s1e)
    lib.chkfile.dump(flo, loname, lo_coeff_sorted)
