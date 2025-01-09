''' This script is for k-point LNO-CC (KLNOCC). For Gamma-point LNO-CCSD/CCSD(T), simply
    follow the script `lnoccsd_mol.py` with a Gamma-point periodic mean-field object.
'''


import os
import sys
import numpy as np

from pyscf.pbc import gto, scf, mp
from pyscf import lo
from pyscf.data.elements import chemcore
from pyscf import lib
from pyscf.lib import logger

from lno.cc import KLNOCCSD
from lno.tools import k2s_scf, kscf_remove_trsymm, k2s_aoint, autofrag_iao

log = logger.Logger(sys.stdout, 6)


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


def cache_ints1e_(kmf, mf):
    kh1e = kmf.get_hcore()
    ks1e = kmf.get_ovlp()

    h1e = k2s_aoint(kmf.cell, kmf.kpts, kh1e)
    s1e = k2s_aoint(kmf.cell, kmf.kpts, ks1e)

    mf.get_hcore = lambda *args: h1e
    mf.get_ovlp = lambda *args: s1e

    return mf


def get_iao_frag(mf, nkpts, frozen):
    s1e = mf.get_ovlp()
    occ_coeff = mf.mo_coeff[:,mf.mo_occ>1e-6][:,frozen:]
    iao_coeff = lo.iao.iao(mf.cell, occ_coeff)
    lo_coeff = lo.orth.vec_lowdin(iao_coeff, s1e)

    # check if IAOs span the entire occ space
    u = np.linalg.multi_dot((lo_coeff.T.conj(), s1e, occ_coeff))
    err = abs(np.dot(u.T.conj(), u) - np.eye(u.shape[1])).max()
    log.warn('IAO span occ error: %.3e', err)

    celliao = lo.iao.reference_mol(mf.cell)
    frag_lolist = autofrag_iao(celliao, 'atom') # this is fragment list for supercell
    nfrag = len(frag_lolist)
    nfrag_per_cell = nfrag//nkpts
    frag_lolist = frag_lolist[:nfrag_per_cell]  # frag list for unit cell only

    return lo_coeff, frag_lolist


def get_pm_frag(mf, nkpts, frozen):
    occ_coeff = mf.mo_coeff[:,mf.mo_occ>1e-6][:,frozen:]
    mlo = lo.PipekMezey(mf.cell, occ_coeff).set(verbose=4)
    lo_coeff = mlo.kernel()
    # while True:
    #     lo_coeff1 = mlo.stability_jacobi()[1]
    #     if lo_coeff1 is lo_coeff:
    #         break
    #     mlo = lo.PipekMezey(mf.cell, lo_coeff1).set(verbose=4)
    #     mlo.init_guess = None
    #     lo_coeff = mlo.kernel()

    nlo = lo_coeff.shape[1]

    from lno.tools import sort_orb_by_cell

    try:
        lo_coeff = sort_orb_by_cell(mf.cell, lo_coeff, nkpts, s=mf.get_ovlp())
        nlo_per_cell = nlo // nkpts
        lo_coeff = lo_coeff[:,:nlo_per_cell]
        frag_lolist = [[i] for i in range(nlo_per_cell)]
    except:
        log.warn('Assigning PM to unit cells failed. This indicates symmetry breaking '
                 'in orbital localization. LNO-CC will be performed for all LOs in the '
                 'supercell.')
        frag_lolist = [[i] for i in range(nlo)]

    return lo_coeff, frag_lolist


if __name__ == '__main__':
    try:
        fchk = sys.argv[1]
        flo = sys.argv[2]
        loname = sys.argv[3]
        frag_lolist_range = sys.argv[4]
        thresh_occ = float(sys.argv[5])
        thresh_vir = float(sys.argv[6])
        mem = float(sys.argv[7]) * 950
        fcderi = sys.argv[8]
        fovL = sys.argv[9]
    except:
        log.warn('Usage: fchk, flo, loname, frag_lolist_range (e.g., "0,80" means 0,1,2,...,79), thresh_occ, thresh_vir, mem_in_GB, fcderi, fovL')
        sys.exit(1)

    frag_lolist_range = list(map(int, frag_lolist_range.split(',')))
    frag_lolist = [[i] for i in range(*frag_lolist_range)]

    log.info('Command line arguments')
    log.info('fchk = %s', fchk)
    log.info('flo = %s', flo)
    log.info('loname = %s', loname)
    log.info('frag_lolist_range = %s', frag_lolist_range)
    log.info('frag_lolist = %s', frag_lolist)
    log.info('thresh_occ = %s', thresh_occ)
    log.info('thresh_vir = %s', thresh_vir)
    log.info('mem = %s', mem)
    log.info('fcderi = %s', fcderi)
    log.info('fovL = %s', fovL)
    log.info('\n')

    cell, scf_res = scf.chkfile.load_scf(fchk)

    kpts = scf_res['kpts']
    nkpts = Ncell = len(kpts)
    kmf = scf.KRHF(cell, kpts=kpts).density_fit()
    for k,v in scf_res.items():
        setattr(kmf, k, v)
    kmf.converged = True

    handle_cderi(kmf, fcderi)

    nocc_per_kpt = [np.count_nonzero(kmf.mo_occ[k] > 1e-6) for k in range(nkpts)]
    Nocc = sum(nocc_per_kpt)

    frozen_per_cell = chemcore(cell)
    frozen = frozen_per_cell * nkpts

    lo_coeff = lib.chkfile.load(flo, loname)
    Nlo = lo_coeff.shape[1]
    assert( Nlo == Nocc - frozen )
    nlo = Nlo // nkpts
    assert( nlo * nkpts == Nlo )

    # ''' Full MP2 for composite correction
    # '''
    # mmp = mp.KMP2(kmf, frozen=frozen_per_cell)
    # mmp.kernel(with_t2=False)
    # efull_mp2 = mmp.e_corr
    efull_mp2 = 0   # KMP2 has been ran separately

    ''' Reference KCCSD
    '''
    # from pyscf.pbc import cc
    # kcc = cc.KCCSD(kmf, frozen=frozen_per_cell)
    # kcc.kernel()
    # efull_ccsd = kcc.e_corr
    efull_ccsd = -0.208286627091736

    ''' K2Gamma
    '''
    mf = k2s_scf(kmf)
    cache_ints1e_(kmf, mf)  # precompute h1e and s1e for fast later access

    mcc = KLNOCCSD(kmf, lo_coeff, frag_lolist, frozen=frozen, mf=mf).set(verbose=5)
    if os.path.isfile(fovL):
        mcc._ovL = fovL
    else:
        mcc._ovL_to_save = fovL
    mcc.lo_proj_thresh_active = 0.1
    mcc.lno_type = ['1h','1h']
    mcc._h1e = mf.get_hcore()
    mcc._s1e = mf.get_ovlp()
    mcc.lno_thresh = [thresh_occ, thresh_vir]
    mcc.kernel()
    # elno_ccsd_uncorr = np.zeros_like(threshs)
    # elno_mp2 = np.zeros_like(threshs)
    # for i,thresh in enumerate(threshs):
    #     mcc.lno_thresh = [thresh*gamma, thresh]
    #     mcc.kernel()
    #     elno_ccsd_uncorr[i] = mcc.e_corr_ccsd
    #     elno_mp2[i] = mcc.e_corr_pt2
    # elno_ccsd = elno_ccsd_uncorr - elno_mp2 + efull_mp2
    #
    # log.info('')
    # log.info('Reference KCCSD E_corr = %.15g', efull_ccsd)
    # for i,thresh in enumerate(threshs):
    #     log.info('thresh = %.3e  E_corr(LNO-CCSD) = %.15g  E_corr(LNO-CCSD+∆PT2) = %.15g', thresh,
    #              elno_ccsd_uncorr[i], elno_ccsd[i])
