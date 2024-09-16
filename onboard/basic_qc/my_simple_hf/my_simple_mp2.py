import sys
import numpy as np

from my_simple_hf import calculate_hf


if __name__ == '__main__':
    ints_path = '../reference/ints'
    h = np.load(f'{ints_path}/h.npy')
    V = np.load(f'{ints_path}/V.npy')
    S = np.load(f'{ints_path}/S.npy')
    e_nuc = np.load(f'{ints_path}/enuc.npy')
    nelectron = 10

    ehf, mo_energy, mo_coeff = calculate_hf(h, V, S, e_nuc, nelectron)

    # define orbital space
    nocc = nelectron // 2
    nao, nmo = mo_coeff.shape  # C(mu,p)
    nvir = nmo - nocc
    Cocc = mo_coeff[:,:nocc]
    Cvir = mo_coeff[:,nocc:]
    eocc = mo_energy[:nocc]
    evir = mo_energy[nocc:]

    # (ov|ov)-type ERIs
    # ovov1 = np.zeros((nocc,nvir,nocc,nvir))
    ## Cost: nocc^2 * nvir^2 * nao^4 ~ O(N^8)
    # for i in range(nocc):
    #     for a in range(nvir):
    #         for j in range(nocc):
    #             for b in range(nvir):
    #                 for p in range(nao):
    #                     for q in range(nao):
    #                         for r in range(nao):
    #                             for s in range(nao):
    #                                 ovov1[i,a,j,b] += V[p,q,r,s] * Cocc[p,i] * Cvir[q,a] *\
    #                                                                Cocc[r,j] * Cvir[s,b]
    ovov = np.einsum('pqrs,pi,qa,rj,sb->iajb', V, Cocc, Cvir, Cocc, Cvir, optimize=True)

    # print( abs(ovov - ovov1).max() )

    # energy denominator
    # denom1 = np.zeros((nocc,nvir,nocc,nvir))
    # for i in range(nocc):
    #     for a in range(nvir):
    #         for j in range(nocc):
    #             for b in range(nvir):
    #                 denom1[i,a,j,b] = (eocc[i] - evir[a]) + (eocc[j] - evir[b])
    eov = (eocc[:,None] - evir).reshape(nocc*nvir)
    denom = (eov[:,None] + eov).reshape(nocc,nvir,nocc,nvir)
    # print(abs(denom - denom1).max())

    # Correlation energy evaluation
    t2 = ovov / denom
    ed = 2 * np.einsum('iajb,iajb->', t2, ovov)
    ex = -   np.einsum('iajb,ibja->', t2, ovov)
    e_corr = ed + ex
    e_tot = ehf + e_corr
    print('MP2 total energy = %.15g  E_corr = %.15g' % (e_tot, e_corr))
