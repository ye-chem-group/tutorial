import sys
import numpy as np
import scipy.linalg


def eig(h, s):
    ''' Diaognalizing 'h' under the metric 's'
    '''
    # h C = S C E
    # X = S^{-1/2} orthogonalize AO
    # (X h X) (X^{-1} C) = (X^{-1} C) epsilon
    # define: h' = X h X, C' = X^{-1} C
    # h' C' = C' epsilon
    # C = X C'
    return scipy.linalg.eigh(h, s)

def make_rdm1(mo_coeff, nocc):
    # D(mu,nu) = 2 * \sum_{i} C(mu,i) * C(nu,i) = \sum_{i} C(mu,i) C.T(i,nu)
    occ_coeff = mo_coeff[:,:nocc]
    dm = 2 * np.dot(occ_coeff, occ_coeff.T)
    return dm

def get_fock(h, V, dm):
    # Coulomb term: J(mu,nu) = \sum_{lambda,sigma} 2*(mu,nu|lambda,sigma) D(sigma,lambda)
    #               J(p,q) = \sum_{r,s} 2*(pq|rs) D(s,r)
    #                        \sum_{r,s} 2*(pq|rs) D(s,r) -> J(p,q)
    #                           J = einsum('pqrs,sr->pq', V, dm)
    # nao = h.shape[0]
    # J = np.zeros((nao,nao))
    # for p in range(nao):
    #     for q in range(nao):
    #         for r in range(nao):
    #             for s in range(nao):
    #                 J[p,q] += 2 * V[p,q,r,s] * dm[s,r]
    #         # np.trace(np.dot(V[pq], dm))
    J = 2 * np.einsum('pqrs,sr->pq', V, dm)
    # Exchange term: K(mu,nu) = -\sum_{lambda,sigma} (mu,sigma|lambda,nu) D(sigma,lambda)
    #                                                (p,s|r,q)
    K = -np.einsum('psrq,sr->pq', V, dm)

    f = h + 0.5 * (J + K)

    return f


def energy_total(h, f, dm, e_nuc):
    return 0.5 * np.trace( np.dot(dm, h+f) ) + e_nuc


def calculate_hf(h, V, S, e_nuc, nelectron):
    nao = h.shape[0]
    nocc = nelectron // 2
    max_cycle = 50
    conv_tol = 1e-8 # energy
    conv_tol_dm = 1e-4 # dm
    print('Number of AO: %d' % (nao))
    print('Number of electrons: %d' % (nelectron))
    print('Number of occupied MO: %d' % (nocc))
    print('Number of max SCF cycle: %d' % (max_cycle))
    print('conv_tol: %.15g' % (conv_tol))
    print('conv_tol_dm: %.15g' % (conv_tol_dm))

    # Initial guess: h C = S C epsilon
    mo_energy, mo_coeff = eig(h, S)
    dm = make_rdm1(mo_coeff, nocc)
    f = get_fock(h, V, dm)
    e_tot = energy_total(h, f, dm, e_nuc)
    print('Initial energy: %.15g' % (e_tot))

    # SCF
    converged = False
    for cycle in range(1,max_cycle+1):
        e_tot_old = e_tot
        dm_old = dm
        # update mo
        mo_energy, mo_coeff = eig(f, S)
        # update dm
        dm = make_rdm1(mo_coeff, nocc)
        # update e_tot
        e_tot = energy_total(h, f, dm, e_nuc)
        # checking convergence
        delta_e = abs(e_tot - e_tot_old)
        delta_dm = abs(dm - dm_old).max()
        # fds - sdf
        comm = np.dot(np.dot(f,dm),S)
        comm -= np.dot(np.dot(S,dm),f)
        delta_comm = abs(comm).max()
        if delta_e < conv_tol and delta_dm < conv_tol_dm and delta_comm < conv_tol_dm:
            converged = True
            break
        print('cycle= %d  e_tot= %.15g  delta_e= %.5g  delta_dm= %.5g  delta_comm= %.5g' %
              (cycle, e_tot, delta_e, delta_dm, delta_comm))
        # update fock
        f = get_fock(h, V, dm)

    if converged:
        print('SCF convergence is reached in %d cycles' % (cycle))
        print('Converged SCF energy is %.15g' % (e_tot))
    else:
        print('Ooooops! Convergence failed after %d cycles' % (max_cycle))
        print('Unconverged SCF energy is %.15g' % (e_tot))

    return e_tot, mo_energy, mo_coeff


if __name__ == '__main__':
    ints_path = '../reference/ints'
    h = np.load(f'{ints_path}/h.npy')
    V = np.load(f'{ints_path}/V.npy')
    S = np.load(f'{ints_path}/S.npy')
    e_nuc = np.load(f'{ints_path}/enuc.npy')
    nelectron = 10

    calculate_hf(h, V, S, e_nuc, nelectron)
