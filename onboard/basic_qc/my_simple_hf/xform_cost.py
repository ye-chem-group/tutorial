import time
import numpy as np

from pyscf import lib


if __name__ == '__main__':
    nao = 80
    nocc = 30
    nvir = nao - nocc

    V = np.random.rand(nao,nao,nao,nao) # ERI
    Cocc = np.random.rand(nao,nocc) # occ MO
    Cvir = np.random.rand(nao,nvir) # vir MO

    # # algorithm 1: O(N^8)
    # tick = time.time()
    # ovov1 = np.einsum('pqrs,pi,qa,rj,sb->iajb', V, Cocc, Cvir, Cocc, Cvir)
    # tock = time.time()
    # dt = tock - tick
    # print('Time elapse for algorithm 1 : %.2f sec' % (dt))

    # algorithm 2: O(N^5)
    tick = time.time()
    tmp1 = np.einsum('pqrs,pi->iqrs', V, Cocc)
    tmp2 = np.einsum('iqrs,rj->iqjs', tmp1, Cocc)
    tmp3 = np.einsum('iqjs,qa->iajs', tmp2, Cvir)
    ovov2 = np.einsum('iajs,sb->iajb', tmp3, Cvir)
    tock = time.time()
    dt = tock - tick
    print('Time elapse for algorithm 2 : %.2f sec' % (dt))

    # algorithm 3: O(N^5)
    tick = time.time()
    ovov1 = np.einsum('pqrs,pi,qa,rj,sb->iajb', V, Cocc, Cvir, Cocc, Cvir, optimize=True)
    tock = time.time()
    dt = tock - tick
    print('Time elapse for algorithm 3 : %.2f sec' % (dt))

    # algorithm 4: O(N^5)
    tick = time.time()
    ovov1 = lib.einsum('pqrs,pi,qa,rj,sb->iajb', V, Cocc, Cvir, Cocc, Cvir)
    tock = time.time()
    dt = tock - tick
    print('Time elapse for algorithm 4 : %.2f sec' % (dt))

    # print(abs(ovov1 - ovov2).max())
