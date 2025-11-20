from pyscf import gto, scf
from pyscf.tools import cubegen


if __name__ == '__main__':
    atom = 'ch3oh.xyz'
    basis = 'cc-pvdz'

    # run a mean-field calculation and calculate density matrix
    mol = gto.M(atom=atom, basis=basis)
    mf = scf.RHF(mol).density_fit()
    mf.kernel()
    dm = mf.make_rdm1()

    margin = 10  # The default margin = 3 Å is usually too small

    # calculate electron density
    fout = 'den.cube'
    den = cubegen.density(mol, fout, dm, margin=margin)

    # calculate MEP
    fout = 'mep.cube'
    mep = cubegen.mep(mol, fout, dm, margin=margin)
