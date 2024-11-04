## Pre-requisites:
1. Read PySCF's User Guide on [PBC](https://pyscf.org/user/pbc.html). Learn how to define a cell, make *k*-points, and run SCF calculations.
2. Download POSCAR file for diamond from [Materials Project](https://next-gen.materialsproject.org/materials/mp-66/). Rename the suffix to `.vasp`, which allows you to open it using VESTA. Study the `.vasp` file and convince yourself that all you need to know to construct such a file is (i) the XYZ coordinates of atoms in a cell and (ii) the lattice cell vectors. These two are also what you need to define a cell in PySCF: `atom` and `a`.
3. Density fitting (DF) is a commonly used technique for PBC calculations. Read David Sherrill's excellent [notes on DF](http://vergil.chemistry.gatech.edu/notes/df.pdf). This is for molecular calculations but the basic idea is the same even for solids.

## Diamond

Use the following primitive cell for diamond at equilibrium geometry for this exercise:
```python
from pyscf.pbc import gto
atom = '''
C    0.0000000000    0.0000000000    0.0000000000
C    0.8917500000    0.8917500000    0.8917500000
'''
a = np.asarray([
    [1.7835000000, 1.7835000000, 0.0000000000],
    [0.0000000000, 1.7835000000, 1.7835000000],
    [1.7835000000, 0.0000000000, 1.7835000000]
])
cell = gto.Cell(atom=atom, a=a, ...)    # fill in other settings
cell.build()
```

### Thermodynamic limit (TDL)

One key difference between solids and molecules is that solids are intrinsically infinite. As a result, practical simulations of solids using periodic boundary conditions (PBCs) require increasing the simulation cell size to approximate the thermodynamic limit (TDL). This can be achieved either in real space by creating **supercells**, or in reciprocal space through **k-point sampling** -- two equivalent methods. In this exercise, you will focus on exploring the k-point sampling approach.

1. Let's first verify that a SCF calculation using k-point sampling is equivalent to a supercell calculation where the size of the supercell matches that of the k-point mesh (e.g., a 2x2x2 kpoint mesh vs a 2x2x2 supercell). An example script is provided in `pbc/script/kpts_vs_gamma` ([link](https://github.com/ye-chem-group/tutorial/tree/main/onboard/pbc/script/kpts_vs_gamma)). Find the computational cost of both calculations by grepping lines containing "CPU time for SCF". Convince yourself that using k-point sampling is much more efficient.

2. Performing k-point RHF calculations for diamond using cc-pVDZ basis set and k-point mesh of increasing size `(n,n,n)` for `n = 1, 2, 3, 4, 5, 6`. In `pbc/script/kpts_rhf` ([link](https://github.com/ye-chem-group/tutorial/tree/main/onboard/pbc/script/kpts_rhf)), an example pyscf script is provided together with a shell script for automatically submitting jobs (Don't forget to set the `pyscf_path` variable before running it).

3. While waiting for the calculations to finish, we can learn more about the details of these calculations. The KRHF calculations use density fitting, which requires storing the three-center (3c) Coulomb integrals `(mu,nu|L)` whose size scales as `N_integral = Nk*(Nk+1)/2 * nao**2 * naux`, where `Nk = n**3` for a `(n,n,n)` k-point mesh, `nao` and `naux` are the number of atomic orbitals (AOs) and auxiliary orbitals (used for density fitting) per unit cell. You can find `nao` and `naux` by grepping "cGTOs" in any of the output file (should be 28 and 140). Now use this formula to estimate the required storage for the KRHF calculations you are performing (to obtain a number in GB, use `Mem_GB = N_integral * 16 / 1e9` where the factor 16 is due to the k-point integrals being complex numbers). Convince yourself that the memory you gave to the calculation (12 core * 4 GB/core = 48 GB) is enough for even the largest calculations with a `6x6x6` k-point mesh.

4. Let's analyze the convergence behavior of HF energy with the k-point mesh size. Theoretical analysis suggests that for sufficiently large k-point meshes, HF energy converges as `1/Nk` to the TDL. Plot HF energy as a function of `1/Nk`. Convine yourself that for diamond, the asymptotic linear regime (linear in `1/Nk`) is not reached until `Nk = 5^3` and `6^3`. Extrapolate using `Etdl = (E1 * Nk1 + E2 * Nk2) / (Nk1 - Nk2)` for `Nk1 = 5^3` and `Nk2 = 6^3` to obtain your final estimate of the HF energy in the TDL limit using a DZ basis set.

5. Now let's repeat the analysis for correlation energy evaluated at MP2 level. In `pbc/script/kpts_mp2` ([link](https://github.com/ye-chem-group/tutorial/tree/main/onboard/pbc/script/kpts_mp2)), an example pyscf script is provided, together with a shell script for automatically submitting jobs (Don't forget to set the `pyscf_path` variable before running it). [Optional] Before you run these calculations, it may be helpful to install a modified KMP2 code I wrote that in general runs 10x faster ([link](https://github.com/hongzhouye/pyscf/blob/pbc_mp2)). Following these steps:
    - copy `pyscf/pbc/mp/kmp2.py` to your pyscf directory
    - copy the entire directory `pyscf/lib/mp` to your pyscf directory
    - add this [line](https://github.com/hongzhouye/pyscf/blob/pbc_mp2/pyscf/lib/CMakeLists.txt#L139) to `pyscf/lib/CMakeLists.txt`
    - recompile pyscf

6. Plot MP2 correlation energy as a function of `1/Nk`. You should see an asymptotic linear convergence rate in `1/Nk` for large `Nk`. However, unlike HF, the asymptotic regime is reached much earlier at around `Nk = 3^3`. Quantify this by performing extrapolation using adjacent `Nk`'s, i.e., `(1^3, 2^3), (2^3, 3^3), (3^3, 4^3)` etc. Plot all extrapolated MP2 correlation energy on a single plot and you should see the curve levels off very quickly. Using `(5^3, 6^3)`-extrapolation as the reference value, what is the minimum `Nk`-duo that gives an extrapolated result that is within **1 kcal/mol per atom** from the reference value (recall that a unit cell has 2 atoms; 1 Ha = 627.5096 kcal/mol).

7. Finally, let's analyze the computational cost. Grepping `CPU time for SCF` for the cost of HF calculations. Grapping CPU time for `KMP2` for that of MP2 calculations (you should see two lines, one for `ao2mo` and the other for the energy evaluation; add them up). Plot both costs (using either CPU time or wall time should be fine) as a function of `Nk` in a log-log scale. Perform a linear fit using the last few points to obtain the scaling factor with respect to `Nk`. Do you see a **quadratic** scaling for the HF cost and a **cubic** scaling for the MP2 cost?

8. We are one step away to obtain a meaningful quantity: the cohesive energy of diamond. For atomic crystal like diamond, the cohesive energy is defined as the energy required to separates all atoms in the crystal into free atoms, i.e., `Ecoh = Ecrys - Eatom`, where `Ecrys` is the crystal energy in the TDL limit normalized to a single atom and `Eatom` is the energy of a free atom. Now calculate `Eatom` at both HF and MP2 level for a free carbon atom using the same cc-pVDZ basis set (NOTE! The ground state of carbon is a spin *triplet* and you need to do *UHF* and *UMP2* calculations). Calculate `Ecoh` for diamond. Referene values are `5.4 eV` for HF and `8.0 eV` for MP2 taken from literature (e.g., Table IV of this [JCP paper](https://doi.org/10.1063/5.0119633)). Your results should be close to the reference. The remaining differences arise from (i) basis set incompleteness error (DZ is not enough) and (ii) basis set superposition error (BSSE) in the atomic calculations.

9. [Optional] Read David Sherrill's excellent [notes on BSSE](http://vergil.chemistry.gatech.edu/notes/cp.pdf). Re-do the calculation of a single carbon atom with ghost atoms from the nearest neighbor to correct for BSSE (for diamond, there are 4 nearest neighbors). How does this change your estimated $E_{\mathrm{coh}}$? (For learning how to use ghost atoms in PySCF, simply google "pyscf ghost atoms".)

10. [Optional] To correct for the finite-basis error, we can use a larger basis set, e.g., cc-pVTZ, and repeat all the calculations above. An alternative to this brute-force approach (which can become very expensive especially for large k-point meshes) is to use some composite correction, `E(TZ,TDL) ≈ E(TZ,Nk_small) + E(DZ,TDL) - E(DZ,Nk_small)`, where `Nk_small` is a small k-point mesh where a TZ calculation is affordable. Obviously, the larger `Nk_small` is, the more accurate the composite correction will be. Perform TZ calculations using `Nk = 1^3, 2^3, 3^3, 4^3` and calculate the corresponding composite correction-estimated `E(TZ,TDL)`. Using `Nk = 4^3` as a reference, what is the minimum `Nk` needed to achieve accuracy better than 1 kcal/mol per atom?
