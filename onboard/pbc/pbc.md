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


1. Perform k-point RHF (KRHF) calculation for diamond using the cc-pVDZ basis set and k-point mesh of size $n \times n \times n$ for $n = 1, 2, 3, 4, 5, 6$. Plot the HF energy as a function of $N_k^{-1}$. Did you see it shows a trend of convergence towards $N_k \to \infty$. (*You can read this [example](https://github.com/pyscf/pyscf/blob/master/examples/pbc/21-k_points_all_electron_scf.py) to learn how to do a KRHF calculation with density fitting.*)

2. The difference between calculations using a finite k-point mesh and the TDL is termed the **finite-size error**. The finite-size error in the HF energy decays as $N_k^{-1}$ for sufficiently large $N_k$, which suggests that $E(N_k) = \frac{A}{N_k} + E(\infty)$, where $A$ is some constant and $E(\infty)$ is the HF energy in the TDL. If you are given $E(N_{k1})$ and $E(N_{k2})$, i.e., the HF energy evaluated using two different k-point meshes, what's your best estimate of $E(\infty)$? Estimate $E(\infty)$ for diamond using the data you generated above.

3. The cohesive energy of diamond is defined as the energy required to break the crystal into free atoms, i.e., $E_{\mathrm{coh}} = E_{\mathrm{atom}} - E_{\mathrm{cell}}/Z$, where $Z$ is the number of atoms per cell. Your calculation in the previous two questions gave you an estimate of $E_{\mathrm{cell}}$. Now calculate $E_{\mathrm{atom}}$ for a free carbon atom (whose ground state is a spin *triplet* and you need to do a *UHF* calculation) using HF and the same basis set. Then calculate $E_{\mathrm{coh}}$ for diamond. A reference value for the HF $E_{\mathrm{coh}}$ can be found in the [literature](https://doi.org/10.1063/5.0119633), which is about 5.4 eV. How does your estimate compare to this reference? The remaining difference arises from (i) basis set incompleteness error and (ii) basis set superposition error (BSSE) in the atomic calculation.

4. Read David Sherrill's excellent [notes on BSSE](http://vergil.chemistry.gatech.edu/notes/cp.pdf). Re-do the calculation of a single carbon atom with ghost atoms from the nearest neighbor to correct for BSSE. How does this change your estimated $E_{\mathrm{coh}}$? (For learning how to use ghost atoms in PySCF, simply google "pyscf ghost atoms".)

5. To correct for the finite-basis error, we can use a larger basis set, e.g., cc-pVTZ, and repeat all the calculations above. An alternative to this brute-force approach is to use some composite corrections. For example, we can evaluate the energy difference between TZ and DZ using some small k-point mesh and use it to correct the DZ/TDL result, i.e.,, $E_{\mathrm{TZ}}(\infty) \approx E_{\mathrm{DZ}}(\infty) + E_{\mathrm{TZ}}(N_k) - E_{\mathrm{DZ}}(N_k)$. Estimate the E_{\mathrm{TZ}}(\infty) result using this approach for e.g., a 3x3x3 k-point mesh. How does this affect your estimated $E_{\mathrm{coh}}$? (To get effect in the TZ basis set, you also need to calculate $E_{\mathrm{atom}}$ using the TZ basis set.)
