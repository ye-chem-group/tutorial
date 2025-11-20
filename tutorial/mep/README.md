# Creating a MEP map with PySCF and Chimera

### Definition

The **molecular electrostatic potential (MEP)** is the Coulomb potential created by the total charge density (electron + nuclei) of a molecule:

$$v(\vec{r}) = \int \mathrm{d}\vec{r}_1~\frac{\rho_{\mathrm{tot}}(\vec{r}_1)}{|\vec{r} - \vec{r}_1|}$$

The value of MEP at a point in space can be understood as the force felt by a positive probe charge. Alternatively, since electrons are what chemists care about, the *negative MEP value* at some point in space is the force felt by an electron.

A **MEP map** is the plot of MEP on a 3D surface (usually chosen to be the electron density isosurface; *vide infra*). Shown below is the MEP map created using this tutorial for a CH3OH molecule. The red, white, and blue regions correspond to MEP < 0, MEP = 0, and MEP > 0. We see that the red and blue regions correspond well to electron rich and electron poor regions of the molecule based on our chemical understanding.

<img src="https://github.com/ye-chem-group/tutorial/blob/main/tutorial/mep/.resources/mep_ch3oh.png" width="400">


### Calculate MEP using PySCF

In PySCF, the MEP can be calculated using the `mep` function from `pyscf.tools.cubegen` ([link](https://github.com/pyscf/pyscf/blob/master/pyscf/tools/cubegen.py#L158)). The following script calculates the MEP for a methanol molecule (CH3OH).
```Python
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
```
This script Gaussian cube files for both the electron density and the MEP. The electron density will be used to create an isosurface, on which the MEP will be plotted.

### Plot MEP using Chimera

Chimera ([link](https://www.cgl.ucsf.edu/chimera/download.html)) is a very good software for visualizing molecular structure and properties, including MEP. The following steps follow primarily the tutorial by Cadillac Chemistry ([link](https://www.youtube.com/watch?v=cE5YY77XXUs)). The first part of the video is about using ORCA to generate the electron density and MEP cube files, which we already did using PySCF. So you can jump directly to 6:19 of the video for the Chimera part.

##### Step 1: Load molecular structure
1. Load xyz file: `File -> Open -> ch3oh.xyz`
2. Change to ball & stick style: `Actions -> Atoms/Bonds -> ball & stick`
3. Change color: `Actions -> Color -> By element`

<img src="https://github.com/ye-chem-group/tutorial/blob/main/tutorial/mep/.resources/mol.png" width="400">

##### Step 2: Create electron density isosurface
1. Load `den.cube`: `Tools -> Volume Data -> Volume Viewer -> File -> Open map -> den.cube`
2. Adjust `Level` to `0.02`, which roughly corresponds to the vdW surface

<img src="https://github.com/ye-chem-group/tutorial/blob/main/tutorial/mep/.resources/den.png" width="500">

#### Step 3: Color the density isosurface with MEP
1. Load `mep.cube`: `Tools -> Surface color -> By electrostatic potential -> browse -> mep.cube`
2. I recommend setting `MEP = -0.03, 0, +0.03` to be red, white, and blue. A negative/positive MEP means electron rich/poor region.

<img src="https://github.com/ye-chem-group/tutorial/blob/main/tutorial/mep/.resources/mep_color.png" width="600">

We see that the MEP for CH3OH is red near the oxygen atom and blue near the hydrogen atom in the OH group, which is consistent with our chemical intuition that oxygen is electron rich and hydrogen is electron poor in an OH group.
