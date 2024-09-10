# INCAR Tags
A full list of INCAR tags can be found [here](https://www.vasp.at/wiki/index.php?title=Category:INCAR_tag&pageuntil=LREAL#mw-pages).
- `SIGMA`: 0.03 if system unknown (metal, cond, semicond), insulators/semiconductors = 0.1, metals = 0.05
- `ISTART`: 0 if new job, 1 for continuation jobs (i.e. band structure and DOS, which obtain CHGCAR from the SCF calculation)
- `IBRION`: difficult relaxation problems it is recommended to use the conjugate gradient algorithm (IBRION=2), which presently possesses the most reliable backup routines. IBRION = -1, The ions are not moved, but NSW outer loops are performed. In each outer loop the electronic degrees of freedom are re-optimized (for NSW>0 this obviously does not make much sense, except for test purposes). If no ionic update is required use NSW=0 instead
- `ISYM`: whether or not to enforce symmetry upon relaxation. If the desired symmetry is known (i.e. IISF = 2 IBRION = 2), then use ISIF 1 | 2 (2 is more efficiently parallelized). If the desired symmetry is not known (i.e. ISIF = 3), then use ISYM = 0.
- `EDIFFG`: when positive, this is the condition to break the ionic relaxation loop when the change in energy is < EDIFF. When negative, this is the break condition when the total net force on an atom is < |EDIFFG|. A good rule of thumb is to select -EDIFFG no larger than 0.05 ev/Å, though in many cases this may be entirely too large, especially for flexible materials. MOF materials often use 0.03 eV/Å, whereas many recommend using 0.01 eV/Å. For vibrational calculations however, smaller forces on the atoms may have to be used instead, though this can greatly increase the computational time. 0.001 eV/Å is more than low enough for phonon calculations, but will certainly take much longer.


# Relaxations
- The Pulay stress arises from the fact that the plane-wave basis set is not complete with respect to changes in the volume. Thus, unless absolute convergence with respect to the basis set has been achieved; the diagonal components of the stress tensor are incorrect. This error is often called "Pulay stress". Pulay stress can be ignored if the relaxation is volume-conserving, but if it is not, then the structure needs to COMPLETELY RELAX at least 3 times to reduce Pulay stress error (i.e. the "reached required accuracy - stopping structural energy minimisation" message appears).

  
# Converging structures
- The relaxation can meet the required EDIFF and EDIFFG energetics convergence criteria, but it is crucial to check the convergence of forces per atom on the structure. A good convergence criteria is 1 meV-- anything larger means the structure is unhappy. This is especially true if forces are large on the inorganic structure, because then it is very unhappy (it's okay if the ligands are a little unhappy because honestly when are they ever happy).
- CONVERGING EXPERIMENTAL STRUCTURES: Sometimes, an experimental structure can be relaxed with just ISIF = 2 and selective dynamics on the hydrogens. However, when the total force in eV/Å > 0.01 eV, especially on the atoms composing the inorganic structure, then further relaxation is required. Often, this can be resolved by maintaining the unit cell but removing selective dynamics such that all of the atoms in the structure are now relaxing. This can help to reduce the force stress per atom.

# Bandstructures
- [Steps to create a bandstructure](https://www.molphys.org/VASP/band_structure.html) (needs CHG* in addition to regular post-relaxed SCF files)
- [Example of creating a cubic diamond Si band structure using VASP and Pymatgen](https://ma.issp.u-tokyo.ac.jp/en/app-post/1146) (includes a .py file to create the reciprocal lattice KPOINTS file in the Brillouin zone)
- [Tess' code to plot a band structure using color codes for the organic and inorganic parts](https://gist.github.com/blondegeek/0941ccb193e2aeae71bf07ca006c5cca)
- NOTE: To create a good band structure, you need the same NBANDS and ENCUT as grepped from the OUTCAR of the SCF calculation of the post-relaxed structure. Be sure to run an SCF calculation with IBRION=-1 and ISMEAR=0 to obtain a good CHGCAR and CHG for the bandstructure calculation. The SCF calculation must be run on the primitive unit cell for the band structure otherwise it will throw an error for the charge density!

# Partial charge densities
- [Band-decomposed charge densities
 (VASP)](https://www.vasp.at/wiki/index.php/Band-decomposed_charge_densities#:~:text=The%20partial%20(band%2Ddecomposed),specific%20region%20in%20real%20space)
- [PARCHG output file](https://www.vasp.at/wiki/index.php/PARCHG)

# DOS
- [Steps to create a DOS calculation](https://molphys.org/gnuplot_tutorial/dos_plot.html) (needs CHG* in addition to regular post-relaxed SCF files)
- [notebook for PDOS by orbital](https://gist.github.com/lan496/ee0bd7a52df99029ac0aacbe69f2bf57)
- NOTE: Make sure that the ENCUT value in the DOS INCAR is set to the same as the ENCUT value in the OUTCAR of the SCF calculation, otherwise you may get an error stating that the "dimensions on CHGCAR file are different" or that the calculation cannot read a charge density of a file greater than ICHARG>10.

# Phonon calculations with Phonopy
- [How to run Phonopy on a terminal](https://materials-lab.io/blog/2021/11/19/phonopy-guide/) (same as [VASP Phonopy documentation](https://phonopy.github.io/phonopy/vasp.html))
- [Python Phonopy](https://phonopy.github.io/phonopy/phonopy-module.html)
- [How to set up VASP phonon calculations](https://rehnd.github.io/tutorials/vasp/phonons)
- in order to generate phonons, make sure to base off a structure that is fully relaxed: EDIFF = 1e-7, EDIFFG <= 1e-4, ENCUT >= 700. ISIF can be either 3 or 2, depending on what kind of relaxation you want to perform. If ISIF = 2, then set ISYM = 1|2 as the desired symmetry of the structure is known. If ISIF = 3. then set ISYM = 0 as the current symmetry of the structure may not be correct. Then, enforce symmetry by using the CONTCAR of the final relaxation and running another relaxation where ISIF = 2, IBRION = 2, and ISYM = 2.
- Finite differences allow computation of VDW interactions, where as Density Functional Perturbation Theory (DFPT) does not. DFPT is where the electric field perturbations enter the hamiltonian and thus is calculated self consistently in the electronic and geometric optimization. This allows for several physical properties to be extracted (dielectric tensor, polarizabilities, Born charges, etc.) besides the phonon frequencies. Computational limitations arise in DFPT to be restricted to small unit cells, while finite differences allow handling much larger unit cells, but normally finite differences and DFPT are quite comparable.
- 

# Effective Mass
- [effmass documentation](https://effmass.readthedocs.io/en/latest/API%20documentation.html)
- [effmass example notebook](https://nbviewer.org/github/lucydot/effmass/blob/master/tutorials/Tutorial.ipynb)
