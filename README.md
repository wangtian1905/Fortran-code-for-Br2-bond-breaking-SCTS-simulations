# Fortran-code-for-Br2-bond-breaking-SCTS-simulations
This program code package is used to simulate photoelectron momentum distributions from tunnel ionization of dissociating Br2 molecule. The algorithm is basing on Semiclassical Two Step (SCTS) model for tunnel ionization. And the simulated system includes the full states of Br2 bond-breaking: the ground Br2 molecule, bond breaking state (BBS), and the ground Br atom.  All codes are run on Intel Fortran compiler with Intel OpenMPI configuration.

## Recommended System requirements: Computer with Intel processors (8 gen Core or higher), Windows11, Visual Studio 2022 or later releases (version: Visual Studio Community)
Installation time of Visual Studio (VS): about 40 mins.
## Environment setup: After installation of VS, then install Intel oneAPI Base Toolkit and Intel oneAPI HPC Toolkit one after one. The links for these two toolkits are below: 
Intel oneAPI Base Toolkit: https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html
Intel oneAPI HPC Toolkit: https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html
Installation time: 30-60 mins for oneAPI Base Toolkit and 45-90 mins for oneAPI HPC Toolkit. The toolkits will automatically be integrated in VS.
## Parameters setup for Intel Fortran compiler: 
(a) create a Fortran project in VS and paste the code. 
(b) go to project properties of the project drop-down manu; under Fortran\Language, set "Process OpenMP Directives" to the value: "Generate Parallel Code (/Qopenmp)".
(c) Under Linker\system, set values of both "Stack Reserve Size" and "Stack Commit Size" to 50000000 or higher.
## Now, after build solution, the code will be compiled and can be run by select a runing mode under Debug drop-down manu. For efficient calculation, please select "Start without debugging".

# Parameters setting guide for SCTS simulation

## Laser pulse parameters: within module "lsr_para" in top. 
Wavelength: L_L (in unit meter)
Intensity: L_I (in unit PW per cm squre)
Pulse shape: half trapezoid defined by plateau cycle number L_PC and attenuation cycle number L_TC.

## Simulation parameters: within module "sim_para"
Ensemble size: gn (recommended to set 500000000 or higher)
Tunnel ionization windows: tnl_t_lb, tnl_t_ub, and tnl_t_lb2, tnl_t_ub2. These two sets of parameters allow to study the interference between any two defined time windows. They are set for the same time window as default. The window length for multicycle simulation is 4 and single cycle simulation is 1.  

## Symplectic numerical integrator parameters: within module "Sym_ctrl"
The radius criterions r_crt1 and r_crt2 devide the real space to three regions, in which three different time steps h1,2,3 are applied to accelerate the integration speed within a fair accuracy.  

## Open MPI parameters: within module "omp_para"
The total number of parallel threads: n_threads (The value should below the maximum cpu threads of the processor)

## Molecular system setting parameters: within module "MO_para"
### The physical property of the simulated system is set by the parameters:
(a) charge on each nuclei: nn_Z1, nn_Z2
(b) Inernuclear distance: n_R_Atm, n_R_Mlc (set to zero for atomic case)
(c) ionization potential: Ip1,Ip2
(d) polarizability: apha_II_Atm, apha_II_Mlc
### The relative contributions of two trajectory groups are set by the parameter: prb_orb1. 
This parameter defines the relative probability of one trajectory group to another. For example, in the bond-breaking state, the bound charge contribution of atomic (Atm) Br+ is equal to the diffuse charge contribution of molecular (Mlc) Br2+, so prb_orb1=0.5. While, when study the effect of ionization potential contribution from a pecific orbital (\pi_u orbital or 4p orbital), prb_orb1 can be set any value within (0,1).
