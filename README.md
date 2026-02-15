# Fortran-code-for-Br2-bond-breaking-SCTS-simulations
This program code package is used to simulate photoelectron momentum distributions from tunnel ionization of dissociating Br2 molecule, as a part of manuscript: "Multi-messenger tracking of coherence loss during bond breaking". The algorithm is basing on Semiclassical Two Step (SCTS) model for tunnel ionization. And the simulated systems are the full states of Br2 bond-breaking: the ground Br2 molecule, bond breaking state (BBS), and the ground Br atom. Since the identical algorithm, all simulation codes have the same structure, and the modeling of different dissociation state is realized by setting the physical system parameters. All codes are run on Intel Fortran compiler with Intel OpenMPI configuration to acheive the best optimized calculation speed.

## Recommended System requirements: 
- Computer with Intel processors (8 gen Core or higher) <br>
- Windows10,11 <br>
- Visual Studio 2022 or later releases (version: Visual Studio Community) <br>
Installation time of Visual Studio (VS): about 40 mins. <br>
## Environment setup 
After installation of VS, then install Intel oneAPI Base Toolkit and Intel oneAPI HPC Toolkit one after one. The links for these two toolkits are below: <br>
- Intel oneAPI Base Toolkit: https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html <br>
- Intel oneAPI HPC Toolkit: https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html <br>
Installation time: 30-60 mins for oneAPI Base Toolkit and 45-90 mins for oneAPI HPC Toolkit. The toolkits will automatically be integrated in VS. <br>
## Parameters setup for Intel Fortran compiler: <br>
(i) create a Fortran project in VS and paste the code. <br>
(ii) go to project properties of the project drop-down manu; under Fortran\Language, set "Process OpenMP Directives" to the value: "Generate Parallel Code (/Qopenmp)". <br>
(iii) Under Linker\system, set values of both "Stack Reserve Size" and "Stack Commit Size" to 50000000 or higher. <br>
Now, after build solution, the code will be compiled and can be run by select a runing mode under Debug drop-down manu. For efficient calculation, please select "Start without debugging". <br>

# Parameters setting guide for SCTS simulation
In the simulation code, the top part shows the parameter modules, which set the simulated physics model and control the program implementation. An example module for laser control is shown: <br>
```fortran
module lsr_para
implicit none
......
end module  
```
Below are the instructions for the key parameter modules:
## Laser pulse parameters (within module "lsr_para" in top) <br>
- Wavelength: L_L (in unit meter) <br>
- Intensity: L_I (in unit PW per cm squre) <br>
- Pulse shape: half trapezoid defined by plateau cycle number L_PC and attenuation cycle number L_TC. <br>

## Simulation parameters (within module "sim_para") <br>
- Ensemble size: gn (recommended to set 500000000 or higher) <br>
- Tunnel ionization windows: tnl_t_lb, tnl_t_ub, and tnl_t_lb2, tnl_t_ub2. These two sets of parameters allow to study the interference between any two defined time windows. They are set for the same time window as default. The window length for multicycle simulation is 4 and single cycle simulation is 1. <br>  

## Symplectic numerical integrator parameters (within module "Sym_ctrl") <br>
The radius criterions r_crt1 and r_crt2 devide the real space to three regions, in which three different time steps h1,2,3 are applied to accelerate the integration speed within a fair accuracy. <br>  

## Open MPI parameters (within module "omp_para") <br>
The total number of parallel threads: n_threads (The value should below the maximum cpu threads of the processor)

## Molecular system setting parameters (within module "MO_para") <br>
### The physical properties of the simulated system is set by the parameters <br>
- charge on each nuclei: nn_Z1, nn_Z2 <br>
- Inernuclear distance: n_R_Atm, n_R_Mlc (set to zero for atomic case) <br>
- ionization potential: Ip1,Ip2 <br>
- polarizability: apha_II_Atm, apha_II_Mlc <br>
### The relative contributions of two trajectory groups <br> 
 The contribution of the first trajectory group is defined by parameter prb_orb1. For example, in the bond-breaking state, the bound charge contribution of atomic (Atm) Br+ is equal to the diffuse charge contribution of molecular (Mlc) Br2+, so prb_orb1=0.5. While, when study the effect of ionization potential contribution from a pecific orbital (\pi_u orbital or 4p orbital), prb_orb1 can be set any value within (0,1). <br>

By setting the molecular system parameters, the tunnel ionization dynamics can be established for the neutral dissociation states: Br2 molecule ground state, Br2 bond-breaking state (BBS), and Br atom ground state. The codes for these states are shown in the repository. <br> 

## Output simulation results <br>
The output electron momentum distribution is written as a matrix in the data file: "Spec_Z.txt". The two files "Spec_gridX.txt" and "Spec_gridY.txt" provide the coordinate axes for the matrix. The demo result datasets are shown in the repository. <br>
Typical run time for a simulation (on laptop computer, Intel i7-8750H, use 10 cpu threads): about 10 hours for the ensemble size of 500000000. <br>

## results analysis by normalized differential (ND) plot
The simulated electron momentum distributions from different dissociation states can be further analyzed by doing normalized differential plot between a signal disttribution and a reference distribution. The results analysis code is: "SCTS_spec_ND_plot_fan_area_analysis.f90". The demo analysis results on Br2 ground molecule ionization (ref.) and Br2 BBS ionization (sig.) are shown in the repository. The ND plot matrix is written within the files named after "ND_spec", and its coordinate axes are the same with that of original electron momentum distributions. The radial momentum distributions of the ND plot in two directions are written in the files named after "diff_r_dis_test".
