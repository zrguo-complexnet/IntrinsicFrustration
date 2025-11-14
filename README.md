# IntrinsicFrustration

Project: Higher-order Interactions Induce Chimera States in Globally Coupled Oscillators

This project reports the discovery of a novel Chimera state in a higher-order Kuramoto model.
We perform both numerical simulations and theoretical analysis of this state.
Details can be found in our paper:
"Higher-order Interactions Induce Chimera States in Globally Coupled Oscillators".

The project consists of 4 parts:

1. Anime
   Run the script to view an animation of the Chimera state we discovered.

2. Simulation
   This corresponds to the numerical simulations shown in Fig. 4 of the paper.

   a. Run the script `simulate_z.m` to obtain the partial order parameters `z_p` and `z_n`,
      representing the oscillators rotating in the positive and negative directions, respectively,
      under different values of K and omega_0.

   b. Simulation results are also included in this folder.
      Files are named like `z1ChimeraSim910D1.mat`, where:
         - `z1` denotes the partial order parameter for positive-frequency oscillators;
           `z2` would represent the negative-frequency ones;
         - `910` indicates the ratio of oscillators in the three synchronous clusters (9:1:0);
         - `D1` means each peak has a standard deviation D = 1.

3. Theory
   This corresponds to the theoretical analysis shown in Fig. 4 of the paper.

   a. `Theory_R31_vs_K_omega0.m` computes the boundary between the chimera state and traveling wave state.
      `Theory_R31_vs_K_omega0_Chimera.m` computes the boundary between the chimera and incoherent states.

   b. `R1_R3.m` is used to numerically determine the functional relationship between
      the partial order parameters R_{1,p} and R_{3,p} in each synchronous group.
      The outputs are saved in files like `R1_910.mat` and `R3_910.mat`.

   c. The results from `Theory_R31_vs_K_omega0.m` are saved in files named like
      `z1_D1_TR910_200_150.mat`, where:
         - z1 denotes the partial order parameter for the positive-frequency group (z2 would correspond to the negative-frequency group);
         - `D1` denotes the standard deviation;
         - `TR` stands for "theoretical result";
         - `910` is the same group ratio as above;
         - `150_200` indicates the range of omega_0 and K values.

4. Supplementary Material
  This folder corresponds to the bifurcation analysis and model extensions presented in the Appendix. It contains the following subfolders:

  a. Bifurcation
  
	- HC_DrawManifold.m: Reproduces Appendix Fig. 4, showing manifold calculations near the homoclinic (HC) bifurcation.
	
	- SN_DrawAndJacob.m: For a chosen parameter set, plots variable trajectories and computes the Jacobian matrix and eigenvalues at the final state.
	
	- SN_Jacob_K.m：Reproduces Appendix Fig. 5, showing eigenvalue evolution across the saddle-node (SN) bifurcation.
	
	- SNLC_FloquetFit.m：Plots fitted curves for the Floquet multipliers.
	
	- SNLC_scanFoquet.m：Reproduces Appendix Fig. 6(a). Floquet multiplier evolution for the SNLC bifurcation of z3-.
	
	- SNLC_scanFoquet2.m：Reproduces Appendix Fig. 6(a). Floquet multiplier evolution for the SNLC bifurcation of z3+.
	
  b. OtherDistribution
  
	- generate_normal_array.m: generates Gaussian samples
	
	- generate_sorted_lorentz_array.m：generates Lorentzian samples
	
	- generate_sorted_uniform_array.m：generates uniformly distributed samples
	
	- simulate_z.m: main script (select frequency distribution inside)
	
	- simulate_oscillators.m：simulation function
	
	- unimodal.m：reproduces Appendix Fig. 2 (dynamical asymmetry with unimodal frequency distributions)
	
  c. OtherInitialCondition
  
	- Running simulate_z.m here produces Appendix Fig. 11:The chimera state induced by intrinsic frustration is robust to initial-condition perturbations.
	
  d. OtherNetwork
  
	- BrainNetwork:
		Using empirical network data and running work.m produces Appendix Fig. 10: Dynamical asymmetry arises on a brain network.
		
	- RandomNetwork:
		build_flag_complex_and_save：constructs a higher-order random network
		
		randomFlagComplex：helper function used for construction
		
		run_highorder_kuramoto_on_flag：runs the higher-order Kuramoto model on the generated network, producing Appendix Fig. 9: Dynamical asymmetry arises even in Erdős–Rényi random networks.
		
  e. OtherOrder
  
	- AddOrder1or2
	
		Asimulate_z_K1: omega0=50; evaluates chimera formation when adding first-order coupling (Appendix Fig. 7, left panel).
		
		Bsimulate_z_K2: omega0=50; evaluates chimera formation when adding second-order coupling (Appendix Fig. 7, right panel).
		
	- Order4
	
		Csimulate_z_O4: evaluates chimera formation in a fourth-order Kuramoto model (Appendix Fig. 8).
		
Environment:
- MATLAB R2021a or later
Optional:
- MATLAB Parallel Computing Toolbox is required if you want to use parpool 
  for large-scale simulations (recommended).
- At least 8 GB RAM is recommended; 16 GB+ for flag complex simulations.


Folder Structure:
- /Anime                   : Chimera animation script
- /Simulation              : Simulation scripts and .mat result files
- /Theory                  : Theoretical analysis scripts and .mat result files
- /Supplementary Material  : Appendix-related analyses and extended numerical results.

How to Run:
1. View animation:
   Go to the `Anime` folder and run the animation script.

2. Run simulation:
   Go to the `Simulation` folder and run `simulate_z.m`.

3. Perform theoretical analysis:
   Go to the `Theory` folder and run `Theory_R31_vs_K_omega0.m` or `Theory_R31_vs_K_omega0_Chimera.m`.

Citation:
If you use this code in your research, please cite our paper:
"Higher-order Interactions Induce Chimera States in Globally Coupled Oscillators".

License:
This project is distributed under the MIT License.
