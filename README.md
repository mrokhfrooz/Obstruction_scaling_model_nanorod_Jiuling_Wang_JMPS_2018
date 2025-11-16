# Obstruction_scaling_model_nanorod_Jiuling_Wang_JMPS_2018
Integration of Obstruction_scaling model developed by Wang et al. 2018

📘 Overview

This repository contains a MATLAB implementation of the model developed by:

Wang, Jiuling, et al. "Diffusion of rod-like nanoparticles in non-adhesive and adhesive porous polymeric gels." Journal of the Mechanics and Physics of Solids 112 (2018): 431-457.

The code computes the expected number of fibers intersecting a rod-like nanoparticle in a polymer network and the corresponding relative diffusivity:
κ=exp(−N)

where

N=N1+N2+N3+N4

Each N_i is computed using geometric integrations derived in the paper (main text + Appendix A).

This repository reproduces the model shown in Figure 7a of Wang et al. (2018).

🧮 What the MATLAB code does

The script:
*Defines physical and geometric parameters
*Fiber density
*Rod length
*Rod radius
*Aspect ratio

Computes 
N_1 ​ from Eq. (9) 
N_2 ​ from Eq. (A2) 
N_3 ​ = N_3a ​ + N_3b ​ (Appendix A) 
N_4 ​ = N_4a ​ + N_4b ​ (Appendix A)

Calculates total obstruction count:
N_total​=N_1​+N_2​+N_3​+N_4

Computes the relative diffusivity: 
𝜅 = exp(-N_total)

▶️ How to Run

Open MATLAB
Place the script in your working directory
Run the code.

The script prints:

Each contribution N_1 ​ ,N_2 ​ ,N_3a, N_3b, N_4a, N_4b 

Their sum 

The diffusivity ratio κ

🧱 Dependencies

MATLAB R2019b or later
No toolboxes required
Script uses only built-in functions (sin, acos, sqrt, trapz, etc.)


📌 Notes

The integration scheme follows the trapezoidal approach described in Appendix A of the paper.

Numerical guards (e.g., small-argument corrections near trig singularities) are included to ensure stability.

Parameters are chosen to match Figure 7a of the original publication.

📧 Contact

For questions or suggestions, feel free to open an issue or pull request.
