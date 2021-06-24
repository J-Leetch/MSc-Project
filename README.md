# MSc-Project

## Progress week 24/06/21:

-   Burgers' equation originally created to mimic turbulence with common advective nonlinearity and diffusion -> allowing defition of Re. However, this was dealt a blow when Hopf/ Cole showed it was explicitly integrable (and therefore not chaotic). This removes the sensitivity to initial conditions so important to N-S in setups with boundaries/ forcings at sufficiently high Re. 

-   So, can we add suitable forcing to restore this property as well as mimicking turbulence? Random forcing brings back the chaotic behaviour, making the process more similar to turbulence.
-   https://arxiv.org/pdf/nlin/0012033.pdf main source for technicalities on burgulence
-   https://doi.org/10.1063/1.4916132 for LES of Burgulence and low-wavenumber forcing. Can be used as a good validation paper and comparison.

Have created a periodic Burgers' solver to replicate the results in LES paper. The forcing adds random energy at low wavenumbers to allow us to mimic the turbulent energy cascade.Image shows the evolution of the equation:

![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/Burgulence.gif)

-   Spent time working out how to match the Re and forcing parameters. Coded using Runge-Kutta 4th order timestepping and 2nd order FD on spatial grid rather than pseudospectral. Now verifying energy spectrum and grid convergence for DNS. After reproduced, I will add the wall back and examine the new energy spectrum/ turbulence characteristics. Then will validate turbulence models using the periodic case (as have the paper as reference) before using similar methods on the BDIM case and observing differences. 

Lots of coding to do this week!

## Update on turbulence mimicking - 17/06/21:

-	The 1D, bounded Burgers’ problem is good from the point of view of:

    - its **simplicity** and cheap computation;
    - convective-diffusive effects similar to N-S;
    - natural wall boundaries.


-	However, the lessons learned from it must be **transferrable** to 3D turbulence. Therefore, there are the following problems with the current test case:

    - the discontinuous forcing in time means that there is no steady turbulent “energy balance” and the fake turbulence is always **decaying**. We want **sustained** turbulence so should add forcing constantly (since the 1D equation can't sustain itself with no dimensions for energy to flow in).
    - additionally, the current forcing is Gaussian and with no spatial correlations - unlike turbulence. Therefore, giving it turbulent-like spatial/temporal statistics would greatly increase the chances of findings which apply to real turbulent flows. The complete roughness of the Gaussian-forced profiles also brings issues with the accuracy of filtering and downsampling the forcing for LES simulations, so this issue could be avoided.


- The goal of the changes would be to:

    - create a new test case for near-wall turbulent flow which **captures the features of the 3D physics** needed for simple but informative investigations into how boundary immersion interacts with turbulence. The foucs is on simplicity, with as little complexity as possible to incorporate the features required.
    
    
- To achieve this I will need to:

    - decide which turbulence statistics to **mimick**. It is probably not necessary to try to accurately replicate 3D turbulence statistics from near-wall regions. However, some basic moments and importantly the energy spectrum should be reproduced - or similar. A suitable Reynolds number should also be chosen to ensure a wide range of scales are present. DNS channel flow/ DNS papers and "Buruglence" papers can be used to inform decision.
    - ensure that the dissipation rate is large enough/ balances the energy added by the forcing. Preliminary tests have shown that the equations can blow up with too aggressive a forcing even in DNS. 
    
    
- After that my objectives are to:

    - calculate accurate DNS simulations when happy with the problem setup;
    - filter the forcing and DNS results to find LES forcing and exact subgrid model. This should be fairly well handled by standard LES subgrid models, and this is a good test to show physical similarity to 3D turbulence problems.
    - identify the additional modelling required and find methods to satisfy it *a posteriori*.
    
    
Other notes:
    - since there is no pressure/continuity in the Burgers' equation, this will always be different to turbulent fluid flow. Therefore, it may be beneficial to look at the modelling for different terms separately to find general principles for modelling these kind of terms rather than looking at the overall LES errors. Most likely, adding the exact LES subgrid stress will not be sufficient in the smoothing region, and other terms will need to be considered in the immersed model. Hopefully, these will be related to immersion of the remaining terms.


## Update on replicating Burgers' problem - 06/21:

See here:
    https://github.com/J-Leetch/MSc-Project/blob/main/1D%20BDIM%20Burgers'/Update%20Burgers'%20BDIM.docx
