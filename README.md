# MSc-Project

## Correlations in ImLES - 29/07/21:

Before diving into the ML, important to work out exactly what we need to model and its characteristics. This requires a detailed summary of the modelling required and related statistics:

### Conventional subgrid stress modelling still works everywhere in ImLES

Can investigate the effectiveness of Smagorinsky modelling of the subgrid stress across the entire domain. Plot the correlation between the Smagorinsky model and the perfect Tij term at different spatial locations. This looks like e.g.:

![Figure_1](https://user-images.githubusercontent.com/46157966/127565972-3edf6483-9b74-44f7-b86c-415836acdde7.png)

![Figure_1](https://user-images.githubusercontent.com/46157966/127566020-c25d5ee3-165d-46a8-9517-94907c15cf8e.png)

with a quality of linear fit which varies like

![Figure_1](https://user-images.githubusercontent.com/46157966/127566123-515b5fdb-fa74-40bd-9a4e-0dd112a8892f.png)

The slope of the linear fit is also not altered significantly:

![Figure_1](https://user-images.githubusercontent.com/46157966/127566276-0a9189ec-f277-4f8a-a5ef-4fec6b569746.png)

demonstrating that subgrid modelling for the convective stresses is working well even in the smoothing region. (note to Gabe: the reason it didn't look so good before was I forgot to apply the BDIM operator to the perfect Tij term when calculating it). This is really nice, since it gives evidence to support the proposal that we shouldn't need to modify the LES modelling we all know and love too significantly.


### Subgrid stress is not enough - there is also a wall model term 

We can see this qualitatively through looking at the snapshot plots of the perfect closure with the LES subgrid stress subtracted off. This is the initial procedure to **isolate the wall model**.
![dudx correlation with wm snapshot 2](https://user-images.githubusercontent.com/46157966/127566560-8579cdc7-b322-426b-954f-ac808983aabd.png)

Now, we can see that this term has a good correlation with the velocity gradient (only in the wall smoothing region of course). This is quantified by looking at density (scatter) plots between the two terms:

![density -0 5](https://user-images.githubusercontent.com/46157966/127566600-6459cf2f-f5e9-41a2-97f5-76412f2cc7a6.png)
![density 0 02](https://user-images.githubusercontent.com/46157966/127566608-32339ec1-27bb-45d6-a2e7-62c12f89404c.png)
![density 0 53](https://user-images.githubusercontent.com/46157966/127566615-2782e4cb-b79b-406b-9f9e-fa73b939fca3.png)
![density 1 04](https://user-images.githubusercontent.com/46157966/127566621-ba53d932-a3b2-4d3d-bb32-b597c3bd4238.png)


and then evaluating the Pearson r value for a linear fit. We could fit a polynomial too if needed, since it seems like the data collapses well onto a slightly nonlinear region. The Pearson r value is high here and the slope of the linear fit varies following a simple form through the smoothing region.

![dudx correlation Pearson r vs yeps](https://user-images.githubusercontent.com/46157966/127566684-8e2ff3b2-d583-421c-9f63-bf36f286e39c.png)

![dudx correlation linear fit slope vs yeps](https://user-images.githubusercontent.com/46157966/127566722-8d040c76-2369-461f-85c2-951c6db4c8a1.png)

This shape of the slope variation fits the kernel shape when seen on this grid resolution as shown below. (Another objective is to downsample the LES to different, staggered grids to fill out the resolution in these plots.)

![kernel](https://user-images.githubusercontent.com/46157966/127566816-ba8319d2-8300-43d2-9ded-f9b6a9faa8ec.png)

Therefore, on first glance it would seem the wall model is very simple to model. **BUT** these above results are found using the DNS timestep (effectively an instantaneous time or limit to small time). After coarsening the timestep by a factor of 10, we see a different picture. Correlations are weaker, although the wall model term is smaller due to the larger timestep, and the nice form of the linear fit slopes seen above is distorted. However, (note to Gabe again), the effect is not hugely dramatic (phew..) as seen in our meeting (due to the same coding error as above), which I guess is a good thing.

This brings up a main point of curiosity for the week; we have terms which we calculate based on our velocity field which do not change at all dependending on the timestep used. However, when we look for a perfect closure, it is obvious that this is dependent on the timestep used. There is clearly modelling to be done which considers the state of the flow relative to the grid resolution **and** the time resolution. I think this is the space to be filled by the machine learning models. Using information about our strong, instantaneous time correlations and the state of the flow field, it seems like a good task to give ML to correct our trajectory through time for the massive timesteps we want to take in LES.

However, whilst talking about the difference between perfect closure and the terms we write down in the equations, I need to mention this observation for context:

This plot shows the difference seen between perfect closure and the analytic Tij term where we take a coarsened timestep: #

![snapshot Cij closure 3](https://user-images.githubusercontent.com/46157966/127568059-fa25906b-6b22-48ef-9c84-dfd29ab91fdf.png)

The timestep must have an effect on the degree of dynamic filtering occuring in the simulation, which is not usually considered in analytics. In "classic" LES equations we have no way of considerating this time component. However, when considering ImLES we can use the immersion stress terms to investigate since I think they consider grid filtering and so can possibly translate this dynamic filtering error into something analysable analytically.

Now, Tij is the subgrid stress term calculated by filtering DNS data and the equation
![image](https://user-images.githubusercontent.com/46157966/127562789-184c70ac-fcdd-410b-9e70-f142c2de7918.png)

and Cij is the term 
![image](https://user-images.githubusercontent.com/46157966/127562680-1bdadce3-c476-49ed-97fe-396d13d00903.png)
also calculated as defined using DNS data and filtering.

In the central channel, without the wall model to confuse us, we can think about the correction required to add to Tij so that it matches the perfect closure for the timestep in question. It turns out this correction is very well correlated with the convective subgrid stress Cij: (this plot is made at roughly mid-channel)

![density Cij3 to Tij correction 3](https://user-images.githubusercontent.com/46157966/127568654-7751e1a5-01c1-494f-b068-8f34ead8ff1a.png)

where the y axis is the correction term and the x axis is Cij. This effect is also shown in the snapshot plot above (#) where Cij is plotted and the collapse is demonstrated qualitatively in a snapshot. This seems too strange to be a coincidence since it is found from stepping forward a simulation, including numerical error etc and looking at the perfect closure (relative to DNS) vs the subgrid stress which we would expect to summarise the modelling required based on any normal version of LES equations.

From ImLES equations, this makes sense since Cij appears as a term; immersion stress just recognises that the LES field is still not well resolved on the grid used. I'm not sure where the factor of ~1/3 comes from though!? Could be something to do with having chosen the two filtering operations to match each other, when in reality they should not.

It turns out this analysis using Cij might also be relevant to the wall model and its isolation from the other modelling

Naturally, we want the wall model to be confined to the wall smoothing region. As shown by the plot below, this is achieved by utilising Cij in the isolation of the wall model. The thinking here is that based on the observations seen above - far from the wall where we assume there is no wall model term present at all - Cij must be combined with Tij to give a perfect closure for the subgrid stress. Therefore, we should remove this perfect subgrid stress from the perfect closure at the wall to leave the wall model term. The plot below shows that this is in fact effective, and the wall model is negligble outside the smoothing region when its isolation is treated this way.

![Figure_1](https://user-images.githubusercontent.com/46157966/127569770-5864d4c0-bc2f-437f-8f77-efe7462ca9c4.png)


### Summary

-   A machine learning model based on top of these already good correlations is the way to go. We know that short-time correlation is excellent for the wall model, but perhaps stabilisation is needed using ML. When we take tiny timsteps, the wall model is a HUGE term, and I can imagine it might be an unstable addition. However, with longer timesteps, the term can be made smaller but the correlations worsen. Maybe there is a compromise between reliance on the ML model to compensate for larger correlation errors with longer timesteps and, on the other end of the spectrum, the need for the ML model to somehow maintain stability if we were to take small LES timesteps.

-   I have found some surprising application of the ImLES ideas in digging deeper into how subgrid modelling is working. I think it is suprising that Cij pops up and correlates so well in fixing the conceptual gap between the Tij term and the real, perfect subgrid model. I would like to think more about the meaning of this term and how it is practically relevant.


## Time to let the computer find a solution? - 22/07/21:

After trying to manually find a closure for the ImLES modelling needed on top of the LES subgrid stress - now more convinced that it is not a simple function. I don't feel too guilty about this anymore after remembering the close similarity between this and standard LES wall modelling and the fact that lots of people are using AI to tackle that problem now.

This brings up something to remember for the thesis; it would be very nice to be able to somehow quantify the "ease of modelling" for the current method compared to LES wall modelling. If the modelling is more easily achieved in the ImLES form, that is a huge win since we get the ability to use complex boundary conditions for free in terms of cost. 

-> It is appropriate to now move onto "breaking out the big guns" with ML. However, need to show this is necessary in the thesis to avoid perpetuating the "throw AI at it to fix it" attitude... Therefore, plan to include statistical plots:

-   Single-point correlations between velocity/velocity gradients in the smoothing region and the additional modelling term.
-   Some attempts at fixing a parameterised form of the wall model term and then looking at the correlation of the parameters with flow field data.
-   Need to give a reasonable guess for the best simple correlation type closure and then show that it is not very good.

Now thinking about the specifics of an ML "wall model". The key part will be determining the input and output data. This will involve some trial and error. Can again use parameterised forms for the model term or let the ML model independently add forcing at each grid point. 

When investigating the additional modelling I noticed that the LES timestep used had a large impact on the modelling required. Combined with this, there is the point made by https://arxiv.org/pdf/2106.11144.pdf on modelling based on training not just with single timesteps; this demonstrates the temporal picture is quite important to not forget and allows the model to take account for compounding errors/other numerical errors so that more confidence can be had in the a posteriori stability. This approach naturally lends itself to reinforcement learning (RL) since we can create an environment where the agent learns in a run-through of a simulation, but is not necessarily restricted to just RL. The arxiv paper gives little justification for why they used RL and we haven't seen any paper comparing the effectiveness of RL methods to deep learning that takes into account multi-timestep dynamics - so that could be an interesting niche to check out.

Plan for this week -

-   Produce final statistical plots for the wall modelling which should be very close to what I will present in the final report.
-   Integrate code for learning a model. For deep learning, the data can just be saved and the ML done separately. However, for RL I will code up an environment for the agent to learn in. Hopefully this shouldn't be too bad since I have fortunately done this before! Initially, will use a basic deep Q-learning setup (because I am confident I know exactly how that works, unlike some of the state of the art RL seen in the arxiv paper!). Best to keep it all as simple as possible since I am researching the CFD, not the ML.

## Progress report - 15/07/21:

-   Part of the process here is to characterise this new test case. This kind of documentation will be done with the main terms/KE and immersion stresses so everything can be visualised nicely since it is not something most readers will be used to at all. That involves non-dimensionalising everything and presenting nicely in the thesis. Example of this:

![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/ImLES/u.png)

In terms of non-dimensionalisation: space x/L (by channel width of course); time can be non-dimensionalised by the forcing timescale (which is just t=1, you can see this timescale in the DNS plot of u in the group chat); still thinking about the velocity scale though I have the forcing coefficient (although non-D I think) and the viscosity to use, so I could use a local Reynolds number u x grid_spacing/viscosity and then divide by the forcing coefficient to give scaling in terms of the forcing.

Can now visualise modelling that is required additional to LES subgrid stress like so:
![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/ImLES/diff%20log.png)

Have found there is this large peak in the smoothing region (for most of the time). These plots below show the correspondence between the velocity field near the wall and the "delta model" required. Here, a log x scale is used to highlight structure in the smoothing region, and a "symlog" scale is used for the colormap to allow us to see the wide range of the additional model term. (yes, the colormap definitely needs work matching between the plots to improve the visualisation..)

![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/ImLES/near%20wall%20u.png)
![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/ImLES/near%20wall%20perfect.png)

This demonstrates qualitatively that the delta model is correlated with the wall gradient in the velocity field, producing a large, discontinuous jump in the diffusive term, and leaving extra modelling to be found as predicted in IP. It seemed that the implicit modelling for BDIM was enough to format the resultant error so that it takes the approximate form of the kernel - hence the "delta model". The same is observed here in this test case on a real simulation. 

Provided this is the actual reason for this phenomenon, the way to model it is to predict the wall gradient in the DNS flow. This is nice, since this is necessary for evaluating real wall stresses from an ImLES simulation anyway. My first instinct in terms of finding a way to do this is by examining the second derivative in the ImLES flow at the outer edge of the smoothing region. This idea is illustrated in the figure below:

![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/ImLES/wall%20deriv%20from%202nd%20deriv%20ImLES.png)


## Progress report - 08/07/21:

### Pt.1

Time to make very clear distinction between different cases:

-   **Wall-resolved**: can be DNS or LES, filtering for aposteriori analysis is easy since we can simply downsample the DNS grid to the LES grid (given that the number of grid points is set up correctly).
-   **BDIM**: aims to compute a solution which matches the DNS *(with $\epsilon$ -> 0)*. Grid resolution can be changed independently of the smoothing filter width if wanted and grid convergence can be achieved independently of $\epsilon$ convergence. However, wide grids compared to flow features lead us into more of an LES simulation. If we are in "LES mode", we can't expect to recover a close approximation to the DNS solution using BDIM, so BDIM aren't the right equations to use.
-   Therefore, we are led to the **immersed LES** perspective. Instead, we try to compute the wall-filtered solution directly, with the introduction of additional modelling to make sure we account for the smoothing of functions (BDIM effects) and the smoothing of the field (LES effect). It is possible that there is a degree of implicit modelling such that BDIM achieves this objective already, but realising this allows us to try to find more accuracy. Here, the methods are combined to give a solution with fundamentally different character. This ties together the two length scales: smoothing region width and LES filter width. The two scales must match (or are the same thing 😉) since the wall-filtered solution must have its function transition occur across the same region as the underlying DNS is extended by the LES filtering. This could be a key realisation for attempting to model this problem as well... 

![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/key%20plots/ImLES%20scale%20equiv.png)

Key points this week with regard to above:

-   To investigate models for ImLES equations, the wall-filtered field is needed. Since we need reasonable resolution in the smoothing region, downsampling is not ideal. Therefore, choose an explicit filtering to apply to DNS field to give a "target solution" for ImLES. Therefore, there is dependance on the filter chosen. Seems natural to choose the cosine filter used by the BDIM kernel so that the target matches the only explicit filter that is involved in a simulation.


Other points:
-   There is the need to consider how body forces are evaluated from an ImLES solution.

### Pt.2

Interesting result examining the effectiveness of Smagorinsky model in BDIM vs wall-resolved simulations:

![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/key%20plots/tuned%20smagorinsky%20BDIM.png)
![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/key%20plots/tuned%20smagorinsky.png)

FRS = fully resovled simulation (DNS); ROM = reduced order model (LES grid but with no model); Smagorinsky optimum (LES simulation with globally optimised - aposteriori - Smagorinsky coefficient)

This showed that the Smagorinsky model can do a better job at correcting sub-grid error in the immersed boundary (IB) context; this is backed up by the higher correlation coefficient between the perfect closure and Smagorinksy terms in the IB case. However, should be taken with a pinch of salt; this was not an ImLES simulation so the "fully resolved BDIM" is a strange thing to try to calculate... Regardless of this though, the result is promising in that it demonstrates how subgrid-like terms are being more effectively modelled with a smoothed boundary than a resolved boundary with all other things equal.


## Energy spectra and channel Burgulence - 01/07/21:

Validated the -2 slope for periodic Burgulence, showing the DNS/BDIM solver is working properly. The solver is nice and stable even using simple 2nd order central finite differencing (although QUICK can be used as well). Also observed the steeper slope (~-16/3) of the channel/bounded Burgulence spectrum in line with expectations.

![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/1D%20Channel%20Burgulence/corrected_spectrum_with_channel_slope.png)

![Alt Text](https://github.com/J-Leetch/MSc-Project/blob/main/1D%20Channel%20Burgulence/Channel%20Burgulence%20BDIM.gif)

Have achieved a continuously forced problem which has a distinct energy spectrum characteristic of the shock dissipation. Similarly to LES, the energy piles up at a cutoff scale without any subgrid modelling; this is promising since it shows close similarity to turbulence LES modelling. Also, the near-wall features seem to be fairly rich in detail (the turbulence hasn't been killed off at the wall) which is important since this is the main area of investigation.

Now, I am in the process of implementing simple subgrid models (starting with Smagorinsky), show them working - as well as they can - on the periodic case and then implement into the channel. Then, I can examine the main objective which is identifying the optimal modelling required at the wall and how this differs from conventional subgrid models.

Notes:

-   Should I be concerned with the "spatio-temporal" LES approach championed in LaBryer? 
-   As well as being concerned with the energy spectrum in the context of the modelling, must examine how the models preserve statistical structure. Well known that these features are key to turbulence production, transition etc in real flows so this extra layer of carefullness could be important for transferability. 


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
