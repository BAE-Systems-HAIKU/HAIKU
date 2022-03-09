# Koopman Assisted Climate Models

The Hybrid AI Koopman-Climate Model (HKCM) aims to augment the ability of current climate models (CICE5 to start) to more accurately model real world measurements by training a Koopman model in series with the climate model to apply missing dynamics at each time step.  But from the programmatic level, we have two other components which should fall under this category as well, the Fast Koopman Proxy Model (FKPM) and our Analysis Toolkit.  These work together leveraging a Koopman proxy model of the full climate system that can be trained and run in minutes instead of days.  This enables our analysis toolkit and allows users to quickly test hypothesis and automatically determine uncertainty on forecasts as well as driving causal factors and potential interventions to avoid tipping points.

Across both the HKCM and the FKPM, we will map from the inputs described in Section 3 into Koopman observables space and back.  We will also specify the mapping between model parameters (or Koopman observables) and variables of interest described in the same section.  By leveraging pregenerated data from the models of interest, we reduce the initial computational cost considerably, while hopefully maintaining sufficient diversity of model parameters to properly train the Koopman operators.  Once we’ve validated this approach, we will revisit running the full climate models with specific parameter values as necessary to explore potential tipping points and causal factors.

The FKPM and mapping to factors of interest are required for us to present an understandable result and allow for expert interaction through the Analysis Toolkit, which consists of a semantic graph of those factors, an explainer that shows pathways that drive specific factors (like rapid sea ice loss), and an identification of which variables will reach a tipping point (along with further statistics of when, how, and why).

More detail into Koopman Operator theory and how the various model components are constructed and serve the climate analyses of interest follow.

##Koopman Operator Theory

<figure>
<img src="/figs/Koopman_diagram.png" alt="Koopman Diagram" style="width:80%">
<figcaption align = "center" style="width:80%"><b>Figure 1:</b> Koopman models the dynamics of a reference system (1). Koopman transforms (2) the state space into an observable space and learns a linear operator (3) in the observable space. Various Koopman publications<sup>1,2,3</sup> contain more details.</figcaption>
</figure>

The Koopman Operator (Figure 1) represents the dynamics of a nonlinear, complex, and uncertain system with an infinite-dimensional space (called “observable space”) that evolves under a linear operator. The lifting function that maps the system into observable space can incorporate prior knowledge, resulting in an abstracted, but understandable representation of the system. The Koopman compressed linear form enables extremely rapid simulations to explore numerous what-if scenarios and generate data for causal discovery of semantic causal factors. Moreover, specific eigenmodes of the system describe its long run behavior and can be used directly for tipping point analysis (see section 5.1.2). Compared to other data-driven approaches (e.g., neural networks), Koopman extrapolates well. Even though it preserves complex non-linear dynamics, it is linear in observable space and is explainable because the operator can be expressed algebraically relating variables of interest, like ice cover. Additionally, Koopman learns on sparse data.

The Koopman evolution equation, **Ψ**(t+1) = K **Ψ**(t), is the closed-form climate dynamics equation we use. In a controlled climate system, we describe the evolution of climate states or observables as **x**(t+1)=F<sub>1</sub> (**x**(t),**u**(t),θ) where t is the time index, **x** is the climate model state vector, **u**, is the vector representing climate forcing, and θ is the vector representing strength of climate interactions. We apply a lifting function to go from the states or observables known to current climate models into a set of collective observables that Koopman can operate: **Ψ**(t+1)=A(θ)**Ψ**(t)+B(θ)**u**(t).

##Hybrid AI Koopman-Climate Model (HKCM)
We plan to train a Koopman operator to learn the dynamics between the predictions from current climate models to actual measured records at each measured time step.  This has the potential to identify physics that may be important to the quantification of tipping points or runaway sea-ice loss (Figure 2).

<figure>
<img src="/figs/HKCM_diagram.png" alt="HKCM Diagram" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 2:</b> Hybrid AI Koopman-Climate Model (HKCM)—The Koopman Operator will predict the dynamics of the system given climate model output for the current step and the estimate of observables from the previous one. HKCM leverages modeled physics in the climate model while accounting for un-modeled physics.</figcaption>
</figure>

By placing the Koopman model in series with the Climate Model, it learns the dynamics of any missing physics and to account for any mismodeled physics. The details of which modes are excited in this Koopman model can be used to trace back to the original climate level variables and help identify which physics properties have been mismodeled.  Our initial goal is to verify that the HKCM can learn the missing dynamics between sea ice volume and ocean temperature and to investigate the relevant modes to potentially determine the missing physics interactions in the CICE5 model.

##Fast Koopman Proxy Model (FKPM)
We will also train a full Koopman model to create a fast proxy model of the full climate simulation.  This will enable a suite of analytics that can extract causality and better characterize tipping points and their associated uncertainty that would take months or years to do with the current best climate models alone. This is shown below (Figure 3).
<figure>
<img src="/figs/FKPM_diagram.png" alt="FKPM Diagram" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 3:</b> Fast Koopman Proxy Model (FKPM)—The FKPM learns the full dynamics of the HKCM or stand-alone climate model, but is able to operate much faster than the either, enabling the analytic toolkit. Analysis of the eigenfunctions and eigenvalues will help identify tipping points and regions of interest for deeper analysis.</figcaption>
</figure>

Initially, we will train a FKPM using the CICE5 model on its own so we can begin developing the Analysis Toolkit, but as the HKCM becomes more capable, we will train an improved FKPM that leverages it as well and compare the analytics for robustness.


<sup>1</sup>Arbabi, H., & Mezic, I. (2017). Ergodic theory, dynamic mode decomposition, and computation of spectral properties of the koopman operator. SIAM Journal on Applied Dynamical Systems, 16(4):2096–2126.
<sup>2</sup>Mezic, I. (2005). Spectral properties of dynamical systems, model reduction and decompositions. Nonlinear Dynamics, 41(1-3), 309-325.
<sup>3</sup>