# Analysis Toolkit
Our extensible toolkit leverages data from the fast counterfactual simulation of the FKPM (Figure 1). 
We highlight several planned analyses below, but envision the set of analyses growing to meet the needs of our key analytic questions. FKPM provides the computational engine enabling the analyses, and we leverage HKCM as a validation tool. 
These components are still in development and will be expanded upon as the implementation is finalized.

<figure>
<img src="/figs/Analysis_toolkit.png" alt="Analysis Toolkit" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 1:</b> HAIKU includes an extensible analysis toolkit that leverages the FKPM to run millions of What-Ifs in support of explanatory, exploratory, and quantitative analyses.</figcaption>
</figure>

##Semantic Graph Builder
The Semantic Graph Builder will build a causal network of key causal factors impacting forecasts, augmenting the set of model factors with user-defined factors computable from the model (e.g., annual variation in sea ice concentration). 
The FKPM will generate many time series of these factors. We plan to leverage Granger Graphs to apply Granger causality to every pair of factors.  
This captures the causality between factors given the generated time series. 
This summary of the climate model will improve understanding by showing predictions and connections in terms of the variables most meaningful.

## Tipping Point Analysis
The Tipping Analysis identifies regions of input space where significant, qualitative differences occur in model prediction under small changes in inputs. 
HAIKU will exploit the linearity of the Koopman Operator to utilize classical results from dynamical systems and control theory, specifically through the eigendecomposition. The eigenvalues identify unstable modes through the sign of the real part. 
Tipping point analysis uses this eigendecomposition in three separate ways to 
1) characterize initial conditions that excite runaway behavior; 
2) locate changes in model parameters that introduce new unstable modes; or 
3) identify points in time that a control, bounded in magnitude, cannot be designed to counteract an unstable mode.

##Explainer
Explainer augments projection models with traceback information to identify key drivers of quantities of interest in the Semantic Graph. 
For example, if a projection shows sea ice concentration decreasing, Explainer identifies the pathways through the graph that accounts for that change. 
Explainer adds interpretability both to ‘baseline’ model projection as well as What-if runs.

##Value of New Data Estimator (VoNDE)
The Value of New Data Estimator will estimate the value of new data and help identify where to focus resources for data collection. 
To improve accuracy of climate forecasts and tipping point estimates, HAIKU assesses how new data will contribute to new, previously unseen model parameters or reduces the uncertainly of target specific downstream effects. 
Either case indicates the new data are valuable.