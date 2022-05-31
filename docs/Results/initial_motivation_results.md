#Preliminary results from Simple Fast Koopman Proxy Model
We are initially investigating the simulation data from the CESM1 Large Ensemble Community Project (LENS) as well as NSIDC data that covers the years 1979-2020. 
When evaluating accuracy of climate forecasting, training and testing occur in separate 4 year windows: 1997-2000,2001-2004 respectively.

All results shown here are of a FKPM that only models the dynamics of the Sea Ice levels provided as output from the full CESM1 simulation. 
This means that we cannot extract causality or execute any what-if analyses outside of self-coupling of sea-ice concentration to itself, but those external forcing functions are represented intrinsically in the output of the CESM1 simulation output.

##Proof of concept for IDing Modes of interest
The lifting function and eigenfunctions of the Koopman Mode Decomposition currently used in the system is straigtforward to visualize. The FKPM learns spatial weights for each element through fitting the dynamical system. By identifying eigenvalues of interest, exponential modes with decadal time-scales and large coefficients, we can find spatial regions that may experience rapid ice loss.

Training a FKPM on NSIDC Sea Ice Concentration (ICEFRAC) data, we are able to identify spatial regions that experience rapid ice loss as seen in Figure 1. A region of interest is highlighted in the Barents and Kara seas as having a large components of rapid exponential sea ice loss.

<figure>
<img src="../../figs/results/NSIDC_1979-2020_exp_mode_both.png" alt="sample NSIDC exponential mode" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 1:</b> Example eigenfunctions for exponentially decaying Koopman modes in North (left) and left (right) hemispheres. This FKPM was trained using NSIDC data (sea ice concentration only) from 1979-2020.  </figcaption>
</figure>  <br/><br/>

The FKPM generated only through Sea Ice Concentration time-series data are insufficient to capture all the relevant climate dynamics, so the above analysis should be considered a proof of concept.   

##Accuracy of FKPMs
<figure>
<img src="../../figs/results/NSIDC_vs_CESM_1979-2020_mean-annual_comparison.png" alt="NSIDC FKPM vs CESM1 FKPM" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 2:</b> The Koopman Mode Decomposition allows us to plot each separately and remove Mean and Annual modes for later analyses at any point during analysis. Two FKPM were trained on NSIDC (left) and CESM1 ensemble (right) data from 1979-2020. The two most prominent NSIDC Mean modes (left-most) can be compared to the two most prominent CESM Mean modes(middle-right) while the two most prominent NSIDC annual modes (middle-left) can be compared to the two most prominent CESM annual modes (right-most).  </figcaption>
</figure>  <br/><br/>

We see that the extent of the eigenfunctions for Mean and Annual modes are qualitatively similar between a CESM-based and NSIDC-based FKPMs. Additionally, the regions of high annual variance are highlighted as expected in the annual mode eigenfunctions, while the regions of highest sea-ice concentration are highlighted in the mean mode eigenfunctions.  This all matches expectation and shows that qualitatively, the FKPM are able to capture the more well-behaved dynamics of the system.

Next we move to a more quantitative analysis looking at the accuracy of FKPM trained on NSIDC data.

<figure>
<img src="../../figs/results/Sea_ice_only_forecasting.png" alt="FKPM sea ice only monthly accuracy" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 3:</b> A FKPM was trained on monthly NSIDC (sea ice concentration only) data from 1997-2000. The Koopman model was then run foreward to forecast monthly predictions from 2001-2004. The pointwise RMSE of the FKPM (blue), CESM1 (orange), and climatalogical mean (green) were then computed from the observed monthly NSIDC data. </figcaption>
</figure>  <br/><br/>

We see that the FKPM trained only on Sea Ice data seems to mirror the climatological mean, while remaining slightly worse. The CESM1 model is generally not used for the purpose of month-month prediction, but it is worth noting that the relative cycles of error for the CESM1 simulation opposes those of the FKPM and climatalogical mean. This hints both the NSIDC based FKPM and the CESM1 may be good candidates for a hybrid model.

It is unsurprising that the Sea Ice Concentration only FKPM would provide similar performance to the climatological mean. In essence, it forecasts based solely on the recent Sea Ice concentration and is unaware of other casual climate variables whose variance will impact Sea Ice levels directly.  In order to improve the accuracy of the HAIKU system, we must include additional climate variables into the FKPM.

##Additional climate variables improve model accuracy

We continue the development of the HAIKU system by including additional climate variables in the FKPM training. In particular, we include: Sea Ice Concentration, Sea Surface Temperature, and Atmoshperic Surface Temperature.  These generally improve the performance of the models and we plan to investigate the inclusion of additional variables (eg. 3D gridded atmospheric temperature, atmospheric C02).

First, we train an FKPM with all three of these inputs for CESM1 data and observed data and compare to make sure the FKPM is providing qualitatively reasonable output. The mean modes are plotted in Figure 4 and match expectation for the typical average values of these variables across the training timespan.

<figure>
<img src="../../figs/results/CESM_NSIDC_compare_multiple_variables_mean_mode.png" alt="CESM" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 4:</b> The Mean eigenfunctions are summed and visualized for the CESM1-based FKPM (top) and the observational-based FKPM (bottom). The results are qualitatively similar and also match physical expectations for annual variance of these variables.  </figcaption>
</figure>  <br/><br/>

The annual modes are shown in Figure 5 and match expectation for the typical annual fluctuation of these variables across the training timespan.

<figure>
<img src="../../figs/results/CESM_NSIDC_compare_multiple_variables_annual_mode.png" alt="CESM" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 5:</b> The annual eigenfunctions are summed and visualized for the CESM1-based FKPM (top) and the observational-based FKPM (bottom). The results are qualitatively similar and also match physical expectations for annual variance of these variables.  </figcaption>
</figure>  <br/><br/>

We have preliminary evidence that the FKPM is accurately modeling the new climate variables and now want to quantitatively assess how much additional prediction power or modeling accuracy we gain by the inclusion of these variables. The RMSE is computed between the monthly predictions of various models and the NSIDC monthly sea ice concentration data in Figure 6. 

<figure>
<img src="../../figs/results/Koopman_accuracy_vs_variabls.png" alt="CESM" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 6:</b> Four FKPM were trained with differing sets of input variables (from observational data) and the pointwise RMSE of monthly Sea Ice concentration is shown (blue) compared to CESM1 data (orange) and the climatalogical mean (green). We see an accuracy gain when including ORAS5 Atmoshperic Temperature (right), but not when including ERA5 Sea Surface Temperature (bottom) compared to Sea Ice alone (top-left).   </figcaption>
</figure>  <br/><br/>

From this result, it is clear that the the FKPM becomes more accurate with the inclusion of Atmospheric temperature. We speculate that the Sea Surface temperature has little positive impact because the ERA5 Sea Surface Temperature data is known to be fairly spotty.  Further analysis will leverage CESM1-based analysis to sea if more accurate Sea Surface data could improve the accuracy of these models.  This will be a nice proof of concept for our Phase II approach to identify new measurements that could improve the accuracy of models built from observational data.

## Next Steps

Inclusion of just one additional variable has brought the FKPM monthly accuracy up to the level of the climatalogical mean. We plan to test the addition of a few new variables to get generate a more accurate model on a monthly timescale.  Additionally, we will present accuracy measurements on longer timescales to more directly compare the climatalogical accuracy and compare directly to the CESM model for identifying long-term Sea Ice Trends.  We expect that some tuning of the Koopman model itself will be involved (such as modifying the lifting functions or adding physical restrictions to the eigenfunctions and eigenvalues representing real physics interaction limitations.)

Early summaries of causal models using the updated FKPM will come as well, allowing for further validation of the approach and the relevance of the produced models.

We also plan to produce the first results of the Hybrid Koopman Proxy Model to verify that the approach described in section X can improve accuracy over the climate model itself (eg. CESM1).