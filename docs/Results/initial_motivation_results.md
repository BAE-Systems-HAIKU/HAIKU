# Preliminary Results

We are initially investigating the simulation data from the CESM1  Large Ensemble Community Project (LENS) as well as observational satellite data from the National Snow and Ice Data Center (NSIDC) from 1979 to 2020. All results shown here are based on a Fast Koopman Proxy Models (FKPM). 

## Proof of Concept: Identifying Regions of Interest 
Unlike most Machine Learning methodologies, the Koopman Operator Theoretical (KOT) framework learns a physically interpretable model. Applying Koopman Mode Decomposition (KMD) reveals underlying low-dimensional dynamics contained in the Koopman modes and eigenvalues. The Koopman modes associated with important eigenvalues can explain the spatial extent of the low-dimensional dynamics.

Training a FKPM on NSIDC sea ice concentrations (ICEFRAC) enables identification of spatial regions that experience rapid ice loss as seen in Figure 1. A region of interest is highlighted in the Barents and Kara seas as having a large components of rapid exponential sea ice loss.

<figure>
<img src="../../figs/results/NSIDC_1979-2020_exp_mode_both.png" alt="sample NSIDC exponential mode" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 1:</b> Koopman modes associated with exponential decay dynamics in the Northern (left) and Southern (right) hemispheres. This FKPM was trained using NSIDC data (sea ice concentration only) from 1979-2020.  </figcaption>
</figure>
&nbsp;  

## Initial FKPM Results

<figure>
<img src="../../figs/results/NSIDC_vs_CESM_1979-2020_mean-annual_comparison.png" alt="NSIDC FKPM vs CESM1 FKPM" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 2:</b> Koopman modes representing the mean and annual variation in sea ice concentration over five year windows for the Northern hemisphere. (Top) 1979–1983 period and (bottom) 2016–2020 period. Two FKPM were trained on observational data (left) and CESM ensemble data (right).   </figcaption>
</figure>
&nbsp;  
We see that the extent of the mean and annual modes is qualitatively similar between  observational-based and CESM-based FKPMs. Additionally, the regions of high annual variance are highlighted as expected in the annual modes, while the regions of highest sea-ice concentration are highlighted in the mean modes. Next, we move to a more quantitative analysis looking at the accuracy of predictions generated by a FKPM trained from observational data.

<figure>
<img src="../../figs/results/Sea_ice_only_forecasting2.png" alt="FKPM sea ice only monthly accuracy" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 3:</b> A FKPM was trained on monthly NSIDC (sea ice concentration only) data from 1997-2000. The Koopman model was then run forward to forecast monthly predictions from 2001-2004. The pointwise RMSE of the FKPM (blue), CESM1 (orange), and climatological mean (green) were then computed from the observed monthly NSIDC data. </figcaption>
</figure>
&nbsp;  

We see that the FKPM trained only on Sea Ice data seems to mirror the climatological mean, while remaining slightly worse. The CESM1 model is generally not used for the purpose of month-month prediction, but it is worth noting that the relative cycles of error for the CESM1 simulation opposes those of the FKPM and climatological mean. This hints both the NSIDC-based FKPM and the CESM1 may be good candidates for a hybrid model.

The model is unaware of the dynamics of other causal climate variables that directly influence the sea ice concentration dynamics. Expanding the dictionary of observables used in the KOT framework by including additional climate variables will improve the modelling power of the HAIKU system.

## Additional Climate Variables as Observables

We continue the development of the HAIKU system by including additional climate variables in the FKPM training. In particular, sea (potential) temperature, atmospheric temperature, and sea ice thickness were added. These generally improve the performance of the models and we plan to investigate the inclusion of additional variables (eg. solar irradiance and greenhouse gases).

We compare observational-based and CESM-based FKPMs trained on all climate variables. The mean modes are shown in Figure 4 and match expectation for the time average of these variables across the training interval.

<figure>
<img src="../../figs/results/CESM_NSIDC_compare_multiple_variables_mean_mode2.PNG" alt="mean modes" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 4:</b> The mean modes are visualized for the CESM-based FKPM (left) and the observational-based FKPM (right). The results are qualitatively similar and also match physical expectations for the time average of these variables.  </figcaption>
</figure>
&nbsp;  

The annual modes are shown in Figure 5 and match expectation for the typical annual fluctuation of these variables across the training interval.

<figure>
<img src="../../figs/results/CESM_NSIDC_compare_multiple_variables_annual_mode2.PNG" alt="annual modes" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 5:</b> The annual modes are visualized for the CESM-based FKPM (left) and the observational-based FKPM (right). The results are qualitatively similar and also match physical expectations for annual variance of these variables.  </figcaption>
</figure>
&nbsp;  

We have preliminary evidence that shows an improvement in modelling the sea ice concentration dynamics with the inclusion of these additional climate variables. We are interested in assessing the quality of the predictions generated from models with additional variables. The Root Mean Squared Error (RMSE) is computed between the monthly predictions of various models and the NSIDC monthly sea ice concentration data in Figure 6.

<figure>
<img src="../../figs/results/Koopman_accuracy_vs_variabls.png" alt="CESM" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 6:</b> Four FKPM were trained with differing sets of climate variables (from observational data) and the pointwise RMSE of monthly Sea Ice concentration is shown (blue) compared to CESM1 data (orange) and the climatological mean (green). We see an accuracy gain when including ORAS5 Atmospheric Temperature (right), but not when including ERA5 Sea Surface Temperature (bottom) compared to Sea Ice Concentrations alone (top-left).   </figcaption>
</figure>
&nbsp;  

From this result, it is clear that the FKPM captures the dynamics more accurately with the inclusion of Atmospheric temperature. We speculate that the Sea Surface temperature has little positive impact we cannot obtain the sea surface temperatures in regions that contain sea ice. Including CESM simulated sea (potential) temperatures a few feet below the sea ice as an additional climate variable could improve the accuracy of these models. This will be a nice proof of concept for our Phase II approach to identify new measurements that could improve the accuracy of models built from observational data.

## Next Steps

Inclusion of just one additional Climate variable has brought the FKPM monthly accuracy up to the level of the climatological mean. We plan to test the addition of a few new climate variable to improve the modelling of dynamics on a monthly timescale. We will present predictions of the September sea ice over decadal time scales to compare long-term trends with CESM simulations and with the observational data. 

Early summaries of causal models using the updated FKPM will come as well, allowing for further validation of the approach and the relevance of the produced models.

We also plan to produce the first results of the Hybrid Koopman Proxy Model to verify that the approach described in section X can improve accuracy over the climate model itself (eg. CESM1).