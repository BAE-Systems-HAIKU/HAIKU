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

To visualize the decadal trends of the data and our model, we looked at the September sea ice area over the last 40 years. Figure 3 shows that while the climatological mean is only able to capture the annual variation, our model is able to learn the long-term trend of the sea ice in addition to the annual variation. 

<figure>
<img src="../../figs/results/Koopman_model_september_sea_ice_beta.png" alt="FKPM sea ice only monthly accuracy" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 3:</b> September sea ice area from 1979 to 2020. (Blue) Predictions from a Koopman model trained from 01/1979 to 12/2009. (Orange) Mean of CESM1 ensemble members’ predictions. (Green) Climatological mean of NSIDC observations. Each calendar month in the climatological mean is the mean of the data values of that calendar month over the entire training interval.  </figcaption>
</figure>
&nbsp;  


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

## Inclusion of forcing terms in the model
Global Climate Models (GCMs) such as the Community Earth System Model (CESM) are driven by forcings such as greenhouse gases and other anthropogenic factors. 
We can apply the same forcings to our model to understand the impact of different scenarios on our climate systems. 
Figure 7 shows our model’s predictions under two different forcing scenarios: one with constant forcing after 2009, and one with historical forcing. 
Importantly, we find that the sea ice dynamics are mostly internally driven with little direct impact from the CO2 forcing.

<figure>
<img src="../../figs/results/FKPM_with_forcing_beta.png" alt="CESM" style="width:80%">
<figcaption align = "center" style="width:90%"><b>Figure 7:</b> – September sea ice area from 1979 to 2027. (Blue) Predictions from a Koopman model trained from 01/1979 to 12/2009. (Orange) Mean of CESM1 ensemble members’ predictions. (Green) CO2 volume mixing ratio. Two scenarios are shown: constant forcing after 2009 (top) and complete historical forcing (bottom). The red star represents a forecast of the sea ice in September 2022 based on the current sea ice: https://www.arcus.org/sipn/sea-ice-outlook/2022/july</figure>
&nbsp;  

## Proof of concept robustness metric
The speed at which the koopman models can be trained allows for a much more extensive analysis of the model robustness to parameter choices and data uncertainty than with a traditional climate model. Climate is generally a chaotic system, but if the same tipping points or decadal trends are present in a cohort of 50 koopman models are trained with slightly varied input data, then we can increase our confidence that those trends or tipping points are real or should be looked into in more depth.

As an initial proof of concept, we train a koopman model on the NSIDC dataset over the range of 1978-2009 and a prediction/validation range from 2009-2021. We additionally generate 50 datasets by assuming gaussian noise of 5% on the underlying measurements as described in the [NOAA paper](https://nsidc.org/sites/default/files/cdrp-atbd-final_7.pdf). Each of these datasets is quickly modeled by using Koopman methods to produce 50 models. We then generate a distribution of predictions using all 50 koopman models to understand how much of an impact the specific values of the measurements in the NSIDC data had as compared to a possible alternative set assuming the 5% relative variance in observation described by NOAA.

The two questions of interest are:
   - do we see the flattening of September sea ice extent across all model predictions?
   - do we see the exponential mode in the Barents Sea across all models?

<figure>
<img src="../../figs/results/Month_3_Sea_Ice_Extent_robustness.png" alt="Robustness Sea Ice Extent March" style="width:48%">
<img src="../../figs/results/Month_9_Sea_Ice_Extent_robustness.png" alt="Robustness Sea Ice Extent March" style="width:48%">
<figcaption align = "center" style="width:96%"><b>Figure 8:</b> Sea Ice Extent as predicted by the Koopman mode
l as compared to NSIDC observation, Climatalogical Mean, and CESM1 Large Ensemble member 002 with the 2sigma uncertainty bands on the model prediction due to measurement uncertainty. This plot is automatically generated during the prediction step. The annual data for the March monthly average (<b>Left</b>) and September monthly average (<b>Right</b>) are shown for comparison.
</figure>
&nbsp;

We plan to extend this analysis to other forms of uncertainty or parameter choice as described in the [robustness section](../../metrics/#robustness-of-haiku-models). And it can be used to estimate the [value of new data](../../analyses/#value-of-new-data-estimator-vonde), although the infrastructure for that is not yet built out.

## Sub-region analysis to improve stability
We can improve the overall robustness with respect to spatially distributed noise by leveraging sub-region analysis. When we investigated the Barents Sea subregion, the Koopman methods reveal an exponentially decaying Koopman mode as seen in Figure 10. The data from the entire hemisphere is subject to more noise than a small subregion, and can cause the model to pick up extraneous dynamics. We enhance the predictions from our global model by updating with information from local models, such as the model trained on the Barents Sea subregion. By feeding the known dynamics back into the global FKPM, we can average nearby Koopman modes to produce global support for the locally learned dynamics.

<figure>
<img src="../../figs/results/region_analysis.jpg" alt="regional_analysis" style="width:80%">
<figcaption align = "center" style="width:90%"><b>Figure 9:</b> <b>Left:</b> Koopman methods are applied to data from the Barents Sea subregion. <b>Right:</b> Koopman methods are applied to the entire northern hemisphere. The red circle on the complex plane is centeredat the exponential decay dynamics learned for the Barents Sea subregion and produces a similar support.
</figure>
&nbsp;  

## Initial HKCM results
In the proposed HKCM architecture, the CESM1 is augmented with a Koopman model to improve predictive ability. This approach is described in more detail in the [Hybrid Koopman-Climate Model](../../koopman/#hybrid-ai-koopman-climate-model-hkcm) section.

The coupled Koopman model was trained on the error between the CESM1 simulation and NSIDC observational data. The results from this model are shown in Figure 10, specifically the error that arises in the mean and annual variations.

<figure>
<img src="../../figs/results/hybrid_difference_modes.jpg" alt="hybrid_differences" style="width:95%">
<figcaption align = "center" style="width:95%"><b>Figure 10:</b> Koopman modes as part of HKCM. The (left) mean and (right) annual variation Koopman modes trained on the error between CESM1 simulations and NSIDC observational data.
</figure>
&nbsp;

By looking at the regions of interest and the dynamics of the error across multiple variables and with forcing terms in place, we can start to track down the sources of the descrepancy and possibly improve the underlying model or at least understand why it makes errors to better contextualize results.

## Next Steps

We intend to extend and publish our preliminary results on Koopman Climate modeling and robustness analysis of climate data system. And hopefully will be able to continue developing this software and analysis.