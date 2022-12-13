#Running HAIKU

##Set up environment
Follow the steps to [set up your environment](../../Tutorials/quickstart/#quickstart-set-up-haiku-and-python-env) if you haven't already.
If your python environment has been set up already, make sure to activate it:

   `source ./haiku-venv/bin/activate`

This ensures that any python packages as well as the full HAIKU library are at your disposal locally.


##Download Sea Ice datasets

HAIKU has been applied to Arctic Sea ice and while we're in the process of generalizing the system to support a wider range of data, right now only certain datasets are fully supported in the system. These are described in detail [here](../../data_models/#data-citations). Steps required to download the data as automatically as possible are described below:

   1. Install CDO: `bash make_install_cdo.sh`
   2. Navigate to [NCAR User API token](https://www.earthsystemgrid.org/ac/user/apiTokenDisplay.html) to copy your API token (you will need to make and account/login)
   3. `bash data_downloader.sh pasted_api_token` to download and preprocess the relevant data.
   4. The resultant files are split into the data/CESM1/ and data/NSIDC/ directories with subdirectories related to the data frequency and variables of interest are ready to be ingested by the HAIKU code
   5. `data_downloader.sh` will also spit out text describing urls and rough instructions for any additional observational/reanalysis data that we've used in our analysis

For this program, we did not automate some of the more complex data download processes which require codes and multiple web interactions, but these are described in the script as mentioned in step 5 above.

##Alternative Tutorial
In the HAIKU repository, a jupyter notebook (haiku/core/HAIKU_training.ipynb) has some additional documentation and commands to call the scripts, etc below. HAIKU was written to be a command line utility and eventually callable via API, and the jupyter notebook implementation is suboptimal as a result. But it can serve as an alternative look into the system (you will still need to set up the environment and download the datasets of interest).


##Data Preprocessing and ClimateData objects

The HAIKU software works with ClimateData objects that are defined [here](../../software_framework/#climate-data) for both plotting and modeling.
To convert the downloaded climate data into ClimateData objects (including formatting, temporal interpolation, etc) execute the following steps for each downloaded dataset of interest:

   1. `python scripts/preprocessing/preprocess.py input_directory target_grid_file output_directory`
   2. Do step 1 for all variables of interest for your downstream plotting or analysis

The preprocessing step first applies a specific grid to any files that we are considering using CDO. We recommend using the built in grid for the northern hemisphere defined in the polar stereoscopic grid by NSIDC. The NSIDC v4 data are missing some needed descriptors of the grid and the CESM1 datasets are saved in lon/lat coordinates. So, the first step inside preprocessing is to load all data from a given dataset, then to modify remap (using CDO remapbil) the dataset to a consistent NSIDC polar stereoscopic grid aligned with v4 Sea Ice data. Then the data are all aligned temporally, any missing timesteps are linearly interpolated between neighboring timesteps for each spatial gridpoint. The Climate Data object is created and saved to disk in the output_directory. You will need to rerun the above script for NSIDC, CESM1 variables of interest. To begin, we recommend starting with observational NOAA/G02202_V4/, ORAS5 SST, CESM1 ICEFRACE, and CESM1 sst datasets.

optional arguments for the preprocessing script allow for a recursive file structure (you can call it once for all your datasets if you use the root data directory as the input directory and apply the -r flag). The --pattern option is required if you have any non-data files in the directories of interest.

After this step is finished (which may take a while due to the remapping step), you should have 4 pickled datasets stored in output_directory/variablenames/ alongside 2-d snapshots pngs for each timestep in the datasets for validation purposes.

<center>
<figure>
<img src="../../figs/results/2-d_snapshot_example.png" alt="example_climate_data_object" style="width:80%">
<figcaption align = "center" style="width:80%"><b>Figure 1:</b> One of these images will be generated for each of the timesteps after a ClimateData object is created and saved for processing. In this case, the POLAR grid is the default option, but N or S will store the data and plot in latitude and longitudinal coordinates if prefered.</figcaption>
</figure>
</center>

##Creating a Mask File

In order to remain consistent over the CESM1 and NSIDC datasets or over multiple models and datasets, you can generate a mask of regions to exclude when training your koopman model. The default script currently masks out any region that is flagged as land and focusses only on the sea, but could be modified to support other sub-regional analyses.

```
python scripts/create_mask.py data_file.nc output_mask_file.pkl 0 100 --variable cdr_seaice_conc_monthly
```
   - <i>data_file.nc</i> should correspond to one of the nc files which has the same grid shape you used in the data preprocessing step (if using the POLAR grid type, this is any of the NSIDC_v4 files).
   - <i>output_mask_file.pkl</i> is the location to save your mask file (you'll want to update your configuration file to point here for the later steps)
   - <i>variable</i> should be set to cdr_seaice_conc_monthly if you're using NSIDC_v4, if using CESM1 as your target file, you can use ICEFRAC as the variable instead.
   - <i>minimum and maximum</i> describe the range of possible values. If using sea ice (cdr_seaice_conc_monthly or ICEFRAC) stick with 0, 100. Since the aim here is to automatically mask the land, one of these two options should be sufficient.

You can check that this code worked correctly by reviewing the image of the mask that is saved alongside the .pkl file (example shown in Figure 2).

<figure>
<img src="../../figs/results/example_mask.png" alt="example_mask_object" style="width:80%">
<figcaption align = "center" style="width:80%"><b>Figure 2:</b> Automatically generated image of the mask. You should verify that it has the same coordinate system as the data produced in the previous step and that the red region on the right corresponds to the region you'd like unmasked.</figcaption>
</figure>
</center>

## Training a Koopman Model

Once a user has downloaded either the CESM or NSIDC datasets, they can use HAIKU to train a Koopman Model. The following steps outline how to do so:

1. Copy the `configs/example_config.yml` file
2. Update the new configuration file appropriately for your environment. Descriptions of each configuration variable are described in place in the example configuration file. 
    - Data directories can contain either CESM or NSIDC dataset files
    - Specifying data directories containing different dataset types (e.g., ICEFRAC and SST) will result in a model combining the two variable types
    - The other variable most likely of interest is to change the start or end time for the training. Just make sure that the training data sources you supply all span the training data window.
    - Make sure you point to the mask file you created in the previous step so that you properly mask the non-sea regions from both koopman modeling and plotting.
3. Run `python scripts/training/train.py path_to_configuration_file`
    - System output will be directed to the log file specified in the configuration file
    - This will produce some plots (in output_directory defined in configuration file) related to the Koopman models that was trained (showing eigenvalues and eigenfunctions)
4. The generated model can then be operated on using the `prediction` and `plotting` modules
    - Methods for plotting the Koopman models and Climate Data Obects are located in `haiku/plotting/plotFunctions.py`
    - `python scripts/prediction/predict.py path_to_configuration_file.yaml trained_koopman_model_file.pkl YYYYMM01 path_to_output_directory`
    - This will produce a set of diagnostic plots related to accuracy of the Koopman model over the training data window extended to the data listed (YYYYMM01: 20201201 for instance). Examples of how to call many of the more relevant plotting functions are found near the end of the `predict.py` script.

## Evaluating a Koopman Model

By default, some diagnostic plots will be produced during the training step alongside the model pickle file. Default plots with brief descriptions are shown below.

###Model characteristics Diagnostic Plots

<figure>
<img src="../../figs/results/eigenvaluesCT_merged_sst.png" alt="auto_eigenvalues" style="width:60%">
<figcaption align = "center" style="width:60%"><b>Figure 3:</b> Eigenvalues plot. This plot is automatically generated during the Koopman model training step. The distribution of eigenvalues and their relative magnitudes can be seen. 
</figure>
&nbsp;

Figure 1 shows the distribution of eigenvalues and their relative magnitudes from the Koopman model. These plots are automatically generated and can be used to verify that the model fit to the data falls within expected parameters. For a stable system, we expect the largest modes to fall along the x=0 axis. We can see that most of the variance falls along the annual frequency and its associate harmonics. The modes which fall on the y=0 axis and have a non-zero x value are exponential modes which can indicate tipping points in the system.

Figure 4 highlights the largest mode, the mean mode and shows the spatial distribution of the sea ice averages in the arctic region as learned by the Koopman model.

<figure>
<img src="../../figs/results/mode1_cdr_seaice.png" alt="auto_mean_mode" style="width:60%">
<figcaption align = "left" style="width:60%"><b>Figure 4:</b> Spatial distribution of Mean mode. This plot is automatically generated during the Koopman model training step. This shows the spatial distribution of the Mean mode for the NOAA sea ice concentration (a separate plot will be produced for each mode and each variable included in Koopman model training).
</figure>
&nbsp;

It should be noted that each individual mode, especially in the case of the annual harmonics, may not have independent physical meaning, there may appear to be some spatial waves present in one annual harmonic mode that when summed with additional annual harmonic modes generates a much smoother annual harmonic picture.

Plots of the most interest to this program are the highest magnitude exponential modes (spatial distribution in Figure 5). By identifying regions of interest, we can cross-check predictions and potentially identify spatial regions of interest.

<figure>
<img src="../../figs/results/exponential_mode_cdr_seaice_conc_monthly.png" alt="auto_exponential_mode" style="width:45%">
<img src="../../figs/results/exponential_mode_sosst.png" alt="auto_exponential_sst_mode" style="width:45%">
<figcaption align = "center" style="width:90%"><b>Figure 5:</b> Spatial distribution of largest exponential mode. A large exponential mode is the sign of a tipping point in the variable of interest. <b>Left:</b> a region of exponential decay in the Sea Ice Coverage variable in the Barents Sea region. <b>Right:</b> some correlation here in the Sea Surface Temperature and the Sea Ice Coverage. In particular, there is some rapid warming predicted by the Koopman model in the same region as the potential tipping point region. 

</figure>
&nbsp;

By comparing the same eigenmodes across different observable variables, we may be able to help infer more about potential tipping points. In the case of the potential Sea Ice Coverage tipping point in the Barents Sea region shown if Figure 5, we can look at the same exponential mode in the Sea Surface Temperature variable. There seems to be general warming in the Arctic on the 20 year timescale and in particular, we can see the region of the Barents Sea is the edge of that drastic warming region and there is likely some causal relationship involved in this interaction.  A more detailed analysis of this sort is planned and described generally in the [Causal Analysis section](../../software_framework/#analytics-toolkit)

###Model Prediction Diagnostic Plots

When the prediction code is run, several default plots will be produced. The first (Figure 6) is a set of 2-d plots for sea ice coverage alonside the measured NSIDC sea ice coverage on the same timestep. By default, these are plotted at each timestep to allow for careful analysis.

<figure>
<img src="../../figs/results/sea_ice_coverage_prediction_shapshot.png" alt="prediction_snapshot" style="width:95%">
<figcaption align = "center" style="width:95%"><b>Figure 6:</b> Spatial distribution of Sea Ice Extent as predicted by a koopman model. One of these plots are automatically generated during the prediction step for each timestep. This shows the spatial distribution of the predicted sea ice coverage at each time step to compare with the NSIDC observation (specifically predicting 11 years into the future for September sea ice estimates).
</figure>
&nbsp;

It is hard to quickly understand how the Koopman model prediction is doing when looking at individual time steps, especially at the decadal time-scale of interest in this program. To cover that portion, the prediction script also generates a series of time series comparisons over composed variables of interest. In this case, the most interesting variable is likely sea ice extent (or the number of square kilometers that are coverage by sea ice in the arctic) seen in Figure 7.

<figure>
<img src="../../figs/results/Month_3_Sea_Ice_Extent.png" alt="Sea Ice Extent March" style="width:48%">
<img src="../../figs/results/Month_9_Sea_Ice_Extent.png" alt="Sea Ice Extent March" style="width:48%">
<figcaption align = "center" style="width:96%"><b>Figure 7:</b> Sea Ice Extent as predicted by the Koopman model as compared to NSIDC observation, Climatalogical Mean, and CESM1 Large Ensemble member 002. This plot is automatically generated during the prediction step. The annual data for the March monthly average (<b>Left</b>) and September monthly average (<b>Right</b>) are shown for comparison.
</figure>
&nbsp;

The Koopman model will generally align perfectly with the training data and diverge somewhat in the prediction window. Additional plots are automatically generated for each month and for spatial correlation as well as RMSE comparing the prediction with the NSIDC data.


## Evaluating Robustness of Koopman models to noise (and other parameters)

Once a Koopman model has been trained, as in the previous step, we'd like to evaluate it's ability to predict more quantitatively. One way to do this with the limited validation data is to measure the robustness of its predictions to various perturbations or assumptions we made in the data process. We currently have one such robustness analysis defined (robustness to measurement noise), but others are in the works and described in more detail [here](../../metrics/#robustness-of-haiku-models).

1. Copy the `configs/example_config.yml` file
2. Update the new configuration file appropriately for your environment
    - Data directories can contain either CESM or NSIDC dataset files
    - Specifying data directories containing different dataset types (e.g., ICEFRAC and SST) will result in a model combining the two variable types
3. Run `python scripts/robustness.py path_to_configuration_file`
    - System output will be directed to the log file specified in the configuration file
    - This will produce n ClimateData objects (saved to disk) with the magnitude of random noise defined in your configuration file applied.
    - This will also produce n Koopman models for each ClimateData object. We can then analysis the distribution of Koopman models and their predictions to understand how robust HIAKU is to the input noise for this dataset. Sub-folders for each Koopman model parameter can be found
4. The generated models can then be operated on using the `robustness_plotting` module
    - `python scripts/robustness_plots.py path_to_configuration_file.yaml path_to_robustness_models/ path_to_default_koopman_model.pkl YYYYMM01 output_directory_for_plots/`
    - This will produce a set of diagnostic plots related to the robustness of the Koopman model over the training data window extended to the data listed (YYYYMM01: 20201201 for instance)


The production of each of the n datasets and koopman models will (by default) have the full set of analytics plots should you wish to analyze any of them more closely, but the main goal is to use the distribution of models and their predictions to understand how the Koopman model predictions are susceptible to the specific selection analyzed. In the default case, this is a 5% uncertainty described in the NOAA NSDIC measurements. When running the robustness_plots script, plots with uncertainty bands around the predictions will be generated to help determine how consistent the predictions of the Koopman model are. In particular, the focus here is on determining if the decadal timescale trends are consistent across the distribution of perturbed datasets.

Each dataset will have its own number from 0-N and the full plots associated with prediction/training will be found in associated subdirectories. Inside the Original subdirectory are the same diagnostic plots for the koopman model trained on the unperturbed dataset. And finally, inside the 'average' directory, one can find the original model with uncertainty bands based on the full distribution of predictions. Figure 8 shows the sea ice extent predicted for March and September.

<figure>
<img src="../../figs/results/Month_3_Sea_Ice_Extent_robustness.png" alt="Robustness Sea Ice Extent March" style="width:48%">
<img src="../../figs/results/Month_9_Sea_Ice_Extent_robustness.png" alt="Robustness Sea Ice Extent March" style="width:48%">
<figcaption align = "center" style="width:96%"><b>Figure 8:</b> Sea Ice Extent as predicted by the Koopman model as compared to NSIDC observation, Climatalogical Mean, and CESM1 Large Ensemble member 002 with the 2sigma uncertainty bands on the model prediction due to measurement uncertainty. This plot is automatically generated during the prediction step. The annual data for the March monthly average (<b>Left</b>) and September monthly average (<b>Right</b>) are shown for comparison.
</figure>
&nbsp;

We typically expect to see that the Koopman model trends are consistent although the smearing will provide potentially large error bands due to the chaotic nature of the climate system. If we see that the Koopman models are predicting different trends on decadal timescales, then we must improve the input data or verify that the value of the varied parameter chosen for the model generation was better motivated than the alternates used in the first step of the robustness analsysis. 