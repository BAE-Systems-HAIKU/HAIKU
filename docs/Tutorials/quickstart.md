#Quickstart

##Quickstart: set up HAIKU and python env
Thesteps below are explained in full detail in [Configure HAIKU](/Tutorial/configure_haiku)

```
git clone https://github.com/BAE-Systems-HAIKU/HAIKU.git HAIKU
cd HAIKU
python3.8 -m venv ./haiku-venv
source ./haiku-venv/bin/activate
pip install -r requirements.txt
#export PYTHONPATH=$PYTHONPATH:/your/absolute/path/to/HAIKU
```

##Quickstart: Data Download

   1. Install CDO: `bash make_install_cdo.sh`
   2. Navigate to [NCAR User API token](https://www.earthsystemgrid.org/ac/user/apiTokenDisplay.html) to copy your API token (you will need to make and account/login)
   3. `bash data_downloader.sh pasted_api_token` to download and preprocess the relevant data.
   4. The resultant files are split into the data/CESM1/ and data/NSIDC/ directories with subdirectories related to the data frequency and variables of interest are ready to be ingested by the HAIKU code
   5. `data_downloader.sh` will also spit out text describing urls and rough instructions for any additional observational/reanalysis data that we've used in our analysis

##Quickstart: Data Preprocessing and ClimateData objects

The HAIKU software works with ClimateData objects that are defined [here](/software_framework/#climate-data) for both plotting and modeling.
To convert the downloaded climate data into ClimateData objects (including formatting, temporal interpolation, etc) execute the following steps for each downloaded dataset of interest:

   1. `python scripts/preprocessing/preprocess.py input_directory target_grid_file output_directory`
   2. Do step 1 for all variables of interest for your downstream plotting or analysis

## Quickstart: Training a Koopman Model

Once a user has downloaded either the CESM or NSIDC datasets, they can use HAIKU to train a Koopman Model. The following steps outline how to do so:

1. Copy the `configs/example_config.yml` file
2. Update the new configuration file appropriately for your environment
    - Data directories can contain either CESM or NSIDC dataset files
    - Specifying data directories containing different dataset types (e.g., ICEFRAC and SST) will result in a model combining the two variable types
3. Run `python scripts/train.py path_to_configuration_file`
    - System output will be directed to the log file specified in the configuration file
    - This will produce some plots related to the Koopman models that was trained (showing eigenvalues and eigenfunctions)
4. The generated model can then be operated on using the `prediction` and `plotting` modules
    - Examples for using the Koopman models are located in `scripts/plotting`
    - `python scripts/predict.py path_to_configuration_file.yaml trained_koopman_model_file.pkl YYYYMM01 path_to_output_directory`
    - This will produce a set of diagnostic plots related to accuracy of the Koopman model over the training data window extended to the data listed (YYYYMM01: 20201201 for instance)

## Quickstart: Evaluating Robustness of Koopman models to noise (and other parameters)

Once a Koopman model has been trained, as in the previous step, we'd like to evaluate it's ability to predict more quantitatively. One way to do this with the limited validation data is to measure the robustness of its predictions to various perturbations or assumptions we made in the data process. We currently have one such robustness analysis defined (robustness to measurement noise), but others are in the works and described in more detail [here](/metrics/#robustness-of-haiku-models).

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
