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

##Quickstart: Data Download and Regridding

   1. Install CDO: `bash make_install_cdo.sh`
   2. Navigate to [NCAR User API token](https://www.earthsystemgrid.org/ac/user/apiTokenDisplay.html) to copy your API token (you will need to make and account/login)
   3. `bash data_downloader.sh pasted_api_token` to download and preprocess the relevant data.
   4. The resultant files are split into the data/CESM1/ and data/NSIDC/ directories with subdirectories related to the data frequency and variables of interest are ready to be ingested by the HAIKU code
   5. `data_downloader.sh` will also spit out text describing urls and rough instructions for any additional observational/reanalysis data that we've used in our analysis


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