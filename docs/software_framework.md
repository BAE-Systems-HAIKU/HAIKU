# Software Framework

This is a summary of the structure of the python-based HAIKU software system currently in development.
As the system is developed, a beta version will be available on our public [github page](https://github.com/BAE-Systems-HAIKU/HAIKU).
This page will include steps required to download and preprocess training data,
and to get HAIKU operational on your system.

 \* **_classes_** and _member functions_ are denoted as such section.

 <font size="5">**Overall structure** </font>

<figure>
<img src="../figs/diagrams/haiku-core-classes-diagram.png" alt="Software capability summary" style="width:95%">
<figcaption align = "center" style="width:80%"><b>Figure 1:</b> Classes representing core functionality in the HAIKU software system.</figcaption>
</figure>
&nbsp;



We generate a **_ClimateData_** object that encapsulates CESM, NSIDC, or Koopman generated data and converts them all to consistent representations (the same coordinate grid and sets of variables) while storing the provenance as a member.



<figure>
<img src="../figs/diagrams/software_capability_summary.png" alt="Software capability summary" style="width:95%">
<figcaption align = "center" style="width:80%"><b>Figure 2:</b> Central objects to the HAIKU system.  The analytics causal analysis and analysis toolkit are still being built out and are not present in the current release.</figcaption>
</figure>
&nbsp;

The **_KoopmanModel_** class contains a trained Koopman model as well as any training parameters associated with the stored data. During training, it takes a single **_ClimateData_** object and learns the dynamics of this system.

The **_Predictor_** object is used to operate on **_KoopmanModel_** objects and generate **_ClimateData_** with a range of parameters.

The **_Plotter_** object takes either a **_KoopmanModel_** or **_ClimateData_** as input, along with runtime parameters, and generates visualizations representing the stored data. As such, this object is also used to visualize predictions.

<figure>
<img src="../figs/diagrams/model_generation_flowchart.png" alt="Model generation software flowchart" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 3:</b> The HAIKU framework ingests data and generates a series of models to enable Tipping Point and other analytics on the climate system.</figcaption>
</figure>
&nbsp;

Finally, the **_KoopmanModel_** itself or the time-series data contained in a **_ClimateData_** object can be passed into the Analytics Toolkit. A **_CausalModel_** object is instantiated and can learn a causal structure from time-series data (**_ClimateData_** object) using its internal methods or from the structure of the **_KoopmanModel_** itself. Similarly, the **_CausalModel_**, **_ClimateData_**, or **_KoopmanModel_** objects are used as input to different tipping-point analyses inside the toolkit.

<figure>
<img src="../figs/diagrams/toolkit_flowchart.png" alt="Analytics Toolkit software flowchart" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 4:</b> Leveraging the generated models and time-series data, several analyses are enabled in the Analytics Toolkit.</figcaption>
</figure>
&nbsp;

Rounding out the system, there are a variety of metrics that are evaluated either as member functions of the systems or as standalone code.



##Climate Data
**_ClimateData_** objects are instantiated by the **_ClimateDataLoader_** class which reads in CESM or NSIDC time-series data and converts it to consistent representations (the same coordinate grid and sets of variables). Similarly, **_ClimateData_** can also be produced from a **_KoopmanModel_** object by running the model through the **_Predictor_** object.
The time-series data is stored in a numpy array and is by default monthly climate variable data.
Internal processing converts between polar and lat-lon coordinates, interpolates missing datapoints, and produces time-series matching the lifted Koopman observables (given a **_KoopmanModel_**).

**_Plotter_** operates on **_ClimateData_** to visually investigate the temporal and spatial behavior of the data.


##Koopman Models
The **_KoopmanModel_** class contains a trained Koopman model. The **_KoopmanModelTrainer_** class contains functions necessary for training a Koopman model based on provided **_ClimateData_**. It is used in conjunction with the **_Predictor_** class to generate prediction time-series data in the form of **_ClimateData_**.

Several model hyperparameters are set at instantiation through the configuration file.

The **_Predictor_** has a member function, **_Predictor_**.predict(**_KoopmanModel**,dt), which returns the predicted **_climate_state_** after the **_KoopmanModel_** has run the original state, x, forward by time dt. This function lifts the original climate state into the Koopman Observables space before propagating the state forward using matrix multiplication, reversing the lifting function, and producing the predicted state in the original **_climate_state_** format.  There is an associated function for bulk processing of the **_KoopmanModel_**._predict_state()_ function which can provide a full **_ClimateData_** object as output.  This is more commonly used in most analytics.  Currently, the lifting function is a relatively straightforward aggregation of the **_ClimateData_**, but we are investigating other approaches as the development continues.

The **_KoopmanModel_** also has external plotting functions to summarize the model structure including plots of eigenfunctions of selected modes and the distribution of eigenvalues for the **_KoopmanModel_**.

The forecasting done by the Koopman Models enables the Analytics Toolkit or can produce stand-alone climate forecasts for public consumption.

##Hybrid Modeling
We're still designing the structure of the Hybrid Koopman-Climate Model (HKPM) implementation.
For the scope of this project we intend to apply a correction on top of pregenerated data from CESM or another climate model rather than running the full CESM climate model locally and applying the correction in place.
This will likely be sufficient to test the HKCM as a proof-of-concept.

The HKPM itself is the correction to apply at each time-step of a climate model.
Input to this system are two **_ClimateData_** time-series with the same variables and over the same time-period.
A **_KoopmanModel_** object is trained on each of the **_ClimateData_** objects constraining them to have the same eigenvalues so that they can be compared directly to one another.  The final result is a **_KoopmanModel_** which is the difference of these two.

This **_KoopmanModel_** can then be used directly to provide a correction factor to **_ClimateData_** used as input through the **_KoopmanModel_**._predict_state()_ function.  Alternatively, it enables analytics (currently done manually) to better understand the causal differences between the two models. It is possible to generate a **_CausalModel_** object from this data, which may further enable understanding of the physical difference between the original datasets, but further study is required.



###Causal Model

This class requires **_ClimateData_** as well as its own member **_parameters_** (which helps define variable transformation from more fine-grained to user oriented causal variables).
The **_CausalModel_**._transform_data(**_ClimateData_**) function generates a user oriented **_ClimateData_** time-series with many fewer variables.
This time-series can then be used as input to train the **_CausalModel_** where it uses pairwise Granger Causality coupled with LASSO to limit number of edges, remove edges explained by other pathways.

This **_CausalModel_** can then be viewed via a **_CausalModel_**._print()_ method.

Other analytics are still being considered that may do things like allow a user to request the variable or pathway with the greatest impact on another variable.


##Analytics Toolkit

The initial implementation of the **_analysis_toolkit_** hinges around the **_CausalModel_** class. The **_analysis_toolkit_** class is envisioned to hold many various analysis method, but not to hold any objects itself.

<figure>
<img src="../figs/diagrams/toolkit_flowchart.png" alt="Analytics Toolkit software flowchart" style="width:90%">
<figcaption align = "center" style="width:90%"><b>Figure 5:</b> Leveraging the generated models and time-series data, several analyses are enabled in the Analytics Toolkit.</figcaption>
</figure>
&nbsp;

The analytics toolkit allows for generation of the causal model class and enables the time-series based what-if analyses that can be conducted on the causal variables themselves.

The analytics toolkit will also enable the metrics described in the [Metrics](../metrics) section.  Specifically, a user will be able to select the metric of interest (eg. robustness to training window bounds) and will automatically run the relevant analysis over the specified parameters and present figures and estimates of the metric in question.

