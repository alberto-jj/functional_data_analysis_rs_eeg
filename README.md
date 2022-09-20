# functional_data_analysis_rs_eeg
Codes for spectral features extraction and Functional Data Analysis (FDA) in resting-state EEG recordings

# Feature Extraction
The extraction of both relative powers (spectral_features.ipynb) and dominant frequency (dominant_freq.ipynb) are implemented in Python.

The user can use preprocessed rs-EEG epoched data in BIDS format as input for feature extraction.

Also, topomaps with the sensor-level results can be generated using the (topoplots.ipynb)

# Statistical Analysis (FDA and sensor-level analysis of spectral features)
These two analyses are implemented in R.

The user can find the conventional approach averaging PSD vectors of all epochs in each channel (boxplots_by_center.R file)
