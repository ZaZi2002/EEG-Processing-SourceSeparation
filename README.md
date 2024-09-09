# EEG Signal Processing Source Separation

## Overview

This project involves processing and analyzing EEG signals using various methods such as GEVD, DSS, PCA, and ICA. The focus is on signal extraction, noise removal, and comparison of different techniques.

## Tasks

### 1. Signal Estimation and Error Analysis

- **Data**: `mat1.Q` (8-channel EEG signal, 100 seconds, 100 Hz sampling rate)
  - **a)** Estimate and plot the triangular signal (`s1`) with a known period of 4 seconds. Calculate and report the estimation error.
  - **b)** Estimate `s1` when the exact period is unknown but within the range of 3 to 7 seconds. Calculate and report the estimation error.
  - **c)** Estimate `s2` using known on/off times (`1T`). Calculate and report the estimation error.
  - **d)** Estimate `s2` with partial on times provided (`2T`). Calculate and report the estimation error.
  - **e)** Estimate `s3` with a known frequency band of 10-15 Hz. Calculate and report the estimation error.
  - **f)** Estimate `s3` with a frequency band range of 5-25 Hz. Calculate and report the estimation error.

### 2. Noise Removal in EEG Signals

- **Data**: Noisy and clean signals from Exercise 2 of Series 2
  - **a)** Determine spike times from the initial clean signal and store them in a vector.
  - **b)** Use GEVD and DSS to extract spike sources. Reconstruct and denoise the signals (`den_X`).
  - **c)** Plot denoised signals for channels 13 and 24 alongside the original and noisy signals.
  - **d)** Compute and report RRMSE (Relative RMSE) for each method, noise type, and SNR.
  - **e)** Compare the results with PCA and ICA from Series 2.

## Requirements

- MATLAB

