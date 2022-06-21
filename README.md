# hmmacc or sleephmm
Hidden Markov Model (HMM) for Sleep/Wake Identification using Actigraphy

## Introduction
This package provides functions on applying HMM to actigraphy/accelerometer data to identify sleep/wake states. The accelerometer data can be either in the summary activity count format, such as activity counts every sixty seconds, or the raw data format. HMM assumes different log (activity count) distribution under sleep and wake states respectively: sleep state has more zeros and low activity counts (zero-inflated truncated Gaussian); wake state has relatively more activity counts (Gaussian). Examine the histogram/density plot of log activity counts in each state to see what distribution assumption is reasonable: i.e. if there are not many zeros but only small numbers in the sleep state, use Gaussian instead.

## Update
A brand-new package has been built to facilitate easy implementation. It can take summary count data directly such as NHANES accelerometer data or raw data that can be pre-processed by GGIR. Examples to be updated soon. 

Reference:
Li X, Zhang Y, Jiang F, Zhao H. A novel machine learning unsupervised algorithm for sleep/wake identification using actigraphy. Chronobiology International. 2020 Apr 30:1-4.

## Circadian analysis
For analysis of periodicities and circadian rhythms as well as interactive data visualization using trelliscope, consider the R package PML.
