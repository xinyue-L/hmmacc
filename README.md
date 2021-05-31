# hmmacc
Hidden Markov Model (HMM) for Sleep/Wake Identification using Actigraphy

## introduction
This package provides functions on applying HMM to actigraphy/accelerometer data to identify sleep/wake states. The accelerometer data should be in the summary activity count format, such as activity counts every thirty seconds. HMM assumes different log (activity count) distribution under sleep and wake states respectively: sleep state has more zeros and low activity counts (zero-inflated truncated Gaussian); wake state has relatively more activity counts (Gaussian). Examine the histogram/density plot of log activity counts in each state to see what distribution assumption is reasonable: i.e. if there are not many zeros but only small numbers in the sleep state, use Gaussian instead.

Reference:
Li X, Zhang Y, Jiang F, Zhao H. A novel machine learning unsupervised algorithm for sleep/wake identification using actigraphy. Chronobiology International. 2020 Apr 30:1-4.

## circadian analysis
For analysis of periodicities and circadian rhythms as well as interactive data visualization using trelliscope, consider the R package PML.
