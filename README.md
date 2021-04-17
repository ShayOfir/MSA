# MSA: Multi-perturbation Shapley Analysis Toolbox
Authors: Shay Ofir-Geva and Isaac Meilijson.

This is Matlab toolbox implements multi-perturbation Shapley-value analysis for lesion-behavior/symptom mapping (LBM).

## Table of Contents
* [General Info](#general-info)
* [Technologies](#technologies)
* [Usage](#usage)
* [Example](#example)


# General Info
Multi-perturbation Shapley-value analysis is a method for multivariate, game-theoretical [lesion-symptom mapping](https://www.nature.com/articles/nrn1521) (finding the neural correlates of certain symptom or behavior using anatomical lesion data of brain-injured patients and their behavioral measurements). A full an detailed descirption of the method could be found in the [here](https://www.mitpressjournals.org/doi/10.1162/0899766041336387) and [here](https://doi.org/10.1162/artl.2006.12.3.333), while an example for usage in lesion-symptom mapping could be found [here](10.1002/hbm.20797), [here](10.1186/s12868-016-0275-6) and [here](10.1002/hbm.23601) . This toolbox was implemented and fully presented in an article under review (a full reference will be included once the article is accepted and published). Shortly, Assuming that a network of brain regions are involved in a certain behavior, each brain region is considered as a player in a coalitional game, and the measured behavioral score when all regions are intact is the game worth. Shapley value is the unique fair division of the game's worth between players, and is usually computed using all the possible player-coalitions.  

The current implementation of MSA uses a uses only coalitions with a certain (user-defined) maximal perturbation depth (i.e., coalitions in which the number of injured / perturbed regions is limited), via a novel formula. This approach is beneficial from both biological and computation point of view. The performance of these coalitions is estimated using a prediction algorithm, related to K-nearest neighbours method, which could be replaced by any predictor the user desires (although we recommend against it to all but most advanced users). The main program recieves a matrix of patients' lesion and behavior data and several parameters (see [Input & Output](#input-output) section for detailes) and outputs the corresponding Shapley values, in addition to its calibrated version (which is used to overcome prediction bias) and inferential statistics (p-values, confidence intervals and others measures related to bootstrap or jacknife techniques).

# Technologies
The toolbox was created and tested using:
* Matlab R2019b (v9.7)

And the following matlab toolboxes:
* Statistics and Machine Learning Toolbox (v11.6)
* Parallel computing toolbox (7.1) 
* MATLAB Parallel Server (7.1)

# Usage

<pre><code>[SV, Calib, coal, Bset, Lset] = PerformMSA (xy, pdepth, nBS, alpha [,alternative_predictor])
</code></pre>

## Input:
* __xy__: a 2D matrix of N rows and M columns, whereas each row is a patients and each column represents a region, excluding the last column which represents the behavioral score. Numbers in region columns are expceted to be in the range [0, 1]. Numbers in the behavior column are not expected to be of certain type. Please avoid NaN, Inf, -Inf etc. as they will cause errors.
* __pdepth__: perturbation depth, i.e. the maximal number of inactive regions in the computed coalition. A rule of thumb, use ~15% of your total number of regions. This number should be an integer between 1 and no. of regions - 1. I used pdepth=4 for 26 regions with relatively short computation time and relatively good accuracy. Please beware: pdepth > 12, espcially if no. of regions is >20, may lead to exhaustion of RAM, and will take very long to compute.
* __nBS__: no. of bootstrap resamples. 0 = No statistical inference. -1 = Jacknife (leave-one-out method; less accurate); any number >0 : bootstrap method. 10,000 provide excellent results, but it seems that 1,000 is almost the same, with much shorter computation time. The number should be integer.
* __alpha__:  Type-I error rate. float in the range [0, 1]. This parameter is ignored if nBS = 0
* __alternative_predictor__: a predictor function which could be used with .predictFcn(XX) method, where XX is the 2D binary matrix of size of N rows (one row for each coalition) and M columns (one column for each region), where 0 means inactive region and 1 an active region (This is an advanced topic and will be addressed in the future. The interested user should refer the code and comments meanwhile).

## Output:
* __SV__ : Shapley Vector. This is actually a matrix with N rows and M columns, when N is pdepth and M is no. of regions. Each row is Shapley vector for a perturbation depth from 1 to pdepth. Hence, SV(end,:) will get you the Shapley-values for the requested pdepth.
* __Calib__ : Calibrated shapley vector: this is a structre, which include the following:
  - __Calib.SV__ - matrix of size (N,M) whereas N is pdepth and M is the number fo regions. The same as SV, but calibrated.
  - __Calib.aver__ - averaged Shapley. 
  - __Calib.mode1__ - the modal Shapley value below the average
  - __Calib.factor1__ - calibration factor
* __Coal__ - this is a structure used for debugging mainly. It contains results of interim computations. Coal.Coal contains explicit lists of all of the coalitions (expressed in bianry working states) for each perturbation depth
* __Bset__ - This cell structure is created when bootstrap option is ON (i.e., nBS > 0). Its length is always equal pdepth. Since some users may want to compare different perturbation depth, it is always a cell-structure in the length of pdepth vector (i.e., if pdepth=1:4, than Bset has four cells, one for each depth). However, this could be very costy in time and memory, so pdepth is a scalar, only Bset{end} will be set. This structure contains numerous fields, some of them are useful only for debugging. The important fields are (all are vectors of length N, where N is the number of regions, unless stated otherwise):
  - __stdestmix__ :estimated standard deviation for each region.
  - __Zscoreestmix__ : Z-score using stdestmix, for difference of each SV(i) from the avergae
  - __pvalestmix__ : p-values using stdestmix, for difference of each SV(i) from the average
  - __pvalestFDRmix__: The same as pvalestmix, but FDR-corrected. You could also use mafdr(pvalestmix,'BHFDR','true') to get the same results.
  - __CImix__: an Nx4 matrix, where N is the no. of regions: 1st coulmn: lower limit of confidence interval; 2nd column - Shapley value (the same as SV); 3rd column - upper limit of confidence interval; 4th column - average shapley (This matrix is useful for quick plotting).
  - __CIcalibmixSHAPL__: The same as CImix, but with calibtrated versions of the CI and Shapley value
  - __Bootstraps__: matrix of N x nBS, where N is no. of regions and nBS is number of boostrap resamples. This one alows visualization and may be more advanced statistics using boostrap distribution
  - __CalibBootstaps__: The same as Bootstraps but calibrated.
 * __Lset__ - This cell structure is created when the Leave-One-Out (Jacknife) option is ON (i.e., nBs == -1). See Bset for explanation regarding the length. This structure contains numerous fields, some of them are useful only for debugging. The important fields are (all are vectors of length N, where N is the number of regions, unless stated otherwise):
  - __stdest__ :estimated standard deviation for each region.
  - __Zscoreest__ : Z-score using stdestmix, for difference of each SV(i) from the avergae
  - __pvalest__ : p-values using stdestmix, for difference of each SV(i) from the average
  - __pvalestFDR__: The same as pvalestmix, but FDR-corrected. You could also use mafdr(pvalestmix,'BHFDR','true') to get the same results.
  - __CI__: an Nx4 matrix, where N is the no. of regions: 1st coulmn: lower limit of confidence interval; 2nd column - Shapley value (the same as SV); 3rd column - upper limit of confidence interval; 4th column - average shapley (This matrix is useful for quick plotting).
  - __CIcalib__: The same as CImix, but with calibtrated versions of the CI and Shapley value
  - __LOO__ - matrix of M x N whereas M is the no. of regions and N is the no. of patients and contains the resultant Shapley vectors for each leave-one-out operation.

## Example
Two real-world [cohorts](https://doi.org/10.1371/journal.pone.0219738) are supplied with the toolbox and could be found in the file `ClinicalCohorts.mat`: `LHD_xy` and `RHD_xy` includes anatomic and behavioral (Fugl-Meyer Score, Upper limb) of 58 left-hemispheric and 49 right-hemispheric first-stroke patients, assesed in the subacute phase (between 1 to 3 months post-stroke). The Fugl-Meyer score of the upper limb range is [0, 66], whereas 66 is normal function. The region names were defined according to the [Automated Anatomic Labeling atlas (AAL)]( http://dx.doi.org/10.1006/nimg.2001.0978) and [White-matter atlas based on DTI](http://dx.doi.org/10.1016/j.neuroimage.2007.12.035). The variable `RegionNames` contains the name of each region, and file `parcellation.docx` describes the exact parcellation.

Please contact many for any question:
Shay Ofir-Geva: shinofir@gmail.com

