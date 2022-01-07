# SpectralAnalysis

A few selected R scripts from my work on spectral analysis of unevenly sampled astronomical time series in the context of the Gaia satellite and other large-scale astronomical surveys (2013-2015). 

The question was how we can find out if a peak in a spectrum of a noisy time series is indeed significant, and we have discovered a periodic signal among the noise. There are several reasons why this is not a trivial issue for astronomical time series. First, astronomical time series are often composed only of about a hundred observations, while the periodogram is computed at hundreds of thousands of frequencies, which means that the periodogram is strongly redundant in its information content, and values at various frequencies are extremely strongly dependent on each other. Second, the times of the observations are a mixture of quasi-regularities and irregularities, which entails the presence of a lot of sidepeaks in the periodogram. 

Significance assessment of a maximum in such a complex object is complex. I proposed extreme-value analysis to solve it, and later extended the method to big, spatially structured data by noticing that the extreme-value parameters of the model appear to depend on the spectral window function (a high-dimensional object) and the number of observations in the time sampling pattern. 

#### Files:

*extrFAP7FirstStar.R*: The codes to compute the False Alarm Probability from a Lomb-Scargle periodogram of a sparse irregular time series.

*betamodelFitting1.R*: Exploration of the structure of the Gaia spectral windows, using mixture-based clustering in a space of strongly reduced dimensions. *cooTransforms.R* contains the coordinate transformations and other auxiliary functions for this. *EMbetaFunctions.R* contains the functions for the EM algorithm to fit the beta mixtures I used.

*gevmodelFitting2.R*: Fitting the generalised extreme-value model to simulated noise sequences with representative Gaia time sampling, and exploring the dependence of the parameters on various features of the spectral window and the time sampling.

*baluev5.R*: Comparison of the extreme-value based methods to other three methods from the literature.

#### Images:

*histLomb_line.pdf*: Perspective view of the histograms of p-values of periodogram maxima of noise sequences with Gaia time samplings from the ecliptic plane to the ecliptic pole.

*map_n_at12.png*: Map of time series lengths (top figure) and spectral window value at 12 c/d (bottom figure) over the whole sky. The red line and patch shows the locations of the representative Gaia time samplings we used in the simulations.

*pardep-gev.png*: The two extreme-value parameters, shown in colour, versus the log length of the time series (x-axis) and the average spectral window peak sizes (y axis).

*maprecovery10sign05-ecl.pdf*: recovery rate of a simulated signal over the small rectangular patch of sky shown in red in image map_n_at12.png, comparing the extreme-value method and the three other methods from the literature. 

#### References:

- M. Süveges. Extreme-value modelling for the significance assessment of periodogram peaks. Monthly Notices of the Royal Astronomical Society, Vol. 440 (2014), 2099–2114.

- M. Süveges, L. P. Guy, L. Eyer, J. Cuypers, B. Holl, I. Lecoeur-Taibi, N. Mowlavi, K. Nienartowicz, D. Ordonez Blanco, L. Rimoldini, I. Ruiz. A comparative study of four significance measures for periodicity detection in astronomical surveys. Monthly Notices of the Royal Astronomical Society, Vol.450 (2015), 2052–2066.



