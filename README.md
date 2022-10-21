# Inferring-Heterogeneous-Causal-Effects-in-Presence-of-Spatial-Confounding
Given a scalar outcome of interest $y$ and a scaler exposure variable $z$, which have been collected over space $s$ so that the observed data in D={y_i,z_i,s_i}_{i=1}^{n}, the code allows to estimate the causal effect $\tau$ of $z$ on $y$ while mitigating for spatial confounding caused by other unobserved spatially varying variable $c$. The effect $\tau$ can itself be spatially varying i.e. heterogeneous. As a motivating example, suppose that the outcome $y$ is crop yield and the exposure $z$ is amount of fertilizer and we are interested in studying whether using more fertilizer increases the crop yield. Here the confounding variable $c$ is the unobserved inherent variability of the soil. Moreover, the effect of fertilizer on crop yield may vary spatially. For details see and cite: Muhammad Osama, Dave Zachariah, Thomas B. Schön.""Inferring Heterogeneous Causal Effects in Presence of Spatial Confounding"", Proceedings of the 36th International Conference on Machine Learning, 2019. [http://proceedings.mlr.press/v97/osama19a.html]

# Functions
The folder 'functions' contain the utility function in MATLAB script while their use is demostrated by the script 1d_cosine_effect.m qnd cont_2d.m file. Running the former script produces the plot below which show the true effect $\tau(s)$, its estimate $\widehat{\tau}(s)$ and the $95\%$ bootstrap confidence interval (CI) (for one dimensional problem only). For an example in two dimensional space, kindly run the other m-file cont_2D.m. For description of input and output of each function, kindly read the comments in the function definitions.

![cosine_gamma_1D_our_model](https://user-images.githubusercontent.com/37805794/58567665-05936700-8233-11e9-8c11-ed357e0feacd.png)

# Real Data 
Outcome y: number of crimes

Exposure z: number of poor families

Space s: state index s = {1,...,50}

Data source: U.S. Census Bureau year 2000

The figures below shows the estimated effect $\widehat{\tau}(s)$ in different states of U.S. and a plot of its significance level at 5% level.

![US_poverty_on_crime](https://user-images.githubusercontent.com/37805794/58567005-d6c8c100-8231-11e9-86e0-3a99732ed502.png)

![US_poverty_on_crime_significance](https://user-images.githubusercontent.com/37805794/58567272-5bb3da80-8232-11e9-9430-a62d4581a4b1.png)
