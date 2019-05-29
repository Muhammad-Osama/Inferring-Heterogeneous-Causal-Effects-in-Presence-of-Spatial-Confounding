# Inferring-Heterogeneous-Causal-Effects-in-Presence-of-Spatial-Confounding
Given a scalar outcome of interest $y$ and a scaler exposure variable $z$, which have been collected over space $s$ so that the observed data in D={y_i,z_i,s_i}_{i=1}^{n}, the code allows to estimate the causal effect $\tau$ of $z$ on $y$ while mitigating for spatial confounding caused by other unobserved spatially varying variable $c$. The effect $\tau$ can itself be spatially varying i.e. heterogeneous. As a motivating example, suppose that the outcome $y$ is crop yield and the exposure $z$ is amount of fertilizer and we are interested in studying whether using more fertilizer increases the crop yield. Here the confounding variable $c$ is the unobserved inherent variability of the soil. Moreover, the effect of fertilizer on crop yield may vary spatially. For details see and cite: Muhammad Osama, Dave Zachariah, Thomas B. Schön.""Inferring Heterogeneous Causal Effects in Presence of Spatial Confounding"", Proceedings of the 36th International Conference on Machine Learning, 2019. [http://proceedings.mlr.press/v97/osama19a.html]

# Functions
The folder 'supporting functions' contain the utility function in MATLAB script while their use is demostrated by the script 1d_cosine_effect.m. Running the script produces the plot below which show the true effect $\tau(s)$, its estimate $\widehat{\tau}(s)$ and the $95\%$ bootstrap confidence interval (CI). For description of input and output of each function, kindly read the comments in the function definitions.

![cosine_gamma_1D_our_model](https://user-images.githubusercontent.com/37805794/58563249-1809a280-822b-11e9-8060-69fdff24b6e6.png)

# Real Data 

