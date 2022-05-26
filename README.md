# CUSP: Cumulative Shrinkage Process

This repository is associated with the article **Bayesian cumulative shrinkage for infinite factorizations** ([Biometrika, 107(3), 745-752](https://academic.oup.com/biomet/advance-article-abstract/doi/10.1093/biomet/asaa008/5847840); also on [arXiv](http://arxiv.org/abs/1902.04349)) and provides detailed materials and codes for posterior inference under the *Cumulative Shrinkage Process (CUSP)* proposed in the article.

The repository contains:

- [`cusp.R`](https://github.com/siriolegramanti/CUSP/blob/master/cusp.R): the source code to perform posterior computation for Gaussian factor models endowed with the CUSP prior via the adaptive Gibbs sampler provided in the paper (Algorithm 2). Gaussian factor models are used here as a notable illustrative example, but the CUSP can be employed in a variety of models (e.g. Poisson factorization). Hence, the parts of `cusp.R` which are not specific to Gaussian factor models can be repurposed for other models. Note that this is a basic `R` implementation meant to provide a clear understanding of the computational routines associated with the proposed model; faster `C++` implementations can be derived. 

- [`tutorial.md`](https://github.com/siriolegramanti/CUSP/blob/master/tutorial.md): a comprehensive tutorial to perform posterior inference for Gaussian factor models endowed with the CUSP prior, leveraging the source code `cusp.R`. As an illustrative example, we reproduce step by step the analysis of the `bfi` dataset from the `R` package `psych` described in Section 5 of the article. The analyses are performed with a laptop equipped with Ubuntu 18.04.5 OS, Intel Core i7-3632QM CPU and 7.7 GB of RAM, using `R` version 3.4.2.

**Acknowledgements** The maintenance of this code has been funded by the LIFEPLAN project (European Union’s Horizon 2020 research and innovation program; grant agreement No 856506).
