# MeanFieldGraph

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jucheval.github.io/MeanFieldGraph.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jucheval.github.io/MeanFieldGraph.jl/dev/)
[![Build Status](https://github.com/jucheval/MeanFieldGraph.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jucheval/MeanFieldGraph.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jucheval/MeanFieldGraph.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jucheval/MeanFieldGraph.jl)

A Julia package devoted to the statistical analysis of a latent (non sparse) high dimensional graphs. The data are assumed to come from a model with mean field interaction driven by the latent graph. Moreover, the model may incorporate a community structure for the nodes of the graph: one excitatory and one inhibitory community. The two objectives are:
1. The estimation of the connectivity parameter of the latent graph (which is independent of the community).
2. The detection of the two communities.

## Ressources

* **Related articles** :
    * [Chevallier, Löcherbach, Ost (2024)](https://arxiv.org/abs/2406.07066) and [Chevallier, Ost (2024)](https://arxiv.org/abs/2411.15627): Discrete time data modeled via a binary Markov Chain. *This model is fully covered by the package.*
    * [Delattre, Fournier (2016)](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-10/issue-1/Statistical-inference-versus-mean-field-limit-for-Hawkes-processes/10.1214/16-EJS1142.full) and [Liu (2020)](https://www.sciencedirect.com/science/article/pii/S0304414919301048) : Continuous time data modeled via a Hawkes process ([Liu (2020)](https://www.sciencedirect.com/science/article/pii/S0304414919301048) extends the first paper by considering partial observation). *This model is not covered by the package for now, but is expected to. There is no community structure here.*

* **Documentation** : <https://jucheval.github.io/MeanFieldGraph.jl/>

## Examples

The folder [examples](examples/) contains all the material used to produce the figures in [Chevallier, Löcherbach, Ost (2024)](https://arxiv.org/abs/2406.07066) (files named ```CLO24_*.jl```, *the plots may differ a little bit because the files are newer than the article*) and [Chevallier, Ost (2024)](https://arxiv.org/abs/2411.15627) (files named ```CO24_*.jl```). The files named:
- ```*_simulation_*.jl``` contain the commands to produce the data used,
- ```*_plot_*.jl``` contain the commands used to produce the plots (assuming that you have already simulated data in the folder ```data```).