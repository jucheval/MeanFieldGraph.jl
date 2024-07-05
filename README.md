# MeanFieldGraph

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jucheval.github.io/MeanFieldGraph.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jucheval.github.io/MeanFieldGraph.jl/dev/)
[![Build Status](https://github.com/jucheval/MeanFieldGraph.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jucheval/MeanFieldGraph.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jucheval/MeanFieldGraph.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jucheval/MeanFieldGraph.jl)

A Julia package devoted to the estimation of the connectivity parameter of a latent (non sparse) high dimensional graph. The data are assumed to come from a model with mean field interaction driven by the latent graph.

## Ressources

* **Related articles** :
    * [Chevallier, Löcherbach, Ost (2024)](https://arxiv.org/abs/2406.07066) : Discrete time data modeled via a binary Markov Chain. *This model is fully covered by the package.*
    * [Delattre, Fournier (2016)](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-10/issue-1/Statistical-inference-versus-mean-field-limit-for-Hawkes-processes/10.1214/16-EJS1142.full) and [Liu (2020)](https://www.sciencedirect.com/science/article/pii/S0304414919301048) : Continuous time data modeled via a Hawkes process ([Liu (2020)](https://www.sciencedirect.com/science/article/pii/S0304414919301048) extends the first paper by considering partial observation). *This model is not covered by the package for now, but is expected to.*

* **Documentation** : <https://jucheval.github.io/MeanFieldGraph.jl/>