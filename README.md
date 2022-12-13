# materialNonLinearity.jl

## Introduction

FEM code to analyze 2D Euler Bernoulli beams with linear and nonlinear constitutive laws. 

## Standard constitutive laws

Implemented isotropic constitutive laws are:

* Linear elastic
* Elastic perfectly plastic
* Linear hardening / softening

## User defined constitutive laws

User can define a constitutive law defined by any kind of function. Examples of a cubic function with softening and a concrete with softening in compression and tension branches are provided.

| **Documentation** | **Tests** | **Coverage** | **License** |
|:-----------------:|:---------------:|:------------:|:------------:|
| [![docs-dev][dev-img]][dev-url] | [![CI][ci-img]][ci-url] | [![codecov][cov-img]][cov-url] | [![license][lic-img]][lic-url] |

[dev-img]: https://img.shields.io/badge/docs-latest-blue.svg
[dev-url]: https://JoaquinViera.github.io/materialNonLinearity.jl/
[ci-img]: https://github.com/JoaquinViera/materialNonLinearity.jl/actions/workflows/ci.yml/badge.svg?branch=main
[ci-url]: https://github.com/JoaquinViera/materialNonLinearity.jl/actions/workflows/ci.yml
[cov-img]: https://codecov.io/gh/JoaquinViera/materialNonLinearity.jl/branch/main/graph/badge.svg?token=PF2QWFHHQ0
[cov-url]: https://codecov.io/gh/JoaquinViera/materialNonLinearity.jl
[lic-img]: https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000
[lic-url]: https://github.com/JoaquinViera/materialNonLinearity.jl/blob/main/LICENSE.md
