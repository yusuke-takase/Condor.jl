# Condor

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yusuke-takase.github.io/Condor.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yusuke-takase.github.io/Condor.jl/dev)
[![Build Status](https://travis-ci.com/yusuke-takase/Condor.jl.svg?branch=master)](https://travis-ci.com/yusuke-takase/Condor.jl)

This tool simulates satellite observations of the CMB.　\
First, a preconvolution file is created according to the HEALPix convention by 2DFFT. \
It then works with Falcons.jl to output the observed universe from the satellite scan information.

## Installation
From the Julia REPL, run

```julia
import Pkg
Pkg.add("Condor")
```
