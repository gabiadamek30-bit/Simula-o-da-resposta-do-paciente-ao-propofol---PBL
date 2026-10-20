# Windkessel model identification

**Tools and example data to identify lumped-parameter models of cardiac afterloads**, specifically 2-element Windkessel, 3-element Windkessel, and parallel 4-element Windkessel models from measured pressure and flow. The identification is done with Newton's method using forward-mode automatic differentiation in Julia.

The data and methods presented here were used to generate the results presented in the article "Identification of cardiac afterload dynamics from data".

The identification code uses [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl), [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl), [LineSearches.jl](https://github.com/JuliaNLSolvers/LineSearches.jl), [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), [GenericLinearAlgebra.jl](https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl/), and [Distributions.jl](https://github.com/JuliaStats/Distributions.jl). Data handling, additional analysis, and plotting is done using [CSV.jl](https://github.com/JuliaData/CSV.jl), [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl), [Plots.jl](https://github.com/JuliaPlots/Plots.jl), [Plotly.jl](https://github.com/plotly/Plotly.jl), [StatsPlots.jl](https://github.com/JuliaPlots/StatsPlots.jl), [FFTW.jl](https://github.com/JuliaMath/FFTW.jl), and [ToeplitzMatrices.jl](https://github.com/JuliaMatrices/ToeplitzMatrices.jl). Thank you to the developers behind those packages.

## Identification tools

[oeid.jl](oeid.jl) contains the functions used to identify Windkessel model parameters (with user-facing functions `identify` for solving with Newton's Method, `identifyglobal` for solving with Nelder-Mead). The models are defined by functions `G2`, `G3`, and `G` for 2, 3, and 4 element Windkessel models, respectively. Estimates of the models' initial state condition are solved using the function `solvex0`, which leverages the periodicity of the input signals. [expm.jl](expm.jl) is used by [oeid.jl](oeid.jl) for natively computing matrix exponentials.

## Data and identification code

Aortic pressure and flow measurements are provided from three sources. In addition to [oeid.jl](oeid.jl), [idanalysis.jl](idanalysis.jl) provides multiple initialization as well as analysis and plotting functions.

[stergiopulos1999_data/](stergiopulos1999_data/) contains previously published human data, from Figure 4A type A beat of [_Total Arterial Inertance as the Fourth Element of the Windkessel Model,_ Stergiopulos et al. (1999)](http://doi.org/10.1152/ajpheart.1999.276.1.H81), digitized using [Webplotdigitizer v4.4]({https://automeris.io/WebPlotDigitizer). Pressure in mmHg, flow in mL/s.

[invivo_data/](invivo_data/) contains data measured in vivo in a 70kg pig. Pressure in mmHg, flow in L/min.

[exvivo_data/](exvivo_data/) contains data measured with a 70kg pig heart beating ex vivo against a synthetic afterload. Pressure in mmHg, flow in L/min.

Windkessel model parameter identification for the three data sources is done in [id_stergiopulos1999total.jl](id_stergiopulos1999total.jl), [id_invivo.jl](id_invivo.jl), and [id_exvivo.jl](id_exvivo.jl).

Persistance of excitation of the three input signals is analysed in [spectrum.jl](spectrum.jl).

## Getting started

1. Download [Julia](https://julialang.org/).
2. Install packages:
   1. Open Julia REPL, e.g. using the command `julia` in terminal,
   2. Enter the Pkg REPL by pressing ]
   3. add packages using command `add DataFrames Dates Plots StatsPlots Distributions LinearAlgebra GenericLinearAlgebra ControlSystems Optim LineSearches ForwardDiff Printf CSV FFTW ToeplitzMatrices`
3. Open and run one of the `id_*.jl` files line-by-line in a text editor with Julia extension (e.g. VS Code, or Atom with Juno)

Note that this code was tested with the following package versions:

````
  [336ed68f] CSV v0.8.3
  [a6e380b2] ControlSystems v0.9.0
  [a93c6f00] DataFrames v0.22.5
  [31c24e10] Distributions v0.24.12
  [7a1cc6ca] FFTW v1.3.2
  [f6369f11] ForwardDiff v0.10.16
  [14197337] GenericLinearAlgebra v0.2.4
  [d3d80556] LineSearches v7.1.1
  [429524aa] Optim v1.2.4
  [91a5bcdd] Plots v1.10.4
  [f3b207a7] StatsPlots v0.14.19
  [c751599d] ToeplitzMatrices v0.6.3
  [37e2e46d] LinearAlgebra
````