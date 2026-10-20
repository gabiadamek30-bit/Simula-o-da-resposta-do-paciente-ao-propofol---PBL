"""
Simulação da resposta do paciente ao propofol.

Junta dois projetos (ver `external/`):
- NeoPKPD (shramish2057/NeoPKPD): motor PK/PD — modelo tricompartimental,
  infusões e compartimento de efeito.
- windkessel-id (hpigot/windkessel-id): modelos Windkessel de 2/3/4 elementos
  e sua identificação a partir de pressão e fluxo aórticos.
"""
module PropofolPBL

using LinearAlgebra, Statistics, Dates
using CSV, DataFrames
using Optim, LineSearches, ForwardDiff
using NeoPKPD

include("windkessel.jl")
include("dados.jl")
include("paciente.jl")
include("propofol.jl")
include("acoplamento.jl")

# Windkessel (windkessel-id)
export WKSistema, G, G2, G3, lsim_wk, estado_final, simulate_wk, solvex0, bode_wk
export IdConditions, identify, identifyglobal, deal, cost, params_stergiopulos1999total
export WINDKESSEL_ID_DIR, importPQ, dados_stergiopulos1999, mLs2Lmin
# Paciente e propofol (NeoPKPD)
export Paciente, imc, massa_magra, superficie_corporal, debito_basal
export parametros_pk, bolus, infusao, esquema_inducao, bis_eleveld, simular_propofol
# Acoplamento
export EfeitoHemodinamico, sensibilidade_hemodinamica, theta_basal
export simular_hemodinamica, simular_paciente

end
