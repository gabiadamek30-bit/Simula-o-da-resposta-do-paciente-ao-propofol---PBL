# PK/PD do propofol sobre o motor do NeoPKPD.
# Tempo em minutos, doses em mg, volumes em L ⇒ concentrações em mg/L = μg/mL.

"""
Parâmetros PK do propofol para o paciente: `(CL, V1, Q2, V2, Q3, V3, ke0)`.
`modelo = :schnider` (Schnider et al. 1998/1999) ou `:marsh` (Marsh et al. 1991, ke0 = 0,26/min).
"""
function parametros_pk(pac::Paciente, modelo::Symbol=:schnider)
    if modelo == :schnider
        lbm = massa_magra(pac)
        V1 = 4.27
        V2 = 18.9 - 0.391 * (pac.idade - 53)
        V3 = 238.0
        CL = 1.89 + 0.0456 * (pac.peso - 77) - 0.0681 * (lbm - 59) + 0.0264 * (pac.altura - 177)
        Q2 = 1.29 - 0.024 * (pac.idade - 53)
        Q3 = 0.836
        ke0 = 0.456
    elseif modelo == :marsh
        V1 = 0.228 * pac.peso
        CL = 0.119 * V1
        Q2 = 0.112 * V1
        V2 = Q2 / 0.055
        Q3 = 0.0419 * V1
        V3 = Q3 / 0.0033
        ke0 = 0.26
    else
        throw(ArgumentError("modelo PK desconhecido: $(modelo)"))
    end
    return (; CL, V1, Q2, V2, Q3, V3, ke0)
end

"Bolus de `dose_mg` em t (min), administrado em `duracao` min (padrão 30 s)."
bolus(t, dose_mg; duracao=0.5) = DoseEvent(Float64(t), Float64(dose_mg), Float64(duracao))

"Infusão contínua de `taxa_mg_min` de t_inicio a t_inicio + duracao (min)."
infusao(t_inicio, duracao, taxa_mg_min) =
    DoseEvent(Float64(t_inicio), Float64(taxa_mg_min * duracao), Float64(duracao))

"""
Esquema típico: bolus de indução (mg/kg) em t = 0 seguido de manutenção (μg/kg/min).
"""
function esquema_inducao(pac::Paciente; bolus_mgkg=2.0, manutencao_mcgkgmin=100.0,
                         inicio_manutencao=1.0, duracao_manutencao=30.0)
    doses = [bolus(0.0, bolus_mgkg * pac.peso)]
    if manutencao_mcgkgmin > 0 && duracao_manutencao > 0
        push!(doses, infusao(inicio_manutencao, duracao_manutencao,
                             manutencao_mcgkgmin * pac.peso / 1000))
    end
    return doses
end

"""
BIS em função da concentração no sítio efetor (modelo sigmoide de Eleveld et al. 2018):
E0 = 93, Ce50 = 3,08·exp(-0,00635·(idade-35)) μg/mL, γ = 1,47 abaixo de Ce50 e 1,89 acima.
Observação: aqui o Ce vem do ke0 do modelo PK escolhido, não do PK de Eleveld.
"""
function bis_eleveld(Ce, idade)
    E0 = 93.0
    Ce50 = 3.08 * exp(-0.00635 * (idade - 35))
    γ = Ce <= Ce50 ? 1.47 : 1.89
    return sigmoid_emax_effect(max(Float64(Ce), 0.0), E0, -E0, Ce50, γ)
end

"""
Simula a PK/PD do propofol com o NeoPKPD.
Retorna `(t, Cp, Ce, bis, pk)` com t em minutos e concentrações em μg/mL.
"""
function simular_propofol(pac::Paciente, doses::Vector{DoseEvent}; modelo=:schnider,
                          tfim=60.0, dt=1 / 60)
    pk = parametros_pk(pac, modelo)
    spec = ModelSpec(ThreeCompIVBolus(), "propofol_$(modelo)",
                     ThreeCompIVBolusParams(pk.CL, pk.V1, pk.Q2, pk.V2, pk.Q3, pk.V3),
                     sort(doses; by=d -> d.time))
    grid = SimGrid(0.0, Float64(tfim), collect(0.0:dt:tfim))
    solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10^7)
    res = NeoPKPD.simulate(spec, grid, solver)
    Cp = res.observations[:conc]
    Ce, _ = simulate_biophase(res.t, Cp, BiophaseEquilibrationParams(pk.ke0, 0.0, 1.0, 1.0), solver)
    bis = bis_eleveld.(Ce, pac.idade)
    return (; t=res.t, Cp, Ce, bis, pk)
end
