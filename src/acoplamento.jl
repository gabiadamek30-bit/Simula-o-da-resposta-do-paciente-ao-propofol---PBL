# Acoplamento PK/PD (NeoPKPD) → hemodinâmica (Windkessel do windkessel-id).
#
# A concentração no sítio efetor Ce(t) reduz, por uma função Emax sigmoide,
# a resistência periférica Rp, o volume sistólico e a FC. A cada batimento o
# Windkessel de 4 elementos é rediscretizado com os parâmetros daquele instante
# e simulado continuando o estado do batimento anterior, gerando a curva de
# pressão arterial e PAS/PAD/PAM por batimento.
#
# ATENÇÃO: os valores padrão de `EfeitoHemodinamico` e de
# `sensibilidade_hemodinamica` são ilustrativos (ordem de grandeza da queda de
# PAM na indução), não ajustados a dados. Devem ser calibrados pelo grupo.

"""
Efeito do propofol sobre a hemodinâmica. Cada variável X é multiplicada por
`1 - Imax_X · s · Ceᵞ/(Ce50ᵞ + Ceᵞ)`, com `s` a sensibilidade do paciente.
"""
Base.@kwdef struct EfeitoHemodinamico
    Ce50::Float64 = 3.0       # μg/mL
    γ::Float64 = 2.0
    Imax_Rp::Float64 = 0.25   # vasodilatação arterial
    Imax_VS::Float64 = 0.15   # depressão miocárdica / redução de pré-carga
    Imax_FC::Float64 = 0.05
end

"""
Fator de sensibilidade hemodinâmica (1 = paciente de referência). Aumenta com
idade > 50 anos, ASA ≥ 3, PVI > 13 % e IVCI > 40 % (sinais de hipovolemia),
fatores associados a hipotensão pós-indução. Coeficientes: hipóteses a calibrar.
"""
function sensibilidade_hemodinamica(pac::Paciente; k_idade=0.01, k_asa=0.15, k_pvi=0.03, k_ivci=0.01)
    return 1.0 + k_idade * max(pac.idade - 50, 0) + k_asa * max(pac.asa - 2, 0) +
           k_pvi * max(pac.pvi - 13, 0) + k_ivci * max(pac.ivci - 40, 0)
end

"""
Parâmetros basais θ = [Rp, C, Rc, L] (unidades mL/s) do Windkessel de 4 elementos.
Rp vem de PAM/DC (em regime, Z(0) = Rp); C, Rc e L de Stergiopulos (1999).
Para usar parâmetros identificados de um registro do próprio paciente, passe
`θbasal` (ex.: resultado de `identify`) para `simular_hemodinamica`.
"""
function theta_basal(pac::Paciente)
    θ, _ = params_stergiopulos1999total()
    Qmedio = debito_basal(pac) * 1000 / 60      # mL/s
    θ[1] = pac.pam / Qmedio
    return θ
end

"""
Forma de onda normalizada do fluxo aórtico (área unitária em τ ∈ [0, 1)),
a partir do batimento de Stergiopulos (1999) distribuído com o windkessel-id.
"""
function forma_fluxo()
    _, _, q = dados_stergiopulos1999()
    return q ./ (sum(q) / length(q))
end

"Reamostra `forma` (um período) para `n` pontos e escala para ejetar `vs` mL em `n·h` s."
function fluxo_batimento(forma, n, vs, h)
    m = length(forma)
    q = map(0:n-1) do k
        s = k * m / n
        i = floor(Int, s)
        f = s - i
        (1 - f) * forma[i+1] + f * forma[mod(i + 1, m)+1]
    end
    return q .* (vs / (sum(q) * h))
end

"Interpolação linear de (tx, yx) em t, constante fora do intervalo."
function interp(tx, yx, t)
    t <= tx[1] && return yx[1]
    t >= tx[end] && return yx[end]
    i = searchsortedlast(tx, t)
    f = (t - tx[i]) / (tx[i+1] - tx[i])
    return (1 - f) * yx[i] + f * yx[i+1]
end

"""
Simula a pressão arterial batimento a batimento durante a PK/PD `pkpd`
(saída de `simular_propofol`).

Retorna `(batimentos, onda)`: `batimentos` é um DataFrame com, por batimento,
t_min, Ce, BIS, FC, VS (mL), DC (L/min), Rp, PAS, PAD e PAM; `onda` é
`(t_s, pressao, fluxo)` se `forma_de_onda=true`, senão `nothing`.
"""
function simular_hemodinamica(pac::Paciente, pkpd; efeito=EfeitoHemodinamico(),
                              θbasal=theta_basal(pac), h=0.005, tfim=pkpd.t[end],
                              forma=forma_fluxo(), forma_de_onda=false)
    s = sensibilidade_hemodinamica(pac)
    fator(Imax, Ce) = max(1 - Imax * s * Ce^efeito.γ / (efeito.Ce50^efeito.γ + Ce^efeito.γ), 0.05)

    Rp0, C, Rc, L = θbasal
    vs0 = debito_basal(pac) * 1000 / pac.fc     # mL por batimento

    # Estado periódico inicial no batimento basal.
    n0 = round(Int, 60 / pac.fc / h)
    x = solvex0(G(θbasal, h), fluxo_batimento(forma, n0, vs0, h))

    colunas = [:t_min, :Ce, :BIS, :FC, :VS, :DC, :Rp, :PAS, :PAD, :PAM]
    batimentos = DataFrame([c => Float64[] for c in colunas])
    ts, ps, qs = Float64[], Float64[], Float64[]

    t = 0.0                                      # s
    while t < tfim * 60
        Ce = interp(pkpd.t, pkpd.Ce, t / 60)
        fc = pac.fc * fator(efeito.Imax_FC, Ce)
        vs = vs0 * fator(efeito.Imax_VS, Ce)
        Rp = Rp0 * fator(efeito.Imax_Rp, Ce)

        n = round(Int, 60 / fc / h)
        u = fluxo_batimento(forma, n, vs, h)
        sys = G([Rp, C, Rc, L], h)
        y, = lsim_wk(sys, u, x)
        x = estado_final(sys, u, x)

        push!(batimentos, (t / 60, Ce, interp(pkpd.t, pkpd.bis, t / 60), fc, vs, vs * fc / 1000,
                           Rp, maximum(y), minimum(y), mean(y)))
        if forma_de_onda
            append!(ts, t .+ h .* (0:n-1))
            append!(ps, y)
            append!(qs, u)
        end
        t += n * h
    end
    return batimentos, forma_de_onda ? (; t_s=ts, pressao=ps, fluxo=qs) : nothing
end

"""
Atalho: PK/PD (NeoPKPD) + hemodinâmica (Windkessel) para um paciente.
Retorna `(pkpd, batimentos, onda)`.
"""
function simular_paciente(pac::Paciente; doses=esquema_inducao(pac), modelo=:schnider,
                          tfim=30.0, kwargs...)
    pkpd = simular_propofol(pac, doses; modelo, tfim)
    batimentos, onda = simular_hemodinamica(pac, pkpd; tfim, kwargs...)
    return pkpd, batimentos, onda
end
