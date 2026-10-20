# Biomarcadores do paciente usados na simulação.

"""
Dados do paciente. Unidades: idade em anos, peso em kg, altura em cm,
pressões em mmHg, FC em bpm, SpO₂ em %, glicemia em mg/dL, PVI e IVCI em %.

Uso no modelo:
- PK (Schnider/Marsh): idade, peso, altura, sexo.
- PD (BIS, Eleveld): idade.
- Hemodinâmica basal (Windkessel): `pam`, `fc`, superfície corporal (`dc` opcional).
- Sensibilidade hemodinâmica ao propofol: idade, ASA, PVI, IVCI
  (ver `sensibilidade_hemodinamica`).
- `spo2` e `glicemia` ainda não entram em nenhuma equação.
"""
Base.@kwdef struct Paciente
    idade::Float64
    peso::Float64
    altura::Float64
    sexo::Symbol = :F           # :M ou :F
    asa::Int = 1
    pam::Float64 = 90.0         # pressão arterial média basal
    fc::Float64 = 75.0          # frequência cardíaca basal
    dc::Float64 = NaN           # débito cardíaco basal [L/min]; NaN ⇒ 3,0 L/min/m² × SC
    spo2::Float64 = 98.0
    glicemia::Float64 = 95.0
    pvi::Float64 = 10.0         # pleth variability index
    ivci::Float64 = 20.0        # índice de colapsabilidade da veia cava inferior
end

imc(p::Paciente) = p.peso / (p.altura / 100)^2

"Massa magra (fórmula de James), usada por Schnider. Perde validade para IMC muito alto."
function massa_magra(p::Paciente)
    if p.sexo == :M
        return 1.1 * p.peso - 128 * (p.peso / p.altura)^2
    elseif p.sexo == :F
        return 1.07 * p.peso - 148 * (p.peso / p.altura)^2
    else
        throw(ArgumentError("sexo deve ser :M ou :F, recebido $(p.sexo)"))
    end
end

"Superfície corporal [m²] (Mosteller)."
superficie_corporal(p::Paciente) = sqrt(p.altura * p.peso / 3600)

"Débito cardíaco basal [L/min]."
debito_basal(p::Paciente) = isnan(p.dc) ? 3.0 * superficie_corporal(p) : p.dc
