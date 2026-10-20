# Simula-o-da-resposta-do-paciente-ao-propofol---PBL
Simulação de anestesia em centro cirúrgico: previsão da resposta do paciente ao propofol a partir de biomarcadores (idade, PA, FC, ASA, oxigenação, glicose, IMC, VCI/PVI). Projeto PBL de Engenharia Biomédica, Einstein.

## Estrutura

O pacote Julia `PropofolPBL` junta dois repositórios, incluídos via `git subtree` em `external/`:

| Origem | O que usamos |
|---|---|
| [shramish2057/NeoPKPD](https://github.com/shramish2057/NeoPKPD) → `external/NeoPKPD` | Motor PK/PD: modelo tricompartimental com infusões e compartimento de efeito (`ke0`) |
| [hpigot/windkessel-id](https://github.com/hpigot/windkessel-id) → `external/windkessel-id` | Modelos Windkessel de 2/3/4 elementos, identificação a partir de pressão/fluxo e dados aórticos |

```
src/
  PropofolPBL.jl   módulo principal
  windkessel.jl    porta de oeid.jl para Julia atual (sem ControlSystems)
  dados.jl         leitura dos dados do windkessel-id
  paciente.jl      struct Paciente com os biomarcadores
  propofol.jl      PK Schnider/Marsh + BIS (Eleveld) sobre o NeoPKPD
  acoplamento.jl   Ce(t) → Rp, volume sistólico, FC → Windkessel → PA batimento a batimento
exemplos/
  simulacao_propofol.jl     indução + manutenção: Cp, Ce, BIS, PAS/PAM/PAD
  identificacao_invivo.jl   identificação do Windkessel nos dados in vivo
test/runtests.jl
```

Fluxo da simulação:

```
Paciente ──► PK (Schnider/Marsh, NeoPKPD) ──► Cp(t) ──► ke0 ──► Ce(t) ──► BIS
                                                                  │
                    PAM, FC, DC basais ──► θ = [Rp, C, Rc, L] ◄───┘ (Emax sigmoide)
                                                  │
                    fluxo aórtico (Stergiopulos) ─┴─► Windkessel 4 elementos ──► pressão arterial
```

## Como rodar

Requer Julia ≥ 1.11.

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"   # primeira vez (demora: DifferentialEquations)
julia --project=. -e "using Pkg; Pkg.test()"
julia --project=. exemplos/simulacao_propofol.jl       # gera resultados/
```

Uso mínimo:

```julia
using PropofolPBL
pac = Paciente(idade=62, peso=80, altura=172, sexo=:M, asa=2, pam=95, fc=72, pvi=15)
pkpd, batimentos, onda = simular_paciente(pac; tfim=30.0, forma_de_onda=true)
```

## Biomarcadores e o que cada um faz no modelo

| Biomarcador | Uso |
|---|---|
| idade, peso, altura, sexo | PK de Schnider (massa magra de James); idade também no Ce50 do BIS |
| PAM, FC (e DC, se medido) | Parâmetros basais do Windkessel (Rp = PAM/DC) |
| ASA, PVI, IVCI, idade | Fator de sensibilidade hemodinâmica ao propofol |
| SpO₂, glicemia | Guardados em `Paciente`, ainda **não** entram em nenhuma equação |

## Limitações

- Os parâmetros de `EfeitoHemodinamico` (Ce50, γ, Imax) e os coeficientes de
  `sensibilidade_hemodinamica` são **ilustrativos**, não ajustados a dados; precisam ser calibrados.
- O Schnider, com V1 fixo de 4,27 L, superestima o pico de concentração após bolus.
- O BIS usa o PD de Eleveld (2018) com o Ce do ke0 do modelo PK escolhido — é uma aproximação.
- A forma de onda do fluxo é o batimento de Stergiopulos (1999) esticado para cada FC.
- A identificação de 4 elementos tem mínimos locais; os scripts originais em
  `external/windkessel-id` usam inicialização múltipla para escolher a solução.

## Licenças

NeoPKPD e windkessel-id são distribuídos sob suas próprias licenças (ver `external/*/LICENSE`).
