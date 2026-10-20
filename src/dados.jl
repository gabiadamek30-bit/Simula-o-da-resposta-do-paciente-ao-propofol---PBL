# Leitura dos dados de pressão e fluxo aórticos distribuídos com o windkessel-id.
# Porta de `importPQ` (external/windkessel-id/idanalysis.jl) e da leitura feita em
# id_stergiopulos1999total.jl.

const WINDKESSEL_ID_DIR = joinpath(@__DIR__, "..", "external", "windkessel-id")

"""
Lê `nfiles` pares de arquivos `<prefixo>P<i>.csv` / `<prefixo>Q<i>.csv`
(colunas `timestamp`, `vb.va`). Pressão em mmHg, fluxo em L/min.
Ex.: `importPQ(joinpath(WINDKESSEL_ID_DIR, "invivo_data", "Ao_baseline1_Ao"), 4)`.
"""
function importPQ(csvprefix::AbstractString, nfiles::Int)
    ler(f) = CSV.read(f, DataFrame; comment="#", types=[Time, Float64], normalizenames=true)
    P = [ler(csvprefix * "P$(i).csv") for i in 1:nfiles]
    Q = [ler(csvprefix * "Q$(i).csv") for i in 1:nfiles]
    for i in 1:nfiles
        Millisecond(-5) < (P[i].timestamp[1] - Q[i].timestamp[1]) < Millisecond(5) ||
            @warn "t0 differs by >5ms between P$(i) and Q$(i)"
        nrow(P[i]) == nrow(Q[i]) || @warn "P Q lengths not equal. Samples may be missing."
    end
    return P, Q
end

"""
Batimento tipo A de Stergiopulos et al. (1999), digitalizado: retorna
`(t, pressao_mmHg, fluxo_mLs)`, amostrado a 5 ms.
"""
function dados_stergiopulos1999()
    pasta = joinpath(WINDKESSEL_ID_DIR, "stergiopulos1999_data")
    ler(f) = CSV.read(joinpath(pasta, f), DataFrame; header=false, types=[Float64, Float64])
    P = ler("apressure.csv")
    Q = ler("aflow.csv")
    return P[:, 1], P[:, 2], Q[:, 2]
end

"Converte θ = [Rp, C, Rc, L] de unidades mL/s para unidades L/min (pressão em mmHg)."
mLs2Lmin(θ) = θ .* [1000 / 60, 60 / 1000, 1000 / 60, 1000 / 60][1:length(θ)]
