# Identificação do Windkessel de 4 elementos nos dados in vivo (porco de 70 kg)
# do windkessel-id — equivalente à primeira parte de id_invivo_2020-12-01.jl.
# Os parâmetros identificados podem ser usados como `θbasal` em `simular_hemodinamica`
# (convertendo para mL/s, já que estes dados estão em L/min).
# Rodar da raiz do repositório:  julia --project=. exemplos/identificacao_invivo.jl

using PropofolPBL, Plots, Statistics

AoP, AoQ = importPQ(joinpath(WINDKESSEL_ID_DIR, "invivo_data", "Ao_baseline1_Ao"), 4)
h = 0.005

u = vcat([q.vb_va for q in AoQ]...)    # fluxo [L/min]
y = vcat([p.vb_va for p in AoP]...)    # pressão [mmHg]

θ0 = [mean(y) / mean(u), 1.0, 1.0, 1.0]   # [Rp, C, Rc, L] em unidades L/min
x00 = solvex0(G(θ0, h), u)
@time phat, res = identify(θ0, x00, u, y, h; x0id=false, x0solve=true)
θhat, _ = deal(phat, IdConditions(θ0, x00, [1, 1, 1, 1], false, true))

yhat = simulate_wk(θhat, solvex0(G(θhat, h), u), u, h)
rmse = sqrt(mean((y .- yhat) .^ 2))
println("θ̂ = [Rp, C, Rc, L] = ", round.(θhat, sigdigits=3), "  (RMSE = $(round(rmse, digits=2)) mmHg)")
println("em mL/s: ", round.(θhat ./ [1000 / 60, 60 / 1000, 1000 / 60, 1000 / 60], sigdigits=3))

t = h .* (0:length(u)-1)
plot(plot(t, [y yhat], label=["medido" "modelo"], ylabel="AoP [mmHg]"),
     plot(t, u, label="AoQ", ylabel="L/min", xlabel="t [s]"), layout=(2, 1), size=(1000, 500))
mkpath("resultados")
savefig("resultados/identificacao_invivo.png")
