# Simulação completa: indução com propofol → BIS e pressão arterial.
# Rodar da raiz do repositório:  julia --project=. exemplos/simulacao_propofol.jl

using PropofolPBL, Plots, CSV

pac = Paciente(idade=62, peso=80, altura=172, sexo=:M, asa=2,
               pam=95, fc=72, spo2=97, glicemia=110, pvi=15, ivci=30)

println("IMC = $(round(imc(pac), digits=1)) kg/m², sensibilidade hemodinâmica = ",
        round(sensibilidade_hemodinamica(pac), digits=2))

# 2 mg/kg em 30 s + 100 μg/kg/min por 30 min (a partir de 1 min)
doses = esquema_inducao(pac; bolus_mgkg=2.0, manutencao_mcgkgmin=100.0)
pkpd, bat, onda = simular_paciente(pac; doses, tfim=30.0, forma_de_onda=true)

i = argmin(bat.PAM)
println("PAM basal $(round(bat.PAM[1], digits=1)) mmHg → mínima $(round(bat.PAM[i], digits=1)) mmHg ",
        "em t = $(round(bat.t_min[i], digits=1)) min (BIS $(round(bat.BIS[i], digits=0)))")

mkpath("resultados")
CSV.write("resultados/batimentos.csv", bat)

p1 = plot(pkpd.t, [pkpd.Cp pkpd.Ce], label=["Cp" "Ce"], ylabel="μg/mL", title="Propofol (Schnider)")
p2 = plot(pkpd.t, pkpd.bis, label="BIS", ylabel="BIS", ylims=(0, 100))
hline!(p2, [40, 60], ls=:dash, c=:gray, label="alvo 40–60")
p3 = plot(bat.t_min, [bat.PAS bat.PAM bat.PAD], label=["PAS" "PAM" "PAD"], ylabel="mmHg", xlabel="t [min]")
p4 = plot(onda.t_s[1:2000], onda.pressao[1:2000], label="pressão aórtica", xlabel="t [s]", ylabel="mmHg")
plot(p1, p2, p3, p4, layout=(2, 2), size=(1100, 700))
savefig("resultados/simulacao_propofol.png")
println("Figura salva em resultados/simulacao_propofol.png")
