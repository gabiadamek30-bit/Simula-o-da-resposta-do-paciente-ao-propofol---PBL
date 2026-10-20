using PropofolPBL
using Test, Statistics

const h = 0.005

@testset "Windkessel (windkessel-id)" begin
    θ, _ = params_stergiopulos1999total()
    _, _, q = dados_stergiopulos1999()
    u = repeat(q, 5)

    @testset "solvex0 dá estado periódico" begin
        sys = G(θ, h)
        x0 = solvex0(sys, q)
        @test estado_final(sys, q, x0) ≈ x0 rtol = 1e-8
    end

    @testset "PAM do modelo ≈ Rp · fluxo médio" begin
        y = simulate_wk(θ, solvex0(G(θ, h), q), q, h)
        @test mean(y) ≈ θ[1] * mean(q) rtol = 0.02
    end

    @testset "bode: |Z(0)| = Rp e |Z(∞)| → Rc" begin
        mag, _ = bode_wk(θ, [1e-6, 1e6])
        @test mag[1] ≈ θ[1] rtol = 1e-4
        @test mag[2] ≈ θ[3] rtol = 1e-3
        mag3, _ = bode_wk(θ[1:3], [1e6]; G_=G3)
        @test mag3[1] ≈ θ[3] rtol = 1e-3
    end

    @testset "identify recupera parâmetros de dados sintéticos" begin
        y = simulate_wk(θ, solvex0(G(θ, h), u), u, h)
        θ0 = θ .* [1.0, 1.4, 0.7, 1.5]
        x00 = solvex0(G(θ0, h), u)
        phat, res = identify(θ0, x00, u, y, h; x0id=false, x0solve=true)
        θhat, _ = deal(phat, IdConditions(θ0, x00, [1, 1, 1, 1], false, true))
        @test θhat ≈ θ rtol = 1e-3
    end

    @testset "importPQ lê os dados in vivo" begin
        P, Q = importPQ(joinpath(WINDKESSEL_ID_DIR, "invivo_data", "Ao_baseline1_Ao"), 4)
        @test length(P) == 4 && length(Q) == 4
        @test 50 < mean(P[1].vb_va) < 150
    end
end

@testset "Propofol (NeoPKPD)" begin
    pac = Paciente(idade=53, peso=77, altura=177, sexo=:M)

    @testset "Schnider no paciente de referência" begin
        pk = parametros_pk(pac, :schnider)
        @test pk.V1 == 4.27
        @test pk.V2 ≈ 18.9
        @test pk.Q2 ≈ 1.29
        @test pk.CL ≈ 1.89 - 0.0681 * (massa_magra(pac) - 59) rtol = 1e-12
    end

    @testset "bolus de 2 mg/kg" begin
        r = simular_propofol(pac, [bolus(0.0, 2 * pac.peso)]; tfim=20.0)
        ipico = argmax(r.Ce)
        @test 1.0 < r.t[ipico] < 4.0          # pico no sítio efetor em ~1,5-2 min (Schnider)
        @test argmax(r.Cp) < ipico            # histerese: Cp atinge o pico antes de Ce
        @test minimum(r.bis) < 60
        @test r.bis[1] ≈ 93.0
    end

    @testset "BIS de Eleveld" begin
        Ce50 = 3.08
        @test bis_eleveld(Ce50, 35) ≈ 93 / 2
        @test bis_eleveld(0.0, 35) == 93.0
        @test bis_eleveld(2.0, 80) < bis_eleveld(2.0, 20)   # idoso mais sensível
    end
end

@testset "Acoplamento PK/PD → Windkessel" begin
    pac = Paciente(idade=45, peso=70, altura=170, sexo=:F, pam=90.0, fc=75.0)
    pkpd, bat, onda = simular_paciente(pac; tfim=10.0, forma_de_onda=true)

    @test bat.PAM[1] ≈ pac.pam rtol = 0.03          # basal reproduz a PAM do paciente
    @test minimum(bat.PAM) < 0.85 * pac.pam          # hipotensão após a indução
    @test all(bat.PAS .> bat.PAM .> bat.PAD)
    @test length(onda.pressao) == length(onda.t_s)

    # Paciente mais sensível (idoso, ASA 3, hipovolêmico) tem queda maior.
    pac2 = Paciente(idade=75, peso=70, altura=170, sexo=:F, pam=90.0, fc=75.0, asa=3, pvi=20.0)
    @test sensibilidade_hemodinamica(pac2) > sensibilidade_hemodinamica(pac)
    _, bat2, _ = simular_paciente(pac2; tfim=10.0)
    @test minimum(bat2.PAM) < minimum(bat.PAM)
end
