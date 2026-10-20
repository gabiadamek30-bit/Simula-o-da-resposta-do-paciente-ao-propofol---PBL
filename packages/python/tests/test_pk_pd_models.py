"""
Comprehensive tests for all Python PK and PD models.
"""

import pytest
from neopkpd import init_julia

# Initialize Julia once for all tests
@pytest.fixture(scope="module", autouse=True)
def setup_julia():
    init_julia()


# =============================================================================
# PK Model Tests
# =============================================================================

class TestPKOneComp:
    """Tests for one-compartment PK models."""

    def test_iv_bolus_basic(self):
        from neopkpd import simulate_pk_iv_bolus
        res = simulate_pk_iv_bolus(
            cl=5.0, v=50.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            t0=0.0, t1=12.0,
            saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0]
        )
        assert res["t"][0] == 0.0
        assert "conc" in res["observations"]
        assert len(res["observations"]["conc"]) == 6
        # Initial concentration should be dose/V = 100/50 = 2
        assert abs(res["observations"]["conc"][0] - 2.0) < 0.01

    def test_iv_infusion(self):
        from neopkpd import simulate_pk_iv_bolus
        res = simulate_pk_iv_bolus(
            cl=5.0, v=50.0,
            doses=[{"time": 0.0, "amount": 100.0, "duration": 1.0}],  # 1-hour infusion
            t0=0.0, t1=12.0,
            saveat=[0.0, 0.5, 1.0, 2.0, 4.0]
        )
        assert "conc" in res["observations"]
        # Concentration at t=0 should be 0 (infusion just started)
        assert res["observations"]["conc"][0] < 0.1

    def test_oral_first_order(self):
        from neopkpd import simulate_pk_oral_first_order
        res = simulate_pk_oral_first_order(
            ka=1.0, cl=5.0, v=50.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            t0=0.0, t1=24.0,
            saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        )
        assert "conc" in res["observations"]
        # Concentration at t=0 should be 0 (not yet absorbed)
        assert res["observations"]["conc"][0] < 0.01
        # Peak should occur around Tmax = ln(ka/ke)/(ka-ke)
        assert max(res["observations"]["conc"]) > 0

    def test_oral_with_alag_bioavailability(self):
        from neopkpd import simulate_pk_oral_first_order
        res = simulate_pk_oral_first_order(
            ka=1.0, cl=5.0, v=50.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            t0=0.0, t1=24.0,
            saveat=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0],
            alag=0.5,  # 30 min lag
            bioavailability=0.8  # 80% F
        )
        assert "conc" in res["observations"]
        # At t=0.5, concentration should still be near 0 due to lag
        assert res["observations"]["conc"][1] < 0.1


class TestPKTwoComp:
    """Tests for two-compartment PK models."""

    def test_twocomp_iv_bolus(self):
        from neopkpd import simulate_pk_twocomp_iv_bolus
        res = simulate_pk_twocomp_iv_bolus(
            cl=10.0, v1=50.0, q=5.0, v2=100.0,
            doses=[{"time": 0.0, "amount": 500.0}],
            t0=0.0, t1=48.0,
            saveat=list(range(49))
        )
        assert "conc" in res["observations"]
        assert "A_central" in res["states"]
        assert "A_peripheral" in res["states"]
        # Initial concentration = dose/V1 = 500/50 = 10
        assert abs(res["observations"]["conc"][0] - 10.0) < 0.1

    def test_twocomp_oral(self):
        from neopkpd import simulate_pk_twocomp_oral
        res = simulate_pk_twocomp_oral(
            ka=1.0, cl=10.0, v1=50.0, q=5.0, v2=100.0,
            doses=[{"time": 0.0, "amount": 500.0}],
            t0=0.0, t1=48.0,
            saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        )
        assert "conc" in res["observations"]
        assert "A_gut" in res["states"]


class TestPKThreeComp:
    """Tests for three-compartment PK models."""

    def test_threecomp_iv_bolus(self):
        from neopkpd import simulate_pk_threecomp_iv_bolus
        res = simulate_pk_threecomp_iv_bolus(
            cl=10.0, v1=50.0, q2=10.0, v2=80.0, q3=2.0, v3=200.0,
            doses=[{"time": 0.0, "amount": 1000.0}],
            t0=0.0, t1=72.0,
            saveat=list(range(73))
        )
        assert "conc" in res["observations"]
        assert "A_central" in res["states"]
        assert "A_periph1" in res["states"]
        assert "A_periph2" in res["states"]


class TestPKAdvanced:
    """Tests for advanced PK models."""

    def test_transit_absorption(self):
        from neopkpd import simulate_pk_transit_absorption
        res = simulate_pk_transit_absorption(
            n=5, ktr=2.0, ka=1.0, cl=10.0, v=50.0,
            doses=[{"time": 0.0, "amount": 500.0}],
            t0=0.0, t1=24.0,
            saveat=[x * 0.5 for x in range(49)]
        )
        assert "conc" in res["observations"]
        # Should have transit compartments in states
        assert len(res["states"]) >= 5

    def test_michaelis_menten(self):
        from neopkpd import simulate_pk_michaelis_menten
        res = simulate_pk_michaelis_menten(
            vmax=100.0, km=5.0, v=50.0,
            doses=[{"time": 0.0, "amount": 500.0}],
            t0=0.0, t1=48.0,
            saveat=list(range(49))
        )
        assert "conc" in res["observations"]
        # Non-linear kinetics: check curve exists
        assert max(res["observations"]["conc"]) > 0


class TestPKCustom:
    """Tests for custom PK models."""

    def test_tmdd(self):
        from neopkpd import simulate_pk_tmdd_custom
        res = simulate_pk_tmdd_custom(
            kel=0.1, kon=0.01, koff=0.001, ksyn=1.0, kdeg=0.1, kint=0.05, v=50.0,
            doses=[{"time": 0.0, "amount": 500.0}],
            t0=0.0, t1=72.0,
            saveat=list(range(73))
        )
        assert "conc" in res["observations"]

    def test_parallel_absorption(self):
        from neopkpd import simulate_pk_parallel_absorption
        res = simulate_pk_parallel_absorption(
            ka1=2.0, ka2=0.5, f1=0.6, cl=10.0, v=50.0,
            doses=[{"time": 0.0, "amount": 500.0}],
            t0=0.0, t1=24.0,
            saveat=list(range(25))
        )
        assert "conc" in res["observations"]

    def test_enterohepatic_recirculation(self):
        from neopkpd import simulate_pk_enterohepatic_recirculation
        res = simulate_pk_enterohepatic_recirculation(
            ka=1.0, cl=10.0, v=50.0, kbile=0.5, kreab=0.3, f_reab=0.7,
            doses=[{"time": 0.0, "amount": 500.0}],
            t0=0.0, t1=48.0,
            saveat=list(range(49))
        )
        assert "conc" in res["observations"]

    def test_autoinduction(self):
        from neopkpd import simulate_pk_autoinduction
        res = simulate_pk_autoinduction(
            cl0=10.0, v=50.0, emax=2.0, ec50=5.0, kenz=0.1,
            doses=[{"time": i * 24.0, "amount": 200.0} for i in range(7)],  # QD for 7 days
            t0=0.0, t1=168.0,
            saveat=list(range(169))
        )
        assert "conc" in res["observations"]


# =============================================================================
# PD Model Tests
# =============================================================================

class TestPDBasic:
    """Tests for basic PD effect models."""

    def test_direct_emax(self):
        from neopkpd import simulate_pkpd_direct_emax
        res = simulate_pkpd_direct_emax(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            e0=0.0, emax=100.0, ec50=5.0,
            t0=0.0, t1=24.0,
            saveat=list(range(25))
        )
        assert "conc" in res["observations"]
        assert "effect" in res["observations"]
        # Effect should be between E0 and E0+Emax
        assert min(res["observations"]["effect"]) >= -1
        assert max(res["observations"]["effect"]) <= 101

    def test_sigmoid_emax(self):
        from neopkpd import simulate_pkpd_sigmoid_emax
        res = simulate_pkpd_sigmoid_emax(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            e0=0.0, emax=100.0, ec50=5.0, gamma=2.0,
            t0=0.0, t1=24.0,
            saveat=list(range(25))
        )
        assert "effect" in res["observations"]

    def test_biophase_equilibration(self):
        from neopkpd import simulate_pkpd_biophase_equilibration
        res = simulate_pkpd_biophase_equilibration(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            ke0=0.5, e0=0.0, emax=100.0, ec50=5.0,
            t0=0.0, t1=24.0,
            saveat=list(range(25))
        )
        assert "effect" in res["observations"]


class TestPDIndirectResponse:
    """Tests for Indirect Response Models (IRM)."""

    def test_irm3_inhibition_kout(self):
        """IRM-III: Inhibition of Kout (original indirect_response)"""
        from neopkpd import simulate_pkpd_indirect_response
        res = simulate_pkpd_indirect_response(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            kin=10.0, kout=0.1, r0=100.0, imax=0.9, ic50=5.0,
            t0=0.0, t1=120.0,
            saveat=list(range(121))
        )
        assert "response" in res["observations"]
        # Response should increase above baseline when Kout is inhibited
        assert max(res["observations"]["response"]) > 100

    def test_irm1_inhibition_kin(self):
        """IRM-I: Inhibition of Kin (production)"""
        from neopkpd import simulate_pkpd_irm1
        res = simulate_pkpd_irm1(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            kin=10.0, kout=0.1, r0=100.0, imax=0.9, ic50=5.0,
            t0=0.0, t1=120.0,
            saveat=list(range(121))
        )
        assert "response" in res["observations"]
        # Response should decrease below baseline when Kin is inhibited
        assert min(res["observations"]["response"]) < 100

    def test_irm2_stimulation_kin(self):
        """IRM-II: Stimulation of Kin (production)"""
        from neopkpd import simulate_pkpd_irm2
        res = simulate_pkpd_irm2(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            kin=10.0, kout=0.1, r0=100.0, smax=2.0, sc50=5.0,
            t0=0.0, t1=120.0,
            saveat=list(range(121))
        )
        assert "response" in res["observations"]
        # Response should increase above baseline when Kin is stimulated
        assert max(res["observations"]["response"]) > 100

    def test_irm4_stimulation_kout(self):
        """IRM-IV: Stimulation of Kout (elimination)"""
        from neopkpd import simulate_pkpd_irm4
        res = simulate_pkpd_irm4(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            kin=10.0, kout=0.1, r0=100.0, smax=2.0, sc50=5.0,
            t0=0.0, t1=120.0,
            saveat=list(range(121))
        )
        assert "response" in res["observations"]
        # Response should decrease below baseline when Kout is stimulated
        assert min(res["observations"]["response"]) < 100

    def test_irm_with_oral_pk(self):
        """Test IRM with oral PK model"""
        from neopkpd import simulate_pkpd_irm1
        res = simulate_pkpd_irm1(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            kin=10.0, kout=0.1, r0=100.0, imax=0.9, ic50=5.0,
            t0=0.0, t1=120.0,
            saveat=list(range(121)),
            pk_kind="OneCompOralFirstOrder",
            ka=1.0
        )
        assert "response" in res["observations"]


class TestPDAdvanced:
    """Tests for advanced PD models."""

    def test_transit_compartment_pd(self):
        """Transit compartment PD with signal transduction delay"""
        from neopkpd import simulate_pkpd_transit_compartment
        res = simulate_pkpd_transit_compartment(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            n_transit=5, ktr=0.5, e0=100.0, emax=-50.0, ec50=5.0, gamma=1.0,
            t0=0.0, t1=168.0,
            saveat=list(range(169))
        )
        assert "effect" in res["observations"]
        # Effect should be delayed compared to concentration
        # Find time of max concentration and min effect
        conc = res["observations"]["conc"]
        effect = res["observations"]["effect"]
        t_max_conc = conc.index(max(conc))
        t_min_effect = effect.index(min(effect))
        # Min effect should occur after max conc due to transit delay
        assert t_min_effect > t_max_conc

    def test_disease_progression_gompertz(self):
        """Disease progression with Gompertz growth"""
        from neopkpd import simulate_pkpd_disease_progression
        res = simulate_pkpd_disease_progression(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            growth_model="gompertz",
            s0=100.0, kgrow=0.05, smax=1000.0, alpha=0.0, kdrug=0.01,
            t0=0.0, t1=168.0,
            saveat=list(range(169))
        )
        assert "tumor_size" in res["observations"]

    def test_disease_progression_logistic(self):
        """Disease progression with Logistic growth"""
        from neopkpd import simulate_pkpd_disease_progression
        res = simulate_pkpd_disease_progression(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            growth_model="logistic",
            s0=100.0, kgrow=0.05, smax=1000.0, alpha=0.0, kdrug=0.01,
            t0=0.0, t1=168.0,
            saveat=list(range(169))
        )
        assert "tumor_size" in res["observations"]

    def test_disease_progression_exponential(self):
        """Disease progression with Exponential growth"""
        from neopkpd import simulate_pkpd_disease_progression
        res = simulate_pkpd_disease_progression(
            cl=1.0, v=10.0,
            doses=[{"time": 0.0, "amount": 100.0}],
            growth_model="exponential",
            s0=100.0, kgrow=0.05, smax=1000.0, alpha=0.0, kdrug=0.01,
            t0=0.0, t1=168.0,
            saveat=list(range(169))
        )
        assert "tumor_size" in res["observations"]


class TestPDTolerance:
    """Tests for tolerance/sensitization models."""

    def test_tolerance_counter_regulation(self):
        """Tolerance with counter-regulatory feedback"""
        from neopkpd import simulate_pkpd_tolerance_counter_regulation
        res = simulate_pkpd_tolerance_counter_regulation(
            cl=1.0, v=10.0,
            doses=[{"time": i * 8.0, "amount": 50.0} for i in range(21)],  # TID for 7 days
            e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
            kin_mod=0.1, kout_mod=0.05, alpha_feedback=1.0,
            t0=0.0, t1=168.0,
            saveat=[i * 1.0 for i in range(169)]
        )
        assert "effect" in res["observations"]
        # With tolerance, later doses should have smaller effect than earlier doses
        # This is hard to test directly, so just check it runs

    def test_receptor_down_regulation(self):
        """Receptor down-regulation tolerance"""
        from neopkpd import simulate_pkpd_receptor_regulation
        res = simulate_pkpd_receptor_regulation(
            cl=1.0, v=10.0,
            doses=[{"time": i * 8.0, "amount": 50.0} for i in range(21)],
            e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
            r_baseline=1.0, kreg=0.1, rmax=2.0, kchange=0.05,
            direction="down",
            t0=0.0, t1=168.0,
            saveat=[i * 1.0 for i in range(169)]
        )
        assert "effect" in res["observations"]

    def test_receptor_up_regulation(self):
        """Receptor up-regulation sensitization"""
        from neopkpd import simulate_pkpd_receptor_regulation
        res = simulate_pkpd_receptor_regulation(
            cl=1.0, v=10.0,
            doses=[{"time": i * 8.0, "amount": 50.0} for i in range(21)],
            e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
            r_baseline=1.0, kreg=0.1, rmax=2.0, kchange=0.05,
            direction="up",
            t0=0.0, t1=168.0,
            saveat=[i * 1.0 for i in range(169)]
        )
        assert "effect" in res["observations"]


# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

class TestErrorHandling:
    """Tests for error handling."""

    def test_invalid_pk_kind(self):
        from neopkpd import simulate_pkpd_irm1
        with pytest.raises(ValueError, match="Unsupported pk_kind"):
            simulate_pkpd_irm1(
                cl=1.0, v=10.0,
                doses=[{"time": 0.0, "amount": 100.0}],
                kin=10.0, kout=0.1, r0=100.0, imax=0.9, ic50=5.0,
                t0=0.0, t1=24.0,
                saveat=list(range(25)),
                pk_kind="InvalidModel"
            )

    def test_missing_ka_for_oral(self):
        from neopkpd import simulate_pkpd_direct_emax
        with pytest.raises(ValueError, match="ka required"):
            simulate_pkpd_direct_emax(
                cl=1.0, v=10.0,
                doses=[{"time": 0.0, "amount": 100.0}],
                e0=0.0, emax=100.0, ec50=5.0,
                t0=0.0, t1=24.0,
                saveat=list(range(25)),
                pk_kind="OneCompOralFirstOrder"
                # ka not provided
            )

    def test_invalid_growth_model(self):
        from neopkpd import simulate_pkpd_disease_progression
        with pytest.raises(ValueError, match="Unsupported growth_model"):
            simulate_pkpd_disease_progression(
                cl=1.0, v=10.0,
                doses=[{"time": 0.0, "amount": 100.0}],
                growth_model="invalid_model",
                s0=100.0, kgrow=0.05, smax=1000.0, alpha=0.0, kdrug=0.01,
                t0=0.0, t1=24.0,
                saveat=list(range(25))
            )

    def test_invalid_receptor_direction(self):
        from neopkpd import simulate_pkpd_receptor_regulation
        with pytest.raises(ValueError, match="direction must be"):
            simulate_pkpd_receptor_regulation(
                cl=1.0, v=10.0,
                doses=[{"time": 0.0, "amount": 100.0}],
                e0=0.0, emax=100.0, ec50=5.0, gamma=1.0,
                r_baseline=1.0, kreg=0.1, rmax=2.0, kchange=0.05,
                direction="invalid",
                t0=0.0, t1=24.0,
                saveat=list(range(25))
            )
