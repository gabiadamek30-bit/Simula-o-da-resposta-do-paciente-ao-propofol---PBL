# Test helper functions for NeoPKPD

using Test
using NeoPKPD

"""
Analytic solution for one-compartment IV bolus concentration.
"""
function analytic_onecomp_ivbolus_conc(
    t::Float64, doses::Vector{DoseEvent}, CL::Float64, V::Float64
)
    k = CL / V
    c = 0.0
    for d in doses
        if t >= d.time
            dt = t - d.time
            c += (d.amount / V) * exp(-k * dt)
        end
    end
    return c
end

"""
Analytic solution for one-compartment oral first-order concentration.
"""
function analytic_onecomp_oral_first_order_conc(
    t::Float64, doses, Ka::Float64, CL::Float64, V::Float64
)
    k = CL / V
    c = 0.0
    for d in doses
        if t >= d.time
            dt = t - d.time
            # Handle Ka close to k for numerical stability
            if abs(Ka - k) < 1e-12
                # Limit as Ka -> k:
                # C(t) = (Dose/V) * Ka * dt * exp(-k*dt)
                c += (d.amount / V) * Ka * dt * exp(-k * dt)
            else
                c += (d.amount / V) * (Ka / (Ka - k)) * (exp(-k * dt) - exp(-Ka * dt))
            end
        end
    end
    return c
end

"""
Direct Emax function for PD testing.
"""
function direct_emax(C::Float64, E0::Float64, Emax::Float64, EC50::Float64)
    return E0 + (Emax * C) / (EC50 + C)
end

"""
Analytic solution for turnover response (no drug effect).
"""
function analytic_turnover_R(t::Float64, Kin::Float64, Kout::Float64, R0::Float64)
    # dR/dt = Kin - Kout*R
    return R0 * exp(-Kout * t) + (Kin / Kout) * (1.0 - exp(-Kout * t))
end

"""
Sigmoid Emax (Hill equation) function.
"""
function sigmoid_emax_ref(C::Float64, E0::Float64, Emax::Float64, EC50::Float64, gamma::Float64)
    if C <= 0.0
        return E0
    end
    return E0 + (Emax * C^gamma) / (EC50^gamma + C^gamma)
end
