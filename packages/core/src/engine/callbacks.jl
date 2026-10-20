export preset_dose_callback

function preset_dose_callback(
    doses::Vector{DoseEvent}, t0::Float64, t1::Float64, target_index::Int
)
    _, dose_times, dose_amounts = normalize_doses_for_sim(doses, t0, t1)

    if isempty(dose_times)
        return nothing
    end

    function affect!(integrator)
        idx = findfirst(==(integrator.t), dose_times)
        if idx === nothing
            error("Internal error: dose time not found for t=$(integrator.t)")
        end
        integrator.u[target_index] += dose_amounts[idx]
    end

    return PresetTimeCallback(dose_times, affect!)
end
