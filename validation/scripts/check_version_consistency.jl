using NeoPKPD

root = normpath(joinpath(@__DIR__, "..", ".."))
version_file = joinpath(root, "VERSION")

v = strip(read(version_file, String))

if v != NeoPKPD.NEOPKPD_VERSION
    error("VERSION file ($(v)) does not match NeoPKPD.NEOPKPD_VERSION ($(NeoPKPD.NEOPKPD_VERSION))")
end

println("Version consistency check passed: " * v)
