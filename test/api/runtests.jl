println("=============== test/api: correlation-function API tests ===============")
for f in ("imagtime", "realtime", "mixedtime", "bcs")
    println(">>> $f")
    include(joinpath(@__DIR__, "$f.jl"))
end
