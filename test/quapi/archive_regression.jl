# 存档回归测试：performance/ 基准输出与本目录存档 CSV 逐字节一致。
#
#   一致 = 数值行为不变的最终防线：任何 branch 约定 / 构造器 / 内核改动只要
#   扰动一个字节，这里立刻暴露。CSV 更新需显式重新拷贝存档（有意识的行为变更）。
#
# 覆盖：bosons 的 firstorder/secondorder（实时/虚时/混合）。
using Test

const ROOT = joinpath(@__DIR__, "..", "..")
const PERF = joinpath(ROOT, "performance")
const ARCHIVES = (
    (dir = joinpath(@__DIR__, "archive_bosons_firstorder"),
     suite = joinpath(PERF, "firstorder", "bosons"), tag = "firstorder bosons"),
    (dir = joinpath(@__DIR__, "archive_bosons_secondorder"),
     suite = joinpath(PERF, "secondorder", "bosons"), tag = "secondorder bosons"),
)

function readcsv(path)
    lines = readlines(path)
    return lines[1], lines[2:end]
end

function run_and_compare(script, archive, live, name)
    @testset "$name" begin
        old = isfile(live) ? read(live, String) : nothing
        try
            run(`$(Base.julia_cmd()) --startup-file=no --project=$ROOT $script`)
            @test isfile(archive)
            @test isfile(live)
            ahead, arows = readcsv(archive)
            lhead, lrows = readcsv(live)
            @test lhead == ahead
            @test length(lrows) == length(arows)
            for (a, b) in zip(arows, lrows)
                @test a == b
            end
        finally
            if old === nothing
                rm(live; force=true)
            else
                write(live, old)
            end
        end
    end
end

@testset "archive regression (performance vs archived CSV)" begin
    for A in ARCHIVES
        for name in ("realtime", "imagtime", "mixedtime")
            run_and_compare(joinpath(A.suite, name * ".jl"),
                            joinpath(A.dir, name * ".csv"),
                            joinpath(A.suite, "results", name * ".csv"),
                            A.tag * " " * name)
        end
    end
end
