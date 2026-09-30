# The scripts under examples/ are documentation that runs. Each one is run
# here in a scratch copy of its directory, so an API change that breaks a
# documented recipe breaks the suite too.

const EXAMPLES = joinpath(dirname(dirname(@__DIR__)), "examples")

function run_example(name)
    mktempdir() do dir
        cp(joinpath(EXAMPLES, name), joinpath(dir, name))
        cd(joinpath(dir, name)) do
            redirect_stdout(devnull) do
                Base.include(Module(), abspath("run.jl"))
            end
            return readdir(".")
        end
    end
end

@testset "examples/ run" begin
    out = run_example("dada2")
    @test "out" in out

    out = run_example("qiime2")
    @test "lineages.tsv" in out

    out = run_example("epa-ng")
    @test "confident.jplace" in out
end
