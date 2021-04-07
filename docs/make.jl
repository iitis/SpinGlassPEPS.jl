using Documenter, SpinGlassPEPS
using DocumenterTools: Themes
using SpinGlassTensors, SpinGlassNetworks, SpinGlassEngine
using MetaGraphs

cd(@__DIR__)
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using Pkg
Pkg.activate(@__DIR__)
CI && Pkg.instantiate()

# %%
# download the themes
#for file in ("juliadynamics-lightdefs.scss", "juliadynamics-darkdefs.scss", "juliadynamics-style.scss")
#    download("https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/$file", joinpath(@__DIR__, file))
#end
# create the themes
for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "juliadynamics-style.scss"), String)
    theme = read(joinpath(@__DIR__, "juliadynamics-$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "juliadynamics-$(w).scss"), header*"\n"*theme)
end
# compile the themes
Themes.compile(joinpath(@__DIR__, "juliadynamics-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "juliadynamics-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

format = Documenter.HTML(
    edit_link = "master",
    prettyurls = !("local" in ARGS),
    mathengine = mathengine = MathJax3(Dict(
        :tex=>Dict(
            "inlineMath"=>[ ["\$","\$"], ["\\(","\\)"] ],
            "processEscapes"=>true,
            "macros"=>Dict(
                "bra"=> ["{\\langle #1|}",1],
                "ket"=> ["{| #1\\rangle}",1],
                "ketbra"=> ["{\\left| #1 \\rangle \\langle #2 \\right|}",2],
                "braket"=> ["{\\langle #1|#2\\rangle}",2],
                "Tr"=> ["\\mathrm{Tr}",0],
                "tr"=> ["\\Tr",0],
                "ee"=> ["\\mathrm{e}"],
                "ii"=> ["\\mathrm{i}"],
                "dd"=> ["\\mathrm{d}"],
                "1"=> ["{\\mathbb{1}}"]
            )
        )
    ))
)

makedocs(
    clean = true,
    format = format,
    modules=[SpinGlassPEPS, SpinGlassTensors, SpinGlassNetworks, SpinGlassEngine],
    sitename = "SpinGlassPEPS.jl",
    authors = "Krzysztof Domino, Bartłomiej Gardas, Konrad Jałowiecki, Łukasz Pawela, Marek Rams, Anna Dziubyna",
    pages = [
        "Home" => "index.md",
        "Contents" => "contents.md",
        "Examples" => "examples.md",
        "Library" => ["Index" => "lib/SpinGlassPEPS.md",
                       "SpinGlassEngine" => "lib/SpinGlassEngine.md",
                       "SpinGlassTensors" => "lib/SpinGlassTensors.md",
                       "SpinGlassNetwork" => "lib/SpinGlassNetwork.md"]
    ]
)

if "deploy" in ARGS
    include("../../faketravis.jl")
end

deploydocs(
    #root = "github.com/iitis/SpinGlassPEPS.jl/tree/ad/docs",
    #repo = "github.com/iitis/SpinGlassPEPS.jl.git",
    repo = "github.com/iitis/SpinGlassPEPS.jl.git",
    branch = "ad/docs",
    devbranch = "ad/docs"
    #devbranch = "lp/docs-example",
    #branch = "ad/docs"
)