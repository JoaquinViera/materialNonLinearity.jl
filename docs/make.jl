using Documenter, materialNonLinearity

makedocs(
  modules=[materialNonLinearity], 
  sitename="materialNonLinearity.jl",
  format = Documenter.HTML(),
  pages=[
        "Home" => "index.md"
        ]
)

deploydocs(
  repo="github.com/JoaquinViera/materialNonLinearity.jl.git",
  push_preview = true
  )
