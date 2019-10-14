using Documenter, AlternatingProjections

makedocs(sitename="Alternating Projections")

# using DocumenterLaTeX
# makedocs(
#     format = LaTeX(),
# #     format = LaTeX(platform = "docker"),
#     sitename="Alternating Projections"
# )

deploydocs(repo = "github.com/olejorik/AlternatingProjections.jl.git")