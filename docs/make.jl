using Documenter, AlternatingProjections

makedocs(sitename="Alternating Projections.jl",
    modules = [lternating Projections],
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://olejorik.github.io/AlternatingProjections.jl/stable/",
        # assets = ["assets/favicon.ico"],
        highlights = ["yaml"],
    ),
    clean = false,
    authors = "O.S.",)

# using DocumenterLaTeX
# makedocs(
# #     format = LaTeX(),
#     format = LaTeX(platform = "docker"),
#     sitename="Alternating Projections"
# )

deploydocs(repo = "github.com/olejorik/AlternatingProjections.jl.git",
    target = "build",)