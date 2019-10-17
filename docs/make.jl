using Documenter, AlternatingProjections

makedocs(sitename="AlternatingProjections.jl",
    modules = [AlternatingProjections],
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = true,
        canonical = "https://olejorik.github.io/AlternatingProjections.jl/stable/",
        # assets = ["assets/favicon.ico"],
        # highlights = ["yaml"],
    ),
    clean = false,
    authors = "O.S.",)

# using DocumenterLaTeX
# makedocs(
# #     format = LaTeX(),
#     format = LaTeX(platform = "docker"),
#     sitename="AlternatingProjections"
# )

deploydocs(repo = "github.com/olejorik/AlternatingProjections.jl.git",
    target = "build",)