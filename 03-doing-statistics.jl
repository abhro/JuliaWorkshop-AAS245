### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 16a16ccf-9ab2-4bb4-a6dd-f2a9991d81b4
using Dates, LinearAlgebra, Downloads

# ╔═╡ 7354f707-51cd-48d5-bdbb-efc2c33b061b
using Plots, PlotThemes

# ╔═╡ 0ad95c30-3f3d-48e5-86ef-7a2b816342e0
using StatsBase, Interpolations, SmoothingSplines, Distributions, LsqFit

# ╔═╡ dce8dbce-c8b4-11ed-3263-65232dc16f8d
begin
    using PlutoUI
    theme(:default) # :dark)
    TableOfContents()
end

# ╔═╡ 9273a0e2-6233-42c2-8419-3190a0b79a18
html"<button onclick='present()'>present</button>"

# ╔═╡ 3c4ee4fc-3daf-4570-889b-88c256ccbb1e
md"""
# Statistics of Supernova Remnants
"""

# ╔═╡ 754dbb34-631a-4aea-8660-443f70f11ea9
md"""
This tutorial is based on the paper, *A Statistical Analysis of Galactic Radio Supernova Remnants,* by S. Ranasinghe and D. Leahy, whom we acknowledge for communication that assisted our understanding.
"""

# ╔═╡ a7f8f920-5ee8-420a-bcb2-fcd30dde7d28

md"""
!!! tip

    The `PlotThemes` package allows for selection of a dark background for plots.


"""

# ╔═╡ 12b1a873-a22d-4fbc-bacb-f37b3767d2a2
md"""
If you desire a dark background, then check this box: $(@bind darkTH CheckBox(default = false))
"""

# ╔═╡ 94062cd4-b278-444d-b93e-694d26d30b50
theme(darkTH ? :dark : :default)

# ╔═╡ c8bacb96-8c31-4f7c-a42c-c94b5df0dd23
module MRTables

export MRTable, Field

struct Field
    range::UnitRange
    type::String
    unit::Union{String, Nothing}
    label::Union{String, Nothing}
    comment::Union{String, Nothing}
end

function Field(r::AbstractString, f::AbstractString, u::AbstractString, l::AbstractString, c::AbstractString)
    Field(range(parse.(Int, split(r, '-'))...),
          f,
          u == "---" ? nothing : u,
          l == "---" ? nothing : l,
          c == "---" ? nothing : c)
end

struct MRTable
    title::AbstractString
    authors::AbstractString
    table::AbstractString
    filename::AbstractString
    notes::Vector{Any}
    format::Vector{Field}
    data::Vector{Tuple}
end

const header_regex = [r"Title:\s*(?<article_title>\N*)\s*",
                      r"Authors:\s*(?<authors>\N*)\s*",
                      r"Table:\s*(?<table_title>\N*)\s*",
                      r"Description of file:\s*(?<filename>\N*)\s*"]
const note_regex   = r"Note\s+\((?<number>\d+)\):\s+(?<label>\N*)"
const field_regex  = r"\s*(?<range>\d+\s*-\s*\d+)\s+(?<format>\S+)\s+(?<unit>\S+)\s+(?<label>\S+)\s+(?<comment>\N+)"

const real_values  = Dict('I'=>Int, 'F'=>Float64)

function parse_data(field, line)
    if field.type[1] in keys(real_values) && length(strip(line[field.range])) > 0
        value = parse(real_values[field.type[1]], line[field.range])
    elseif length(strip(line[field.range])) == 0
        value = missing
    else
        value = strip(line[field.range])
    end
    value
end

function MRTable(filename::String)

    lines = eachline(filename)

    #   Read header lines
    header = Dict()
    for line in lines
        if match(r"^-+$", line) != nothing
            break
        end
        for regex in header_regex
            if (matched = match(regex, line)) != nothing
                header[keys(matched)[1]] = String(matched[1])
            end
        end
    end

    #   Read format lines
    format = []
    for line in lines
        if match(r"^-+$", line) != nothing && length(format) > 0
            break
        elseif (matched = match(field_regex, line)) != nothing
            push!(format, Field(matched.captures...))
        end
    end

    #   Read note lines
    notes = []
    note::Vector{Union{Tuple, String}} = []
    j = 0
    for line in lines
        if match(r"^-+$", line) != nothing
            break
        elseif (matched = match(note_regex, line)) != nothing
            j = parse(Int, (matched["number"]))
            note = [(j, matched["label"])]
            push!(notes, note)
        elseif j > 0
            push!(note, line)
        end
    end

    #   Read data lines
    data = []
    for line in lines
        if match(r"^-+$", line) != nothing
            break
        end
        push!(data, Tuple([parse_data(field, line) for field in format]))
    end

    MRTable(
        header["article_title"], header["authors"],
        header["table_title"], header["filename"],
        notes, format, data)
end

end

# ╔═╡ 01bbc40e-94be-4153-8a93-62f103121cae
md"""
## Read the machine readable table (MRT) file
"""

# ╔═╡ d10e1824-f773-41cf-9b69-5b296f245fd6
md"""
The data for this tutorial is in a standard "Machine Readable Table" from _The Astrophysical Journal_.
 - [Here](https://journals.aas.org/mrt-overview/) is the link to the standard.
 - [Here](https://iopscience.iop.org/article/10.3847/1538-4365/acc1de#apjsacc1det1) is the link to table.

However, the tutorial website already has a copy, which we will use.
"""

# ╔═╡ f07f26a4-54be-4561-8565-0aaace111969
allInfo = MRTables.MRTable(Downloads.download("https://github.com/barrettp/JuliaAASWorkshop/raw/main/apjsacc1det1_mrt.txt"))

# ╔═╡ 81265717-bcbf-4a70-a7a6-9f770f47ce42
md"""

We find $(length(allInfo.data)) rows of SNR data and clarifying notations.

!!! note
    Note that the module that reads the MRT data takes advantage of Julia's "missing" singleton, which is a recommended Julia standard. Where data are not provided, one need not worry about any particular form or special symbol(s), and in what follows, we shall benefit from this uniformity.
"""

# ╔═╡ bb8b0bd9-1b04-4ef0-b064-359b16ea59c9
md"""
## Producing (and then fitting) a histogram of SNR diameters.
"""

# ╔═╡ 84c802b9-69a6-4524-8e5d-766043fb98dc
md"""
Because we'll be looking at data in columns, to save a little typing we define a (one-line) function, `goDownColumn`, that does the marching down the provided column and also implements a `skipmissing`.

The `skipmissing` takes effect (typically with `collect`) later to produce arrays of the values.  Without `skipmissing`, Julia creates a "Union" of {Missing, Float64}, and the `missing`s would then be included in the final data array.

```julia
goDownCol(col::Int64) = skipmissing([row[col] for row in allInfo.data])
typeof(([row[colRadius] for row in allInfo.data]))
typeof(skipmissing([row[colRadius] for row = allInfo.data]))
```
"""

# ╔═╡ 77406430-6687-461e-bb4e-b6dabb873266


# ╔═╡ 1e7de711-6ddb-479b-8b97-cd42e0bb4421


# ╔═╡ cbaa2844-7054-4e8a-8d92-9706a6a49cef


# ╔═╡ 39044a6f-c1e4-48a9-988a-50a8194d0224
md"""
Columns we will be using throughout are defined here by hand, and we can confirm the correct number has been chosen by looking at `format`. The "comment" for `colRadius` lets us know that for additional information, there is a note, (4), to read.

```julia
begin
    colRadius = 41
    colDistanceLimit = 40
    colSNRname = 3
    colFluxDensity = 44
    colSNRType = 6
    colAge = 32
    colAgeLimit = 31
    xAxisLabel_Dhist = "SNR Diameter, D (pc)"
    yAxisLabelCount = "Number"
    yAxisLabelPDF = "Probability Density (pc\$\\textrm{^{-1}}\$)"
end
```
"""

# ╔═╡ 15439e8c-e40a-44fe-a881-f54eafcef952


# ╔═╡ e9629629-7870-4b7a-a457-cc4890245d99
md"""
Show the format for `colRadius`.

```julia
allInfo.format[colRadius]
```
"""

# ╔═╡ 4d375940-79b4-4764-8d06-10982b37c809


# ╔═╡ 62b0759b-d78d-4e03-95fc-0d87e7935857
md"""
The format for `colRadius` references Note 4. Let's take a look at it.

```julia
allInfo.notes[4]
```
"""

# ╔═╡ 2c9baa7f-1949-4ba7-bbf0-33ea49877b19


# ╔═╡ b48e97b0-8e8e-4710-87fe-7e3a7bee007d
md"""
As an aside, with `findall`, which yields the indices satisfying a given condition, I start the variable with `ind_` to remind myself that these are the indices of some variable, not the variable itself.

We first "go down" the column that contains the radii information, taking only the ones that are not missing (`!ismissing(x)`).  As noted above, the `goDownCol` function contains the `skipmissing` command.

```julia
ind_nonMissingRadii = findall(!ismissing, goDownCol(colRadius))
```
"""

# ╔═╡ 0618b8ae-17fd-432c-a84b-d079c9038658


# ╔═╡ 1948458d-31ba-46d2-915b-b3aba9dabf61
md"""
Collect the non-missing diameter indices; initially, we find 219.

We then need to exclude one specific SNR and because some distances are given as upper or lower limits, another 4.
"""

# ╔═╡ 3cfd5406-a4e3-4ecf-86ea-e30e67e18bf2
md"""
We need the index to the problematic SNR that is too close to the Galactic Center (G1.4 -0.1).

```julia
ind_G1p4M0p1 = findall(x -> x == "G1.4 -0.1", goDownCol(colSNRname))
```
"""

# ╔═╡ 48d4e255-df17-40ab-9edb-a8a770a116bb


# ╔═╡ f3822529-62fa-4275-bd0e-ea410b22fa59
md"""
By looking at the two variables above, we can see immediately that the index for that SNR, index 5, is already not included in the "nonMissingRadii," but in a more general case, we can check  with `in(item, collection)`, and the "`false`" response confirms its absence.

```julia
in(ind_G1p4M0p1, ind_nonMissingRadii)
```
"""

# ╔═╡ 619b64b4-fe02-4954-a02e-d7f8d979546d


# ╔═╡ 2e3cf549-0616-4809-ba38-73f7c6a94c34
md"""
Now we find the indices of the SNRs to be excluded because their distances are given only as limits, and then we use `symdiff` (if you will, the opposiste of `intersect`: "set of non-overlapping items between two sets") to remove them from the array of desired indices.

```julia
ind_UpperOrLowerLimits_OnDistance = findall(in([">", "<"]), goDownCol(colDistanceLimit))
```
"""

# ╔═╡ bde10148-2a5f-4597-8e81-f54d71f05728


# ╔═╡ 1b13dc19-adcb-493b-a571-aa2811f79398
md"""
Get indices of good radii.

```julia
ind_goodRadii = symdiff(ind_nonMissingRadii, ind_UpperOrLowerLimits_OnDistance)
```
"""

# ╔═╡ e11428ed-8814-417c-8746-6d98b14d28f8


# ╔═╡ 1734423f-9bf5-4d0e-b135-c550405f8783
md"""
We go down the column again, but this time selecting only the indices that we have whittled down to the good ones, and we multiply those radii by 2 to get the SNR diameters. The yield is 215 SNRs.

```julia
diameters = 2*goDownCol(colRadius)[ind_goodRadii]
nDiam = length(diameters)
```
"""

# ╔═╡ 06e9ea84-bb8e-4c29-93d2-9e8db8fa7dd3


# ╔═╡ 41d6c20f-eb56-4b63-ba30-39cd4c1e070b
md"""
## Histogram of SNR diameters
"""

# ╔═╡ c7939c4f-c446-4099-8f3b-9ab3f410129f
md"""
We will begin by fitting a histogram to the diameters. This unadorned statistical function performs simple counting with automatic binning into bins of 20 pc. The bins' `edges` are given, as are the resulting `weights` (counts) in each bin. We see, too, that the left side of each bin is included but the right side is not. (If there was a diameter exactly at the right edge, an additional bin would have been needed.)

```julia
yFit = fit(Histogram, diameters)
```
"""

# ╔═╡ aded4bdf-3511-49ee-9072-3f6f7f55e8a9


# ╔═╡ 7da41643-f486-4e2f-9ee2-06af0d7395d4
md"""
A check on the sum of the weights shows we have not missed any SNRs.

```julia
sum(yFit.weights)
```
"""

# ╔═╡ 2ec1a630-cde0-4d8f-b9b5-170fa1d122c6


# ╔═╡ 14a44aba-00c1-4b0a-abfd-cfd04147d2e2
md"""
Plot the histogram.

```julia
plot(yFit, label = "automatic binning", xlabel = xAxisLabel_Dhist, ylabel = yAxisLabelCount)
```
"""

# ╔═╡ 029336a2-e4db-47f9-8fff-6d751bcf4c4d


# ╔═╡ b54e32a0-b0b7-4fae-826d-2c61cf85dcc1
md"""
Let's change the bin size to 10 pc and create `yFit2`.

We retain the first and last edges from the automatic procedure.

```julia
begin
    binSize = 10
    firstEdge = first(yFit.edges[1]) # is zero
    lastEdge =  last(yFit.edges[1])  # is 180
    edgeRange = firstEdge:binSize:lastEdge
end
```

```julia
yFit2 = fit(Histogram, diameters, edgeRange)
```
"""

# ╔═╡ 53fcf13c-c5d1-425e-bc1a-241aaea7f9c5


# ╔═╡ bc29b8ce-4d25-4faa-b53a-17d2c8f81705


# ╔═╡ 0c2ed5dd-1085-4ee7-8969-2961b319978d
md"""
Replot histogram.

```julia
plot(yFit2,
     label = "bin size specified to 10 pc",
     xlabel = xAxisLabel_Dhist, ylabel = yAxisLabelCount)
```
"""

# ╔═╡ 940ee3ff-0e76-4806-8345-d3acd351fbff


# ╔═╡ 22a4c491-37f0-4e10-82e1-a9512ab1ef17
md"""
Now change the y axis to a probability density function: `isdensity` becomes `true` when we "normalize" the fit: `h = normalize(yFit2, mode= :pdf)`.

Make the x axis more informative.

```julia
h = normalize(yFit2, mode= :pdf)
```
"""

# ╔═╡ edd772bd-1152-4aca-ac19-58df089fae40


# ╔═╡ 4ec97f31-089e-4e89-bdd3-d8f03892ad31


# ╔═╡ 3432a558-5db8-4b65-8c4b-613586d8bb5e
md"""
## Cubic Spline interpolation
"""

# ╔═╡ 23bebd3d-c744-44ad-b01d-bff674a36b07
md"""
Can do a cubic spline interpolation on the pdf and then turn it into a function.

First, the interpolation, which goes through every point.

```julia
cubicInterp = Interpolations.interpolate(h.weights, BSpline(Cubic(Line(OnCell()))))
```
"""

# ╔═╡ 754a4149-8d0e-40bf-9f9e-a18b193dfe71


# ╔═╡ a3092e77-6888-4d58-aca1-37944a4896fb
md"""
Plot spline interpolation.
```julia
begin
    scatter(yFit2.edges[1][1:18, 1], h.weights,
            xlabel = xAxisLabel_Dhist, ylabel = yAxisLabelPDF,
            xlims = (firstEdge, lastEdge),
            xticks = (firstEdge:(2*binSize):lastEdge), xminorticks = 2,
            widen = true, label = "")

   plot!(yFit2.edges[1][1:18, 1], cubicInterp,
         label = "", title = "Interpolation")
end
```
"""

# ╔═╡ 10476ffe-2b37-4e15-9ef8-543c5c146e77


# ╔═╡ c1c15ecc-1f72-402c-bcc5-d1adc4814f38
md"""
The "scaling" turns the interpolation into a function that can be evaluated at intermediate values to prodce a smooth line.

Because the orignal has 18 points (18 bins, 19 edges), the scaling function requires 18; hence we start it at 1 (second edge divided by binSize:

```julia
yFit2.edges[1][2])/binSize = $(yFit2.edges[1][2]/binSize))
```
and not at `firstEdge`.

```julia
cubicScaled = Interpolations.scale(cubicInterp, 1:binSize:lastEdge)
```

```julia
typeof(cubicScaled)
```
"""

# ╔═╡ 635ec2cd-a3b3-4629-ac6d-9fc4853a3749


# ╔═╡ 4dda961e-a2ff-403f-9087-48254dac1d77


# ╔═╡ 2409c827-9d9d-44b4-91fb-a1516dcdc632
md"""
As an example of how the new variable can be evaluated at any (internal) point, `cubicScaled(166)`.

`cubicScaled(176)` fails because `1:binSize:lastEdge` corresponds to the 18 points 1.0:10.0:171.0.

For plotting, we could create a new variable, `cubicFiner = [cubicScaled(i) for i=0:176]`, with values at about ten times as many points that could then be plotted.

We would have to be careful about the x variable, which corresponds to the index, which starts at 1 and would correspond to the first index of `cubicFiner` unless explicitly directed.

Better, we can simply plot the scaled variable function directly, and the limits being plotted are clear: `plot(cubicScaled, lowerLimit, upperLimit)`.

```julia
cubicFiner = [cubicScaled(i) for i=0:176]
```
"""

# ╔═╡ 7c5fb45c-a687-4ed2-8c5f-e93f1d48fac7


# ╔═╡ b3607a9a-2472-4e15-8e90-d2da6675949f
md"""
Plot both the histogram and the spline on the same plot.

```julia
begin
    plot(h, label = "",
         xlabel = xAxisLabel_Dhist, ylabel = yAxisLabelPDF,
         xlims = (firstEdge, lastEdge),
         xticks = (firstEdge:2*binSize:lastEdge), xminorticks = 2,
         widen = true)

    plot!(cubicScaled, 0, 176,
          linewidth=2, color = :black, label = "interpolated")
end
```
"""

# ╔═╡ a5bc399a-3f49-4e17-a136-df02f5bd7e4a


# ╔═╡ eeba065c-942e-46f0-a3eb-d89e90bc514a
md"""
However, suppose instead we want a smoothing (as in the figure in the paper) that does not hit the points exactly but, instead, well, smooths them.

Again we need 18 points on the x axis, but now we want to start at zero, so we can use:

```julia
yFit2.edges[1][1:18, 1]
```

or
```julia
firstEdge:binSize:lastEdge - 1
```
"""

# ╔═╡ afe53e77-d73b-439f-8746-dc27c399dc59


# ╔═╡ ae34a53f-da42-42ac-b283-eefa09f3c86d
md"""
```julia
smooth = fit(SmoothingSpline, yFit2.edges[1][1:18, 1], h.weights, 250.0) # λ=250.0
```

!!! note

     This smoothing function is based on:

     - Reinsch, Christian H. "Smoothing by spline functions." Numerische mathematik 10.3 (1967): 177-183.

     - Green, Peter J., and Bernard W. Silverman. Nonparametric regression and generalized linear models: a roughness penalty approach. CRC Press, 1993.
"""

# ╔═╡ e2f76f80-4f4a-4304-ad6f-fea7b4042019


# ╔═╡ cdfc3aeb-169d-438f-8482-a78282b6c674
md"""
And again we have the choice between plotting a new variable with a finer distribution or plotting the function directly, but with this package, the transition to a functional form of `smooth` is:

```julia
smoothFiner = [predict(smooth, i) for i = 0:176.]
```

The `176.` ensures that we are using floating point numbers and not integers.
"""

# ╔═╡ 54fb1116-4887-4690-8c14-e70e206ce033


# ╔═╡ 539d14a0-627e-482a-91c7-27df784ca323
md"""
We now add the means and the standard error to the plot as vertical lines with `vline!`.

    begin
        ArithMean = mean(diameters)
        GeoMean = geomean(diameters)
        DStdError = std(diameters)/sqrt(length(diameters))
    end
"""

# ╔═╡ 09f7e727-c45a-424c-8a12-62c952a85814


# ╔═╡ 758a6437-156f-4b95-b522-341950203175
md"""
Recreate histogram-spline plot marking the means.

    begin
        annGeoMean = (65, 0.025, text("$(round(GeoMean, digits = 1)) = geometric mean", :red, 10))
        annArithMean = (72, 0.022, text("$(round(ArithMean, digits = 1)) \$\\pm\$ $(round(DStdError, digits = 1)) = arithmetic mean", :red, 10))
        plot(h, xlabel = xAxisLabel_Dhist, ylabel = yAxisLabelPDF, label = "", xlims = (firstEdge, lastEdge), xticks = (firstEdge:2*binSize:lastEdge), xminorticks = 2, widen = true)
       plot!(x -> cubicScaled(x), 0, 176, linewidth=2, color = :black, label = "interpolated")
       plot!(x -> predict(smooth, x), 0, 176, linewidth = 3, color = :orange, label = "smoothed")
       vline!([ArithMean], linewidth = 3, label = "", ann = annArithMean)
       vline!([GeoMean], linewidth = 3, label = "", ann = annGeoMean)
    end
"""

# ╔═╡ 66cc97b4-77f7-4cb3-8c7d-fc58990e6cc3


# ╔═╡ a9861db2-4382-4ed1-a94a-34c14c4d1cfc
md"""
## Empirical Cumulative Distribution Function
"""

# ╔═╡ 7ae4956c-b5e0-41fa-8b34-2e71e3cfd6f5
md"""
With standard functional distributions, e.g., Normal, use the Distributions package. Here, the function `ecdf` could be used, resulting in a type "ECDF{Vector{Float64}, Weights{Float64, Float64, Vector{Float64}}}"

Instead, we need only two lines to contruct the CDF manually (and the fitting that then can follow is more direct). We create a range from 1 to `nDiam` and broadcast a division by `nDiam`.  Then we sort the diameters and make the plot.

    begin
        empiricalCDF = (1:nDiam)./nDiam
        sortedD = sort(diameters)
    end
"""

# ╔═╡ 2f902640-65c2-4e97-8cca-ebc58c3bd5b5


# ╔═╡ 417dd70f-6895-4117-b7d5-8a3159a2c1fc
md"""
Plot the cumulative distribution.

    scatter(sortedD, empiricalCDF, xlabel = xAxisLabel_Dhist, ylabel = "Fraction of SNRs", title = "Empirical Cumuluative Distribution", label = "", xlims = (firstEdge, lastEdge), xticks = (firstEdge:2*binSize:lastEdge), xminorticks = 2, widen = true)
"""

# ╔═╡ 08c95754-1919-48d5-9610-5e9667864d91


# ╔═╡ f6938470-42a8-471d-99a2-fe22d81aefa7
md"""
We can also do an interpolated ECDF.

    interpECDF = Interpolations.interpolate(empiricalCDF, BSpline(Cubic(Line(OnCell()))))
"""

# ╔═╡ 1362c5e8-a8df-4daa-ae02-3b9c3e779e80


# ╔═╡ 1209bcc5-1fc4-4c31-8031-836c191f2e15
md"""
As before, an interpolated version would allow determination at intermediate values, but scaling would not be needed because the interprolation is already scaled to the proper x-axis limits.

    scatter(sortedD, interpECDF)

would reproduce the above plot.
"""

# ╔═╡ e77886fa-4145-402e-8674-ec0146608073


# ╔═╡ 010944f3-d959-4c7e-94ea-a645a170d0d3
md"""
## The SNR Σ -- D Relation
"""

# ╔═╡ a58f3563-ec5d-459d-89ea-892a964d53a6
md"""
All the SNRs with diameter estimations do not  have 1 GHz flux densities. A further subset of those that remain  are shell-type, and those are the ones we wish to plot.

Again we use `findall` and `intersect` to find the indices of the ones with flux densities that are shell type, plotting those with "good Sigma."  Additionally, we will use a best-fit line to Σ--D and a color bar of the ages of the remnants.

```julia
begin
    ind_FluxDensities = findall(x -> !ismissing(x), goDownCol(colFluxDensity))
    ind_shellType = findall(x -> x=="S", goDownCol(colSNRType))
    ind_RadiiWFluxDen = intersect(ind_goodRadii, ind_FluxDensities)
    ind_goodSigma = intersect(ind_shellType, ind_RadiiWFluxDen)
    length(ind_RadiiWFluxDen) # we get 185; the paper says 187
    diametersΣShells = 2*goDownCol(colRadius)[ind_goodSigma]
    ΣShells = goDownCol(colFluxDensity)[ind_goodSigma]
end
```
"""

# ╔═╡ 21a191f5-4b90-48e3-ab46-57352552040b


# ╔═╡ f709e059-2ff0-4320-ba01-038d9d643945
md"""
### Age in kilo years
"""

# ╔═╡ b8725478-19c6-4c62-9c8f-274cf768bfc9
md"""
Here we're ignoring scientific niceties such as limits.

And we note the function `extrema`, which delivers the minimum and maximum of a variable at one shot.

```julia
ind_AgeLimits = findall(in([">", "<"]), goDownCol(colAgeLimit))
```
"""

# ╔═╡ 575922c3-2007-45ad-9ced-a68f7c56f789


# ╔═╡ 6ca08891-0490-489c-8d8c-0dcccd399794
md"""
Find all ages of SNR in the table.

```julia
begin
    ind_MissingAge = findall(ismissing, [row[colAge] for row in allInfo.data])
    ind_goodSigmaMissingAge = intersect(ind_MissingAge, ind_goodSigma)
    ind_goodSigmaWithAge = symdiff(ind_goodSigma, ind_goodSigmaMissingAge)
    diametersToPlot = 2*goDownCol(colRadius)[ind_goodSigmaWithAge]
    ΣShellsToPlot = goDownCol(colFluxDensity)[ind_goodSigmaWithAge]
    agesToPlot = goDownCol(colAge)[ind_goodSigmaWithAge]
    extrema(agesToPlot)
end
```
"""

# ╔═╡ f2b86512-40a3-4f99-a0c3-937cfe684a14


# ╔═╡ 14daec81-e47e-43d7-9a5c-1f5ed5844eee
md"""
The extremes let us know immediately that we'll need to use the logarithms of the ages also.

```julia
begin
    z = agesToPlot
    ticks = 10.0 .^ collect(range(-2, length = 7))
    Ylabel = "\$\\textrm{\\Sigma_{1GHz} (10^{-21} W m^{-2} Hz^{-1} Sr^{-1}}\$)"

    scatter(diametersToPlot, ΣShellsToPlot, zcolor = log10.(z),
            yaxis=:log, xaxis=:log,
            label="", xlabel=xAxisLabel_Dhist, ylabel=Ylabel,
            xlims=(0.1, 200), xticks=ticks, xminorticks=9,
            ylims=(0.01, 10.0^(4.5)), yticks=ticks, yminorticks=9,
            widen=true, colorbar_title="Log age")
end
```
"""

# ╔═╡ c5615313-31ac-4993-8e88-987b4b6e0851


# ╔═╡ a239a50b-db01-4342-9c8d-821e63290a7a
md"""
We do a best fit on Σ--D and here plot the logarithms directly (instead of letting the plotting function take them).

```julia
begin
    XX, YY = log10.(diametersToPlot), log10.(ΣShellsToPlot)
    minXX, maxXX = extrema(XX)[1], extrema(XX)[2]
end
```
"""

# ╔═╡ ade1768e-c5ee-4cd5-ad9f-19948265b9c3


# ╔═╡ 9ae6b1ec-b446-4652-a1e1-f9629c47f1c3
md"""
The straight-line model to be fit is here.

```julia
begin
    @. model(t, p) = p[1] * t + p[2]
    p0 = [-1.5, 10]
end
```
"""

# ╔═╡ c39c4144-de4a-407d-beda-2f481bdc8f49


# ╔═╡ 295d30f4-4a03-45d5-8dec-86d836da13ae
md"""
Use the Levenberg-Marquart least squares algorithm to fit a straight line through the data.

```julia
LSfit = LsqFit.curve_fit(model, XX, YY, p0)
```
"""

# ╔═╡ 263833e1-3105-4c56-a3dc-30ae6934e1f8


# ╔═╡ 877ccfff-de5d-46d5-aff8-0bd71d14e996
md"""
We use the results of the model to get the Σ values at the minimum and maximum diameters, which we then plot.

```julia
begin
    m = LSfit.param[1]
    b = LSfit.param[2]
    @. results(x) = m * x + b
    twoYs = results([minXX, maxXX])
end
```
"""

# ╔═╡ d7b126b1-f5b6-4e38-a40e-2d6918fda741


# ╔═╡ 2bd1f2be-9318-44e8-9791-c8d1451b4522
md"""
Create a scatter plot with best fit line.

```julia
begin
    scatter(XX, YY,
            zcolor=log10.(z),
            label="",
            xlabel="Log "*xAxisLabel_Dhist, ylabel="Log "*Ylabel,
            xminorticks=9, widen=true, colorbar_title="Log age ")

    scatter!([minXX, maxXX], twoYs, label="")

    plot!([minXX, maxXX], twoYs, label="")
end
```
"""

# ╔═╡ 35405245-31b4-400b-9c42-b957f606188c


# ╔═╡ 57e60f44-342a-4023-8011-cfc425908744
md"""
!!! note

    What have we learned?

    - With the newly created MRTs module, machine-readable ApJ tables can be read.
      - It made major use of structs to create composite types.
      - To disambiguate the text file, Regular Expressions operating on Strings were needed throughout.
      - The `eachline` command executed on the filename allowed for implicit iteration over all the file lines.
      - Perhaps most important, by replacing blank entries with the Julia standard, "missing," subsequent processing was standardized.
    - Histogram and Cumulative Distribution Creation
      - A one-line function to march down the MRT's columns used `skipmissing`, which accepts the missing values initially, but when "collected" into an array, they are not included.
      - To pull out the desired values, `findall` generates indices meeting specified conditions.
      - Then `intersect` and `symdiff` were used to produce arrays containing what was desired.
      - The SNR diameters were `fit` to a Histogram function, first with automatic binning, and then with a specified bin size.
      - Instead of counts, the Histogram output was changed to a density.
      - he Interpolations package was used to create a function that reproduced the points.
      - The SmoothingSplines packaged was used to produce a smooth function over the points.
      - Plotting the function could then be accomplished directly in the plotting command.
      - Two different means were generated and placed as vertical lines.
      - Broadcasting was used to write a one-line empirical CDF.
    - The Σ -- D Relationship
      - After further selection of diameters, densities, and ages, the LsqFit package was used to find a best fit.
      - Two different ways of plotting logarithms were used.

    - Finally, we reiterate that as long as someone has the appropriate MRT file to read (and Julia, of course!), this notebook will work. It can also be exported from Pluto as an html file and viewed without Julia.
"""

# ╔═╡ 8237e174-f6ee-4b4a-9928-dfe7b239b467
html"<button onclick='present()'>present</button>"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SmoothingSplines = "102930c3-cf33-599f-b3b1-9a29a5acab30"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
Distributions = "~0.25.109"
Interpolations = "~0.15.1"
LsqFit = "~0.15.0"
PlotThemes = "~3.2.0"
Plots = "~1.40.5"
PlutoUI = "~0.7.59"
SmoothingSplines = "~0.3.2"
StatsBase = "~0.34.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "e9304ff0616944237591131bfe8e9e541d43e779"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "f54c23a5d304fb87110de62bace7777d59088c34"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.15.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "b8fe8546d52ca154ac556809e10c75e6e7430ac8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.5"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d8a9c0b6ac2d9081bf76324b39c78ca3ce4f0c98"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.6"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9c405847cc7ecda2dc921ccf18b47ca150d7317e"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.109"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0653c0a2396a6da5bc4766c43041ef5fd3efbe57"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.11.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "2de436b72c3422940cbe1367611d137008af7ec3"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.23.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "5b0d630f3020b82c0775a51d05895852f8506f50"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.4"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "Printf", "StatsAPI"]
git-tree-sha1 = "40acc20cfb253cf061c1a2a2ea28de85235eeee1"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.15.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "082f0c4b70c202c37784ce4bfbc33c9f437685bf"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.5"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PtrArrays]]
git-tree-sha1 = "f011fbb92c4d401059b2212c05c0601b70f8b759"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "e237232771fdafbae3db5c31275303e056afaa9f"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.10.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d483cd324ce5cf5d61b77930f0bbd6cb61927d21"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.2+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SmoothingSplines]]
deps = ["LinearAlgebra", "Random", "Reexport", "StatsBase"]
git-tree-sha1 = "3a68e878003f7d6ea0be9e3bafcabfb79f5a70ee"
uuid = "102930c3-cf33-599f-b3b1-9a29a5acab30"
version = "0.3.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "96612ac5365777520c3c5396314c8cf7408f436a"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.1"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "d9717ce3518dc68a99e6b96300813760d887a01d"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.1+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─9273a0e2-6233-42c2-8419-3190a0b79a18
# ╟─16a16ccf-9ab2-4bb4-a6dd-f2a9991d81b4
# ╟─7354f707-51cd-48d5-bdbb-efc2c33b061b
# ╟─0ad95c30-3f3d-48e5-86ef-7a2b816342e0
# ╟─dce8dbce-c8b4-11ed-3263-65232dc16f8d
# ╟─3c4ee4fc-3daf-4570-889b-88c256ccbb1e
# ╟─754dbb34-631a-4aea-8660-443f70f11ea9
# ╟─a7f8f920-5ee8-420a-bcb2-fcd30dde7d28
# ╟─12b1a873-a22d-4fbc-bacb-f37b3767d2a2
# ╟─94062cd4-b278-444d-b93e-694d26d30b50
# ╟─c8bacb96-8c31-4f7c-a42c-c94b5df0dd23
# ╟─01bbc40e-94be-4153-8a93-62f103121cae
# ╟─d10e1824-f773-41cf-9b69-5b296f245fd6
# ╠═f07f26a4-54be-4561-8565-0aaace111969
# ╟─81265717-bcbf-4a70-a7a6-9f770f47ce42
# ╟─bb8b0bd9-1b04-4ef0-b064-359b16ea59c9
# ╟─84c802b9-69a6-4524-8e5d-766043fb98dc
# ╠═77406430-6687-461e-bb4e-b6dabb873266
# ╠═1e7de711-6ddb-479b-8b97-cd42e0bb4421
# ╠═cbaa2844-7054-4e8a-8d92-9706a6a49cef
# ╟─39044a6f-c1e4-48a9-988a-50a8194d0224
# ╠═15439e8c-e40a-44fe-a881-f54eafcef952
# ╟─e9629629-7870-4b7a-a457-cc4890245d99
# ╠═4d375940-79b4-4764-8d06-10982b37c809
# ╟─62b0759b-d78d-4e03-95fc-0d87e7935857
# ╠═2c9baa7f-1949-4ba7-bbf0-33ea49877b19
# ╟─b48e97b0-8e8e-4710-87fe-7e3a7bee007d
# ╠═0618b8ae-17fd-432c-a84b-d079c9038658
# ╟─1948458d-31ba-46d2-915b-b3aba9dabf61
# ╟─3cfd5406-a4e3-4ecf-86ea-e30e67e18bf2
# ╠═48d4e255-df17-40ab-9edb-a8a770a116bb
# ╟─f3822529-62fa-4275-bd0e-ea410b22fa59
# ╠═619b64b4-fe02-4954-a02e-d7f8d979546d
# ╟─2e3cf549-0616-4809-ba38-73f7c6a94c34
# ╠═bde10148-2a5f-4597-8e81-f54d71f05728
# ╟─1b13dc19-adcb-493b-a571-aa2811f79398
# ╠═e11428ed-8814-417c-8746-6d98b14d28f8
# ╟─1734423f-9bf5-4d0e-b135-c550405f8783
# ╠═06e9ea84-bb8e-4c29-93d2-9e8db8fa7dd3
# ╟─41d6c20f-eb56-4b63-ba30-39cd4c1e070b
# ╟─c7939c4f-c446-4099-8f3b-9ab3f410129f
# ╠═aded4bdf-3511-49ee-9072-3f6f7f55e8a9
# ╟─7da41643-f486-4e2f-9ee2-06af0d7395d4
# ╠═2ec1a630-cde0-4d8f-b9b5-170fa1d122c6
# ╟─14a44aba-00c1-4b0a-abfd-cfd04147d2e2
# ╠═029336a2-e4db-47f9-8fff-6d751bcf4c4d
# ╟─b54e32a0-b0b7-4fae-826d-2c61cf85dcc1
# ╠═53fcf13c-c5d1-425e-bc1a-241aaea7f9c5
# ╠═bc29b8ce-4d25-4faa-b53a-17d2c8f81705
# ╟─0c2ed5dd-1085-4ee7-8969-2961b319978d
# ╠═940ee3ff-0e76-4806-8345-d3acd351fbff
# ╟─22a4c491-37f0-4e10-82e1-a9512ab1ef17
# ╠═edd772bd-1152-4aca-ac19-58df089fae40
# ╠═4ec97f31-089e-4e89-bdd3-d8f03892ad31
# ╟─3432a558-5db8-4b65-8c4b-613586d8bb5e
# ╟─23bebd3d-c744-44ad-b01d-bff674a36b07
# ╠═754a4149-8d0e-40bf-9f9e-a18b193dfe71
# ╟─a3092e77-6888-4d58-aca1-37944a4896fb
# ╠═10476ffe-2b37-4e15-9ef8-543c5c146e77
# ╟─c1c15ecc-1f72-402c-bcc5-d1adc4814f38
# ╠═635ec2cd-a3b3-4629-ac6d-9fc4853a3749
# ╠═4dda961e-a2ff-403f-9087-48254dac1d77
# ╟─2409c827-9d9d-44b4-91fb-a1516dcdc632
# ╠═7c5fb45c-a687-4ed2-8c5f-e93f1d48fac7
# ╟─b3607a9a-2472-4e15-8e90-d2da6675949f
# ╠═a5bc399a-3f49-4e17-a136-df02f5bd7e4a
# ╟─eeba065c-942e-46f0-a3eb-d89e90bc514a
# ╠═afe53e77-d73b-439f-8746-dc27c399dc59
# ╟─ae34a53f-da42-42ac-b283-eefa09f3c86d
# ╠═e2f76f80-4f4a-4304-ad6f-fea7b4042019
# ╟─cdfc3aeb-169d-438f-8482-a78282b6c674
# ╠═54fb1116-4887-4690-8c14-e70e206ce033
# ╟─539d14a0-627e-482a-91c7-27df784ca323
# ╠═09f7e727-c45a-424c-8a12-62c952a85814
# ╟─758a6437-156f-4b95-b522-341950203175
# ╠═66cc97b4-77f7-4cb3-8c7d-fc58990e6cc3
# ╟─a9861db2-4382-4ed1-a94a-34c14c4d1cfc
# ╟─7ae4956c-b5e0-41fa-8b34-2e71e3cfd6f5
# ╠═2f902640-65c2-4e97-8cca-ebc58c3bd5b5
# ╟─417dd70f-6895-4117-b7d5-8a3159a2c1fc
# ╠═08c95754-1919-48d5-9610-5e9667864d91
# ╟─f6938470-42a8-471d-99a2-fe22d81aefa7
# ╠═1362c5e8-a8df-4daa-ae02-3b9c3e779e80
# ╟─1209bcc5-1fc4-4c31-8031-836c191f2e15
# ╠═e77886fa-4145-402e-8674-ec0146608073
# ╟─010944f3-d959-4c7e-94ea-a645a170d0d3
# ╟─a58f3563-ec5d-459d-89ea-892a964d53a6
# ╠═21a191f5-4b90-48e3-ab46-57352552040b
# ╟─f709e059-2ff0-4320-ba01-038d9d643945
# ╟─b8725478-19c6-4c62-9c8f-274cf768bfc9
# ╠═575922c3-2007-45ad-9ced-a68f7c56f789
# ╟─6ca08891-0490-489c-8d8c-0dcccd399794
# ╠═f2b86512-40a3-4f99-a0c3-937cfe684a14
# ╠═14daec81-e47e-43d7-9a5c-1f5ed5844eee
# ╠═c5615313-31ac-4993-8e88-987b4b6e0851
# ╟─a239a50b-db01-4342-9c8d-821e63290a7a
# ╠═ade1768e-c5ee-4cd5-ad9f-19948265b9c3
# ╟─9ae6b1ec-b446-4652-a1e1-f9629c47f1c3
# ╠═c39c4144-de4a-407d-beda-2f481bdc8f49
# ╟─295d30f4-4a03-45d5-8dec-86d836da13ae
# ╠═263833e1-3105-4c56-a3dc-30ae6934e1f8
# ╟─877ccfff-de5d-46d5-aff8-0bd71d14e996
# ╠═d7b126b1-f5b6-4e38-a40e-2d6918fda741
# ╟─2bd1f2be-9318-44e8-9791-c8d1451b4522
# ╠═35405245-31b4-400b-9c42-b957f606188c
# ╟─57e60f44-342a-4023-8011-cfc425908744
# ╟─8237e174-f6ee-4b4a-9928-dfe7b239b467
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
