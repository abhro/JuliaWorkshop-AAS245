### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 3ba5672b-18ff-4320-a4d2-954e0b873d47
using PlutoUI; TableOfContents()

# ╔═╡ a833decc-ef58-46a7-923d-877204e716ff
#=╠═╡
using PyCall
  ╠═╡ =#

# ╔═╡ 6899990c-00b2-40e5-941c-df8f7073183c
md"""
# Calling Python from Julia
"""

# ╔═╡ 2ec6992f-38de-4ff5-a2f3-21f9643657e9
md"""
Julia has two packages for calling Python code from Julia: PyCall & PythonCall. PyCall has an easier learning curve, but has a few limitations. Whereas, PythonCall has better performance, because it is more complete and allows for more control.

Both packages have symmetric interfaces, so the user can call Python from Julia and Julia from Python. This session will only focus on calling Python from Julia.
"""

# ╔═╡ 14c3864e-b094-4c46-b38d-1abb375dc44a
md"""
## PyCall & PyJulia
"""

# ╔═╡ 41103ed5-6145-4772-84d6-8858ebb560ac


# ╔═╡ 306378d0-91a3-4f07-94fb-45b8cc474274
md"""
A simple example, calculate 1/√2 using Python. First, import the math module.

```julia
begin
    math = pyimport("math")
    math.sin(math.pi / 4)
end
```
Try it.
"""

# ╔═╡ 667e42f2-b6b5-4eb9-817f-57aee042d7ca


# ╔═╡ 06d63e41-f9b9-4c75-87d4-241beb454dc5
md"""
Here is another example using matplotlib.pyplot. The example uses the `pyimport_conda` function to load `matplotlib.pyplot`. If `matplotlib.pyplot` fails because `matplotlib` hasn't been installed, then it will automatically install `matplotlib`, which may take some time, or retry `pyimport`.

```julia
let
    plt = pyimport_conda("matplotlib.pyplot", "matplotlib")
    x = range(0; stop=2pi, length=1000)
    y = sin.(3x + 4 .* cos.(2x))
    plt.plot(x, y, color="red", linewidth=2.0, linestyle="--")
    plt.savefig("plt_example.png")
    LocalResource("plt_example.png")
end
```

!!! note
    `x` and `y` are calculated using Julia and the plotting uses Python.
"""

# ╔═╡ 20c226ca-44e6-4565-b601-07f5df7828ac


# ╔═╡ 7728f44a-93e1-4514-8964-6351a5ebff07
md"""
You can also just wrap Python code in `py"..."` or `py\"""...\"""` strings, which are equivalent to Python's `eval` and `exec` commands, respectively.

First install `astropy` using `pyimport_conda`.

```julia
pyimport_conda("astropy.io.fits", "astropy")
```
"""

# ╔═╡ 8648cd40-d805-49b3-b05a-93a4b6ea0a20


# ╔═╡ f790ac60-d874-4e95-8d34-252f44ff32b1
md"""
For example, try the following block of code to access the `time` and `rate` of some *Newton-XMM* X-ray data from a FITS file.

```julia
begin
    py\"""
    import astropy.io.fits as fits
    mos = fits.open("https://github.com/sefffal/AASJuliaWorkshop/blob/main/P0801780301M1S001SRCTSR8001.FIT?raw=true")
    time, rate = mos[1].data["TIME"], mos[1].data["RATE"]
    \"""
    (py"time"[1:3], py"rate"[1:3])
end
```
"""

# ╔═╡ 90c8b0db-f94b-440c-b05a-72e3f6483051


# ╔═╡ d8e3430e-5ad7-44c8-8a08-be067babac72
md"""
### Specifying the Python version
"""

# ╔═╡ 5d883d2b-515f-4286-a2ab-15127ac6b5ea
md"""
In Julia,

```julia
ENV["PYTHON"] = path_of_python_executable
# ENV["PYTHON"] = "/usr/bin/python3.10"
Pkg.build("PyCall")
```

To use Julia's Python distribution, set the path to an empty string, i.e., `""`

!!! note
    Usually, the necessary libraries are installed along with Python, but pyenv on  MacOS requires you to install it with env `PYTHON_CONFIGURE_OPTS="--enable-framework"` pyenv install 3.4.3. The Enthought Canopy Python distribution is currently not supported. As a general rule, we tend to recommend the Anaconda Python distribution on MacOS and Windows, or using the Julia Conda package, in order to minimize headaches.
"""

# ╔═╡ a836055d-9d11-46c7-bfca-a6759259e145
md"""
## PythonCall & JuliaCall
"""

# ╔═╡ a7a6b13e-cda7-4629-9b54-950759e6bab5
md"""
### Getting Started
"""

# ╔═╡ 09aa5c19-7d18-43b1-9623-1609f8792070
md"""
PyCall and PythonCall can be used in the same Julia session on Unix (Linux, OS X, etc.) as long the same interpreter is used for both. On Windows, it appears separate interpreters can be used. Let's ensure the same interpreter is used for both.

```julia
ENV["JULIA_PYTHONCALL_EXE"] = "@PyCall"
```
"""

# ╔═╡ 0db12101-fcf9-48a3-a337-ad53a4713f6d


# ╔═╡ 93e945de-7ceb-4124-b3db-f6b64e2e46bb
md"""
Now import `PythonCall`:

```julia
using PythonCall
```
"""

# ╔═╡ fa825783-6f48-4f0d-b346-220ae9e2fb11


# ╔═╡ 99e16357-0113-4a12-8a22-f627c6fdef84
md"""
By default importing the module will initialize a conda environment in your Julia environment, install Python into it, load the corresponding Python library, and initialize an interpreter.

Here is an example using Python's "re" module. Because `PyCall` and `PythonCall` both define pyimport, this example must qualify which one we are using, namely PythonCall's `pyimport`.

```julia
begin
    re = PythonCall.pyimport("re")
    words = re.findall("[a-zA-Z]+", "PythonCall.jl is great")
    sentence = Py(" ").join(words)
    pyconvert(String, sentence)  # convert Python string to Julia string
end
```

Try it.
"""

# ╔═╡ 5805328f-cb15-41c7-8b4c-9990cd0df319


# ╔═╡ 42a42775-e835-49ca-b0dd-819d1f601a51
md"""
### Wrapper Types
"""

# ╔═╡ d9fb6972-80de-4ee8-9026-615eb136d69d
md"""
A wrapper is a Julia type that wraps a Python object, but gives it Julia semantics. For example, the PyList type wraps Python's list object.

```julia
let
    x = pylist([3,4,5])
    y = PyList{Union{Nothing, Int64}}(x)
    push!(y, nothing)
    append!(y, 1:2)
    x
end
```

Try it.
"""

# ╔═╡ 966da8a4-6cee-47c3-bca1-a90460218f3f


# ╔═╡ 4adcd9bb-40dc-47dc-941d-66d34465b218
md"""
There are wrappers for other container types, such as PyDict and PySet.

```julia
let
    x = PythonCall.pyimport("array").array("i", [3,4,5])
    y = PythonCall.PyArray(x)
    println(sum(y))
    y[1] = 0
    x
end
```

Try this example too.
"""

# ╔═╡ adb77392-c109-4e90-917e-75a23e7f21a9


# ╔═╡ cec37a44-bfa8-40cb-af87-583f8fe47d16
md"""
`PyArray` directly wraps the underlying data buffer, so array operations such as indexing are about as fast as an ordinary Array.
"""

# ╔═╡ e76f12c2-d3ab-4914-9505-5a312793f95b
md"""
### Configuration
"""

# ╔═╡ bcbbef7d-ef13-4050-9f41-7f2ea7997a93
md"""
By default, PythonCall uses CondaPkg.jl to manage its dependencies. This will install Conda and use it to create a Conda environment specific to your current Julia project.
"""

# ╔═╡ b31af97d-0338-4270-a371-a79730cdf0f2
md"""
### Using your current Python installation
"""

# ╔═╡ 28d52fcd-a6ec-40f7-b709-74edeb40ec89
md"""

```julia
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/path/to/python"   # optional
ENV["JULIA_PYTHONCALL_EXE"] = "@PyCall"   # optional
```

By setting the CondaPkg backend to Null, no Conda packages will be installed. PythonCall will use your current Python installation.

If Python is not in your `PATH`, you will need to set the `JULIA_PYHTHONCALL_EXE` environment variable to include it in your path.

If you also use PyCall, you can set the `JULIA_PYTHONCALL_EXE` environment variable to use the same interpreter.
"""

# ╔═╡ 0988a074-8d2f-4f3c-a66e-759d9b749b44
md"""
### Using your current Conda environment
"""

# ╔═╡ 4ac0b58b-583d-48b9-95d1-125c565569e3
md"""
```julia
ENV["JULIA_CONDAPKG_BACKEND"] = "Current"
ENV["JULIA_CONDAPKG_EXE"] = "/path/to/conda"   # optional
```

Note that this configuration will still install any required Conda packages into your Conda envirnment.

If `conda`, `mamba`, and `micromamba` are not in your `PATH` you will need to set `JULIA_CONDAPKG_EXE` to include them.
"""

# ╔═╡ ce032143-258e-4dc8-aee5-0c2f3839520b
md"""
### Using your current Conda, Mamba, and MicroMamba environment
"""

# ╔═╡ e32b228b-50d5-4024-9281-d387445830bc
md"""
```julia
ENV["JULIA_CONDAPKG_BACKEND"] = "System"
ENV["JULIA_CONDAPKG_EXE"] = "/path/to/conda"   # optional
```

The System backend will use your preinstalled Conda environment.
"""

# ╔═╡ 69dfef23-dd0d-427b-9ef5-a7561d303311
md"""
### Installing Python packages
"""

# ╔═╡ 7ef3bb3a-6341-4c45-874d-131b45877ece
md"""
Assuming you are using `CondaPkg.jl`, PythonCall uses it to automatically install any Python packages. For example,

```julia
using CondaPkg
# enter package manager
conda add some_package
```
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "6e7bcec4be6e95d1f85627422d78f10c0391f199"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

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

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─3ba5672b-18ff-4320-a4d2-954e0b873d47
# ╟─6899990c-00b2-40e5-941c-df8f7073183c
# ╟─2ec6992f-38de-4ff5-a2f3-21f9643657e9
# ╟─14c3864e-b094-4c46-b38d-1abb375dc44a
# ╠═a833decc-ef58-46a7-923d-877204e716ff
# ╠═41103ed5-6145-4772-84d6-8858ebb560ac
# ╟─306378d0-91a3-4f07-94fb-45b8cc474274
# ╠═667e42f2-b6b5-4eb9-817f-57aee042d7ca
# ╟─06d63e41-f9b9-4c75-87d4-241beb454dc5
# ╠═20c226ca-44e6-4565-b601-07f5df7828ac
# ╟─7728f44a-93e1-4514-8964-6351a5ebff07
# ╠═8648cd40-d805-49b3-b05a-93a4b6ea0a20
# ╟─f790ac60-d874-4e95-8d34-252f44ff32b1
# ╠═90c8b0db-f94b-440c-b05a-72e3f6483051
# ╟─d8e3430e-5ad7-44c8-8a08-be067babac72
# ╟─5d883d2b-515f-4286-a2ab-15127ac6b5ea
# ╟─a836055d-9d11-46c7-bfca-a6759259e145
# ╟─a7a6b13e-cda7-4629-9b54-950759e6bab5
# ╟─09aa5c19-7d18-43b1-9623-1609f8792070
# ╠═0db12101-fcf9-48a3-a337-ad53a4713f6d
# ╟─93e945de-7ceb-4124-b3db-f6b64e2e46bb
# ╠═fa825783-6f48-4f0d-b346-220ae9e2fb11
# ╟─99e16357-0113-4a12-8a22-f627c6fdef84
# ╠═5805328f-cb15-41c7-8b4c-9990cd0df319
# ╟─42a42775-e835-49ca-b0dd-819d1f601a51
# ╟─d9fb6972-80de-4ee8-9026-615eb136d69d
# ╠═966da8a4-6cee-47c3-bca1-a90460218f3f
# ╟─4adcd9bb-40dc-47dc-941d-66d34465b218
# ╠═adb77392-c109-4e90-917e-75a23e7f21a9
# ╟─cec37a44-bfa8-40cb-af87-583f8fe47d16
# ╟─e76f12c2-d3ab-4914-9505-5a312793f95b
# ╟─bcbbef7d-ef13-4050-9f41-7f2ea7997a93
# ╟─b31af97d-0338-4270-a371-a79730cdf0f2
# ╟─28d52fcd-a6ec-40f7-b709-74edeb40ec89
# ╟─0988a074-8d2f-4f3c-a66e-759d9b749b44
# ╟─4ac0b58b-583d-48b9-95d1-125c565569e3
# ╟─ce032143-258e-4dc8-aee5-0c2f3839520b
# ╟─e32b228b-50d5-4024-9281-d387445830bc
# ╟─69dfef23-dd0d-427b-9ef5-a7561d303311
# ╟─7ef3bb3a-6341-4c45-874d-131b45877ece
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
