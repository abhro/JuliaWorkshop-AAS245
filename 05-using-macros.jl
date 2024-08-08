### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 53265062-78c7-434a-8d96-089cdd758bf0
using PlutoUI; TableOfContents()

# ╔═╡ 21ea940f-d90d-40da-aafe-959bb25041f2
md"""
# Using Macros
"""

# ╔═╡ 21ad70d1-4d8f-40f7-adc2-744430040e07
md"""
Unlike Matlab, Python, and R, Julia has macros. They are a very powerful programming construct. Macros change existing source code or generate entirely new code. In essence, they simplify programming by automating mundane coding tasks.

Preprocessor "macro" systems, like that of C/C++, work at the source code level by perform textual manipulation and substitution before any parsing or interpretation is done. Whereas, Julia macros work at the level of the abstract syntax tree after parsing, but before compilation is done.
"""

# ╔═╡ a7f62fc0-e877-4496-bc8f-347f0af1f730
md"""
## The Abstract Syntax Tree
"""

# ╔═╡ aa6a3246-c767-4e7a-8e5f-f6ff81a2d8e1
md"""
```julia
begin
    str1 = "3 * 4 + 5"
    ex1 = Meta.parse(str1)
    println(typeof(ex1))
end
```

Try it.
"""

# ╔═╡ 9c45b534-cd6c-4b3d-9ff4-d035347c87e5


# ╔═╡ 6bc277f1-3ab6-4e12-92fc-65e6288af2c1
md"""
Note that the parsed expression is an `Expr` type.

```julia
dump(ex1)
```

Show or dump the abstract syntax tree (AST) of the expression.
"""

# ╔═╡ f09d458f-7846-46a5-a814-b51771bf91b7


# ╔═╡ 41a9007c-e9a1-4228-be45-55becfffa944
md"""
The `Expr` type has two fields, a `head` and an `args`. The `head` declares the item is a function `call` of type `Symbol` and the `args` is an `Array` of function arguments, namely the `+` operator of type `Symbol` and the two values to be added. The first value is another `Expr` type and the second is a 64-bit integer. The second nested `Expr` is also composed of a `Symbol` call. Its arguments are the multiplication `*` operation of type `Symbol` and two 64-bit integers. Note that the multiplication operation is nested in the addition operation, since multiplication has precedence over addition.

Note that the expression `ex1` can be written directly in AST syntax.

```julia
begin
    ex2 = Expr(:call, :+, Expr(:call, :*, 3, 4), 5)
    dump(ex2)
    println(ex1 == ex2)
end
```

Try it.
"""

# ╔═╡ 06f3afb2-869a-4e2b-af58-4366b27d0a65


# ╔═╡ 631362df-1aba-4769-b1cb-bf6d252799e7
md"""

!!! note
    If you prefer writing code as an AST, you can do so in Julia.

    However, you may not be popular with other coders.
"""

# ╔═╡ df05c665-111a-461c-be17-801c69b931c4
md"""
## The `:` character
"""

# ╔═╡ cd172701-d49f-4ae5-bfae-509c861438af
md"""
The `:` character has two syntatic purposes in Julia. It signifies a `Symbol` or a `Quote`.
"""

# ╔═╡ 64e8a372-b7ca-497e-aa02-d5be96f3616c
md"""
### Symbols
"""

# ╔═╡ f768336e-10d5-47b6-b640-9a2876981f2a
md"""
A `Symbol` is an interned string used as a building-block of an expression.

```julia
begin
    sym = :foo
    println(typeof(sym))
    println(sym == Symbol("foo"))
end
```

Try it.
"""

# ╔═╡ bc5709bd-8c44-4d58-82cc-e14eef20cfe6


# ╔═╡ c1e446d9-63e9-4b9d-a0da-c1321c7f0a38
md"""
The `Symbol` constructor takes a variable number of arguments and concatenates them together to create a `Symbol` string.

```julia
Symbol("func", 10)
```

Try this.
"""

# ╔═╡ 8863cbcc-4ff6-4714-b7d8-c39373af384e


# ╔═╡ eaf210c7-5605-48ab-b025-05a71d792267
md"""
In the context of expressions, symbols are used to indicate access to variables. When an expression is evaluated, a symbol is replaced with the value bound to that symbol in the appropriate scope.
"""

# ╔═╡ 2599447f-ad96-4832-ba3c-980e67c8ee57
md"""
### Quoting
"""

# ╔═╡ 638b51c7-e928-4708-b0b8-dd75cc8e9dbc
md"""
The second purpose of the `:` character is to create expression objects without using the explicit `Expr` constructor. This is referred to as *quoting*.

```julia
begin
    ex3 = :(3 * 4 + 5)
    dump(ex3)
    println(ex1 == ex2 == ex3)
end
```

Try this example.
"""

# ╔═╡ 611e0640-dc39-4846-89c9-172bc001e994


# ╔═╡ e7934544-0344-48d3-8046-c3796874e5b3
md"""
So, there are three ways to construct an expression allowing the programmer to use whichever one is most convenient.

There is a second syntatic form of quoting, called a `quote` block (i.e. `quote ... end`), that is commonly used for multiple expressions.

Try the following example.

```julia
begin
    ex4 = quote
        x = 1
        y = 2
        x + y
    end
    println(typeof(ex4))
end
```
"""

# ╔═╡ c47eb25e-ff42-47a0-9949-11228e3c0535


# ╔═╡ ee0bf6a4-969b-4f58-ad8e-d175498fa6bf
md"""
### Interpolation
"""

# ╔═╡ a7c09c1a-71da-41cb-8ab9-ec1b57e8e77b
md"""
Julia allows *interpolation* of literals or expressions into quoted expressions by prefixing a variable with the `$` character. For example.

```julia
begin
    a = 1;
    ex5 = :($a + b)
    dump(ex5)
end
```

Try it.
"""

# ╔═╡ 24368ee3-43e1-44c0-8293-75cc1e8f4f6c


# ╔═╡ 12b97064-d8da-49cc-92aa-f5155b5c1b52
md"""
Splatting is also possible.

```julia
begin
    args = [:x, :y, :x];
    dump(:(f(1, $(args...))))
end
```

And this too.
"""

# ╔═╡ 9de33d97-b08d-4a5a-b53a-862c0bfd2cb0


# ╔═╡ 3f9a3e17-a75c-46d1-8d15-167bf00d3008
md"""
### Expression evaluation
"""

# ╔═╡ 7075a117-9f0f-46e6-ba77-d356e1352300
md"""
Julia will evalution an expression in the global scope using `eval`.
```julia
println(eval(ex1))
```
"""

# ╔═╡ ade3dfbb-23bb-4ca6-8a33-3ec16dbed523


# ╔═╡ 6d8b869a-ecdc-4002-8ed5-1ab99749616e
md"""
## Macros
"""

# ╔═╡ e0d48a51-bf75-4688-b0f9-250872eef3d2
md"""
### Basics
"""

# ╔═╡ 7748a6ed-946e-4116-8c9c-7ee57f2c8cf4
md"""
Defining a macro is like defining a function. For example,

```julia
macro sayhello(name)
    return :( println("Hello", $name))
end
```
"""

# ╔═╡ f1ba0108-fbb7-4128-98f7-89c4f1e25790


# ╔═╡ 8bbe63e8-6a6e-4f96-8025-bbae25983a7e
md"""
Now invoke the `@sayhello` macro.

```julia
@sayhello("world")
```
"""

# ╔═╡ 6e6829b6-1d2a-46e8-bee2-dd4e116d4b3d


# ╔═╡ 2494a756-0b40-4b6e-b79f-deb616183384
md"""
Parentheses around the arguments are optional.

```julia
@sayhello "world"
```
"""

# ╔═╡ 3f9fcacd-05e6-4dd7-ba01-491750f3940d


# ╔═╡ f7fffc76-3607-4d90-ba75-d521366f7fb1
md"""
!!! note

    When using the parenthesized form, there should be no space between the macro name and the parentheses. Otherwise the argument list will be interpreted as a single argument containing a tuple.
"""

# ╔═╡ 0ec77681-e0ed-4d80-bf40-ddc1b92fbcc5
md"""
### Usage
"""

# ╔═╡ 437a78f9-b2fd-432c-b3c6-e87b55d81b05
md"""
Macros can be used for many tasks such as performing operations on blocks of code, e.g. timing expressions (`@time`), threading `for` loops (`@threads`), and converting expressions to broadcast form (`@.`). Another use is to generate code. When a significant amount of repetitive boilerplate code is required, it is common to generate it programmatically to avoid redundancy.

Consider the following example.

```julia
struct MyNumber
    x::Float64
end
```
"""

# ╔═╡ bcce9741-7e63-4609-9529-0e9323b8eb21


# ╔═╡ a278e94d-1f4c-422c-9ac7-c950750a6ef9
md"""
We want to add various methods to it. This can be done programmatically using a loop.

```julia
begin
    println(length(methods(sin)))
    for op = (:sin, :tan, :log, :exp)
        eval(quote
            Base.$op(a::MyNumber) = MyNumber($op(a.x))
        end)
    end
    println(length(methods(sin)))
end
```

Try this example.
"""

# ╔═╡ 386c7942-700d-4486-bede-853a6c2ef811


# ╔═╡ 801548cf-98f6-4be8-9bfc-1fbe48a72f21
md"""
A slightly shorter version of the above code generator that use the `:` prefix is:

```julia
for op = (:sin, :cos, :tan, :log, :exp)
    eval(:(Base.$op(a::MyNumber) = MyNumber($op(a.x))))
end
```

An even shorter version of the code uses the `@eval` macro.

```julia
for op = (:sin, :cos, :tan, :log, :exp)
    @eval Base.$op(a::MyNumber) = MyNumber($op(a.x))
end
```

For longer blocks of generated code, the `@eval` macro can proceed a code block:

```julia
@eval begin
    # multiple lines
end
```
"""

# ╔═╡ 80fc7c35-ce72-43a0-af5a-e5ace3b9104b
md"""
## Generated Functions
"""

# ╔═╡ c45c3f31-08b8-4eea-b664-cfcc7a05caa7
md"""
A special macro is `@generated`. It defines so-called *generated functions*. They have the capability to generate specialized code depending on the types of their arguments with more flexibility and/or less code than what can be done with multiple dispatch. Macros work with expressions at the AST level and cannot acces the input types. Whereas, generated functions are expanded at the time when the argument types are known, but before the function is compiled.

Here is an example.

```julia
@generated function foo(x)
    Core.println(x)
    return :(x * x)
end
```
"""

# ╔═╡ 1fbdea1e-ef3a-4ace-ae76-1ccda2273087


# ╔═╡ 1522314b-4982-439e-ab97-933869a5f307
md"""
```julia
begin
    x = foo(2)
    x
end
```

!!! note
    The result of the function call is the return type, not the return value.
"""

# ╔═╡ 63ed8db8-2104-4370-8a75-e9348c50cdf8


# ╔═╡ 3de9644d-b66a-45b8-945d-6e000319f2b9
md"""
```julia
begin
    y = foo("bar")
    y
end
```
"""

# ╔═╡ 28487c59-7df4-41b9-859a-79438b1d3b06


# ╔═╡ 78bb1e78-70ba-4f58-abb0-c0de3272eb50
md"""
Generated functions differ from regular functions in five ways. They:

1. are prefixed by `@generated`. This provides the AST with some additional information.
2. can only access the types of the arguments, not their values.
3. are nonexecuting. They return a quoted expression, not a value.
4. can only call previously defined functions. (Otherwise, `MethodErrors` may result.)
5. must not *mutate* or *observe* any non-constant global state. They can only read global constants, and cannot have side effects. In other words, they must be completely pure.
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
# ╟─53265062-78c7-434a-8d96-089cdd758bf0
# ╟─21ea940f-d90d-40da-aafe-959bb25041f2
# ╟─21ad70d1-4d8f-40f7-adc2-744430040e07
# ╟─a7f62fc0-e877-4496-bc8f-347f0af1f730
# ╟─aa6a3246-c767-4e7a-8e5f-f6ff81a2d8e1
# ╠═9c45b534-cd6c-4b3d-9ff4-d035347c87e5
# ╟─6bc277f1-3ab6-4e12-92fc-65e6288af2c1
# ╠═f09d458f-7846-46a5-a814-b51771bf91b7
# ╟─41a9007c-e9a1-4228-be45-55becfffa944
# ╠═06f3afb2-869a-4e2b-af58-4366b27d0a65
# ╟─631362df-1aba-4769-b1cb-bf6d252799e7
# ╟─df05c665-111a-461c-be17-801c69b931c4
# ╟─cd172701-d49f-4ae5-bfae-509c861438af
# ╟─64e8a372-b7ca-497e-aa02-d5be96f3616c
# ╟─f768336e-10d5-47b6-b640-9a2876981f2a
# ╠═bc5709bd-8c44-4d58-82cc-e14eef20cfe6
# ╟─c1e446d9-63e9-4b9d-a0da-c1321c7f0a38
# ╠═8863cbcc-4ff6-4714-b7d8-c39373af384e
# ╟─eaf210c7-5605-48ab-b025-05a71d792267
# ╟─2599447f-ad96-4832-ba3c-980e67c8ee57
# ╟─638b51c7-e928-4708-b0b8-dd75cc8e9dbc
# ╠═611e0640-dc39-4846-89c9-172bc001e994
# ╟─e7934544-0344-48d3-8046-c3796874e5b3
# ╠═c47eb25e-ff42-47a0-9949-11228e3c0535
# ╟─ee0bf6a4-969b-4f58-ad8e-d175498fa6bf
# ╟─a7c09c1a-71da-41cb-8ab9-ec1b57e8e77b
# ╠═24368ee3-43e1-44c0-8293-75cc1e8f4f6c
# ╟─12b97064-d8da-49cc-92aa-f5155b5c1b52
# ╠═9de33d97-b08d-4a5a-b53a-862c0bfd2cb0
# ╟─3f9a3e17-a75c-46d1-8d15-167bf00d3008
# ╟─7075a117-9f0f-46e6-ba77-d356e1352300
# ╠═ade3dfbb-23bb-4ca6-8a33-3ec16dbed523
# ╟─6d8b869a-ecdc-4002-8ed5-1ab99749616e
# ╟─e0d48a51-bf75-4688-b0f9-250872eef3d2
# ╟─7748a6ed-946e-4116-8c9c-7ee57f2c8cf4
# ╠═f1ba0108-fbb7-4128-98f7-89c4f1e25790
# ╟─8bbe63e8-6a6e-4f96-8025-bbae25983a7e
# ╠═6e6829b6-1d2a-46e8-bee2-dd4e116d4b3d
# ╟─2494a756-0b40-4b6e-b79f-deb616183384
# ╠═3f9fcacd-05e6-4dd7-ba01-491750f3940d
# ╟─f7fffc76-3607-4d90-ba75-d521366f7fb1
# ╟─0ec77681-e0ed-4d80-bf40-ddc1b92fbcc5
# ╟─437a78f9-b2fd-432c-b3c6-e87b55d81b05
# ╠═bcce9741-7e63-4609-9529-0e9323b8eb21
# ╟─a278e94d-1f4c-422c-9ac7-c950750a6ef9
# ╠═386c7942-700d-4486-bede-853a6c2ef811
# ╟─801548cf-98f6-4be8-9bfc-1fbe48a72f21
# ╟─80fc7c35-ce72-43a0-af5a-e5ace3b9104b
# ╟─c45c3f31-08b8-4eea-b664-cfcc7a05caa7
# ╠═1fbdea1e-ef3a-4ace-ae76-1ccda2273087
# ╟─1522314b-4982-439e-ab97-933869a5f307
# ╠═63ed8db8-2104-4370-8a75-e9348c50cdf8
# ╟─3de9644d-b66a-45b8-945d-6e000319f2b9
# ╠═28487c59-7df4-41b9-859a-79438b1d3b06
# ╟─78bb1e78-70ba-4f58-abb0-c0de3272eb50
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
