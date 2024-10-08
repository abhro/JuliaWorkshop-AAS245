### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 35f025bb-7ee1-482b-bdc2-f5448b7ec027
using Pkg

# ╔═╡ 1e47e383-b735-4c57-a300-2afe8491b49a
using PlutoUI; TableOfContents()

# ╔═╡ 40abc83f-b4bd-479f-8671-189cc712d792
using Measurements

# ╔═╡ 88ca2a73-6203-447c-afcc-9e370a82076b
using Unitful

# ╔═╡ d1366a55-b4fc-4ddb-b5c2-5f3381c48b49
html"<button onclick='present()'>present</button>"

# ╔═╡ 09193424-25b9-45ce-840f-f24bbcc46c9d
md"""
## Introduction
"""

# ╔═╡ 54aa12fb-6a56-404d-86e3-4f2a0478319f
md"""
### Historical Context

Twenty-six years ago at ADASS VI, Harrington and Barrett hosted a Birds-of-a-Feather session entitled "Interactive Data Analysis Environments". Based on their review of over a dozen interpreted programming languages such as Glish, GUILE, IDL, IRAF, Matlab, Perl, Python, and Tcl; they recommended that Python be considered the primary language for astronomical data analysis. The primary reasons were that the language was simple to learn, yet powerful; well supported by the programming community; and had FORTRAN-like arrays. However, for good performance, the multi-dimensional arrays needed to be written in a compiled language, namely C. So Numerical Python suffered from the "two language problem".
"""

# ╔═╡ 11643fb1-2898-47cc-a2c7-71a049a5ebd5
md"""
### Why Julia?

In about 2009, four faculty members at MIT, who were not satisfied with the state of scientific computing, decided to develop a high performance, scientific programming language. After ten years of development, they release Julia Version 1.0 on August 8, 2018. Their aims were to create an open-source interpreted language that was concise, extensible, and high performance.
"""

# ╔═╡ 24e6db88-8aaa-4f51-9d67-9e10c9737e39
md"""
### What Differentiates Julia From Other Languages?

* Julia is **composable**.
* Julia is **concise**.
* Julia is **high performance**.
* Julia is **productive**.
* Julia is **easy to maintain**.
* Julia is **free and open-source**.
"""

# ╔═╡ a2142e30-f585-4e7a-951e-c906dcee175c
md"""
### Why Have I migrated to Julia?

Although an early advocate and developer of Numerical Python, I knew its limitations, namely, the two language problem. Therefore, once a better scientific programming language came along, I was prepared to migrate to it. Julia is that language.
"""

# ╔═╡ b1ed2c4e-f5fa-4e5e-87d8-7af6f80a83ca

md"""## Getting Started"""




# ╔═╡ 2ed5cce9-989a-4972-a27c-e3ec9cd4257e
md"""
---
"""

# ╔═╡ 1eff0f30-698b-4e70-848b-c39c8087cbe1
md"""
### Starting Julia
"""

# ╔═╡ 7f3357bc-4103-4a35-af21-9c86f5a0ec2f
md"""
Enter `julia` at the terminal prompt. Set the number of threads to `auto`. Threads will be discussed later in Parallel Computing.

```
$ julia --threads=auto

               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.10.0 (2023-12-25)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```

!!! tip

    The command line option "-q" can be used to remove the start-up banner.
"""

# ╔═╡ cb01ca1b-838a-4137-90d6-fc9ecf9b2468
md"""
---
"""

# ╔═╡ 53bc6b61-e69c-4c62-9eca-0b649a2486bd
md"""
### Stopping Julia
"""

# ╔═╡ 7475c896-d1b1-4429-9ba8-8e78de41e0b0
md"""
To exit Julia, enter `<Ctrl-D>` or `exit()`

```julia-repl
julia> <Ctrl-D>
```

!!! tip
    Don't do this now!

"""

# ╔═╡ 5df8264e-6e37-4674-abdf-2b05c530787f
md"""
---
"""

# ╔═╡ c7ce5c61-b358-4feb-808e-f092e43351a5
md"""
### The command line  or  REPL (Read-Eval-Print-Loop)
"""

# ╔═╡ f646ca14-c01e-47ee-8e2b-052d9db0985b
md"""
Our first command:
```julia
println("Hello World")
```
"""

# ╔═╡ 4a404280-2845-4deb-8eee-2dcdcb9aed27


# ╔═╡ 7813824a-cae9-4b97-ac90-e542fbd630d5
md"""
!!! note
    Unlike Jupyter and the REPL, Pluto prints the result above the line, not below.

Our first calculation
```julia
a = 4
```
"""

# ╔═╡ 6ac51e87-87a2-4ccc-9f08-0028700b3cda


# ╔═╡ 27208179-35c3-43c1-9548-3620c8aa7680
md"""
```julia
b = 2
```
"""

# ╔═╡ 40d8d18c-3713-4e77-812d-9d77a4e1ac50


# ╔═╡ aa3e9db7-49d1-40f8-b745-6c4faa2197e1
md"""
```julia
a + b
```
"""

# ╔═╡ 756d3f69-a5b2-4ee6-bae6-94c513baae6c


# ╔═╡ d882d5f8-6898-4810-aa4a-7b064cff80b0
md"""
---
"""

# ╔═╡ 98fbf4c1-fae3-4cab-abea-9b7746d807a9
md"""
### Other REPL Modes
"""

# ╔═╡ 7f41ba9b-e650-45a6-bc90-b6e119dc60b2
md"""
#### Help, '?'
"""

# ╔═╡ 419a6dec-1db0-477f-911f-049223b5674f
md"""
For help mode,

```
julia> ?
help?> println
search: println printstyled print sprint isprint

println([io::IO], xs...)

Print (using print) xs to io followed by a newline. If io is not supplied, prints to the default output stream stdout.

See also printstyled to add colors etc

Examples
≡≡≡≡≡≡≡≡≡≡

julia> println("Hello, world")
Hello, world

julia> io = IOBuffer();

julia> println(io, "Hello", ',', " world.")

julia> String(take!(io))
"Hello, world.\n"
```
Enter 'delete' or 'backspace' to exit help"""

# ╔═╡ 98340265-f51e-47a0-95d2-df179b87f54b


# ╔═╡ b4f40ffa-138a-4451-8f3c-34f491ea6b13
md"""
#### Shell, ';'
"""

# ╔═╡ 8ee7f43d-bf75-4975-ac64-54c2d5a0174a
md"""
For shell mode,
```julia-repl
julia> ;
shell> pwd
/Users/myhomedir
```
Enter 'delete' or 'backspace' to exit shell
"""

# ╔═╡ d4368e22-60c6-456a-94a5-56e6dfdb26d7


# ╔═╡ 62786c37-744e-4bd9-8fd9-68c35289f5f5
md"""
#### Package Manager, `]`
"""

# ╔═╡ d1e9c51c-efb9-4dcb-9d28-8c54a235fbb4
md"""
There are two ways of using the package manager. Either using the module

```julia-repl
julia> using Pkg
```

or through the pkg-mode in the repl:

```julia-repl
julia> ]
pkg>
```

For package manager help,
```julia-repl
pkg> ? `return`
```

Returns a brief summary of package commands

To add a package,
```julia-repl
pkg> add <package>
pkg> add <package1>, <package2>
```
When adding a package, the Julia on-line repository will be searched. The package and its dependencies will then be downloaded, compiled, and installed. This may take anywhere from a few seconds to a few minutes depending on the size of the package and its dependencies.

To use or load a package (after it has been added),
```julia-repl
julia> using <package>
julia> using <package1>, <package2>
```

A feature of Pluto is running the `using` command will add the package, if it hasn't alaredy been added.
"""

# ╔═╡ 948b4044-aee8-4b21-b9a5-4949adb9e01e


# ╔═╡ b27578b2-f5f5-4e46-82e6-0007be187ba6
md"""
To check the manifest, run either command in pkg mode:

```julia-repl
pkg> status
pkg> st
```

or in the REPL, run
```julia-repl
julia> Pkg.status()
```
"""

# ╔═╡ ea23f70b-b9f2-42fa-96de-9eb5ab6d1a0d
Pkg.status()

# ╔═╡ 065265a5-c9ad-4a39-b14d-f4e2e49d3f7a
md"""
To update a package in the manifest (use any of the following commands):
```julia-repl
pkg> update <package>

pkg> up <package>

julia> Pkg.update("package")
```

To update all packages in the manifest,
```julia-repl
pkg> up

pkg> update

julia> Pkg.update()
```

To garbage collect packages not used for a significant time,
```julia-repl
pkg> gc

julia> Pkg.gc()
```
"""

# ╔═╡ a4392002-6b0c-48bb-beed-1360a54bbafa
md"""
---
"""

# ╔═╡ 563f07ad-6aed-495e-85fb-bae4a1755ac2
md"""
The Measurements package enables variables to have both values and errors.
Let's add the Measurements package using the `using` statement.
"""

# ╔═╡ 297cd86c-5e9d-4f70-b11a-cbae8fa96d1e
md"""
Let's do some more calculations.
```julia
m1 = measurement(4.5, 0.1)
```
"""

# ╔═╡ 8f016c75-7768-4418-8c57-100db3073c85


# ╔═╡ 094b6f30-cbd6-46b1-8e0c-3fdb1ef18261
md"""Typing '`measurements`' is rather awkward. There must be a better way. How about the following?
```julia
m2 = 15 ± 0.3
```
where the plus-minus character is entered using LaTeX syntax followed by tab, i.e., \pm<tab>.
"""

# ╔═╡ 7ba8dc19-e0ca-40de-a778-7583ca70978d


# ╔═╡ 668abc35-fdc3-430f-8c90-de3c2c2cd77b
md"""
One of the features of Julia is that it understands unicode. For example, expressions in a printed document that contain greek characters can be entered as greek characters in your code. Let's calculate the following expression.

```julia
α = m1 + m2
```
"""

# ╔═╡ 232cc444-03b7-442a-8737-8b7725b43421


# ╔═╡ d2a2d0bc-e883-439f-8e34-166e2369caef
md"""
!!! note

    Notice that the error of the result α has been propagated correctly.
"""

# ╔═╡ ce6545c9-15eb-44d9-a999-770c222e249a
md"""
Let's add another package called Unitful, which enables attaching units to variables.
"""

# ╔═╡ c24f1ddd-5e31-4073-a627-86cedb1d44c2
md"""
Now let's create two new values `m3` and `m4` with units attached, and then multiply them together to create a third variable `β`.

```julia
m3 = (32 ± 0.1)u"m/s"
m4 = (9.8 ± 0.3)u"s"
β = m3 * m4
```
"""

# ╔═╡ 63a4b27a-5361-4d95-8787-ae31ca7987fe


# ╔═╡ 15674bb0-2fe1-40b1-a6c0-3a5a64a6a5c3


# ╔═╡ 70f08712-002c-4adc-84b1-73a8655d8a44


# ╔═╡ 3fa06ac2-20f3-4df4-8ce3-b0e1c21ca842
2pi

# ╔═╡ cf4a0e8f-9210-4f1e-84d4-ee7ff09aaf61
md"""
The variable `β`'s value now has an associated error and unit.

Let's see if this works with one dimensional arrays or vectors.

```julia
γ = [10 ± 0.1, 20 ± 0.2, 30 ± 0.3]u"m/s" .* [15 ± 0.01, 25 ± 0.02, 25 ± 0.03]u"s"
```

Note the dot '`.`' before the multiplication character '`*`'.  This means element-wise multiplication. Whereas the multiplication character '`*`' by itself means matrix multiplication. If you are coming from Python, this difference may take a little time.
"""

# ╔═╡ e00b826d-1bbb-4413-a907-eb181369526b


# ╔═╡ b56255c6-9d3b-4e2f-a9a0-c6fe69990f3d
md"""
!!! note

    What have we learned about the Julia command line and features?

    * Julia has four command line modes: **REPL**, **help**, **shell**, and **package manager**.

    * Julia understands **unicode**.

    * Julia packages are **composable**. It means that independent packages are compatible and work together without modification, as demonstrated by the Measurements and Unitful packages.
"""

# ╔═╡ aaf91290-98db-4508-9519-18f7f25439dd
md"""
---
"""

# ╔═╡ 6809e5b1-86da-460e-bcba-ca379cc54dfe
md"""
## Language Basics
"""

# ╔═╡ 5cd072cb-5d71-4a08-8e41-4eaaa7faaa5c
md"""
Because of Julia's multiple dispatch, types and functions are loosely bound. Thus, it is not a true object-oriented language, where functions are selected for execution using single dispatch. Multi-dispatch will be explained later when we dicsuss functions.
"""

# ╔═╡ 90755b03-8a48-40ea-a43a-24bf01248d4b
md"""
---
"""

# ╔═╡ 9f61328e-59df-49b2-8c61-c3f339e4254a
md"""
### Comments
"""

# ╔═╡ f37bc13e-fa91-4166-983b-fd13a8493435
md"""
A comment string begins with a "`#`" and extends to the end of the line.

A comment block begins and ends with "###".
"""

# ╔═╡ 915b1ea3-4563-4daf-bf99-25e1de6ef6ca
md"""
---
"""

# ╔═╡ b3fad295-b319-4eba-b250-b83ca7e2d278
md"""
### Types
"""

# ╔═╡ 0d0c11c0-d39f-462c-9fb6-ab90ca98d230
md"""
The optional type operator "`::`" is used to specify a type to expressions and variables, especially when writing functions. If no type is specified, Julia will infer the type based on context.

There are two primary reasons for type annotation:

1. As an assertion to confirm that the code is working properly, and
2. To provide extra type information to the compiler to improve performance.
"""

# ╔═╡ a02bbbbb-6b3f-47ef-a11f-1db9b802db6f
md"""
```julia
(1+2)::Float32
(1+2)::Int
```

Let's see how this works. Try the above examples.
"""

# ╔═╡ 2262c860-c06c-4293-8e6d-b616228cb301
(1+2)::Float32

# ╔═╡ 68e64f74-8a6b-403e-a404-52fb9cdea54b
(1+2)::Int

# ╔═╡ 0887eca0-6760-4d9b-b44e-d1a14059aede
md"""Julia has various categories of types within a type-hierarchy. The following are some of the more common types.

!!! note
    Types should be capitalized.
"""

# ╔═╡ 86ce3011-3e52-43e4-83ec-1a4422073d45
md"""
#### Abstract Types
"""

# ╔═╡ 0ad9aa76-f6c7-4368-8ae4-58daa548e065
md"""
"`abstract type`" declares a type that cannot be instantiated, and serves only as a node in the type graph, thereby describing sets of related concrete types.

Let's create an abstract type.
```julia
abstract type Widget end
```
"""

# ╔═╡ 1bc3da9e-143c-489c-b8de-a29dc48f17cb


# ╔═╡ f00dd72a-8705-426b-9eb4-b91cf1ea95d4
md"""
And some Widget subtypes using the subtype operator "`<:`".

```julia
abstract type Round <: Widget end
abstract type Square <: Widget end
```
"""

# ╔═╡ d308df6b-14ec-49ec-8270-a3b9efd88517


# ╔═╡ 01805f02-f9f6-4e3e-8e93-a0628753130f


# ╔═╡ a90b9011-714e-41d1-b7a3-fb3eb9dc56da
md"""
The subtype and supertype of a type can be shown using the functions "`subtype`" and "`supertype`".
"""

# ╔═╡ 60e8eeac-0d2d-4148-a564-8c6c3dbc63ab
md"""
Show the supertype and subtypes of `Widget`.
```julia
supertype(Round)
subtypes(Widget)
```
"""

# ╔═╡ b8325403-9744-4a9d-ae64-be88671da89b


# ╔═╡ 4879dae5-442e-4dc6-90c9-366ff76912bb


# ╔═╡ e2e57f49-f848-468a-a6f5-482b6e1ad4ba
typeof(1)

# ╔═╡ 4c278c5a-3324-4245-8ddf-f5390167168f
md"""
!!! note
    The "`Any`" type is at the top of the hierarchy. It is the union of all types. In other words, it is the root node.

    When the type of an expression or variable cannot be inferred from the context, the type defaults to "`Any`".
"""

# ╔═╡ aace56e5-c11a-4306-bcb1-d08fda9dcd84
md"""
#### Primitive Types
"""

# ╔═╡ 3772a828-561d-4600-8e67-49a28cc6cf09
md"""
A primitive type is a concrete type whose data consists of plain old bits. Classic examples of primitive types are integers and floating-point values. Unlike most languages, Julia lets you declare your own primitive types, rather than providing only a fixed set of built-in ones.

Let's see what primitive types `Integer` and `AbstractFloat` contain.

```julia
subtypes(Integer)
```
"""

# ╔═╡ aa4a7ec0-a270-482b-abeb-7168de767938


# ╔═╡ b8e3b72a-e501-4164-b06c-cbb3282d9d11
md"""
```julia
subtypes(Signed)
```
"""

# ╔═╡ d9aa9f5e-31b6-49a3-bae8-a9b149e6ab91


# ╔═╡ 15b0159b-9c8c-4327-b73d-d7e19decde2a
md"""
```julia
subtypes(AbstractFloat)
```
"""

# ╔═╡ 5d5b1283-043b-437a-afda-75801808acc9


# ╔═╡ 6a6b2a0a-6bb6-4a67-b4c1-46631503918d
md"""Theoretically, a primitive type can have any number of bits, e.g., 5 or 17. Practically, the number of bits is constrained to multiples of 8. This is a limitation of the LLVM compiler, not Julia. So the `Bool` type is 8 bits, not 1 bit.
"""

# ╔═╡ 80cfbd64-61f9-46b9-8e55-ffb92e223604
md"""
#### Characters (' ') vs. Strings (" ")
"""

# ╔═╡ 877faa74-7490-44a3-9e97-b36b36050796
md"""
Unlike Python, single and double quotes have different meanings. Single quotes create characters. Double quotes create strings. The reason for this is Unicode.

```julia
'j'
```
"""

# ╔═╡ bba18435-d355-4fca-a6f5-10dacde17413


# ╔═╡ d9e911a8-13f9-41e5-ac36-4aee3ec24c59
md"""
```julia
Char(167)
```

or

```julia
'\u00A7'
```
"""

# ╔═╡ 5f72777b-a174-453c-8b18-ebf1f4bebe0d


# ╔═╡ 734a4185-4001-410f-affc-71b33e339339


# ╔═╡ c349f7b8-bdf0-4b94-b412-06c5e7f3cbc5
md"""    "This is a string" """

# ╔═╡ d8be9383-fb60-4938-9376-f91d59f21559


# ╔═╡ 31dfb05b-ed87-48f9-a74c-0055e46de160
md"""
Triple quotes work the same as in Python.
```julia
\"""
This is line 1.
This is line 2.
\"""
```
Try it.
"""

# ╔═╡ d2ada743-b82d-47c8-9b1d-4bd56de76e62


# ╔═╡ 32085179-7b01-4a3e-b810-d79b5b7aadc1
md"""
#### Composite Types
"""

# ╔═╡ ea15815e-0ae3-4f22-9dce-a17cb3a0560b
md"""
Composite types are called records, structs, or objects in various languages. A composite type is a collection of named fields, an instance of which can be treated as a single value.

In mainstream object oriented languages, such as C++, Java, Python and Ruby, composite types also have named functions associated with them, and the combination is called an "object". In Julia, all types are objects, but the objects have no bound functions. This is necessary because Julia selects the function/method using multiple dispatch, meaning that all argument types of a function are used to select the method, not just the first argument type.

Composite types are defined using the "`struct`" keyword followed by a block of field names. They are immutable (for performance reasons), unless modified by the "`mutable`" keyword.

```julia
struct Longday
    day::Int64
    frac::Float64
end
```

An instance of Longday is created as follows.
```julia
day1 = Longday(1, 0.5)
```
Let's create a Longday type and an instance of it.

"""

# ╔═╡ be09f5d0-daea-4f47-8dc8-33c875fca843


# ╔═╡ 10ec3b0d-1add-4f92-8f4c-b594ab3f0e68


# ╔═╡ 6ee4665d-c5b9-4881-ad65-15c6a8229f3f
md"""
The field can be access using "dot" notation as follows:
```julia
day1.day
day1.frac
```
"""

# ╔═╡ f5596a05-04de-4955-9575-4c035e0f1495


# ╔═╡ a1b4f7bb-8238-40d6-81cb-6d5e6c737134


# ╔═╡ 582fda28-df9c-4d4c-b376-78463a1e41bf
md"""
#### Type Union
"""

# ╔═╡ 3b8e773f-df6e-4b59-9f5d-e14366d02754
md"""
A type union is an abstract type that includes all instances of any of its argument types. The empty union `Union{}` is the leaf node of all Julia types.

```julia
Union{Int, Nothing}
```

The variable "`nothing`" is the singleton instance of the type "`Nothing`".

Try it.
"""

# ╔═╡ 91f35db2-6a17-42aa-8580-1dea220b8c11


# ╔═╡ 6d76f5fc-e59e-4cf3-a2ff-aa2c839d4661
md"""
#### Symbol Type
"""

# ╔═╡ a631464d-e08a-4a89-8c47-fd5a7b2dee16
md"""
A type used to represent identifiers in parsed Julia code, namely the Abstract Syntax Trees (ASTs). Also often used as a name or label to identify an entity (e.g., as a dictionary key). Symbols are created using the colon prefix operator "`:`".

Symbols can be confusing when you first meet them in Julia code.

```julia
:symbol
typeof(:symbol)
```
"""

# ╔═╡ 7a8faa02-34b1-4416-beab-2909fb56c767


# ╔═╡ 6e1a3b46-05f0-487d-933a-6ff0d9d43a2b


# ╔═╡ 9064310f-d82d-44dd-bd76-90d194ee84ac
md"""
#### Using Types
"""

# ╔═╡ 05adfd23-c809-4706-9bf2-1a0a2445748b
md"""
The type hierarchy allows variables and functions to be constrained to a particular set of types. Let's try a simple example.
"""

# ╔═╡ 8c59d1ec-606a-4665-8e9e-262d34415e85
md"""
Enter the following expressions:
```julia
arg1::Float32 = 12.3
typeof(arg1) <: Integer
typeof(arg1) <: AbstractFloat
```
"""

# ╔═╡ 67a4ff9f-c75f-444c-9091-e9b5c17ee773


# ╔═╡ cfd93268-174f-4a7e-9f98-3d5787c9392c


# ╔═╡ a96dd069-09aa-4add-baba-99ffae36bfe8


# ╔═╡ 8a3aa0d3-1ade-4961-975d-b39899731ffe
md"""
!!! note

    What new things have we learned about Julia?

    * Julia has a type hierarchy with the type "Any" at the top.

    * Julia defines characters and strings using single and double quotes, respectively.

    * Julia defines composite types using the "struct" keyword.

    * Julia allows a set of types to be defined using the "Union" type.
"""

# ╔═╡ e4cd8bbc-d090-4be3-9de2-b97fe096724f
md"""
---
"""

# ╔═╡ 1cb05c5b-1cb4-4b88-ba86-b0675c1228d4
md"""
### Functions
"""

# ╔═╡ 62edc512-89e6-4b29-b96e-f43b253654b9
md"""
In Julia, a function is an object that maps a tuple of argument values to a return value.

There are three syntaxes for defining a function. The first two are named functions and the third is an anonymous function. If the return value is the last statement, then the "`return`" keyword is optional.

Standard function definition:

```julia
function myadd(x::Int, y::Int)
    x + y
end
```

One-line function definition:

```julia
myadd(x::Float64, y::Float64) = x + y
```

Anonymous function definition:

```julia
x, y -> x + y
```

Anonymous functions are often used when a function argument expects a function, e.g., the `filter` method that expects a Boolean comparison function.

Let's define the above three functions.
"""

# ╔═╡ 771dee9c-1615-435a-884f-7d274172191c


# ╔═╡ c0c8fde0-1526-4e8a-896a-67c226b0badf


# ╔═╡ d67ba2a6-6fd1-4fae-86a9-d7ce1e8ef41f


# ╔═╡ c3b1713c-1207-427f-bc2b-7ff973f5e35e
md"""Notice that the function "`myadd`" now has two methods; one for `Int`s and one for `Float64`s.

Try adding an `Int` and `Float64` using the "`myadd`" function.
"""

# ╔═╡ c7b43469-232a-46a0-8bb6-c7a928e6d2f2


# ╔═╡ cc19d021-1f25-4469-8239-9924cc01f883
md"""
The compiler returns a MethodError because their is no method that adds a `Int` and `Float64`. We can fix this by defining a generic "`myadd`" function.

```julia
myadd(x, y) = x + y
```
"""

# ╔═╡ 43b6afe5-8c9d-412a-ae68-c190b93c74e6


# ╔═╡ 02296dd4-ddca-4acb-929f-61ef5d9f755f
md"""
!!! note
    Now look at the result above of adding an `Int` and a `Float64` using "`myadd`".

    In many cases, a function with generic arguments is sufficiently performant. But in those cases where extreme performance is needed, defining methods with specific argument types may be necessary.

!!! note
    One-line functions are usually inlined by the compiler. So, there is usually no performance penalty for using them. Multi-lined functions may also be inlined.
"""

# ╔═╡ 92200db6-62a3-4865-b02c-0ff7b41421bd
md"""
#### Optional Arguments
"""

# ╔═╡ 197727b0-f566-4953-94fd-9062f8d4e828
md"""
Functions can often take sensible default values. Julia allows the default values to be defined in the function definition.
```julia
optargs(y::Int, m::Int=1, d::Int=1) = "$y-$m-$d"
```
Define the above function and execute it with a variable number of arguments.

Note how many methods are created when the function is defined.
"""

# ╔═╡ f5ce11b3-c3d5-4464-84d7-a9d4cfc87953


# ╔═╡ c463427e-1584-4eb7-aefe-0eb24a9c01ba
optargs(2, 2, 2)

# ╔═╡ 47e23641-84e5-423f-9783-626ab426b381
md"""
#### Keyword Arguments
"""

# ╔═╡ 3ddf7fd7-9ebd-4f63-a4ac-c6cea8973478
md"""
Some functions have a large number of arguments or a large number of behaviors. Remembering how to call such functions can be annoying. Keyword arguments can make these functions easier to use and extend by allowing arguments to be identified by name instead of only by position.

Keyword arguments are listed after the required and optional arguments. They are delimited by a semicolon in the argument list.

```julia
kwfunc(arg1, arg2=1; kwd1="blue", kwd2="red")
```

!!! note
    Don't confuse keyword arguments and optional arguments. Optional arguments are positional arguments with default values. Keyword arguments are positionless arguments with default values.
"""

# ╔═╡ 246d64cb-edfb-40a7-aeea-baa62217389f
md"""
#### Functors
"""

# ╔═╡ f997567b-b403-4e21-a87f-063b59dcc5a6
md"""
Functors are anonymous functions that are defined only by their argument signature. They are synonymous with callable objects in Python.

```julia
struct Polynomial{R}
    coeffs::Vector{R}
end

function (p::Polynomial)(x)
    v = p.coeffs[end]
    for i = (length(p.coeffs)-1):-1:1
       v = v*x + p.coeffs[i]
    end
    return v
end
```

```julia
p = Polynomial([1,10,100])

p(5)
```
"""

# ╔═╡ 8a1c438b-185c-422b-b689-3304c41f5b6f
md"""
Define the Polynomial type and the functor by placing the struct and function in a `begin ... end` block.
"""

# ╔═╡ 802d9fbf-8a1c-4bb3-aa2d-cd9bab659115


# ╔═╡ 679a571e-d866-4005-a047-028c426fb167
md"""Create a polynomial"""

# ╔═╡ 1e8b04e8-ea02-41d1-94e1-42b02bbafdcc


# ╔═╡ 3ffc37d1-8fd2-4436-bb8d-4bd82291c174
md"""Evaluate the polynomial"""

# ╔═╡ fad1263d-6a0a-435e-a6b5-2e2d394307be


# ╔═╡ 7a35a96c-be9e-4e6e-ba70-7fb9b84a609f
md"""
!!! note
    What have we learned about functions?

    * Julia uses the argument signature, called multiple dispatch, to select the executable function.
    * Julia has two syntaxes for defining functions: one is for many-line functions and the other for one-line functions.
    * Julia has named functions and anonymous functions.
    * Julia function signatures have arguments and keywords. Arguments are required and listed first, but can have optional default values. Whereas, keywords are listed last and are optional.
    * Julia has anonymous functions called "functors" that are defined by their argument signature.

"""

# ╔═╡ 9bbe3e6c-3b86-4724-8d5a-8fd32228ba76
md"""
---
"""

# ╔═╡ 2010e9a0-4196-4db4-b30d-ba93cda98095
md"""
## Multi-dimensional Arrays
"""

# ╔═╡ 33105044-e651-40a5-b928-592032c68e42
md"""
The array library is implemented almost completely in Julia itself, and derives its performance from the compiler. All arguments to functions are passed by sharing (i.e. by pointers). By convention, a function name ending with a "`!`" indicates that it will mutate or destroy the value of one or more of its arguments (compare, for example, "`sort`" and "`sort!`").

Two characteristics of Julia arrays are:

* Column-major indexing
* One-based indexing

Both column-major indexing and one-base indexing follow the matrix convention of vectors being column arrays and the first index being 1. This is the same as FORTRAN and Matlab, and, of course, unlike Python.

!!! tip
    Just remember that the first index varies fastest.
"""

# ╔═╡ 9f7a9b73-1cdf-4e2a-8dd5-1b2b024a499f
md"""
---
"""

# ╔═╡ b3c2831f-1de1-47f4-ba4a-1cc30c30d510
md"""
### Array Construction and Initialization
"""

# ╔═╡ 86cb7030-f054-4a6f-99ef-3c7bb03686d9
md"""

There several ways to create and initialize a new array:

```julia
Array{T}(undef, dims...)    # an unitialized dense array

ones(T, dims...)            # an array of zeros
```

where `T` signifies the array type, and `dims...` is a list of array dimensions.

```julia
[1, 2, 3]                   # an array literal

[2*i + j for i=1:3, j=4:6]  # array comprehension

(2*i + j for i=1:3, j=4:6)  # generator expression
```

!!! note

    A generator expression doesn't create an array, it produces a value on demand.

"""

# ╔═╡ 911167cc-8009-425a-946b-fac4d1e87a91
md"""
Let's create some arrays. Create:
```julia
zeros(Int8, 2, 3)
```
"""

# ╔═╡ 579259ef-3b67-4497-a8a3-5e6bed5b2ce0


# ╔═╡ 76afc0a5-5da0-446d-afbd-1f202d84cf9a
md"""Create
```julia
zeros(Int8, (2,3))
```
"""

# ╔═╡ c92272d7-8729-468d-8bc5-f80f12a53856


# ╔═╡ e87c1b53-da8e-4747-92ea-b8299b9107b7
md"""

The array dimensons can be either a list or tuple.

Now create an array without the type argument.
```julia
zeros(2, 3)
```
"""

# ╔═╡ 13e6db9b-8b75-4f30-b174-ce3623148169


# ╔═╡ cd46d32e-84e0-4d29-892f-b30db3fdcf8a
md"""The type defaults to `Float64`"""

# ╔═╡ 133702ad-1d34-4974-a612-231564d7806a
collect(1:10)

# ╔═╡ 3d1691fb-12ce-4e9f-9643-d691df69ba99
[i for i=1:10]

# ╔═╡ 5f449433-a95c-46fa-9894-f372420506ac
md"""
#### Indexing
"""

# ╔═╡ d41bcf68-f472-48d0-ad82-1883f1d8d8ae
md"""
Indexes may be a scalar integer, an array of integers, or any other supported index. This includes Colon (`:`) to select all indices within the entire dimension, ranges of the form `begin:end` or `begin:step:end` to select contiguous or strided subsections, and arrays of booleans to select elements at their true indices. Slices in Julia are inclusive, meaning the beginning and ending indices are included in the slice.

`begin` and `end` can be used to indicate the first and last index of a slice. So, `end-1` is the penultimate index.

!!! note
    Julia allows the beginning and ending indices to be any value. That is they can be positive, negative, or zero. For example, the indices can be `-3:3`. This feature requires the OffsetArrays package.
"""

# ╔═╡ f44cd2f8-bc47-45bd-92fd-9b455ee580d5
md"""
One supported index that is commonly used is the "`CartesianIndex`". It is an index that represents a single multi-dimensional index.

```julia
A = reshape(1:32, 4, 4, 2)
A[3, 2, 1]
A[CartesianIndex(3, 2, 1)] == A[3, 2, 1] == 7
```

Try the above example.
"""

# ╔═╡ 8aa1bee5-c3f3-425c-8c33-5fed56866342


# ╔═╡ b0ef4e74-1a13-4877-846f-f49155f3297c
CartesianIndex(3, 2, 1)

# ╔═╡ c76b138f-feb1-41af-9bb2-ad045a3675ac
md"""
An array of `CartesianIndex` is also supported. They help simplify manipulating arrays. For example, it enables accessing the diagonal elements from the first "page" of `A` from above:

```julia
page = A[:,:,1]
```

```julia
page[CartesianIndex(1,1),
     CartesianIndex(2,2),
     CartesianIndex(3,3),
     CartesianIndex(4,4)]
```

Try it.
"""

# ╔═╡ d98ad311-6bf5-4f39-8e92-167fb4eea9a5


# ╔═╡ c5e2c91c-3e4f-434e-9c35-d3e8933d28f4


# ╔═╡ cc9eae6f-4cef-4160-9d1d-08f53e0681f6
md"""
This is expressed more simply using dot broadcasting and combining it with a normal integer index (instead of extracting the first page from `A` as a separate step).

```julia
A[CartesianIndex.(axes(A, 1), axes(A, 2)), 1]
```

Try this too.
"""

# ╔═╡ 2793ca45-024c-4289-8075-c48c02acb971


# ╔═╡ bfe40210-0183-4358-959d-3b2e8e3355e8
md"""
#### Iteration
"""

# ╔═╡ a15b5f47-1be5-42ae-91c0-b868382d7e9b
md"""
The preferred way of iterating over an array is:

```julia
for a in A
    # Do something with the element a
end

for i in eachindex(A)
    # Do something with i and/or A[i]
end
```

The first example returns the value and the second returns the index. These methods work with both dense and sparse arrays.
"""

# ╔═╡ 173531b1-4347-4cd7-97dd-213b449087bb


# ╔═╡ 8a080dae-1a60-4384-ad4a-07cad30485bd
md"""
#### Vectors and Matrices

A vector and matrix are just aliases for one and two dimensional arrays. To perform matrix multiplication, use the matrix multiply operator `*`.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Measurements = "~2.11.0"
PlutoUI = "~0.7.59"
Unitful = "~1.21.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "ec642d6e206cef819e1f92268a248252c17b0393"

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

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

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

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "Requires"]
git-tree-sha1 = "bdcde8ec04ca84aef5b124a17684bf3b302de00e"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.11.0"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

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

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

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
# ╟─1e47e383-b735-4c57-a300-2afe8491b49a
# ╟─d1366a55-b4fc-4ddb-b5c2-5f3381c48b49
# ╟─09193424-25b9-45ce-840f-f24bbcc46c9d
# ╟─54aa12fb-6a56-404d-86e3-4f2a0478319f
# ╟─11643fb1-2898-47cc-a2c7-71a049a5ebd5
# ╟─24e6db88-8aaa-4f51-9d67-9e10c9737e39
# ╟─a2142e30-f585-4e7a-951e-c906dcee175c
# ╟─b1ed2c4e-f5fa-4e5e-87d8-7af6f80a83ca
# ╟─2ed5cce9-989a-4972-a27c-e3ec9cd4257e
# ╟─1eff0f30-698b-4e70-848b-c39c8087cbe1
# ╟─7f3357bc-4103-4a35-af21-9c86f5a0ec2f
# ╟─cb01ca1b-838a-4137-90d6-fc9ecf9b2468
# ╟─53bc6b61-e69c-4c62-9eca-0b649a2486bd
# ╟─7475c896-d1b1-4429-9ba8-8e78de41e0b0
# ╟─5df8264e-6e37-4674-abdf-2b05c530787f
# ╟─c7ce5c61-b358-4feb-808e-f092e43351a5
# ╟─f646ca14-c01e-47ee-8e2b-052d9db0985b
# ╠═4a404280-2845-4deb-8eee-2dcdcb9aed27
# ╟─7813824a-cae9-4b97-ac90-e542fbd630d5
# ╠═6ac51e87-87a2-4ccc-9f08-0028700b3cda
# ╟─27208179-35c3-43c1-9548-3620c8aa7680
# ╠═40d8d18c-3713-4e77-812d-9d77a4e1ac50
# ╟─aa3e9db7-49d1-40f8-b745-6c4faa2197e1
# ╠═756d3f69-a5b2-4ee6-bae6-94c513baae6c
# ╟─d882d5f8-6898-4810-aa4a-7b064cff80b0
# ╟─98fbf4c1-fae3-4cab-abea-9b7746d807a9
# ╟─7f41ba9b-e650-45a6-bc90-b6e119dc60b2
# ╟─419a6dec-1db0-477f-911f-049223b5674f
# ╠═98340265-f51e-47a0-95d2-df179b87f54b
# ╟─b4f40ffa-138a-4451-8f3c-34f491ea6b13
# ╟─8ee7f43d-bf75-4975-ac64-54c2d5a0174a
# ╠═d4368e22-60c6-456a-94a5-56e6dfdb26d7
# ╟─62786c37-744e-4bd9-8fd9-68c35289f5f5
# ╟─d1e9c51c-efb9-4dcb-9d28-8c54a235fbb4
# ╠═35f025bb-7ee1-482b-bdc2-f5448b7ec027
# ╠═948b4044-aee8-4b21-b9a5-4949adb9e01e
# ╟─b27578b2-f5f5-4e46-82e6-0007be187ba6
# ╠═ea23f70b-b9f2-42fa-96de-9eb5ab6d1a0d
# ╟─065265a5-c9ad-4a39-b14d-f4e2e49d3f7a
# ╟─a4392002-6b0c-48bb-beed-1360a54bbafa
# ╟─563f07ad-6aed-495e-85fb-bae4a1755ac2
# ╠═40abc83f-b4bd-479f-8671-189cc712d792
# ╟─297cd86c-5e9d-4f70-b11a-cbae8fa96d1e
# ╠═8f016c75-7768-4418-8c57-100db3073c85
# ╟─094b6f30-cbd6-46b1-8e0c-3fdb1ef18261
# ╠═7ba8dc19-e0ca-40de-a778-7583ca70978d
# ╟─668abc35-fdc3-430f-8c90-de3c2c2cd77b
# ╠═232cc444-03b7-442a-8737-8b7725b43421
# ╟─d2a2d0bc-e883-439f-8e34-166e2369caef
# ╟─ce6545c9-15eb-44d9-a999-770c222e249a
# ╠═88ca2a73-6203-447c-afcc-9e370a82076b
# ╟─c24f1ddd-5e31-4073-a627-86cedb1d44c2
# ╠═63a4b27a-5361-4d95-8787-ae31ca7987fe
# ╠═15674bb0-2fe1-40b1-a6c0-3a5a64a6a5c3
# ╠═70f08712-002c-4adc-84b1-73a8655d8a44
# ╠═3fa06ac2-20f3-4df4-8ce3-b0e1c21ca842
# ╟─cf4a0e8f-9210-4f1e-84d4-ee7ff09aaf61
# ╠═e00b826d-1bbb-4413-a907-eb181369526b
# ╟─b56255c6-9d3b-4e2f-a9a0-c6fe69990f3d
# ╟─aaf91290-98db-4508-9519-18f7f25439dd
# ╟─6809e5b1-86da-460e-bcba-ca379cc54dfe
# ╟─5cd072cb-5d71-4a08-8e41-4eaaa7faaa5c
# ╟─90755b03-8a48-40ea-a43a-24bf01248d4b
# ╟─9f61328e-59df-49b2-8c61-c3f339e4254a
# ╟─f37bc13e-fa91-4166-983b-fd13a8493435
# ╟─915b1ea3-4563-4daf-bf99-25e1de6ef6ca
# ╟─b3fad295-b319-4eba-b250-b83ca7e2d278
# ╟─0d0c11c0-d39f-462c-9fb6-ab90ca98d230
# ╟─a02bbbbb-6b3f-47ef-a11f-1db9b802db6f
# ╠═2262c860-c06c-4293-8e6d-b616228cb301
# ╠═68e64f74-8a6b-403e-a404-52fb9cdea54b
# ╟─0887eca0-6760-4d9b-b44e-d1a14059aede
# ╟─86ce3011-3e52-43e4-83ec-1a4422073d45
# ╟─0ad9aa76-f6c7-4368-8ae4-58daa548e065
# ╠═1bc3da9e-143c-489c-b8de-a29dc48f17cb
# ╟─f00dd72a-8705-426b-9eb4-b91cf1ea95d4
# ╠═d308df6b-14ec-49ec-8270-a3b9efd88517
# ╠═01805f02-f9f6-4e3e-8e93-a0628753130f
# ╟─a90b9011-714e-41d1-b7a3-fb3eb9dc56da
# ╟─60e8eeac-0d2d-4148-a564-8c6c3dbc63ab
# ╠═b8325403-9744-4a9d-ae64-be88671da89b
# ╠═4879dae5-442e-4dc6-90c9-366ff76912bb
# ╠═e2e57f49-f848-468a-a6f5-482b6e1ad4ba
# ╟─4c278c5a-3324-4245-8ddf-f5390167168f
# ╟─aace56e5-c11a-4306-bcb1-d08fda9dcd84
# ╟─3772a828-561d-4600-8e67-49a28cc6cf09
# ╠═aa4a7ec0-a270-482b-abeb-7168de767938
# ╟─b8e3b72a-e501-4164-b06c-cbb3282d9d11
# ╠═d9aa9f5e-31b6-49a3-bae8-a9b149e6ab91
# ╟─15b0159b-9c8c-4327-b73d-d7e19decde2a
# ╠═5d5b1283-043b-437a-afda-75801808acc9
# ╟─6a6b2a0a-6bb6-4a67-b4c1-46631503918d
# ╟─80cfbd64-61f9-46b9-8e55-ffb92e223604
# ╟─877faa74-7490-44a3-9e97-b36b36050796
# ╠═bba18435-d355-4fca-a6f5-10dacde17413
# ╟─d9e911a8-13f9-41e5-ac36-4aee3ec24c59
# ╠═5f72777b-a174-453c-8b18-ebf1f4bebe0d
# ╠═734a4185-4001-410f-affc-71b33e339339
# ╟─c349f7b8-bdf0-4b94-b412-06c5e7f3cbc5
# ╠═d8be9383-fb60-4938-9376-f91d59f21559
# ╟─31dfb05b-ed87-48f9-a74c-0055e46de160
# ╠═d2ada743-b82d-47c8-9b1d-4bd56de76e62
# ╟─32085179-7b01-4a3e-b810-d79b5b7aadc1
# ╟─ea15815e-0ae3-4f22-9dce-a17cb3a0560b
# ╠═be09f5d0-daea-4f47-8dc8-33c875fca843
# ╠═10ec3b0d-1add-4f92-8f4c-b594ab3f0e68
# ╟─6ee4665d-c5b9-4881-ad65-15c6a8229f3f
# ╠═f5596a05-04de-4955-9575-4c035e0f1495
# ╠═a1b4f7bb-8238-40d6-81cb-6d5e6c737134
# ╟─582fda28-df9c-4d4c-b376-78463a1e41bf
# ╟─3b8e773f-df6e-4b59-9f5d-e14366d02754
# ╠═91f35db2-6a17-42aa-8580-1dea220b8c11
# ╟─6d76f5fc-e59e-4cf3-a2ff-aa2c839d4661
# ╟─a631464d-e08a-4a89-8c47-fd5a7b2dee16
# ╠═7a8faa02-34b1-4416-beab-2909fb56c767
# ╠═6e1a3b46-05f0-487d-933a-6ff0d9d43a2b
# ╟─9064310f-d82d-44dd-bd76-90d194ee84ac
# ╟─05adfd23-c809-4706-9bf2-1a0a2445748b
# ╟─8c59d1ec-606a-4665-8e9e-262d34415e85
# ╠═67a4ff9f-c75f-444c-9091-e9b5c17ee773
# ╠═cfd93268-174f-4a7e-9f98-3d5787c9392c
# ╠═a96dd069-09aa-4add-baba-99ffae36bfe8
# ╟─8a3aa0d3-1ade-4961-975d-b39899731ffe
# ╟─e4cd8bbc-d090-4be3-9de2-b97fe096724f
# ╟─1cb05c5b-1cb4-4b88-ba86-b0675c1228d4
# ╟─62edc512-89e6-4b29-b96e-f43b253654b9
# ╠═771dee9c-1615-435a-884f-7d274172191c
# ╠═c0c8fde0-1526-4e8a-896a-67c226b0badf
# ╠═d67ba2a6-6fd1-4fae-86a9-d7ce1e8ef41f
# ╟─c3b1713c-1207-427f-bc2b-7ff973f5e35e
# ╠═c7b43469-232a-46a0-8bb6-c7a928e6d2f2
# ╟─cc19d021-1f25-4469-8239-9924cc01f883
# ╠═43b6afe5-8c9d-412a-ae68-c190b93c74e6
# ╟─02296dd4-ddca-4acb-929f-61ef5d9f755f
# ╟─92200db6-62a3-4865-b02c-0ff7b41421bd
# ╟─197727b0-f566-4953-94fd-9062f8d4e828
# ╠═f5ce11b3-c3d5-4464-84d7-a9d4cfc87953
# ╠═c463427e-1584-4eb7-aefe-0eb24a9c01ba
# ╟─47e23641-84e5-423f-9783-626ab426b381
# ╟─3ddf7fd7-9ebd-4f63-a4ac-c6cea8973478
# ╟─246d64cb-edfb-40a7-aeea-baa62217389f
# ╟─f997567b-b403-4e21-a87f-063b59dcc5a6
# ╟─8a1c438b-185c-422b-b689-3304c41f5b6f
# ╠═802d9fbf-8a1c-4bb3-aa2d-cd9bab659115
# ╟─679a571e-d866-4005-a047-028c426fb167
# ╠═1e8b04e8-ea02-41d1-94e1-42b02bbafdcc
# ╟─3ffc37d1-8fd2-4436-bb8d-4bd82291c174
# ╠═fad1263d-6a0a-435e-a6b5-2e2d394307be
# ╟─7a35a96c-be9e-4e6e-ba70-7fb9b84a609f
# ╟─9bbe3e6c-3b86-4724-8d5a-8fd32228ba76
# ╟─2010e9a0-4196-4db4-b30d-ba93cda98095
# ╟─33105044-e651-40a5-b928-592032c68e42
# ╟─9f7a9b73-1cdf-4e2a-8dd5-1b2b024a499f
# ╟─b3c2831f-1de1-47f4-ba4a-1cc30c30d510
# ╟─86cb7030-f054-4a6f-99ef-3c7bb03686d9
# ╟─911167cc-8009-425a-946b-fac4d1e87a91
# ╠═579259ef-3b67-4497-a8a3-5e6bed5b2ce0
# ╟─76afc0a5-5da0-446d-afbd-1f202d84cf9a
# ╠═c92272d7-8729-468d-8bc5-f80f12a53856
# ╟─e87c1b53-da8e-4747-92ea-b8299b9107b7
# ╠═13e6db9b-8b75-4f30-b174-ce3623148169
# ╟─cd46d32e-84e0-4d29-892f-b30db3fdcf8a
# ╠═133702ad-1d34-4974-a612-231564d7806a
# ╠═3d1691fb-12ce-4e9f-9643-d691df69ba99
# ╟─5f449433-a95c-46fa-9894-f372420506ac
# ╟─d41bcf68-f472-48d0-ad82-1883f1d8d8ae
# ╟─f44cd2f8-bc47-45bd-92fd-9b455ee580d5
# ╠═8aa1bee5-c3f3-425c-8c33-5fed56866342
# ╠═b0ef4e74-1a13-4877-846f-f49155f3297c
# ╟─c76b138f-feb1-41af-9bb2-ad045a3675ac
# ╠═d98ad311-6bf5-4f39-8e92-167fb4eea9a5
# ╠═c5e2c91c-3e4f-434e-9c35-d3e8933d28f4
# ╟─cc9eae6f-4cef-4160-9d1d-08f53e0681f6
# ╠═2793ca45-024c-4289-8075-c48c02acb971
# ╟─bfe40210-0183-4358-959d-3b2e8e3355e8
# ╟─a15b5f47-1be5-42ae-91c0-b868382d7e9b
# ╠═173531b1-4347-4cd7-97dd-213b449087bb
# ╟─8a080dae-1a60-4384-ad4a-07cad30485bd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
