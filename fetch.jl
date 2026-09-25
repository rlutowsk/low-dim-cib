using Redis
using LinearAlgebra
using Oscar

import Oscar: generators, denominator
import Base: sort, show

function parse_matrices(s)
    v = eval(Meta.parse(s))
    [matrix(ZZ,permutedims(hcat(M...))) for M in v]
end

struct ZClass
    name::String
    dim::Int
    exp
    group
    normalizer
    cb
end

dim(z::ZClass) = z.dim
exp(z::ZClass) = z.exp
name(z::ZClass) = z.name
group(z::ZClass) = z.group
normalizer(z::ZClass) = z.normalizer
coboundary_basis_int(z::ZClass) = z.cb

function generators(z::ZClass)
    map(matrix, generators(z.group))
end

show(io::IO, z::ZClass) = print(io, name(z))

show(io::IO, ::MIME"text/plain", z::ZClass) = print(io, "ZClass(", name(z), ")")

function coboundary_basis_int(gg, d::Int)
    id    = identity_matrix(ZZ, d)
    eq    = hcat([g-id for g in gg]...)
    S,T,U = snf_with_transform(eq)
    rank  = sum(!iszero(S[i,i]) for i in 1:min(size(S)...))
    long  = inv(U)[1:rank, :]
    len   = length(gg)-1
    [
        [
            matrix(permutedims(long[j, i*d+1:(i+1)*d])) for i in 0:len
        ]
        for j in 1:rank
    ]
end

function fetch_exponent(conn, name::AbstractString)
    p1 = findnext('.', name, firstindex(name))
    if p1 == nothing
        p1 = 1
    end
    p2 = findnext('.', name, nextind(name, p1))
    if p2 == nothing
        p2 = length(name)+1
    end
    exp = hget(conn, "exponents", SubString(name, firstindex(name), prevind(name, p2)))
    if exp == nothing
        return [0]
    end
    return eval(Meta.parse(exp))
end

function fetch_z_class(conn, name::AbstractString)

    dim_str  = hget(conn, name, "dim")
    gens_str = hget(conn, name, "generators")
    norm_str = hget(conn, name, "normalizer")

    dim = parse(Int,dim_str)
    gens = parse_matrices(gens_str)
    normgens = parse_matrices(norm_str)

    grp = matrix_group(gens...)
    N = matrix_group(normgens...)

    return ZClass(name, dim, fetch_exponent(conn, name), grp, N, coboundary_basis_int(gens, dim))
end

function aff_names(conn, z::ZClass)
    sort(collect(parse(Int,x) for x in smembers(conn, name(z)*":anames")))
end

function flat(mat::MatGroupElem{ZZRingElem, ZZMatrix})
    M = matrix(mat)
    Tuple( M[i,j] for j in 1:ncols(M), i in 1:nrows(M) )
end

function flat(mat::ZZMatrix)
    Tuple( mat[i,j] for j in 1:ncols(mat), i in 1:nrows(mat) )
end

struct AffClass
    name
    zclass
    gens
    lcm_den
end

denominator(a::AffClass) = a.lcm_den

show(io::IO, a::AffClass) = print(io, a.name)

show(io::IO, ::MIME"text/plain", a::AffClass) = print(io, "AffClass(", a.name, ")")

function vector_system(zcl::ZClass, vs_gens, ring=QQ)
    n = dim(zcl)

    I = identity_matrix(ZZ, n)
    z = zero_matrix(ring, 1, n)

    key(A) = flat(A)

    sys = Dict(key(I) => (z, I))
    queue = [(z, I)]

    aff_gens = zip(vs_gens,generators(zcl))

    while !isempty(queue)
        (v, A) = popfirst!(queue)

        for (w, B) in aff_gens
            C = A * B
            u = v * B + w

            k = key(C)

            if !haskey(sys, k)
                sys[k] = (u, C)
                push!(queue, (u, C))
            end
        end
    end

    sys
end

function fetch_aff_class(conn, z::ZClass, no::Int)
    aff_name= name(z)*"."*string(no)
    vs_str  = hget(conn, aff_name, "vs")
    vs_str  = replace(vs_str, "/" => "//")
    vs      = eval(Meta.parse(vs_str))

    gens = [matrix(QQ,permutedims(v)) for v in vs]
    ld   = lcm([denominator(e) for row in gens for e in row])
    AffClass( aff_name, z, gens, ld )
end

function cb_aff_class(z::ZClass, no::Int)
    AffClass(name(z)*".cb"*string(no), z, z.cb[no])
end

struct ExpandAffClass
    denom
    vs
    cb
    rk
end

function expand_aff_class(a::AffClass)

    vs_dict = vector_system(a.zclass, a.gens)
    cb_dict = [ vector_system(a.zclass, x, ZZ) for x in a.zclass.cb ]

    pg_keys = keys( vs_dict )

    ExpandAffClass( denominator(a), [vs_dict[k][1] for k in pg_keys], [ [ c[k][1] for k in pg_keys] for c in cb_dict ], length(cb_dict) )
end

# function get_zmodn_group(list::ZZMatrix, g::ZZRingElem)
#     # l = [ mod(r,g) for r in list ]
#     l = [ mod(list[i:i,:],g) for i in 1:nrows(list)]
#     if ! (0*l[1] in l)
#         return nothing
#     end
#     s = length(l)
#     t = ZZ(2)
#     for i in 1:s
#         if !( mod(-l[i],g) in l )
#             return nothing
#         end
#         if !( mod(t*l[i],g) in l)
#             return nothing
#         end
#         for j in i+1:s
#             if l[i] == l[j]
#                 return nothing
#             end
#             if !( mod(l[i]+l[j],g) in l)
#                 return nothing
#             end
#         end
#     end
#     return l
# end

# function cib_by_exp(e::ExpandAffClass, exp::Int)
#     g = lcm(exp, e.denom)

#     vs = [matrix(ZZ,row*g) for row in e.vs]

#     curr_l = reduce(vcat, [ mod(matrix(r),g) for r in vs])

#     b = [ reduce(vcat,x) for x in e.cb ]

#     c = [ ZZ(0) for i in 1:e.rk]

#     cib = QQMatrix[]
#     while true
#         res = get_zmodn_group( curr_l, g )
#         if res != nothing
#             push!(cib, matrix(QQ,curr_l)/g)
#         end
#         pos = e.rk
#         while pos > 0 && c[pos]==g-1
#             curr_l = curr_l - (g-1)*b[pos]
#             c[pos] = 0
#             pos -= 1
#         end
#         if pos == 0
#             break
#         end

#         c[pos] += 1
#         curr_l = curr_l + b[pos]
#     end
#     return cib
# end

#
# Z mod n approach
#
# function get_zmodn_group(l::Vector{ZZModMatrix})

#     if ! (0*l[1] in l)
#         return nothing
#     end

#     s = length(l)

#     for i in 1:s
#         if !( -l[i] in l )
#             return nothing
#         end

#         if !( l[i]+l[i] in l)
#             return nothing
#         end

#         for j in i+1:s
#             if l[i] == l[j]
#                 return nothing
#             end

#             if !( l[i]+l[j] in l)
#                 return nothing
#             end
#         end
#     end
#     return l
# end

# function cib_by_exp(e::ExpandAffClass, exp::Int)

#     g = lcm(exp, e.denom)
#     R = ZZModRing(g)

#     vs = [matrix(R, matrix(ZZ,row*g)) for row in e.vs]

#     b = [ [ matrix(R,x) for x in y ] for y in e.cb ]

#     c = [ ZZ(0) for i in 1:e.rk]

#     cib = Vector{ZZModMatrix}[]

#     while true
#         res = get_zmodn_group( vs )
#         if res != nothing
#             push!(cib, vs)
#         end

#         pos = e.rk
#         while pos > 0 && c[pos]==g-1
#             vs += b[pos]
#             c[pos] = 0
#             pos -= 1
#         end

#         if pos == 0
#             break
#         end
#         c[pos] += 1
#         vs += b[pos]
#     end
#     return cib
# end 


# działa, ale dla dużych eksponentów - długie

# using Nemo

# # Funkcja testująca – używa bufora tmp, nie alokuje nic na stercie
# function is_zmodn_group!(l::Vector{ZZModMatrix}, tmp::ZZModMatrix)
#     # 1. Sprawdzanie obecności macierzy zerowej bez alokacji
#     any(iszero, l) || return false

#     s = length(l)
#     for i in 1:s
#         # -l[i] w miejscu do bufora tmp
#         Nemo.neg!(tmp, l[i])
#         (tmp in l) || return false

#         # l[i] + l[i] w miejscu
#         Nemo.add!(tmp, l[i], l[i])
#         (tmp in l) || return false

#         for j in (i+1):s
#             # Sprawdzanie unikalności element po elemencie
#             l[i] == l[j] && return false

#             # l[i] + l[j] w miejscu
#             Nemo.add!(tmp, l[i], l[j])
#             (tmp in l) || return false
#         end
#     end
#     return true
# end

# # Pomocnicza funkcja dodająca wektory macierzy w miejscu
# function add_vec!(vs::Vector{ZZModMatrix}, b_pos::Vector{ZZModMatrix})
#     @inbounds for k in 1:length(vs)
#         Nemo.add!(vs[k], vs[k], b_pos[k])
#     end
# end

# function cib_by_exp(e::ExpandAffClass, exp::Int)
#     g = lcm(exp, e.denom)

#     R = ZZModRing(g)
#     vs = [matrix(R, matrix(ZZ, row*g)) for row in e.vs]
#     b = [[matrix(R, x) for x in y] for y in e.cb]

#     c = zeros(Int, e.rk)
#     tmp = zero(vs[1])

#     cib = Vector{ZZModMatrix}[]
#     while true
#         if is_zmodn_group!(vs, tmp)
#             # Używamy deepcopy zamiast copy dla obiektów Nemo
#             push!(cib, [deepcopy(m) for m in vs])
#         end

#         pos = e.rk
#         while pos > 0 && c[pos] == g - 1
#             add_vec!(vs, b[pos])
#             c[pos] = 0
#             pos -= 1
#         end

#         if pos == 0
#             break
#         end

#         c[pos] += 1
#         add_vec!(vs, b[pos])
#     end
#     return cib
# end


# kolejna optymalizacja

# using Nemo

# # Szybkie dodawanie macierzy modulo g w miejscu
# @inline function add_mod!(res::Matrix{Int}, a::Matrix{Int}, b::Matrix{Int}, g::Int)
#     @inbounds for i in 1:length(res)
#         v = a[i] + b[i]
#         res[i] = v >= g ? v - g : v
#     end
# end

# # Szybkie negowanie macierzy modulo g w miejscu
# @inline function neg_mod!(res::Matrix{Int}, a::Matrix{Int}, g::Int)
#     @inbounds for i in 1:length(res)
#         v = a[i]
#         res[i] = v == 0 ? 0 : g - v
#     end
# end

# # Sprawdzanie czy macierz jest zerowa
# @inline function is_zero_matrix(m::Matrix{Int})
#     @inbounds for i in 1:length(m)
#         m[i] == 0 || return false
#     end
#     return true
# end

# # Dodawanie wektora macierzy w miejscu
# @inline function add_vec_native!(vs::Vector{Matrix{Int}}, b_pos::Vector{Matrix{Int}}, g::Int)
#     @inbounds for k in 1:length(vs)
#         add_mod!(vs[k], vs[k], b_pos[k], g)
#     end
# end

# # Weryfikacja grupy na natywnych macierzach
# function is_zmodn_group_native!(l::Vector{Matrix{Int}}, tmp::Matrix{Int}, g::Int)
#     any(is_zero_matrix, l) || return false

#     s = length(l)
#     for i in 1:s
#         neg_mod!(tmp, l[i], g)
#         (tmp in l) || return false

#         add_mod!(tmp, l[i], l[i], g)
#         (tmp in l) || return false

#         for j in (i+1):s
#             l[i] == l[j] && return false

#             add_mod!(tmp, l[i], l[j], g)
#             (tmp in l) || return false
#         end
#     end
#     return true
# end

# # Helper do konwersji ZZModMatrix -> Matrix{Int}
# function to_native_int_matrix(m::ZZModMatrix)
#     r, c = nrows(m), ncols(m)
#     res = Matrix{Int}(undef, r, c)
#     @inbounds for i in 1:r, j in 1:c
#         res[i, j] = Int(lift(m[i, j]))
#     end
#     return res
# end

# function cib_by_exp(e::ExpandAffClass, exp::Int)
#     # Konwersja g na natywny Int
#     g = lcm(exp, e.denom)
#     R = ZZModRing(g)

#     # 1. Konwersja wstępna do natywnych typów Julii
#     vs_nemo = [matrix(R, matrix(ZZ, row * g)) for row in e.vs]
#     b_nemo = [[matrix(R, x) for x in y] for y in e.cb]

#     vs = [to_native_int_matrix(m) for m in vs_nemo]
#     b = [[to_native_int_matrix(x) for x in y] for y in b_nemo]

#     c = zeros(Int, e.rk)
#     tmp = zeros(Int, size(vs[1]))

#     cib = Vector{ZZModMatrix}[]

#     g = Int(g)
#     # 2. Główna pętla w natywnym kodzie maszynowym
#     while true
#         if is_zmodn_group_native!(vs, tmp, g)
#             push!(cib, [matrix(R, m) for m in vs])
#         end

#         pos = e.rk
#         while pos > 0 && c[pos] == g - 1
#             add_vec_native!(vs, b[pos], g)
#             c[pos] = 0
#             pos -= 1
#         end

#         if pos == 0
#             break
#         end

#         c[pos] += 1
#         add_vec_native!(vs, b[pos], g)
#     end

#     return cib
# end


using Nemo

# 1. Czysta pętla Julii (inlined) zamiast ccall(:memcmp) – kompilator wkleja to wprost w kod maszynowy
@inline function fast_in_flat(tmp::Vector{Int}, vs_flat::Vector{Int}, N::Int, S::Int)
    @inbounds for k in 0:(S-1)
        offset = k * N
        match = true
        for i in 1:N
            if tmp[i] != vs_flat[offset + i]
                match = false
                break
            end
        end
        match && return true
    end
    return false
end

# 2. Porównywanie dwóch wektorów wewnątrz ciągłego bloku pamięci vs_flat
@inline function vec_eq_flat(vs_flat::Vector{Int}, i::Int, j::Int, N::Int)
    offset_i = (i - 1) * N
    offset_j = (j - 1) * N
    @inbounds for k in 1:N
        vs_flat[offset_i + k] == vs_flat[offset_j + k] || return false
    end
    return true
end

# 3. Sprawdzanie obecności wektora zerowego
@inline function has_zero_vec_flat(vs_flat::Vector{Int}, N::Int, S::Int)
    @inbounds for k in 0:(S-1)
        offset = k * N
        is_zero = true
        for i in 1:N
            if vs_flat[offset + i] != 0
                is_zero = false
                break
            end
        end
        is_zero && return true
    end
    return false
end

# 4. Operacje arytmetyczne na wycinkach ciągłej pamięci
@inline function neg_mod_slice!(tmp::Vector{Int}, vs_flat::Vector{Int}, i::Int, N::Int, g::Int)
    offset = (i - 1) * N
    @inbounds for k in 1:N
        v = vs_flat[offset + k]
        tmp[k] = v == 0 ? 0 : g - v
    end
end

@inline function add_mod_slice!(tmp::Vector{Int}, vs_flat::Vector{Int}, i::Int, j::Int, N::Int, g::Int)
    offset_i = (i - 1) * N
    offset_j = (j - 1) * N
    @inbounds for k in 1:N
        v = vs_flat[offset_i + k] + vs_flat[offset_j + k]
        tmp[k] = v >= g ? v - g : v
    end
end

# 5. Aktualizacja CAŁEGO zestawu wektorów w pojedynczej pętli SIMD (100% ciągły blok pamięci)
@inline function add_vec_flat!(vs_flat::Vector{Int}, b_pos_flat::Vector{Int}, g::Int)
    @inbounds @simd for i in 1:length(vs_flat)
        v = vs_flat[i] + b_pos_flat[i]
        vs_flat[i] = v >= g ? v - g : v
    end
end

# Szybkie skalarne dodawanie współczynnika (dla inicjalizacji stanu w paczkach)
@inline function add_scaled_vec_flat!(vs_flat::Vector{Int}, b_pos_flat::Vector{Int}, coeff::Int, g::Int)
    @inbounds @simd for i in 1:length(vs_flat)
        v = (vs_flat[i] + coeff * b_pos_flat[i]) % g
        vs_flat[i] = v
    end
end

# Test grupy na całkowicie płaskiej tablicy
function is_zmodn_group_flat!(vs_flat::Vector{Int}, tmp::Vector{Int}, N::Int, S::Int, g::Int)
    has_zero_vec_flat(vs_flat, N, S) || return false

    for i in 1:S
        neg_mod_slice!(tmp, vs_flat, i, N, g)
        fast_in_flat(tmp, vs_flat, N, S) || return false

        add_mod_slice!(tmp, vs_flat, i, i, N, g)
        fast_in_flat(tmp, vs_flat, N, S) || return false

        for j in (i+1):S
            vec_eq_flat(vs_flat, i, j, N) && return false

            add_mod_slice!(tmp, vs_flat, i, j, N, g)
            fast_in_flat(tmp, vs_flat, N, S) || return false
        end
    end
    return true
end

# Funkcje pomocnicze do płaskiej konwersji
function to_native_vec(m::ZZModMatrix)
    r, c = nrows(m), ncols(m)
    res = Vector{Int}(undef, r * c)
    idx = 1
    @inbounds for i in 1:r, j in 1:c
        res[idx] = Int(lift(m[i, j]))
        idx += 1
    end
    return res
end

function to_flat_vector(vs_nemo::Vector{ZZModMatrix})
    vecs = [to_native_vec(m) for m in vs_nemo]
    N = length(vecs[1])
    S = length(vecs)
    res = Vector{Int}(undef, N * S)
    for k in 1:S
        offset = (k - 1) * N
        res[(offset + 1):(offset + N)] .= vecs[k]
    end
    return res, N, S
end

# Funkcja przekształcająca płaski wektor Int na wektor macierzy wymiernych QQMatrix (m / g)
function flat_to_rational_matrices(vs_flat::Vector{Int}, N::Int, S::Int, r_dim::Int, c_dim::Int, g::Int)
    res = Vector{QQMatrix}(undef, S)
    for k in 1:S
        offset = (k - 1) * N
        # Każdy element x konwertujemy na ułamek x / g w ciele QQ
        entries = [QQ(vs_flat[offset + j], g) for j in 1:N]
        res[k] = matrix(QQ, r_dim, c_dim, entries)
    end
    return res
end

function cib_by_exp(e::ExpandAffClass, exp::Int)
    g = lcm(exp, e.denom)
    R = ZZModRing(g)

    sample_m = matrix(R, matrix(ZZ, e.vs[1] * g))
    r_dim, c_dim = nrows(sample_m), ncols(sample_m)

    vs_nemo = [matrix(R, matrix(ZZ, row * g)) for row in e.vs]
    b_nemo = [[matrix(R, x) for x in y] for y in e.cb]

    # Spłaszczenie danych do jednolitych buforów 1D
    vs_flat, N, S = to_flat_vector(vs_nemo)
    b_flat = [to_flat_vector(y)[1] for y in b_nemo]

    c = zeros(Int, e.rk)
    tmp = zeros(Int, N)

    # Inicjalizacja kontenera na wektory macierzy wymiernych (QQMatrix)
    cib = Vector{Vector{QQMatrix}}()

    g = Int(g)
    while true
        if is_zmodn_group_flat!(vs_flat, tmp, N, S, g)
            push!(cib, flat_to_rational_matrices(vs_flat, N, S, r_dim, c_dim, g))
        end

        pos = e.rk
        while pos > 0 && c[pos] == g - 1
            add_vec_flat!(vs_flat, b_flat[pos], g)
            c[pos] = 0
            pos -= 1
        end

        if pos == 0
            break
        end

        c[pos] += 1
        add_vec_flat!(vs_flat, b_flat[pos], g)
    end

    return cib
end

using Base.Threads

function run_chunk(k_start::Int, chunk_steps::Int, e_rk::Int, g::Int,
                   vs_flat_init::Vector{Int}, b_flat::Vector{Vector{Int}},
                   N::Int, S::Int, r_dim::Int, c_dim::Int)
    
    # 1. Odzyskanie stanu licznika c dla indekse k_start
    c = zeros(Int, e_rk)
    rem_k = k_start
    for pos in e_rk:-1:1
        c[pos] = rem_k % g
        rem_k = div(rem_k, g)
    end

    # 2. Inicjalizacja vs_flat dla stanu k_start
    vs_flat = copy(vs_flat_init)
    for pos in 1:e_rk
        if c[pos] > 0
            add_scaled_vec_flat!(vs_flat, b_flat[pos], c[pos], g)
        end
    end

    tmp = zeros(Int, N)
    res = Vector{Vector{QQMatrix}}()

    # 3. Pętla wykonująca dokładnie 'chunk_steps' kroków z natywnym Int
    for step in 1:chunk_steps
        if is_zmodn_group_flat!(vs_flat, tmp, N, S, g)
            push!(res, flat_to_rational_matrices(vs_flat, N, S, r_dim, c_dim, g))
        end

        step == chunk_steps && break

        # Krokowy cykl licznika
        pos = e_rk
        while pos > 0 && c[pos] == g - 1
            add_vec_flat!(vs_flat, b_flat[pos], g)
            c[pos] = 0
            pos -= 1
        end

        if pos > 0
            c[pos] += 1
            add_vec_flat!(vs_flat, b_flat[pos], g)
        end
    end

    return res
end

using Nemo
using Base.Threads

# --- DYNAMICZNIE ZRÓWNOLEGLONA FUNKCJA GŁÓWNA ---

function cib_by_exp_par(e::ExpandAffClass, exp::Int)
    g_zz = lcm(exp, e.denom)
    R = ZZModRing(g_zz)
    g = Int(g_zz)

    sample_m = matrix(R, matrix(ZZ, e.vs[1] * g_zz))
    r_dim, c_dim = nrows(sample_m), ncols(sample_m)

    vs_nemo = [matrix(R, matrix(ZZ, row * g_zz)) for row in e.vs]
    b_nemo = [[matrix(R, x) for x in y] for y in e.cb]

    vs_flat_init, N, S = to_flat_vector(vs_nemo)
    b_flat = [to_flat_vector(y)[1] for y in b_nemo]

    total_combs = g^e.rk

    # Rozmiar mikro-paczki (np. 10 000 - 20 000 kombinacji)
    chunk_size = 20_000
    n_chunks = Int(cld(total_combs, chunk_size))

    # Tablica wyników ma dokładnie n_chunks elementów -> zero problemów z indeksowaniem
    chunk_results = Vector{Vector{Vector{QQMatrix}}}(undef, n_chunks)

    # Dynamiczne rozdzielanie paczek na wątki przez harmonogram Julii
    @threads :dynamic for chunk_idx in 1:n_chunks
        k_start = (chunk_idx - 1) * chunk_size
        k_end = min(chunk_idx * chunk_size - 1, total_combs - 1)
        chunk_steps = k_end - k_start + 1

        chunk_results[chunk_idx] = run_chunk(
            k_start, chunk_steps, e.rk, g,
            vs_flat_init, b_flat, N, S, r_dim, c_dim
        )
    end

    return reduce(vcat, chunk_results; init=Vector{Vector{QQMatrix}}())
end