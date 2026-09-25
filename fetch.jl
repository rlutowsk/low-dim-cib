using Redis
using LinearAlgebra
using Oscar
using Nemo
using Base.Threads

# Wymagane importy dla poprawnego rozszerzania metod bez konfliktów
import Oscar: generators, denominator, dim
import Base: sort, show, exp

# --- STRUKTURY DANYCH ---

struct ZClass
    name::String
    dim::Int
    exp::Vector{Int}
    group::MatGroup
    normalizer::MatGroup
    cb::Vector{Vector{ZZMatrix}}
end

dim(z::ZClass) = z.dim
exp(z::ZClass) = z.exp
name(z::ZClass) = z.name
group(z::ZClass) = z.group
normalizer(z::ZClass) = z.normalizer
coboundary_basis_int(z::ZClass) = z.cb

function generators(z::ZClass)
    return map(matrix, generators(z.group))
end

show(io::IO, z::ZClass) = print(io, name(z))
show(io::IO, ::MIME"text/plain", z::ZClass) = print(io, "ZClass(", name(z), ")")

struct AffClass
    name::String
    zclass::ZClass
    gens::Vector{QQMatrix}
    lcm_den::ZZRingElem
end

denominator(a::AffClass) = a.lcm_den
show(io::IO, a::AffClass) = print(io, a.name)
show(io::IO, ::MIME"text/plain", a::AffClass) = print(io, "AffClass(", a.name, ")")

struct ExpandAffClass
    denom::ZZRingElem
    vs::Vector{QQMatrix}
    cb::Vector{Vector{ZZMatrix}}
    rk::Int
    plist
    zclass::ZClass
end

# --- PARSOWANIE I INTEGRACJA Z REDIS ---

function parse_matrices(s::AbstractString)
    v = eval(Meta.parse(s))
    return [matrix(ZZ, permutedims(hcat(M...))) for M in v]
end

function coboundary_basis_int(gg, d::Int)
    id    = identity_matrix(ZZ, d)
    eq    = hcat([g - id for g in gg]...)
    S, T, U = snf_with_transform(eq)
    
    rank  = count(!iszero(S[i, i]) for i in 1:min(nrows(S), ncols(S)))
    long  = inv(U)[1:rank, :]
    len   = length(gg) - 1
    
    return [
        [
            matrix(permutedims(long[j, i*d+1:(i+1)*d])) for i in 0:len
        ]
        for j in 1:rank
    ]
end

function fetch_exponent(conn, name::AbstractString)
    parts = split(name, '.')
    prefix = length(parts) >= 2 ? join(parts[1:2], '.') : name
    
    exp_str = hget(conn, "exponents", prefix)
    isnothing(exp_str) && return [0]
    
    return eval(Meta.parse(exp_str))
end

function fetch_z_class(conn, name::AbstractString)
    dim_str  = hget(conn, name, "dim")
    gens_str = hget(conn, name, "generators")
    norm_str = hget(conn, name, "normalizer")

    dim_val  = parse(Int, dim_str)
    gens     = parse_matrices(gens_str)
    normgens = parse_matrices(norm_str)

    grp = matrix_group(gens...)
    N   = matrix_group(normgens...)

    return ZClass(name, dim_val, fetch_exponent(conn, name), grp, N, coboundary_basis_int(gens, dim_val))
end

function aff_names(conn, z::ZClass)
    members = smembers(conn, name(z) * ":anames")
    return sort(parse.(Int, members))
end

function flat(mat::ZZMatrix)
    r, c = nrows(mat), ncols(mat)
    return Tuple(mat[i, j] for i in 1:r for j in 1:c)
end

function flat(mat::MatGroupElem)
    return flat(matrix(mat))
end

function vector_system(zcl::ZClass, vs_gens, ring=QQ)
    n = dim(zcl)
    I = identity_matrix(ZZ, n)
    z = zero_matrix(ring, 1, n)

    key(A) = flat(A)

    sys = Dict(key(I) => (z, I))
    queue = [(z, I)]

    aff_gens = collect(zip(vs_gens, generators(zcl)))

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

    return sys
end

function fetch_aff_class(conn, z::ZClass, no::Int)
    aff_name = name(z) * "." * string(no)
    vs_str   = hget(conn, aff_name, "vs")
    vs_str   = replace(vs_str, "/" => "//")
    vs       = eval(Meta.parse(vs_str))

    gens = [matrix(QQ, permutedims(v)) for v in vs]
    ld   = reduce(lcm, (denominator(e) for row in gens for e in row); init=ZZ(1))
    return AffClass(aff_name, z, gens, ld)
end

function cb_aff_class(z::ZClass, no::Int)
    return AffClass(name(z) * ".cb" * string(no), z, z.cb[no], ZZ(1))
end

function expand_aff_class(a::AffClass)
    vs_dict = vector_system(a.zclass, a.gens)
    cb_dict = [vector_system(a.zclass, x, ZZ) for x in a.zclass.cb]

    # Sortowanie kluczy dla 100% determinizmu kolejności wektorów
    pg_keys = sort(collect(keys(vs_dict)))

    return ExpandAffClass(
        denominator(a),
        [vs_dict[k][1] for k in pg_keys],
        [[c[k][1] for k in pg_keys] for c in cb_dict],
        length(cb_dict),
        [vs_dict[k][2] for k in pg_keys],
        a.zclass
    )
end

# --- OPTYMALIZOWANE FUNKCJE POMOCNICZE DLA CIB ---

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

@inline function vec_eq_flat(vs_flat::Vector{Int}, i::Int, j::Int, N::Int)
    offset_i = (i - 1) * N
    offset_j = (j - 1) * N
    @inbounds for k in 1:N
        vs_flat[offset_i + k] == vs_flat[offset_j + k] || return false
    end
    return true
end

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

@inline function add_vec_flat!(vs_flat::Vector{Int}, b_pos_flat::Vector{Int}, g::Int)
    @inbounds @simd for i in 1:length(vs_flat)
        v = vs_flat[i] + b_pos_flat[i]
        vs_flat[i] = v >= g ? v - g : v
    end
end

@inline function add_scaled_vec_flat!(vs_flat::Vector{Int}, b_pos_flat::Vector{Int}, coeff::Int, g::Int)
    @inbounds @simd for i in 1:length(vs_flat)
        v = (vs_flat[i] + coeff * b_pos_flat[i]) % g
        vs_flat[i] = v
    end
end

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

function flat_to_rational_matrices(vs_flat::Vector{Int}, N::Int, S::Int, r_dim::Int, c_dim::Int, g::Int)
    res = Vector{QQMatrix}(undef, S)
    for k in 1:S
        offset = (k - 1) * N
        entries = [QQ(vs_flat[offset + j], g) for j in 1:N]
        res[k] = matrix(QQ, r_dim, c_dim, entries)
    end
    return res
end

# --- SEKWENCYJNE I ZRÓWNOLEGLONE OBLICZANIE CIB ---

function cib_by_exp(e::ExpandAffClass, exp::Int)
    g_zz = lcm(exp, e.denom)
    R = ZZModRing(g_zz)
    g = Int(g_zz)

    sample_m = matrix(R, matrix(ZZ, e.vs[1] * g_zz))
    r_dim, c_dim = nrows(sample_m), ncols(sample_m)

    vs_nemo = [matrix(R, matrix(ZZ, row * g_zz)) for row in e.vs]
    b_nemo = [[matrix(R, x) for x in y] for y in e.cb]

    vs_flat, N, S = to_flat_vector(vs_nemo)
    b_flat = [to_flat_vector(y)[1] for y in b_nemo]

    c = zeros(Int, e.rk)
    tmp = zeros(Int, N)
    cib = Vector{Vector{QQMatrix}}()

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

        pos == 0 && break

        c[pos] += 1
        add_vec_flat!(vs_flat, b_flat[pos], g)
    end

    return cib
end

# function run_chunk(k_start::Int, chunk_steps::Int, e_rk::Int, g::Int,
#                    vs_flat_init::Vector{Int}, b_flat::Vector{Vector{Int}},
#                    N::Int, S::Int, r_dim::Int, c_dim::Int)
    
#     c = zeros(Int, e_rk)
#     rem_k = k_start
#     for pos in e_rk:-1:1
#         c[pos] = rem_k % g
#         rem_k = div(rem_k, g)
#     end

#     vs_flat = copy(vs_flat_init)
#     for pos in 1:e_rk
#         if c[pos] > 0
#             add_scaled_vec_flat!(vs_flat, b_flat[pos], c[pos], g)
#         end
#     end

#     tmp = zeros(Int, N)
#     res = Vector{Vector{QQMatrix}}()

#     for step in 1:chunk_steps
#         if is_zmodn_group_flat!(vs_flat, tmp, N, S, g)
#             push!(res, flat_to_rational_matrices(vs_flat, N, S, r_dim, c_dim, g))
#         end

#         step == chunk_steps && break

#         pos = e_rk
#         while pos > 0 && c[pos] == g - 1
#             add_vec_flat!(vs_flat, b_flat[pos], g)
#             c[pos] = 0
#             pos -= 1
#         end

#         if pos > 0
#             c[pos] += 1
#             add_vec_flat!(vs_flat, b_flat[pos], g)
#         end
#     end

#     return res
# end

# function cib_by_exp_par(e::ExpandAffClass, exp::Int)
#     g_zz = lcm(exp, e.denom)
#     R = ZZModRing(g_zz)
#     g = Int(g_zz)

#     sample_m = matrix(R, matrix(ZZ, e.vs[1] * g_zz))
#     r_dim, c_dim = nrows(sample_m), ncols(sample_m)

#     vs_nemo = [matrix(R, matrix(ZZ, row * g_zz)) for row in e.vs]
#     b_nemo = [[matrix(R, x) for x in y] for y in e.cb]

#     vs_flat_init, N, S = to_flat_vector(vs_nemo)
#     b_flat = [to_flat_vector(y)[1] for y in b_nemo]

#     total_combs = g^e.rk

#     chunk_size = 20_000
#     n_chunks = Int(cld(total_combs, chunk_size))

#     chunk_results = Vector{Vector{Vector{QQMatrix}}}(undef, n_chunks)

#     @threads :dynamic for chunk_idx in 1:n_chunks
#         k_start = (chunk_idx - 1) * chunk_size
#         k_end = min(chunk_idx * chunk_size - 1, total_combs - 1)
#         chunk_steps = k_end - k_start + 1

#         chunk_results[chunk_idx] = run_chunk(
#             k_start, chunk_steps, e.rk, g,
#             vs_flat_init, b_flat, N, S, r_dim, c_dim
#         )
#     end

#     res = Vector{Vector{QQMatrix}}()
#     for chunk_res in chunk_results
#         append!(res, chunk_res)
#     end

#     return res
# end

# Zamiast Vector{Vector{QQMatrix}} wracamy po prostu wektory całkowite Vector{Int}
function run_chunk(k_start::Int, chunk_steps::Int, e_rk::Int, g::Int,
                   vs_flat_init::Vector{Int}, b_flat::Vector{Vector{Int}},
                   N::Int, S::Int)
    
    c = zeros(Int, e_rk)
    rem_k = k_start
    for pos in e_rk:-1:1
        c[pos] = rem_k % g
        rem_k = div(rem_k, g)
    end

    vs_flat = copy(vs_flat_init)
    for pos in 1:e_rk
        if c[pos] > 0
            add_scaled_vec_flat!(vs_flat, b_flat[pos], c[pos], g)
        end
    end

    tmp = zeros(Int, N)
    res = Vector{Vector{Int}}()

    for step in 1:chunk_steps
        if is_zmodn_group_flat!(vs_flat, tmp, N, S, g)
            # Zapisujemy czysty wektor Int zamiast konwertować na QQMatrix!
            push!(res, copy(vs_flat))
        end

        step == chunk_steps && break

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

function cib_by_exp_par(e::ExpandAffClass, exp::Int)
    g_zz = lcm(exp, e.denom)
    g = Int(g_zz)

    vs_nemo = [matrix(ZZModRing(g_zz), matrix(ZZ, row * g_zz)) for row in e.vs]
    b_nemo = [[matrix(ZZModRing(g_zz), x) for x in y] for y in e.cb]

    vs_flat_init, N, S = to_flat_vector(vs_nemo)
    b_flat = [to_flat_vector(y)[1] for y in b_nemo]

    total_combs = g^e.rk

    chunk_size = 20_000
    n_chunks = Int(cld(total_combs, chunk_size))

    chunk_results = Vector{Vector{Vector{Int}}}(undef, n_chunks)

    @threads :dynamic for chunk_idx in 1:n_chunks
        k_start = (chunk_idx - 1) * chunk_size
        k_end = min(chunk_idx * chunk_size - 1, total_combs - 1)
        chunk_steps = k_end - k_start + 1

        chunk_results[chunk_idx] = run_chunk(
            k_start, chunk_steps, e.rk, g,
            vs_flat_init, b_flat, N, S
        )
    end

    res = Vector{Vector{Int}}()
    for chunk_res in chunk_results
        append!(res, chunk_res)
    end

    # Zwracamy tuple: (wyniki, S, wymiar_wierszy, wymiar_kolumn, g)
    r_dim, c_dim = 1, N
    return (cib = res, S = S, r_dim = r_dim, c_dim = c_dim, g = g)
end

using JSON

function cib_flat_string(x::NamedTuple)
    d = hasproperty(x, :d) ? x.d : x.c_dim
    n_cib = length(x.cib)

    if isempty(x.cib)
        header = [0, 0, 0, 1]
        return JSON.json(header)
    end

    # 1. Spłaszczenie wszystkich wektorów reprezentantów
    payload = reduce(vcat, x.cib)

    # 2. Obliczenie NWD mianownika g oraz wszystkich elementów z payloadu
    gcd_val = foldl(gcd, payload; init=x.g)

    # 3. Uproszczenie mianownika g oraz wartości w wektorze
    g_reduced = div(x.g, gcd_val)
    payload_reduced = div.(payload, gcd_val)

    # 4. Nagłówek ze zmniejszonym g
    header = [n_cib, x.S, d, g_reduced]

    # 5. Składanie pełnego ciągu
    full_vec = vcat(header, payload_reduced)

    return JSON.json(full_vec)
end

function extract_cib_header(json_str::String)
    # Znajduje pierwsze 4 liczby całkowite w napisie
    matches = collect(first(eachmatch(r"-?\d+", json_str), 4))
    
    length(matches) < 4 && error("Ciąg nie zawiera pełnego 4-elementowego nagłówka")
    
    return [parse(Int, m.match) for m in matches]
end

"""
Transformuje macierz systemów CIB w miejscu dla generatora normalizatora N_mat:
out[j, :] = (in_mat[src_idx[j], :] * N_mat) % g
"""
@inline function transform_system!(out::Matrix{Int}, in_mat::Matrix{Int},
                                   src_idx::Vector{Int}, N_mat::Matrix{Int}, g::Int)
    m, d = size(in_mat)
    @inbounds for j in 1:m
        idx = src_idx[j]
        for c in 1:d
            s = 0
            for k in 1:d
                s += in_mat[idx, k] * N_mat[k, c]
            end
            out[j, c] = mod(s, g)
        end
    end
end

"""
Filtruje zbiór CIB do reprezentantów orbit pod działaniem normalizatora N,
a następnie obcina wynik wyłącznie do wierszy odpowiadających generatorom grupy P.
"""
function cib_reps_generators(e::ExpandAffClass, cib_data::NamedTuple)
    cib_raw = cib_data.cib  # Vector{Vector{Int}} (pełne płaskie wektory)
    g       = cib_data.g    # Int
    m       = cib_data.S    # Int (liczba wszystkich elementów P)
    d       = dim(e.zclass) # Int (wymiar macierzy)

    # 1. Pobieramy macierze wszystkich elementów P w kolejności zgodnej z e.plist
    Plist_mats = [matrix(QQ, matrix(mat)) for mat in e.plist]

    # 2. Rekonstrukcja macierzy m x d do przeszukiwania BFS
    cib_mats = [permutedims(reshape(flat_sys, d, m)) for flat_sys in cib_raw]

    # 3. Wyznaczenie generatorów N oraz ich odwrotności
    N_gens_raw = generators(e.zclass.normalizer)
    N_mats = Matrix{Int}[]
    for n_mat in N_gens_raw
        n_qq = matrix(QQ, matrix(n_mat))
        push!(N_mats, Matrix{Int}(matrix(ZZ, n_qq)))
        push!(N_mats, Matrix{Int}(matrix(ZZ, inv(n_qq))))
    end

    # 4. Precomputacja indeksów sprzężeń: e_i = n * e_j * n^-1
    src_indices = [Vector{Int}(undef, m) for _ in 1:length(N_mats)]
    for (k, n_mat) in enumerate(N_mats)
        n_qq  = matrix(QQ, n_mat)
        n_inv = inv(n_qq)

        for j in 1:m
            ej = Plist_mats[j]
            ei = n_qq * ej * n_inv
            i  = findfirst(==(ei), Plist_mats)
            src_indices[k][j] = i
        end
    end

    # 5. Wyznaczenie indeksów generatorów grupy P na liście e.plist
    P_gens = generators(e.zclass.group)
    P_gen_mats = [matrix(QQ, matrix(g_p)) for g_p in P_gens]
    
    gen_indices = [findfirst(==(gm), Plist_mats) for gm in P_gen_mats]
    filter!(!isnothing, gen_indices)
    if isempty(gen_indices)
        gen_indices = [1] # Awaryjnie dla grupy trywialnej
    end

    # 6. Podział na orbity za pomocą BFS
    cib_set   = Set(cib_mats)
    visited   = Set{Matrix{Int}}()
    reps_full = Vector{Vector{Int}}()

    tmp_sys = Matrix{Int}(undef, m, d)
    for (idx, sys) in enumerate(cib_mats)
        sys in visited && continue

        # Zapamiętujemy reprezentanta pełnej orbity
        push!(reps_full, cib_raw[idx])
        push!(visited, sys)

        queue = [sys]
        while !isempty(queue)
            curr = popfirst!(queue)

            for k in 1:length(N_mats)
                transform_system!(tmp_sys, curr, src_indices[k], N_mats[k], g)

                if tmp_sys in cib_set && !(tmp_sys in visited)
                    nxt = copy(tmp_sys)
                    push!(visited, nxt)
                    push!(queue, nxt)
                end
            end
        end
    end

    # 7. Wyciągnięcie WYŁĄCZNIE generatorów dla każdego reprezentanta
    reps_gens_only = Vector{Vector{Int}}()
    num_gens = length(gen_indices)

    for full_flat in reps_full
        # Długość wynikowego wektora = num_gens * d
        gen_flat = Vector{Int}(undef, num_gens * d)
        
        for (g_idx, row_i) in enumerate(gen_indices)
            # Wycinamy wiersz odpowiadający generatorowi o indeksie row_i w e.plist
            src_start = (row_i - 1) * d + 1
            src_end   = row_i * d
            
            dest_start = (g_idx - 1) * d + 1
            dest_end   = g_idx * d
            
            gen_flat[dest_start:dest_end] .= @view full_flat[src_start:src_end]
        end
        
        push!(reps_gens_only, gen_flat)
    end

    # 8. Zwracamy podsumowanie – 'S' reprezentuje teraz liczbę generatorów
    return (
        cib   = reps_gens_only,
        S     = num_gens,       # Liczba generatorów grupy P
        r_dim = cib_data.r_dim,
        c_dim = cib_data.c_dim, # Wymiar d
        g     = cib_data.g
    )
end