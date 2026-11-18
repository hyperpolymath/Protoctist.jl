# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Empirical checks for the proof obligations in docs/catalogue (PO-nn).
# These are TESTS on randomly generated cases, not proofs. Each line of
# output names the obligation, the number of cases and the number that held.
#
#   julia --project=. docs/catalogue/evidence/obligation_checks.jl
#
# Standard library and Protoctist only.

using Random
using Protoctist
using Protoctist: Formats

const RNG = MersenneTwister(20261118)
report(id, what, held, total) =
    println(rpad(id, 7), rpad(what, 66), held, " / ", total, held == total ? "" : "   <-- NOT ALL")

# ---------------------------------------------------------------- NJ ----

"Symmetric matrix of p-distances k//L from random sequences of length L."
function pdist(n, L, T)
    seqs = [rand(RNG, 1:4, L) for _ = 1:n]
    # make the data tree-ish: each sequence mutates its predecessor's parent
    for i = 2:n
        seqs[i] = copy(seqs[rand(RNG, 1:i-1)])
        for _ = 1:rand(RNG, 0:max(1, L ÷ 12))
            seqs[i][rand(RNG, 1:L)] = rand(RNG, 1:4)
        end
    end
    D = zeros(T, n, n)
    for i = 1:n, j = 1:n
        D[i, j] = T(count(seqs[i] .!= seqs[j])) / T(L)
    end
    return D
end

"Full scan: least (i, j), j < i, among the minimisers of Q."
function argmin_full(D, live)
    r = length(live)
    t = Dict(i => sum(D[i, k] for k in live) for i in live)
    best = nothing
    bi = bj = 0
    for i in live, j in live
        j < i || continue
        q = (r - 2) * D[i, j] - t[i] - t[j]
        if best === nothing || q < best
            best, bi, bj = q, i, j
        end
    end
    return bi, bj, best, t
end

"RapidNJ row rule on sorted rows; returns the same triple and the entries skipped."
function argmin_pruned(D, live)
    r = length(live)
    t = Dict(i => sum(D[i, k] for k in live) for i in live)
    tmax = maximum(values(t))
    best = nothing
    bi = bj = 0
    skipped = 0
    violations = 0
    for i in live
        row = sort([j for j in live if j != i]; by = j -> (D[i, j], j))
        for (pos, j) in enumerate(row)
            qb = (r - 2) * D[i, j] - t[i] - tmax
            if best !== nothing && qb > best          # strict: a tie is never pruned
                skipped += length(row) - pos + 1
                for jj in row[pos:end]                 # PO-08a: bound holds on the suffix
                    (r - 2) * D[i, jj] - t[i] - t[jj] >= qb || (violations += 1)
                end
                break
            end
            q = (r - 2) * D[i, j] - t[i] - t[j]
            a, b = max(i, j), min(i, j)
            if best === nothing || q < best || (q == best && (a, b) < (bi, bj))
                best, bi, bj = q, a, b
            end
        end
    end
    return bi, bj, best, skipped, violations
end

"Neighbour joining; returns the set of nontrivial splits and the tie count."
function nj(D0; pruned = false)
    n = size(D0, 1)
    D = copy(D0)
    live = collect(1:n)
    members = Dict(i => BitSet(i) for i = 1:n)
    splits = Set{BitSet}()
    ties = 0
    skipped = 0
    viol = 0
    agree = true
    while length(live) > 3
        i, j, q, t = argmin_full(D, live)
        r = length(live)
        nt = count(a -> any(b -> b < a && (r - 2) * D[a, b] - t[a] - t[b] == q, live), live)
        nt > 1 && (ties += 1)
        if pruned
            pi, pj, pq, sk, vi = argmin_pruned(D, live)
            skipped += sk
            viol += vi
            agree &= (pi, pj) == (i, j) && pq == q
        end
        for k in live
            (k == i || k == j) && continue
            d = (D[i, k] + D[j, k] - D[i, j]) / 2
            D[i, k] = d
            D[k, i] = d
        end
        union!(members[i], members[j])
        filter!(!=(j), live)
        s = members[i]
        2 <= length(s) <= n - 2 && push!(splits, 1 in s ? setdiff(BitSet(1:n), s) : copy(s))
    end
    return splits, ties, skipped, viol, agree
end

let cases = 0, agree = 0, viol = 0, skipped = 0, tied = 0, joins = 0
    for _ = 1:300
        n = rand(RNG, 6:24)
        D = pdist(n, 240, Rational{BigInt})
        _, ties, sk, vi, ag = nj(D; pruned = true)
        cases += 1
        agree += ag
        viol += vi
        skipped += sk
        tied += ties > 0
        joins += n - 3
    end
    report("PO-08", "bound holds on every skipped entry (violations = $viol), runs", cases - (viol > 0), cases)
    report("PO-08", "pruned search = full scan at every join, exact rationals", agree, cases)
    println("       entries skipped by the row rule: $skipped over $joins joins")
    report("PO-08", "runs on 240-site p-distances with at least one tied join", tied, cases)
end

let cases = 0, same = 0, permsame = 0
    for _ = 1:300
        n = rand(RNG, 8:24)
        D = pdist(n, 240, Rational{BigInt})
        exact, _ = nj(D)
        float, _ = nj(Float64.(D))
        cases += 1
        same += exact == float
        p = randperm(RNG, n)
        ip = invperm(p)
        permuted, _ = nj(D[p, p])
        back = Set(let b = BitSet(p[k] for k in s)
            1 in b ? setdiff(BitSet(1:n), b) : b
        end for s in permuted)
        permsame += back == exact
    end
    report("PO-08", "Float64 tree = exact rational tree (same input order)", same, cases)
    report("PO-08", "exact tree unchanged by permuting the input order", permsame, cases)
end

# ---------------------------------------------------- community formulae ----

hill(p, q) = q == 1 ? exp(-sum(x -> x > 0 ? x * log(x) : 0.0, p)) : sum(x -> x > 0 ? x^q : 0.0, p)^(1 / (1 - q))

let cases = 2000, mono = 0, ends = 0
    for _ = 1:cases
        n = rand(RNG, 1:500, rand(RNG, 2:40))
        p = n ./ sum(n)
        d = [hill(p, q) for q in (0, 0.5, 1, 2, 3)]
        mono += all(d[k] >= d[k+1] - 1e-9 for k = 1:4)
        ends += isapprox(d[1], length(n)) && d[end] >= 1 - 1e-9
    end
    report("PO-10", "Hill numbers non-increasing in q over (0, .5, 1, 2, 3)", mono, cases)
    report("PO-10", "order 0 equals richness, and every order is at least 1", ends, cases)
end

let cases = 2000, ok = 0, mono = 0
    for _ = 1:cases
        n = rand(RNG, 1:30, rand(RNG, 2:25))
        N = sum(n)
        f1, f2 = count(==(1), n), count(==(2), n)
        chao = length(n) + (f2 > 0 ? f1^2 // (2 * f2) : f1 * (f1 - 1) // 2)
        ok += chao >= length(n) && 0 <= 1 - f1 // N <= 1
        es(m) = sum(1 - binomial(big(N - x), big(m)) // binomial(big(N), big(m)) for x in n)
        ms = sort(rand(RNG, 1:N, 4))
        mono += issorted([es(m) for m in ms]) && es(N) == length(n)
    end
    report("PO-21", "Chao1 >= observed richness; coverage in [0, 1]", ok, cases)
    report("PO-23", "expected rarefied richness non-decreasing; equals S at full depth", mono, cases)
end

clr(x) = (l = log.(x); l .- sum(l) / length(l))

let cases = 2000, zero = 0, scale = 0, sub = 0
    for _ = 1:cases
        S = rand(RNG, 3:30)
        x, y = rand(RNG, S) .+ 0.01, rand(RNG, S) .+ 0.01
        zero += abs(sum(clr(x))) < 1e-9
        d = sqrt(sum(abs2, clr(x) .- clr(y)))
        scale += isapprox(d, sqrt(sum(abs2, clr(3.7 .* x) .- clr(0.2 .* y))); atol = 1e-9)
        k = rand(RNG, 2:S-1)
        sub += sqrt(sum(abs2, clr(x[1:k]) .- clr(y[1:k]))) <= d + 1e-9
    end
    report("PO-13", "centred log-ratio sums to zero (tolerance 1e-9)", zero, cases)
    report("PO-14", "Aitchison distance invariant under rescaling each sample", scale, cases)
    report("PO-14", "Aitchison distance of a subcomposition never exceeds the full", sub, cases)
end

function bh(p)
    m = length(p)
    o = sortperm(p)
    q = similar(p)
    run = one(eltype(p))
    for k = m:-1:1
        run = min(run, min(one(eltype(p)), p[o[k]] * m // k))
        q[o[k]] = run
    end
    return q
end

let cases = 2000, ok = 0
    for _ = 1:cases
        p = [rand(RNG, 0:1000) // 1000 for _ = 1:rand(RNG, 1:40)]
        q = bh(p)
        o = sortperm(p)
        ok += all(q .>= p) && all(q .<= 1) && issorted(q[o])
    end
    report("PO-15", "Benjamini-Hochberg: q >= p, q <= 1, order kept (exact rationals)", ok, cases)
end

# ------------------------------------------------------ Protoctist itself ----

const NAMES = ["Eukaryota", "TSAR", "Alveolata", "Rhizaria", "Ciliophora", "Cercozoa", "Dinoflagellata", "Excavata", "Discoba"]
randpath() = join([rand(RNG, NAMES) for _ = 1:rand(RNG, 1:6)], ";")

let cases = 1000, cons = 0, part = 0, mono = 0
    for _ = 1:cases
        tree = build_taxonomy_tree([randpath() for _ = 1:rand(RNG, 1:12)])
        rows = [
            (
                path = TaxPath(randpath()),     # often NOT in the tree: counted at an ancestor
                count = rand(RNG, 0:50),
                status = rand(RNG, (FACTIVE, BELIEF, COLLAPSED, SANS_FIBRE)),
            ) for _ = 1:rand(RNG, 1:30)
        ]
        a = rollup_counts(tree, rows)
        total = sum(r.count for r in rows)
        root = findfirst(n -> n.parent == 0, tree.nodes)
        cons += sum(a.counts) == total && a.cum[root] == total
        part += all(a.f[i] + a.b[i] <= a.cum[i] for i in eachindex(a.cum)) &&
                a.f[root] == sum(r.count for r in rows if r.status == FACTIVE; init = 0)
        mono += all(n.parent == 0 || a.cum[i] <= a.cum[n.parent] for (i, n) in enumerate(tree.nodes))
    end
    report("PO-01", "rollup conserves the total, including rows counted at an ancestor", cons, cases)
    report("PO-01", "a child's cumulative count never exceeds its parent's", mono, cases)
    report("PO-19", "factive + belief <= cumulative; root factive = sum of factive rows", part, cases)
end

let cases = 500, ok = 0, mass = 0
    dir = mktempdir()
    for c = 1:cases
        qs = ["q$(i)" for i = 1:rand(RNG, 1:6)]
        ents = String[]
        for q in qs
            w = rand(RNG, rand(RNG, 1:5))
            w ./= sum(w) * (1 + rand(RNG))            # reported candidates need not sum to 1
            ps = join(["[$(rand(RNG, 0:20)), -$(1000 + 10k), $(w[k]), 0.01, 0.02]" for k in eachindex(w)], ", ")
            push!(ents, "{\"p\": [$ps], \"n\": [\"$q\"]}")
        end
        js = "{\"tree\": \"((A:1{0},B:1{1}):1{2},C:1{3});\", \"placements\": [$(join(ents, ", "))], " *
             "\"fields\": [\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"distal_length\", \"pendant_length\"], \"version\": 3}"
        path = joinpath(dir, "c$c.jplace")
        write(path, js)
        b = Formats.read_jplace(path)
        best = Formats.best_placements(b)
        shuffled = Formats.PlacementBatch(b.tree, shuffle(RNG, b.placements), b.fields, b.version)
        key(v) = sort([(p.name, p.like_weight_ratio) for p in v])
        ok += key(best) == key(Formats.best_placements(shuffled)) &&
              all(p.like_weight_ratio == maximum(x.like_weight_ratio for x in b.placements if x.name == p.name) for p in best)
        mass += all(sum(x.like_weight_ratio for x in b.placements if x.name == q) <= 1 + 1e-9 for q in qs)
    end
    report("PO-06", "best placement is the maximum weight, whatever the listing order", ok, cases)
    report("PO-05", "weights of a query's candidates sum to at most 1 (as generated)", mass, cases)
end

let cases = 1000, sym = 0, refl = 0
    refs = ("PR2", "SILVA")
    vers = ("5.0.0", "5.1.1", "144")
    man() = build_manifest(
        (reference = rand(RNG, refs), ladder = "pr2", version = rand(RNG, vers)),
        (factive = rand(RNG, (0.95, 0.99)), belief = rand(RNG, (0.7, 0.8))),
        Dict("in" => string(rand(RNG, UInt32))),
    )
    for _ = 1:cases
        a, b = man(), man()
        ab, ba = validate_merge(a, b), validate_merge(b, a)
        sym += ab.ok == ba.ok && length(ab.reasons) == length(ba.reasons)
        refl += validate_merge(a, a).ok
    end
    report("PO-07", "merge verdict symmetric in outcome and number of reasons", sym, cases)
    report("PO-07", "a manifest always merges with itself", refl, cases)
end

println("\nJulia ", VERSION, "; seed 20261118; these are tests on generated cases, not proofs.")
