#!/usr/bin/env python3
"""
Computational tools for the idempotent total graph Gamma(Z_n).

Vertices : the non-zero zero-divisors of Z_n.
Edges    : x ~ y  (x != y)  iff  x + y is an idempotent of Z_n.

The module
  * builds Gamma(Z_n) in O(|V| * 2^m) time (m = number of distinct prime
    divisors of n), using the fact that the neighbours of x are among the
    2^m elements u - x with u idempotent;
  * evaluates the closed formulas of the paper (degree formula, degree
    bounds, edge-count formula) and checks them against the graph;
  * produces the tables and figures used in the paper.

Usage
  python idempotent_total_graph.py info 30          # vertices, degrees, structure
  python idempotent_total_graph.py draw 77 -o g.pdf # draw Gamma(Z_77)
  python idempotent_total_graph.py verify 3000      # check all formulas for n <= 3000
  python idempotent_total_graph.py paper            # regenerate all tables/figures

Running without the console ("Run Code" / Run button in VS Code, PyCharm,
Spyder, IDLE ...): edit the SETTINGS block below and run the file.

Requirements: Python >= 3.8, networkx, matplotlib
  (install once with:  pip install networkx matplotlib)
"""
from __future__ import annotations

# ==========================================================================
#  SETTINGS  --  change these values and press "Run Code" / Run.
#  (Only used when the file is run without command-line arguments.)
# ==========================================================================

TASK = "info_and_draw"
#   "info_and_draw" : print vertices/degrees of Gamma(Z_N) and draw it
#   "info"          : only print vertices/degrees
#   "draw"          : only draw
#   "verify"        : check the degree/edge formulas for all 2 <= n <= VERIFY_UP_TO
#   "paper"         : regenerate all figures and tables of the paper

N = 30                  # the ring Z_N for "info" / "draw"
SAVE_FIGURE_AS = None   # e.g. "gamma_Z30.pdf" -> saved next to this file;
                        # None -> the figure opens in a window
VERIFY_UP_TO = 3000     # upper bound for "verify" and "paper"

# Limits for large N (the graph itself is built quickly even for N = 100000):
MAX_DRAW_VERTICES = 300   # larger graphs: a degree histogram is drawn instead
MAX_TABLE_ROWS = 60       # "info" prints at most this many vertex rows
MAX_DIAMETER_VERTICES = 3000  # the diameter (one BFS per vertex) only up to this size

# ==========================================================================

import argparse
import sys
from collections import defaultdict
from math import gcd, prod
from pathlib import Path

try:
    import networkx as nx
except ImportError:  # friendly message instead of a traceback
    sys.exit("networkx is not installed.  Install it once with:  pip install networkx matplotlib")

HERE = Path(__file__).resolve().parent  # outputs go next to this file

# --------------------------------------------------------------------------
# Elementary number theory
# --------------------------------------------------------------------------


def factorize(n: int) -> dict[int, int]:
    """Prime factorisation {p: alpha} of n >= 1 (trial division)."""
    f: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def euler_phi(n: int) -> int:
    result = n
    for p in factorize(n):
        result = result // p * (p - 1)
    return result


def idempotents(n: int) -> list[int]:
    """All idempotents of Z_n, obtained from the Chinese Remainder Theorem.

    For every subset S of the prime-power factors q_i = p_i^alpha_i, the
    element e_S = sum_{i in S} (n/q_i) * ((n/q_i)^{-1} mod q_i) satisfies
    e_S = 1 (mod q_i) for i in S and e_S = 0 (mod q_i) otherwise.
    This gives exactly 2^m idempotents, without scanning all of Z_n.
    """
    q = [p**a for p, a in factorize(n).items()]
    basis = [(n // qi) * pow(n // qi, -1, qi) % n for qi in q]
    result = []
    for mask in range(1 << len(q)):
        result.append(sum(b for i, b in enumerate(basis) if mask >> i & 1) % n)
    return sorted(result)


def is_zero_divisor(x: int, n: int) -> bool:
    """x is a zero-divisor of Z_n  iff  gcd(x, n) > 1  (0 counts as one)."""
    return gcd(x, n) > 1


def zero_divisors(n: int) -> list[int]:
    """Non-zero zero-divisors of Z_n, i.e. the vertex set of Gamma(Z_n)."""
    return [x for x in range(1, n) if gcd(x, n) > 1]


# --------------------------------------------------------------------------
# The graph
# --------------------------------------------------------------------------


def neighbours(x: int, n: int, idem: list[int] | None = None) -> set[int]:
    """Neighbourhood of the vertex x: {u - x : u idempotent} minus 0, x, units."""
    idem = idempotents(n) if idem is None else idem
    out = set()
    for u in idem:
        y = (u - x) % n
        if y != 0 and y != x and gcd(y, n) > 1:
            out.add(y)
    return out


def idempotent_total_graph(n: int) -> nx.Graph:
    """Gamma(Z_n) as a networkx graph."""
    idem = idempotents(n)
    g = nx.Graph(n=n)
    V = zero_divisors(n)
    g.add_nodes_from(V)
    for x in V:
        for y in neighbours(x, n, idem):
            if x < y:
                g.add_edge(x, y)
    return g


# --------------------------------------------------------------------------
# Closed formulas from the paper
# --------------------------------------------------------------------------


def r_value(x: int, n: int) -> int:
    """r(x) = #{ i : p_i divides x(x-1) }."""
    return sum(1 for p in factorize(n) if (x * (x - 1)) % p == 0)


def degree_formula(x: int, n: int) -> int:
    """deg(x) = 2^m - 2^(m - r(x)) - [x in I(Z_n)] - [2x in I(Z_n)]  (Theorem 3.1)."""
    f = factorize(n)
    m = len(f)
    r = sum(1 for p in f if (x * (x - 1)) % p == 0)
    I = set(idempotents(n))
    return 2**m - 2 ** (m - r) - (x in I) - ((2 * x) % n in I)


def edge_count_formula(n: int) -> int:
    """|E(Gamma(Z_n))| (Theorem 3.5)."""
    f = factorize(n)
    m = len(f)
    U = prod(p ** (a - 1) * (2 * p - 3) for p, a in f.items())
    twice = 2**m * (n - 2 * euler_phi(n) - 3) + U + 4 + (1 if n % 4 else 0)
    assert twice % 2 == 0
    return twice // 2


def predicted_degree_range(n: int) -> tuple[int, int] | None:
    """(min degree, max degree) predicted by Corollary 3.2; None if Gamma is empty."""
    f = factorize(n)
    m = len(f)
    if n == 4:
        return (0, 0)
    squarefree = all(a == 1 for a in f.values())
    if m == 1:
        if squarefree:  # n = p: no vertices
            return None
        return (0 if 2 in f else 1, 1)
    return (2 ** (m - 1) - 1, 2**m - 2 if squarefree else 2**m - 1)


# --------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------


def verify(N: int, verbose: bool = True) -> dict:
    """Check Theorem 3.1, Corollary 3.2 and Theorem 3.5 for every 2 <= n <= N.

    Degrees are computed directly from the definition (neighbour sets),
    independently of the formulas being tested.
    """
    stats = defaultdict(lambda: {"n": 0, "vertices": 0, "min": None, "max": None,
                                 "pred_min": None, "pred_max": None})
    failures = []
    for n in range(2, N + 1):
        f = factorize(n)
        m = len(f)
        idem = idempotents(n)
        I = set(idem)
        V = zero_divisors(n)
        degs = []
        for x in V:
            d = len(neighbours(x, n, idem))
            degs.append(d)
            r = sum(1 for p in f if (x * (x - 1)) % p == 0)
            pred = 2**m - 2 ** (m - r) - (x in I) - ((2 * x) % n in I)
            if d != pred:
                failures.append(("degree", n, x, d, pred))
        if sum(degs) // 2 != edge_count_formula(n):
            failures.append(("edges", n, sum(degs) // 2, edge_count_formula(n)))
        if V:
            if (min(degs), max(degs)) != predicted_degree_range(n):
                failures.append(("range", n, (min(degs), max(degs)), predicted_degree_range(n)))
            s = stats[m]
            s["n"] += 1
            s["vertices"] += len(V)
            s["min"] = min(degs) if s["min"] is None else min(s["min"], min(degs))
            s["max"] = max(degs) if s["max"] is None else max(s["max"], max(degs))
            lo, hi = predicted_degree_range(n)
            s["pred_min"] = lo if s["pred_min"] is None else min(s["pred_min"], lo)
            s["pred_max"] = hi if s["pred_max"] is None else max(s["pred_max"], hi)
    if verbose:
        print(f"Checked 2 <= n <= {N}")
        for m in sorted(stats):
            s = stats[m]
            print(f"  m={m}: {s['n']:5d} rings, {s['vertices']:8d} vertices, "
                  f"degrees in [{s['min']}, {s['max']}]  "
                  f"(predicted by Cor. 3.2: [{s['pred_min']}, {s['pred_max']}])")
        print("  failures:", len(failures))
        for fl in failures[:10]:
            print("   ", fl)
    return {"stats": dict(stats), "failures": failures}


def connectivity_check(N: int) -> tuple[int, list[int]]:
    """Return (#n checked, list of n <= N with m >= 2 whose graph is disconnected)."""
    checked, bad = 0, []
    for n in range(6, N + 1):
        if len(factorize(n)) < 2:
            continue
        checked += 1
        if not nx.is_connected(idempotent_total_graph(n)):
            bad.append(n)
    return checked, bad


# --------------------------------------------------------------------------
# Structural summary
# --------------------------------------------------------------------------


def factor_tex(n: int) -> str:
    return r" \cdot ".join(f"{p}^{{{a}}}" if a > 1 else f"{p}" for p, a in factorize(n).items())


def summary(n: int, g: nx.Graph | None = None) -> dict:
    g = idempotent_total_graph(n) if g is None else g
    degs = [d for _, d in g.degree()]
    comps = nx.number_connected_components(g)
    small = g.number_of_nodes() <= MAX_DIAMETER_VERTICES
    return {
        "n": n,
        "factorization": factor_tex(n),
        "m": len(factorize(n)),
        "V": g.number_of_nodes(),
        "E": g.number_of_edges(),
        "E_formula": edge_count_formula(n),
        "min_deg": min(degs),
        "max_deg": max(degs),
        "components": comps,
        "diameter": nx.diameter(g) if comps == 1 and small else None,
        "diameter_skipped": comps == 1 and not small,
        "is_path": comps == 1 and g.number_of_edges() == g.number_of_nodes() - 1 and max(degs) <= 2,
    }


def print_info(n: int) -> None:
    g = idempotent_total_graph(n)
    I = idempotents(n)
    Iset = set(I)
    f = factorize(n)
    m = len(f)
    print(f"Gamma(Z_{n}),  n = {factor_tex(n).replace(' ', '').replace(chr(92) + 'cdot', '*').replace('{', '').replace('}', '')}")
    print(f"  I(Z_{n}) = {I}")
    print(f"  |V| = {g.number_of_nodes()},  |E| = {g.number_of_edges()} "
          f"(formula: {edge_count_formula(n)})")

    def formula(x):  # Theorem 3.1 (factorization and idempotents computed once)
        r = sum(1 for p in f if (x * (x - 1)) % p == 0)
        return 2**m - 2 ** (m - r) - (x in Iset) - ((2 * x) % n in Iset)

    nodes = sorted(g.nodes)
    print(f"  {'x':>6} {'deg':>4} {'formula':>7}  neighbours")
    for x in nodes[:MAX_TABLE_ROWS]:
        print(f"  {x:>6} {g.degree(x):>4} {formula(x):>7}  {sorted(g[x])}")
    if len(nodes) > MAX_TABLE_ROWS:
        print(f"  ... ({len(nodes) - MAX_TABLE_ROWS} more vertices not printed)")
    wrong = [x for x in nodes if g.degree(x) != formula(x)]
    print(f"  Theorem 3.1 agrees with the graph for {len(nodes) - len(wrong)} of {len(nodes)} vertices")
    hist = defaultdict(int)
    for x in nodes:
        hist[g.degree(x)] += 1
    print("  degree distribution: " + ", ".join(f"deg {d}: {c}" for d, c in sorted(hist.items())))
    s = summary(n, g)
    diam = "not computed (graph too large)" if s["diameter_skipped"] else s["diameter"]
    print(f"  components = {s['components']}, diameter = {diam}, path = {s['is_path']}")


# --------------------------------------------------------------------------
# Drawing
# --------------------------------------------------------------------------


def _path_layout(g: nx.Graph, per_row: int = 8) -> dict:
    """Serpentine layout for a path graph."""
    ends = [v for v in g if g.degree(v) == 1]
    order = nx.shortest_path(g, ends[0], ends[1]) if ends else list(g)
    pos = {}
    for i, v in enumerate(order):
        row, col = divmod(i, per_row)
        if row % 2:
            col = per_row - 1 - col
        pos[v] = (col, -row * 0.9)
    return pos


def _component_layout(g: nx.Graph) -> dict:
    """Lay out components side by side (small components in a row)."""
    pos, x0 = {}, 0.0
    comps = sorted(nx.connected_components(g), key=lambda c: (-len(c), min(c)))
    for comp in comps:
        h = g.subgraph(comp)
        if len(comp) == 1:
            p = {next(iter(comp)): (0.0, 0.4)}
        elif len(comp) == 2:
            a, b = sorted(comp)
            p = {a: (0.0, 0.8), b: (0.0, 0.0)}
        else:
            p = nx.kamada_kawai_layout(h)
        xs = [c[0] for c in p.values()]
        shift = x0 - min(xs)
        for v, (x, y) in p.items():
            pos[v] = (x + shift, y)
        x0 = max(x + shift for x, _ in p.values()) + 0.8
    return pos


def layout(g: nx.Graph) -> dict:
    degs = [d for _, d in g.degree()] or [0]
    if nx.is_connected(g) and max(degs) <= 2 and g.number_of_edges() == g.number_of_nodes() - 1:
        return _path_layout(g, per_row=min(8, g.number_of_nodes()))
    if not nx.is_connected(g):
        return _component_layout(g)
    g.graph["layout"] = "kk"
    return _separate(nx.kamada_kawai_layout(g), min_dist=0.16)


def _separate(pos: dict, min_dist: float, rounds: int = 200) -> dict:
    """Push apart vertices closer than min_dist (keeps labels legible)."""
    import numpy as np
    keys = list(pos)
    P = np.array([pos[k] for k in keys], dtype=float)
    P = (P - P.min(0)) / max(np.ptp(P, 0).max(), 1e-9)
    for _ in range(rounds):
        moved = False
        for i in range(len(P)):
            d = P - P[i]
            dist = np.hypot(d[:, 0], d[:, 1])
            close = (dist < min_dist) & (dist > 0)
            if close.any():
                moved = True
                push = (d[close] / dist[close, None]) * (min_dist - dist[close, None]) / 2
                P[close] += push
                P[i] -= push.sum(0)
        if not moved:
            break
    return {k: tuple(p) for k, p in zip(keys, P)}


def draw(n: int, filename: str | Path | None = None, ax=None, title: bool = True,
         node_size: int = 480, font_size: int = 9, legend: bool = True,
         crt_labels: bool = False):
    """Draw Gamma(Z_n) in a print-friendly (greyscale) style.

    Non-trivial idempotents are shaded grey; vertices x with 2x idempotent
    (and x not idempotent) are drawn as squares - these are the vertices that
    lose a neighbour in Theorem 3.1.  Everything else is a white circle.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    g = idempotent_total_graph(n)
    if g.number_of_nodes() > MAX_DRAW_VERTICES:
        return draw_degree_histogram(n, g, filename, ax)
    I = set(idempotents(n))
    pos = layout(g)
    own = ax is None
    if own:
        xs = [p[0] for p in pos.values()]
        ys = [p[1] for p in pos.values()]
        span_x, span_y = max(xs) - min(xs), max(ys) - min(ys)
        if g.graph.get("layout") == "kk":  # generic 2-D layout: size by vertex count
            w = max(4.0, min(6.3, 1.1 * len(g) ** 0.5 + 0.8))
            hgt = w * max(0.6, min(1.0, span_y / max(span_x, 1e-9)))
        else:  # grid-like layouts (paths, disjoint components)
            w = max(3.2, min(6.3, 0.6 * (span_x + 1.2)))
            hgt = max(1.6, min(4.0, 0.75 * (span_y + 1.2)))
        fig, ax = plt.subplots(figsize=(w, hgt))

    def fill(v):
        return "#bdbdbd" if v in I else "white"

    def halving(v):  # 2v idempotent but v not idempotent
        return (2 * v) % n in I and v not in I

    nx.draw_networkx_edges(g, pos, ax=ax, width=0.9, edge_color="#222222")
    for shape, sel in (("o", lambda v: not halving(v)), ("s", halving)):
        nodes = [v for v in g.nodes if sel(v)]
        if nodes:
            nx.draw_networkx_nodes(g, pos, nodelist=nodes, ax=ax, node_shape=shape,
                                   node_size=node_size if shape == "o" else int(node_size * 0.85),
                                   node_color=[fill(v) for v in nodes],
                                   edgecolors="black", linewidths=0.9).set_clip_on(False)
    labels = None
    if crt_labels:  # label x by its CRT coordinates (x mod q_1, ..., x mod q_m)
        qs = [p**a for p, a in factorize(n).items()]
        labels = {v: "(" + ",".join(str(v % q) for q in qs) + ")" for v in g.nodes}
    for t in nx.draw_networkx_labels(g, pos, labels=labels, ax=ax, font_size=font_size,
                                     font_family="serif").values():
        t.set_clip_on(False)
    ax.set_axis_off()
    ax.margins(0.12)
    if title:
        ax.set_title(rf"$\Gamma(\mathbb{{Z}}_{{{n}}})$", fontsize=11)

    handles = []
    if any(v in I for v in g):
        handles.append(Line2D([], [], marker="o", ls="", ms=8, mfc="#bdbdbd", mec="black",
                              label=r"$x\in I(\mathbb{Z}_n)$"))
    if any(halving(v) for v in g):
        handles.append(Line2D([], [], marker="s", ls="", ms=7, mfc="white", mec="black",
                              label=r"$x\notin I(\mathbb{Z}_n),\ 2x\in I(\mathbb{Z}_n)$"))
    if handles and legend:
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.02),
                  ncol=len(handles), frameon=False, fontsize=9, handletextpad=0.3)
    if own and filename:
        fig.tight_layout()
        fig.savefig(filename, bbox_inches="tight")
        plt.close(fig)
    return ax


def draw_degree_histogram(n: int, g: nx.Graph | None = None,
                          filename: str | Path | None = None, ax=None):
    """For graphs too large to draw: bar chart of the vertex degrees."""
    import matplotlib.pyplot as plt

    g = idempotent_total_graph(n) if g is None else g
    print(f"\nGamma(Z_{n}) has {g.number_of_nodes()} vertices: too many to draw legibly "
          f"(limit MAX_DRAW_VERTICES = {MAX_DRAW_VERTICES}).\n"
          f"The degree distribution is drawn instead.")
    hist = defaultdict(int)
    for _, d in g.degree():
        hist[d] += 1
    degs = sorted(hist)
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(5.5, 3.4))
    bars = ax.bar([str(d) for d in degs], [hist[d] for d in degs],
                  color="#bdbdbd", edgecolor="black", linewidth=0.8)
    ax.bar_label(bars, fontsize=9)
    ax.set_xlabel("degree")
    ax.set_ylabel("number of vertices")
    lo, hi = predicted_degree_range(n) or (None, None)
    bounds = f",  Cor. 3.2: [{lo}, {hi}]" if lo is not None else ""
    ax.set_title(rf"Degrees in $\Gamma(\mathbb{{Z}}_{{{n}}})$:  $|V|={g.number_of_nodes()}$, "
                 rf"$|E|={g.number_of_edges()}${bounds}", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    if own and filename:
        fig.tight_layout()
        fig.savefig(filename, bbox_inches="tight")
        plt.close(fig)
    return ax


# --------------------------------------------------------------------------
# LaTeX tables for the paper
# --------------------------------------------------------------------------


def latex_structure_table(ns: list[int]) -> str:
    rows = []
    for n in ns:
        s = summary(n)
        diam = "--" if s["diameter"] is None else str(s["diameter"])
        rows.append(
            rf"{n} & ${s['factorization']}$ & {s['m']} & {s['V']} & {s['E']} & {s['E_formula']} & "
            rf"{s['min_deg']} & {s['max_deg']} & {s['components']} & {diam} \\"
        )
    return "\n".join([
        r"\begin{tabular}{rlcrrrcccc}",
        r"\toprule",
        r"$n$ & factorization & $m$ & $|V|$ & $|E|$ & $|E|$ (Thm.~\ref{thm:edges}) & $\delta$ & $\Delta$ & comp. & diam. \\",
        r"\midrule",
        *rows,
        r"\bottomrule",
        r"\end{tabular}",
    ])


def latex_verification_table(result: dict, N: int) -> str:
    rows = []
    for m in sorted(result["stats"]):
        s = result["stats"][m]
        lo = 0 if m == 1 else 2 ** (m - 1) - 1
        rows.append(rf"{m} & {s['n']} & {s['vertices']} & $[{s['min']},{s['max']}]$ & $[{lo},{2**m - 1}]$ \\")
    total = sum(s["vertices"] for s in result["stats"].values())
    return "\n".join([
        r"\begin{tabular}{crrcc}",
        r"\toprule",
        r"$m$ & rings $\mathbb{Z}_n$ & vertices checked & observed degrees & bounds of Cor.~\ref{cor:bounds} \\",
        r"\midrule",
        *rows,
        r"\midrule",
        rf"\multicolumn{{5}}{{l}}{{All $2\le n\le {N}$: {total} vertices, "
        rf"{len(result['failures'])} disagreements with Theorems~\ref{{thm:degree}} and \ref{{thm:edges}}.}} \\",
        r"\bottomrule",
        r"\end{tabular}",
    ])


def paper(outdir: str | Path = ".", N: int = 3000) -> None:
    """Regenerate every figure and table used in the revised manuscript."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})

    out = Path(outdir)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)

    for ext in ("pdf", "png"):  # Gamma(Z_20) with CRT labels (Z_4 x Z_5)
        draw(20, out / "figures" / f"gamma_Z20_crt.{ext}", title=False, node_size=900,
             font_size=10, legend=False, crt_labels=True)
    for n in (32, 36, 77, 20, 30):
        # Z_20 and Z_36 are printed side by side: larger labels, legend in the text
        big = {"node_size": 620, "font_size": 11, "legend": False} if n in (20, 36) else {}
        for ext in ("pdf", "png"):
            draw(n, out / "figures" / f"gamma_Z{n}.{ext}", title=False, **big)

    res = verify(N)
    (out / "tables" / "verification.tex").write_text(latex_verification_table(res, N) + "\n")
    ns = [16, 27, 32, 81, 12, 20, 36, 72, 100, 15, 35, 77, 30, 60, 105, 180, 210]
    (out / "tables" / "structure.tex").write_text(latex_structure_table(ns) + "\n")
    checked, bad = connectivity_check(min(N, 2000))
    (out / "tables" / "connectivity.txt").write_text(
        f"n <= {min(N, 2000)} with m >= 2: {checked} graphs checked, disconnected: {bad}\n")
    print(f"connectivity: {checked} graphs with m>=2 checked, disconnected: {bad}")


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------


def run_from_settings() -> None:
    """Entry point for "Run Code": uses the SETTINGS block at the top."""
    task = TASK.strip().lower()
    if task in ("info", "info_and_draw"):
        print_info(N)
    if task in ("draw", "info_and_draw"):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            sys.exit("matplotlib is not installed.  Install it once with:  pip install matplotlib")
        plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})
        if SAVE_FIGURE_AS:
            target = HERE / SAVE_FIGURE_AS
            draw(N, target)
            print(f"\nFigure saved to: {target}")
        else:
            draw(N)
            plt.show()
    elif task == "verify":
        verify(VERIFY_UP_TO)
    elif task == "paper":
        paper(HERE, VERIFY_UP_TO)
        print(f"\nFigures and tables written to: {HERE / 'figures'} and {HERE / 'tables'}")
    elif task not in ("info", "info_and_draw"):
        sys.exit(f"Unknown TASK = {TASK!r}.  Use one of: info_and_draw, info, draw, verify, paper")


def main() -> None:
    if len(sys.argv) == 1:  # started with "Run Code" / Run button
        run_from_settings()
        return
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("info", help="print vertices, degrees and structure of Gamma(Z_n)")
    a.add_argument("n", type=int)
    a = sub.add_parser("draw", help="draw Gamma(Z_n)")
    a.add_argument("n", type=int)
    a.add_argument("-o", "--output", default=None, help="output file (pdf/png); shows window if omitted")
    a = sub.add_parser("verify", help="check the closed formulas for all n <= N")
    a.add_argument("N", type=int)
    a = sub.add_parser("paper", help="regenerate all figures and tables")
    a.add_argument("--out", default=".")
    a.add_argument("--N", type=int, default=3000)
    args = ap.parse_args()

    if args.cmd == "info":
        print_info(args.n)
    elif args.cmd == "draw":
        import matplotlib.pyplot as plt
        if args.output:
            draw(args.n, args.output)
        else:
            draw(args.n)
            plt.show()
    elif args.cmd == "verify":
        verify(args.N)
    elif args.cmd == "paper":
        paper(args.out, args.N)


if __name__ == "__main__":
    main()
