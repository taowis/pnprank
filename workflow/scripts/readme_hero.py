"""Animated README banner (docs/assets/readme_hero.svg).

A myelinated axon doubles as the time axis (nodes of Ranvier at the cutoff years). An
action potential jumps from node to node; as it passes, each tier-1 gene appears at the
position PNPrank gave it genome-wide in the year it was reported
(results/backtest/interval_simulation.tsv, column pct_full). The terminal card shows
real `pnprank` commands with output copied from the released ranking (2026-11-03).

Pure SVG + CSS, no script, so GitHub renders it in an <img>. Light and dark colour
schemes; with prefers-reduced-motion the final frame is shown.

Usage: python workflow/scripts/readme_hero.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from xml.sax.saxutils import escape

OUT = Path("docs/assets/readme_hero.svg")
rows = list(csv.DictReader(open("results/backtest/interval_simulation.tsv"), delimiter="\t"))
genes = [(r["hgnc_symbol"], int(r["year"]), float(r["pct_full"])) for r in rows]
n_top1 = sum(p >= 0.99 for _, _, p in genes)

W, H = 1280, 480
T = 16.0                       # seconds per chart cycle
X0, X1 = 664, 1224             # plot x range
YT, YB = 84, 352               # plot y range (top of genome at YT)
AXON_Y = 404
Y_MIN, Y_MAX = 2004.4, 2026.6
NODES = [2005, 2008, 2011, 2014, 2017, 2020, 2023, 2026]
LMIN = math.log10(0.0005)      # top 0.05% at the very top


def xs(year: float) -> float:
    return X0 + (year - Y_MIN) / (Y_MAX - Y_MIN) * (X1 - X0)


def ys(pct: float) -> float:
    top = max(1 - pct, 0.0005)
    return YT + (math.log10(top) - LMIN) / (0 - LMIN) * (YB - YT)


def pc(t: float) -> str:
    return f"{100 * t / T:.2f}%"


# pulse timing: arrive at node k at A[k] seconds, dwell, then jump
STEP, JUMP = 1.45, 0.22
A = [0.6 + k * STEP for k in range(len(NODES))]
END = A[-1] + 0.6              # sweep finished
FADE0, FADE1 = T - 1.4, T - 0.5

css = []
el = []

# ---------------------------------------------------------------- style tokens
css.append("""
svg{--bg:#fbfaf7;--ink:#1f2328;--muted:#59636e;--grid:#d8dee4;--axon:#8c959f;
--myelin:#dbe4f5;--myelin-edge:#8aa4d6;--dot:#a3abb5;--top10:#1a7f86;--top1:#c96a14;
--glow:#f0a050;--term:#1f2328;--term-ink:#e6edf3;--term-muted:#9198a1;--term-ok:#7ee2b8;
--chip:#eef1f4}
@media (prefers-color-scheme: dark){svg{--bg:#0d1117;--ink:#e6edf3;--muted:#9198a1;
--grid:#30363d;--axon:#6e7681;--myelin:#1c2b45;--myelin-edge:#4b6aa6;--dot:#5b636d;
--top10:#39c5cf;--top1:#f0a050;--glow:#ffb86b;--term:#161b22;--term-ink:#e6edf3;
--term-muted:#8b949e;--term-ok:#7ee2b8;--chip:#161b22}}
.bg{fill:var(--bg)}.ink{fill:var(--ink)}.muted{fill:var(--muted)}
.grid{stroke:var(--grid);stroke-width:1}.gridd{stroke:var(--grid);stroke-width:1;stroke-dasharray:3 4}
.axon{stroke:var(--axon);stroke-width:2.5}.myelin{fill:var(--myelin);stroke:var(--myelin-edge);stroke-width:1.2}
.d0{fill:var(--dot)}.d10{fill:var(--top10)}.d1{fill:var(--top1)}.lab{fill:var(--top1);font-weight:600}
.glow{fill:var(--glow)}.ring{fill:none;stroke:var(--glow);stroke-width:2}
.term{fill:var(--term)}.tink{fill:var(--term-ink)}.tmut{fill:var(--term-muted)}.tok{fill:var(--term-ok)}
.chip{fill:var(--chip);stroke:var(--grid)}
text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.pop{transform-box:fill-box;transform-origin:center}
@media (prefers-reduced-motion: reduce){*{animation:none!important}}
""")

# ---------------------------------------------------------------- background + left text
el.append(f'<rect class="bg" x="0" y="0" width="{W}" height="{H}" rx="18"/>')
el.append('<text class="muted" x="56" y="92" font-size="13" letter-spacing="2.2" font-weight="600">'
          'INHERITED PERIPHERAL NEUROPATHY</text>')
el.append('<text class="ink" x="54" y="152" font-size="62" font-weight="700" letter-spacing="-1">PNPrank</text>')
el.append('<text class="muted" x="56" y="192" font-size="18">Ranking candidate genes with public single-cell atlases,</text>')
el.append('<text class="muted" x="56" y="217" font-size="18">tested on the genes reported after each cutoff year.</text>')
x = 56
for chip in ("19,229 genes ranked", "11 atlases", "back-test AUROC 0.75–0.79"):
    w = 7.3 * len(chip) + 24
    el.append(f'<rect class="chip" x="{x}" y="240" width="{w:.0f}" height="28" rx="14"/>')
    el.append(f'<text class="ink" x="{x + 12}" y="259" font-size="13">{escape(chip)}</text>')
    x += w + 10

# ---------------------------------------------------------------- terminal card
TX, TY, TW, TH = 56, 292, 520, 150
el.append(f'<rect class="term" x="{TX}" y="{TY}" width="{TW}" height="{TH}" rx="12"/>')
el.append(f'<defs><clipPath id="tclip"><rect x="{TX}" y="{TY}" width="{TW}" height="{TH}" rx="12"/></clipPath></defs>')
el.append('<g clip-path="url(#tclip)">')
for i, c in enumerate(("#ff5f57", "#febc2e", "#28c840")):
    el.append(f'<circle cx="{TX + 20 + 16 * i}" cy="{TY + 18}" r="5" fill="{c}"/>')
SCENES = [
    ("pnprank top -n 3", ["  1  TUBA1A   unlabelled           1.000",
                          "  2  L1CAM    listed_not_in_truth  1.000",
                          "  3  MBP      unlabelled           1.000"]),
    ("pnprank query KIF5A", ["KIF5A   rank 35   percentile 0.9986",
                             "        known tier 2 (never used in training)"]),
    ("pnprank interval DHMN1=chr7:138,500,001-155,200,000",
     ["DHMN1  145 genes scored",
      "  1  EPHB6   0.981     2  KCNH2   0.981"]),
]
TS = 6.0                       # seconds per terminal scene
TT = TS * len(SCENES)
CW = 7.8                       # monospace advance at 13px
for k, (cmd, out) in enumerate(SCENES):
    t0 = k * TS
    typed = 0.045 * len(cmd)
    def tp(t):
        return f"{100 * t / TT:.2f}%"
    css.append(f".s{k}{{animation:s{k} {TT}s linear infinite}}"
               f"@keyframes s{k}{{0%,{tp(t0)}{{opacity:0}}{tp(t0 + 0.01)},{tp(t0 + TS - 0.3)}{{opacity:1}}"
               f"{tp(t0 + TS - 0.05)},100%{{opacity:0}}}}")
    css.append(f".c{k}{{animation:c{k} {TT}s infinite}}"
               f"@keyframes c{k}{{0%,{tp(t0 + 0.35)}{{transform:translateX(0);animation-timing-function:steps({len(cmd)},end)}}"
               f"{tp(t0 + 0.35 + typed)},100%{{transform:translateX({CW * len(cmd) + 4:.1f}px)}}}}")
    css.append(f".o{k}{{animation:o{k} {TT}s infinite}}"
               f"@keyframes o{k}{{0%,{tp(t0 + 0.5 + typed)}{{opacity:0}}{tp(t0 + 0.75 + typed)},100%{{opacity:1}}}}")
    vis = "1" if k == 0 else "0"
    g = [f'<g class="s{k}" opacity="{vis}">']
    g.append(f'<text class="mono tok" x="{TX + 20}" y="{TY + 52}" font-size="13">$</text>')
    g.append(f'<text class="mono tink" x="{TX + 36}" y="{TY + 52}" font-size="13" '
             f'textLength="{CW * len(cmd):.1f}" lengthAdjust="spacingAndGlyphs">{escape(cmd)}</text>')
    g.append(f'<rect class="term c{k}" x="{TX + 34}" y="{TY + 38}" width="{TW - 50}" height="20" '
             f'transform="translate({CW * len(cmd) + 4:.1f} 0)"/>')
    g.append(f'<g class="o{k}">')
    for j, line in enumerate(out):
        g.append(f'<text class="mono tmut" x="{TX + 20}" y="{TY + 78 + 21 * j}" font-size="13" '
                 f'xml:space="preserve">{escape(line)}</text>')
    g.append("</g></g>")
    el.extend(g)

el.append("</g>")

# ---------------------------------------------------------------- chart frame
el.append(f'<text class="ink" x="{X0 - 14}" y="44" font-size="15" font-weight="600">'
          'Each tier-1 gene, ranked genome-wide in the year it was reported</text>')
el.append(f'<text class="muted" x="{X0 - 14}" y="64" font-size="12.5">'
          'Model trained only on the genes known before that year (2005–2026, n = 70)</text>')
for top, lab in ((0.001, "top 0.1%"), (0.01, "top 1%"), (0.1, "top 10%"), (1.0, "100%")):
    y = ys(1 - top)
    cls = "gridd" if top < 1 else "grid"
    el.append(f'<line class="{cls}" x1="{X0}" x2="{X1}" y1="{y:.1f}" y2="{y:.1f}"/>')
    el.append(f'<text class="muted" x="{X0 - 8}" y="{y + 4:.1f}" font-size="11" text-anchor="end">{lab}</text>')

# ---------------------------------------------------------------- axon with myelin
el.append(f'<line class="axon" x1="{X0 - 20}" x2="{X1 + 12}" y1="{AXON_Y}" y2="{AXON_Y}"/>')
for a, b in zip(NODES, NODES[1:]):
    xa, xb = xs(a) + 6, xs(b) - 6
    el.append(f'<rect class="myelin" x="{xa:.1f}" y="{AXON_Y - 10}" width="{xb - xa:.1f}" height="20" rx="10"/>')
for nd in NODES:
    el.append(f'<text class="muted" x="{xs(nd):.1f}" y="{AXON_Y + 32}" font-size="11.5" text-anchor="middle">{nd}</text>')
el.append(f'<text class="muted" x="{X1 + 12}" y="{AXON_Y + 52}" font-size="11" text-anchor="end">'
          'myelinated axon as the time axis; nodes of Ranvier at the cutoff years</text>')

# ---------------------------------------------------------------- dots, appearing as the pulse passes
def node_index(year: int) -> int:
    for k, nd in enumerate(NODES):
        if year <= nd:
            return k
    return len(NODES) - 1


labels = []
for i, (sym, yr, p) in enumerate(sorted(genes, key=lambda g: (g[1], -g[2]))):
    k = node_index(yr)
    t = A[k] + 0.04 * (i % 6)
    cls = "d1" if p >= 0.99 else ("d10" if p >= 0.9 else "d0")
    r = 5.5 if p >= 0.99 else 4.5
    name = f"g{i}"
    css.append(f".{name}{{animation:{name} {T}s infinite}}"
               f"@keyframes {name}{{0%,{pc(t)}{{opacity:0;transform:scale(0)}}"
               f"{pc(t + 0.18)}{{opacity:1;transform:scale(1.6)}}{pc(t + 0.4)},{pc(FADE0)}{{opacity:1;transform:scale(1)}}"
               f"{pc(FADE1)},100%{{opacity:0;transform:scale(1)}}}}")
    cx, cy = xs(yr), ys(p)
    el.append(f'<circle class="{cls} pop {name}" cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" opacity="0.92"/>')
    if p >= 0.99:
        labels.append((sym, cx, cy, name))
for sym, cx, cy, name in labels:
    el.append(f'<text class="lab {name}" x="{cx + 8:.1f}" y="{cy + 4:.1f}" font-size="11.5">{sym}</text>')

# ---------------------------------------------------------------- pulse (action potential) + cursor
kf = ["0%{transform:translateX(0);opacity:0}", f"{pc(A[0] - 0.3)}{{transform:translateX(0);opacity:0}}",
      f"{pc(A[0])}{{transform:translateX(0);opacity:1}}"]
for k in range(1, len(NODES)):
    dx = xs(NODES[k]) - xs(NODES[0])
    prev = xs(NODES[k - 1]) - xs(NODES[0])
    kf.append(f"{pc(A[k] - JUMP)}{{transform:translateX({prev:.1f}px);opacity:1}}")
    kf.append(f"{pc(A[k])}{{transform:translateX({dx:.1f}px);opacity:1}}")
last = xs(NODES[-1]) - xs(NODES[0])
kf.append(f"{pc(END + 0.4)}{{transform:translateX({last:.1f}px);opacity:0}}")
kf.append(f"100%{{transform:translateX({last:.1f}px);opacity:0}}")
css.append(f".pulse{{animation:pulse {T}s infinite}}@keyframes pulse{{{''.join(kf)}}}")
x0n = xs(NODES[0])
el.append(f'<g class="pulse" opacity="0" transform="translate({last:.1f} 0)">'
          f'<line class="gridd" x1="{x0n:.1f}" x2="{x0n:.1f}" y1="{YT - 6}" y2="{AXON_Y - 12}" opacity="0.9"/>'
          f'<circle class="glow" cx="{x0n:.1f}" cy="{AXON_Y}" r="13" opacity="0.25"/>'
          f'<circle class="glow" cx="{x0n:.1f}" cy="{AXON_Y}" r="6.5"/></g>')
for k, nd in enumerate(NODES):
    name = f"r{k}"
    css.append(f".{name}{{animation:{name} {T}s ease-out infinite}}"
               f"@keyframes {name}{{0%,{pc(A[k] - 0.01)}{{opacity:0;transform:scale(.4)}}{pc(A[k] + 0.05)}{{opacity:.9}}"
               f"{pc(A[k] + 0.9)},100%{{opacity:0;transform:scale(2.4)}}}}")
    el.append(f'<circle class="ring pop {name}" cx="{xs(nd):.1f}" cy="{AXON_Y}" r="8" opacity="0"/>')

# ---------------------------------------------------------------- closing caption
css.append(f".cap{{animation:cap {T}s infinite}}@keyframes cap{{0%,{pc(END)}{{opacity:0}}"
           f"{pc(END + 0.5)},{pc(FADE0)}{{opacity:1}}{pc(FADE1)},100%{{opacity:0}}}}")
el.append(f'<g class="cap"><text class="ink" x="{X1}" y="{YB + 24}" font-size="13" text-anchor="end">'
          f'<tspan class="lab">{n_top1} of {len(genes)}</tspan> genes were in the top 1% of the genome '
          'in the year they were reported</text></g>')

title = "PNPrank"
desc = (f"Animated banner. A myelinated axon serves as a time axis from 2005 to 2026. An action "
        f"potential jumps between nodes of Ranvier; as it passes, each of {len(genes)} neuropathy genes "
        f"appears at the position PNPrank gave it genome-wide in the year it was reported, and "
        f"{n_top1} of them land in the top 1%. A terminal shows pnprank top, query and interval commands.")
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
       f'role="img" aria-labelledby="t d">\n<title id="t">{title}</title>\n<desc id="d">{escape(desc)}</desc>\n'
       f'<style>{"".join(css)}</style>\n' + "\n".join(el) + "\n</svg>\n")
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(svg)
print(f"wrote {OUT} ({len(svg) / 1024:.0f} KB, {len(genes)} genes, {n_top1} in top 1%)")
