"""Look-and-feel helpers for the Streamlit app: CSS, hero banner, cards, helix graphic.

Everything user-supplied (record names, descriptions) is HTML-escaped before it is
drawn, because FASTA headers can contain anything.
"""

from __future__ import annotations

import html
import math
import random

import streamlit as st

from .insights import Insight
from .plots import BASE_STRONG

TEAL = "#0F9D8F"
BLUE = BASE_STRONG["C"]
ORANGE = BASE_STRONG["G"]
GREEN = BASE_STRONG["A"]
VERMILLION = BASE_STRONG["T"]
PURPLE = "#CC79A7"

CSS = """
:root{--ss-surface:rgba(127,127,127,.09);--ss-line:rgba(127,127,127,.28)}
.block-container{padding-top:2rem;max-width:1280px}
.stTabs [data-baseweb="tab"]{font-weight:600}
.ss-hero{position:relative;overflow:hidden;display:flex;align-items:center;justify-content:space-between;gap:1.5rem;flex-wrap:wrap;border-radius:28px;padding:2.8rem 3rem;color:#fff;background:radial-gradient(900px 380px at 92% -15%,rgba(45,212,191,.55),transparent 62%),radial-gradient(600px 320px at -5% 110%,rgba(124,92,191,.35),transparent 60%),linear-gradient(135deg,#0B2545 0%,#12355F 50%,#0E6B67 100%);box-shadow:0 24px 60px rgba(11,37,69,.35);margin-bottom:1.6rem}
.ss-hero-text{flex:1 1 420px;min-width:260px}
.ss-eyebrow{font-size:.8rem;letter-spacing:.22em;text-transform:uppercase;font-weight:600;color:#7DE3D2;margin-bottom:.6rem}
.ss-title{font-family:SpaceGrotesk,Inter,sans-serif;font-weight:700;font-size:4rem;line-height:1;letter-spacing:-.03em;margin:0 0 1rem}
.ss-title span{background:linear-gradient(90deg,#7DE3D2,#F5CD7A);-webkit-background-clip:text;background-clip:text;color:transparent}
.ss-sub{font-size:1.12rem;line-height:1.55;max-width:34rem;color:rgba(255,255,255,.88);margin:0 0 1.3rem}
.ss-pills{display:flex;flex-wrap:wrap;gap:.5rem}
.ss-pill{font-size:.78rem;font-weight:600;padding:.3rem .8rem;border-radius:999px;border:1px solid rgba(255,255,255,.28);background:rgba(255,255,255,.1);color:#fff}
.ss-helix{flex:0 1 440px;max-width:100%;animation:ss-float 7s ease-in-out infinite}
.ss-helix svg{width:100%;height:auto;display:block}
@keyframes ss-float{0%,100%{transform:translateY(0)}50%{transform:translateY(-9px)}}
@media (prefers-reduced-motion:reduce){.ss-helix{animation:none}}
@media (max-width:760px){.ss-hero{padding:1.7rem 1.4rem;border-radius:20px}.ss-title{font-size:2.9rem}.ss-helix{display:none}.ss-banner-id{font-size:1.5rem}}
.ss-banner{border-radius:22px;padding:1.4rem 1.8rem;color:#fff;background:linear-gradient(135deg,#0B2545 0%,#12355F 55%,#0E6B67 100%);margin-bottom:1.2rem;display:flex;flex-wrap:wrap;align-items:center;justify-content:space-between;gap:1rem;box-shadow:0 14px 34px rgba(11,37,69,.28)}
.ss-banner-id{font-family:SpaceGrotesk,Inter,sans-serif;font-weight:700;font-size:2rem;line-height:1.1;letter-spacing:-.02em;word-break:break-all}
.ss-banner-desc{font-size:.95rem;opacity:.86;margin-top:.3rem;max-width:46rem}
.ss-chips{display:flex;flex-wrap:wrap;gap:.5rem}
.ss-chip{font-size:.8rem;font-weight:600;padding:.35rem .8rem;border-radius:999px;background:rgba(255,255,255,.14);border:1px solid rgba(255,255,255,.25);color:#fff}
.ss-stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin:.2rem 0 1.2rem}
.ss-stat{position:relative;border:1px solid var(--ss-line);background:var(--ss-surface);border-radius:18px;padding:1rem 1.2rem 1rem 1.5rem;overflow:hidden}
.ss-stat:before{content:"";position:absolute;left:0;top:0;bottom:0;width:7px;background:var(--c)}
.ss-stat-label{font-size:.74rem;text-transform:uppercase;letter-spacing:.09em;opacity:.7;font-weight:600}
.ss-stat-value{font-family:SpaceGrotesk,Inter,sans-serif;font-size:2.1rem;font-weight:700;line-height:1.15;letter-spacing:-.02em}
.ss-stat-hint{font-size:.82rem;opacity:.65}
.ss-section{margin:1.1rem 0 .7rem}
.ss-section-title{font-family:SpaceGrotesk,Inter,sans-serif;font-weight:700;font-size:1.35rem;letter-spacing:-.01em}
.ss-section-sub{font-size:.88rem;opacity:.7;margin-top:.1rem}
.ss-section:after{content:"";display:block;width:46px;height:4px;border-radius:4px;margin-top:.5rem;background:linear-gradient(90deg,#0F9D8F,#E69F00)}
.ss-insight{display:flex;gap:.8rem;border:1px solid var(--ss-line);border-left:6px solid var(--c);border-radius:14px;background:var(--ss-surface);padding:.85rem 1rem;margin-bottom:.7rem}
.ss-dot{flex:0 0 1.5rem;height:1.5rem;border-radius:50%;background:var(--c);color:#fff;font-weight:700;font-size:.85rem;display:flex;align-items:center;justify-content:center;margin-top:.1rem}
.ss-insight b{display:block;margin-bottom:.15rem}
.ss-insight p{margin:0;font-size:.9rem;opacity:.9;line-height:1.45}
.ss-tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px;margin:1.4rem 0 .6rem}
.ss-tile{border:1px solid var(--ss-line);background:var(--ss-surface);border-radius:20px;padding:1.2rem 1.25rem 1.3rem;border-top:5px solid var(--c)}
.ss-tile svg,.ss-action svg{width:30px;height:30px;color:var(--c);margin-bottom:.5rem}
.ss-tile b,.ss-action b{display:block;font-family:SpaceGrotesk,Inter,sans-serif;font-size:1.12rem;margin-bottom:.25rem}
.ss-tile p,.ss-action p{margin:0;font-size:.9rem;opacity:.8;line-height:1.45}
.ss-action{padding:.3rem .2rem .6rem}
.ss-wait{border:2px dashed var(--ss-line);border-radius:22px;padding:2.2rem 1.5rem;text-align:center;margin:1rem 0}
.ss-wait b{font-family:SpaceGrotesk,Inter,sans-serif;font-size:1.3rem}
.ss-brand{display:flex;align-items:center;gap:.6rem;margin-bottom:.4rem}
.ss-brand svg{width:34px;height:34px}
.ss-brand-name{font-family:SpaceGrotesk,Inter,sans-serif;font-weight:700;font-size:1.55rem;letter-spacing:-.02em;line-height:1}
.ss-brand small{display:block;opacity:.6;font-size:.75rem;margin-top:.15rem}
.ss-footer{text-align:center;opacity:.6;font-size:.82rem;margin:2.2rem 0 .5rem}
"""


def _compact(markup: str) -> str:
    """Strip indentation and blank lines: Markdown treats those as the end of an HTML block."""
    return "\n".join(line.strip() for line in markup.splitlines() if line.strip())


def render(markup: str) -> None:
    st.markdown(_compact(markup), unsafe_allow_html=True)


def inject_css() -> None:
    st.markdown("<style>" + _compact(CSS) + "</style>", unsafe_allow_html=True)


# --------------------------------------------------------------------------- #
# Graphics
# --------------------------------------------------------------------------- #
def helix_svg(width: int = 440, height: int = 260, rungs: int = 28, turns: float = 2.0, seed: int = 7) -> str:
    """A double helix whose rungs are coloured like the bases in the sequence view."""
    rnd = random.Random(seed)
    pairs = [("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")]
    mid, amp = height / 2, height * 0.30
    left, right = 24, width - 24
    parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="DNA double helix">']

    for sign in (1, -1):  # the two backbones
        d = []
        for k in range(241):
            t = k / 240
            x = left + t * (right - left)
            y = mid + sign * amp * math.sin(2 * math.pi * turns * t)
            d.append(f"{'M' if k == 0 else 'L'}{x:.1f},{y:.1f}")
        parts.append(
            f'<path d="{" ".join(d)}" fill="none" stroke="rgba(255,255,255,.55)" stroke-width="3" stroke-linecap="round"/>'
        )

    for i in range(rungs):
        t = i / (rungs - 1)
        x = left + t * (right - left)
        s = math.sin(2 * math.pi * turns * t)
        y1, y2 = mid + amp * s, mid - amp * s
        a, b = rnd.choice(pairs)
        parts.append(f'<line x1="{x:.1f}" y1="{y1:.1f}" x2="{x:.1f}" y2="{mid:.1f}" stroke="{BASE_STRONG[a]}" stroke-width="3.4" stroke-linecap="round"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{mid:.1f}" x2="{x:.1f}" y2="{y2:.1f}" stroke="{BASE_STRONG[b]}" stroke-width="3.4" stroke-linecap="round"/>')
        parts.append(f'<circle cx="{x:.1f}" cy="{y1:.1f}" r="5.2" fill="{BASE_STRONG[a]}" stroke="#fff" stroke-width="1.2"/>')
        parts.append(f'<circle cx="{x:.1f}" cy="{y2:.1f}" r="5.2" fill="{BASE_STRONG[b]}" stroke="#fff" stroke-width="1.2"/>')

    parts.append("</svg>")
    return "".join(parts)


def _logo_svg() -> str:
    pts_a = " ".join(f"{x * 2:.1f},{16 + 9 * math.sin(x / 34 * 2 * math.pi * 1.5):.1f}" for x in range(0, 18))
    pts_b = " ".join(f"{x * 2:.1f},{16 - 9 * math.sin(x / 34 * 2 * math.pi * 1.5):.1f}" for x in range(0, 18))
    return (
        '<svg viewBox="0 0 34 32" xmlns="http://www.w3.org/2000/svg" aria-hidden="true">'
        f'<polyline points="{pts_a}" fill="none" stroke="{TEAL}" stroke-width="3" stroke-linecap="round"/>'
        f'<polyline points="{pts_b}" fill="none" stroke="{ORANGE}" stroke-width="3" stroke-linecap="round"/></svg>'
    )


def _icon(name: str) -> str:
    """Small line icons drawn for this app."""
    paths = {
        "spark": '<path d="M12 3l2.2 5.8L20 11l-5.8 2.2L12 19l-2.2-5.8L4 11l5.8-2.2z"/>',
        "upload": '<path d="M12 16V4m0 0l-4 4m4-4l4 4M4 16v3a1 1 0 001 1h14a1 1 0 001-1v-3"/>',
        "paste": '<rect x="6" y="5" width="12" height="16" rx="2"/><path d="M9 5V4h6v1M9 10h6M9 14h6"/>',
        "ring": '<circle cx="12" cy="12" r="8.5"/><circle cx="12" cy="12" r="4.5"/><path d="M12 3.5V7"/>',
        "wave": '<path d="M2.5 12c3-9 6-9 9.5 0s6.500 9 9.500 0"/>',
        "orf": '<path d="M3 6.500h10M8 12h13M5 17.500h9" stroke-width="3.2"/>',
        "bases": (
            f'<rect x="3" y="3" width="8" height="8" rx="2" fill="{GREEN}" stroke="none"/>'
            f'<rect x="13" y="3" width="8" height="8" rx="2" fill="{BLUE}" stroke="none"/>'
            f'<rect x="3" y="13" width="8" height="8" rx="2" fill="{ORANGE}" stroke="none"/>'
            f'<rect x="13" y="13" width="8" height="8" rx="2" fill="{VERMILLION}" stroke="none"/>'
        ),
    }
    return (
        '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" '
        f'stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">{paths[name]}</svg>'
    )


# --------------------------------------------------------------------------- #
# Components
# --------------------------------------------------------------------------- #
def hero() -> None:
    render(
        f"""
        <div class="ss-hero">
        <div class="ss-hero-text">
        <div class="ss-eyebrow">DNA sequence explorer</div>
        <div class="ss-title">Sidero<span>Scope</span></div>
        <p class="ss-sub">See a bacterial genome at a glance: coloured bases, GC content and skew,
        six-frame ORFs and a circular genome map. Made for microbial genomics, with siderophore
        gene screening on the roadmap.</p>
        <div class="ss-pills"><span class="ss-pill">FASTA &amp; GenBank</span>
        <span class="ss-pill">Runs on your computer</span><span class="ss-pill">Open source (MIT)</span></div>
        </div>
        <div class="ss-helix">{helix_svg()}</div>
        </div>
        """
    )


def action_card(icon: str, title: str, text: str, color: str) -> None:
    render(
        f'<div class="ss-action" style="--c:{color}">{_icon(icon)}<b>{html.escape(title)}</b>'
        f"<p>{html.escape(text)}</p></div>"
    )


def feature_tiles() -> None:
    tiles = [
        ("bases", "#009E73", "Coloured sequence view", "Read the bases as a colour-blind-friendly mosaic instead of a wall of letters."),
        ("wave", "#E69F00", "GC content and skew", "Look for GC dips, which can mark genomic islands, and for skew shifts near the origin and terminus."),
        ("orf", "#0072B2", "Six-frame ORFs", "Candidate genes on both strands, with protein translations and downloads."),
        ("ring", "#CC79A7", "Circular genome map", "ORFs, GC and skew on a single ring that is ready to screenshot."),
    ]
    cells = "".join(
        f'<div class="ss-tile" style="--c:{c}">{_icon(i)}<b>{html.escape(t)}</b><p>{html.escape(d)}</p></div>'
        for i, c, t, d in tiles
    )
    render(f'<div class="ss-tiles">{cells}</div>')


def sidebar_brand(version: str) -> None:
    render(
        f'<div class="ss-brand">{_logo_svg()}<div><div class="ss-brand-name">SideroScope</div>'
        f"<small>v{html.escape(version)}</small></div></div>"
    )


def banner(record_id: str, description: str, chips: list[str]) -> None:
    chip_html = "".join(f'<span class="ss-chip">{html.escape(c)}</span>' for c in chips)
    desc = f'<div class="ss-banner-desc">{html.escape(description)}</div>' if description else ""
    render(
        f'<div class="ss-banner"><div><div class="ss-banner-id">{html.escape(record_id)}</div>{desc}</div>'
        f'<div class="ss-chips">{chip_html}</div></div>'
    )


def stat_cards(cards: list[tuple[str, str, str, str]]) -> None:
    """cards: (label, value, hint, colour)"""
    cells = "".join(
        f'<div class="ss-stat" style="--c:{color}"><div class="ss-stat-label">{html.escape(label)}</div>'
        f'<div class="ss-stat-value">{html.escape(value)}</div><div class="ss-stat-hint">{html.escape(hint)}</div></div>'
        for label, value, hint, color in cards
    )
    render(f'<div class="ss-stats">{cells}</div>')


def section(title: str, subtitle: str | None = None) -> None:
    sub = f'<div class="ss-section-sub">{html.escape(subtitle)}</div>' if subtitle else ""
    render(f'<div class="ss-section"><div class="ss-section-title">{html.escape(title)}</div>{sub}</div>')


def insight_cards(insights: list[Insight]) -> None:
    colors = {"info": BLUE, "good": TEAL, "warning": ORANGE}
    glyphs = {"info": "i", "good": "&#10003;", "warning": "!"}
    for item in insights:
        render(
            f'<div class="ss-insight" style="--c:{colors[item.level]}"><div class="ss-dot">{glyphs[item.level]}</div>'
            f"<div><b>{html.escape(item.title)}</b><p>{html.escape(item.text)}</p></div></div>"
        )


def waiting_card(mode: str) -> None:
    text = {
        "upload": ("Drop a FASTA or GenBank file above", "Files stay on your computer. Genomes of several Mb are fine."),
        "paste": ("Paste a sequence above", "FASTA with a header line, or just the raw bases."),
    }.get(mode, ("Choose a sequence to begin", ""))
    render(f'<div class="ss-wait"><b>{html.escape(text[0])}</b><p>{html.escape(text[1])}</p></div>')


def footer(version: str) -> None:
    render(f'<div class="ss-footer">SideroScope v{html.escape(version)} &middot; MIT licence &middot; fonts: Inter, Space Grotesk, JetBrains Mono (SIL OFL)</div>')
