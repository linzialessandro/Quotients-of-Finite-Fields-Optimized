"""
Atlas of finite-field quotient hyperfields by order.

Groups all F_q / G_r (r fixed so |K|=r+1) into isomorphism classes using
structure fingerprints, and records char / C-char / sample prime powers.

Used for:
- Massouros order-7 comparison (r=6)
- Ameri-style H_r vs Q_r^{fin} tables
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from quotient_hyperfields.experiments import (
    baker_jin_residue,
    prime_powers_with_index,
    resolve_scan_limit,
    structure_fingerprint,
)
from quotient_hyperfields.hyperfield import QuotientHyperfield
from quotient_hyperfields.isomorphism import baker_jin_bound_sharp


@dataclass
class QuotientClass:
    """One isomorphism class of finite-field quotients of fixed order."""

    class_id: int
    r: int
    order: int
    fingerprint: tuple
    sample_qs: list[int]
    char: float
    c_char: float
    baker_jin_residues: set[str] = field(default_factory=set)
    is_stable_large: bool = False  # has a large q witness

    def summary_line(self) -> str:
        ch = "∞" if math.isinf(self.char) else str(int(self.char))
        cc = "∞" if math.isinf(self.c_char) else str(int(self.c_char))
        qs = ",".join(str(q) for q in self.sample_qs[:6])
        if len(self.sample_qs) > 6:
            qs += ",..."
        res = ",".join(sorted(self.baker_jin_residues))
        flag = "stable" if self.is_stable_large else "sporadic"
        return (
            f"  class {self.class_id:3d}  [{flag:8s}]  char={ch:>3} C-char={cc:>3}  "
            f"res={{{res}}}  #q={len(self.sample_qs):3d}  samples=[{qs}]"
        )


@dataclass
class OrderAtlas:
    """All finite-field quotient classes of a given hyperfield order."""

    order: int
    r: int
    max_q: int
    classes: list[QuotientClass]
    theoretically_complete: bool

    @property
    def q_fin(self) -> int:
        return len(self.classes)

    def summary(self) -> str:
        n12 = baker_jin_bound_sharp(self.r) if self.r >= 2 else 0
        lines = [
            f"Finite-field quotient atlas: order n={self.order} (r={self.r})",
            f"  max_q={self.max_q}  Remark1.2 N_r={n12}  "
            f"complete={self.theoretically_complete}",
            f"  Q_r^{{fin}} = {self.q_fin} isomorphism classes",
            "",
        ]
        # stable first, then sporadic by min sample q
        ordered = sorted(
            self.classes,
            key=lambda c: (not c.is_stable_large, min(c.sample_qs), c.class_id),
        )
        for c in ordered:
            lines.append(c.summary_line())
        return "\n".join(lines)


def build_order_atlas(
    order: int,
    max_q: int | None = None,
    *,
    max_q_cap: int = 4000,
    progress: bool = False,
) -> OrderAtlas:
    """
    Enumerate iso classes of F_q/G with |F_q/G| = order.

    Parameters
    ----------
    order :
        Hyperfield order n = r+1 >= 2.
    max_q :
        Optional scan ceiling; default resolves past Remark 1.2 N_r.
    """
    if order < 2:
        raise ValueError("order must be >= 2")
    r = order - 1
    if r >= 2:
        limit = resolve_scan_limit(r, max_q, max_q_cap=max_q_cap)
        n12 = baker_jin_bound_sharp(r)
        complete = limit >= n12
    else:
        limit = max_q if max_q is not None else 64
        limit = min(limit, max_q_cap)
        n12 = 0
        complete = True

    qs = prime_powers_with_index(r, limit)
    # fingerprint -> aggregation
    buckets: dict[tuple, dict] = {}
    stable_cache: dict[tuple[int, str], tuple] = {}
    for i, q in enumerate(qs):
        if progress:
            print(f"  [atlas n={order}] {i + 1}/{len(qs)} q={q}", flush=True)
        h = QuotientHyperfield.from_q_r(q, r)
        # Large-q fingerprints reuse one Aut-canonical table per Baker–Jin
        # residue (same classes as full computation; v0.2 scan optimisation).
        fp = structure_fingerprint(
            h, use_stable=True, _stable_cache=stable_cache
        )
        if fp not in buckets:
            buckets[fp] = {
                "qs": [],
                "char": h.char,
                "c_char": h.c_char,
                "residues": set(),
            }
        buckets[fp]["qs"].append(q)
        buckets[fp]["residues"].add(baker_jin_residue(q, r))
        # invariants must match within class
        if (buckets[fp]["char"], buckets[fp]["c_char"]) != (h.char, h.c_char):
            # keep first; mark inconsistency via extreme values
            buckets[fp]["char"] = float("nan")

    classes: list[QuotientClass] = []
    for cid, (fp, data) in enumerate(sorted(buckets.items(), key=lambda kv: min(kv[1]["qs"]))):
        sample = sorted(data["qs"])
        is_large = any(q >= n12 for q in sample) if r >= 2 else True
        classes.append(
            QuotientClass(
                class_id=cid,
                r=r,
                order=order,
                fingerprint=fp,
                sample_qs=sample,
                char=data["char"],
                c_char=data["c_char"],
                baker_jin_residues=set(data["residues"]),
                is_stable_large=is_large,
            )
        )

    return OrderAtlas(
        order=order,
        r=r,
        max_q=limit,
        classes=classes,
        theoretically_complete=complete,
    )
