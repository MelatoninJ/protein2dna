"""Host profiles.

Host-specific behaviour lives in YAML, never in ``if host == "ecoli"``
branches.  Adding a new expression host should be a pull request that adds one
file to ``p2d/hosts/`` and nothing else.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

import yaml

HOSTS_PACKAGE = "p2d.hosts"


@dataclass
class HostProfile:
    slug: str
    name: str
    raw: dict

    codon_usage: dict[str, dict[str, float]] = field(default_factory=dict)
    stop_usage: dict[str, float] = field(default_factory=dict)
    rare_codons: set[str] = field(default_factory=set)
    forbidden_motifs: list[dict] = field(default_factory=list)

    @property
    def codon_table_verified(self) -> bool:
        return bool(self.raw.get("codon_usage", {}).get("verified", False))

    def best_codon(self, aa: str) -> str:
        """Most-used codon for ``aa``.

        v0.1 placeholder strategy.  Always-best-codon is known to be a poor
        choice in practice (tRNA depletion, ribosome jamming); the real
        optimiser lands in v0.2 and this becomes one strategy among several.
        """
        if aa == "*":
            return max(self.stop_usage, key=self.stop_usage.get)
        try:
            codons = self.codon_usage[aa]
        except KeyError as exc:
            raise ValueError(f"no codons for residue {aa!r} in host {self.slug!r}") from exc
        return max(codons, key=codons.get)

    def codons_for(self, aa: str) -> dict[str, float]:
        return dict(self.codon_usage[aa])


def _parse(data: dict) -> HostProfile:
    usage = data.get("codon_usage", {})
    table = usage.get("table", {})
    profile = HostProfile(
        slug=data["slug"],
        name=data["name"],
        raw=data,
        codon_usage={aa: dict(c) for aa, c in table.items()},
        stop_usage=dict(usage.get("stop", {"TAA": 1.0})),
        rare_codons=set(data.get("rare_codons", [])),
        forbidden_motifs=list(data.get("forbidden_motifs", [])),
    )
    if not profile.codon_table_verified:
        warnings.warn(
            f"host {profile.slug!r} ships an unverified placeholder codon table "
            "-- replace it with a cited usage table before relying on output",
            UserWarning,
            stacklevel=3,
        )
    return profile


def load_host(slug_or_path: str | Path) -> HostProfile:
    """Load a bundled profile by slug, or any YAML file by path."""
    path = Path(slug_or_path)
    if path.suffix in {".yaml", ".yml"} and path.exists():
        return _parse(yaml.safe_load(path.read_text()))
    try:
        text = resources.files(HOSTS_PACKAGE).joinpath(f"{slug_or_path}.yaml").read_text()
    except FileNotFoundError as exc:
        raise ValueError(
            f"unknown host {slug_or_path!r}; available: {', '.join(available_hosts())}"
        ) from exc
    return _parse(yaml.safe_load(text))


def available_hosts() -> list[str]:
    return sorted(
        p.name.removesuffix(".yaml")
        for p in resources.files(HOSTS_PACKAGE).iterdir()
        if p.name.endswith(".yaml")
    )
