"""Step 1b: find anchor proteins by HMM search.

Only protein models are permitted here; see configs/blind/anchors.json. The anchor set
defines where we look, and nothing about what we expect to find there.
"""

from __future__ import annotations

import dataclasses
import json
import subprocess
import tempfile
from pathlib import Path

from .orfs import ORF
from .tools import which


@dataclasses.dataclass(frozen=True)
class Anchor:
    orf: ORF
    family: str
    model: str  # best-scoring model within the assigned family
    evalue: float
    score: float
    all_models: tuple[str, ...] = ()  # every model that hit this ORF, best first
    all_families: tuple[str, ...] = ()  # every family represented, highest priority first

    @property
    def is_ambiguous(self) -> bool:
        """The ORF hit more than one anchor family — e.g. a fused or chimeric protein."""
        return len(self.all_families) > 1


def load_config(path: str | Path) -> dict:
    with open(path) as fh:
        cfg = json.load(fh)
    return validate_config(cfg)


def validate_config(cfg: dict) -> dict:
    """Enforce the blind and the model-to-family invariant.

    Shared with the neighbourhood config, which carries protein marker groups
    under the same rules: no RNA model may enter either file.
    """
    seen: dict[str, str] = {}
    for fam, spec in cfg["families"].items():
        for m in spec["models"]:
            if m["acc"].startswith("RF"):
                raise ValueError(f"{fam}/{m['name']}: Rfam RNA model in an anchor config")
            # model -> family is a single mapping downstream, so a model listed
            # under two families would silently belong to whichever was parsed last.
            if m["name"] in seen:
                raise ValueError(
                    f"model {m['name']} listed under both {seen[m['name']]} and {fam}; "
                    f"each model must belong to exactly one family")
            seen[m["name"]] = fam
    return cfg


def family_priority(cfg: dict) -> dict[str, int]:
    """Family -> priority. Higher wins when one ORF hits several families."""
    return {fam: spec.get("priority", 0) for fam, spec in cfg["families"].items()}


def build_anchor_hmm(cfg: dict, out_hmm: str | Path) -> dict[str, str]:
    """Fetch the configured models out of Pfam-A into one small HMM file.

    Returns model-name -> family.
    """
    out_hmm = Path(out_hmm)
    db = cfg["hmm_db"]
    # hmmfetch writes <db>.ssi for a plain HMM file but <db>.h3m.ssi once the
    # database has been hmmpress'ed, and refuses to overwrite either. Pfam-A
    # here is pressed, so check both before building the index.
    if not any(Path(db + suffix).exists() for suffix in (".ssi", ".h3m.ssi")):
        subprocess.run([which("hmmfetch"), "--index", db], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    model_family: dict[str, str] = {}
    names = []
    for fam, spec in cfg["families"].items():
        for m in spec["models"]:
            model_family[m["name"]] = fam
            names.append(m["name"])

    with tempfile.NamedTemporaryFile("w", suffix=".keys", delete=False) as kf:
        kf.write("\n".join(names) + "\n")
        keyfile = kf.name
    with open(out_hmm, "wb") as out:
        subprocess.run([which("hmmfetch"), "-f", db, keyfile], stdout=out, check=True,
                       stderr=subprocess.DEVNULL)
    Path(keyfile).unlink(missing_ok=True)
    return model_family


def _write_faa(orfs: list[ORF], path: Path) -> None:
    with open(path, "w") as fh:
        for o in orfs:
            fh.write(f">{o.index}\n{o.aa}\n")


def find_anchors(orfs: list[ORF], anchor_hmm: str | Path, model_family: dict[str, str],
                 cfg: dict, workdir: str | Path | None = None) -> list[Anchor]:
    """hmmsearch the anchor models against one genome's ORFs.

    One ORF yields at most one anchor. Family assignment is by **priority first,
    score second**: a group II intron maturase hits GIIM and RVT_1 both, and
    picking the higher-scoring domain would scatter the same biological family
    between GroupII_RT and RT depending on sequence. The specific family wins,
    and every model and family that hit is kept on the record so the ambiguity
    stays visible rather than being silently resolved.
    """
    th = cfg["thresholds"]
    orfs = [o for o in orfs if len(o.aa) >= th["min_aa_len"]]
    if not orfs:
        return []

    tmp = tempfile.TemporaryDirectory(dir=workdir) if workdir else tempfile.TemporaryDirectory()
    with tmp as td:
        td = Path(td)
        faa, dom = td / "q.faa", td / "hits.domtbl"
        _write_faa(orfs, faa)

        cmd = [which("hmmsearch"), "--domtblout", str(dom), "--noali", "--cpu", "1"]
        if th.get("use_cut_ga", True):
            cmd.append("--cut_ga")
        cmd += [str(anchor_hmm), str(faa)]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        by_index = {o.index: o for o in orfs}
        prio = family_priority(cfg)

        # Collect every passing domain hit per ORF, then resolve once.
        hits: dict[int, list[tuple[str, str, float, float]]] = {}
        with open(dom) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                idx, model = int(f[0]), f[3]
                ievalue, score = float(f[12]), float(f[13])
                if ievalue > th["max_domain_evalue"]:
                    continue
                fam = model_family.get(model)
                if fam is None:
                    continue
                hits.setdefault(idx, []).append((fam, model, ievalue, score))

        anchors = []
        for idx, hs in hits.items():
            # Best model per ORF overall, for the record.
            hs.sort(key=lambda h: -h[3])
            # Winning family: highest priority, then highest score within it.
            fam = max(hs, key=lambda h: (prio.get(h[0], 0), h[3]))[0]
            in_fam = [h for h in hs if h[0] == fam]
            _, model, ievalue, score = max(in_fam, key=lambda h: h[3])
            # Length floor is per family and applied *after* assignment: RVT_3
            # otherwise admits a 138 aa host protein as a reverse transcriptase,
            # and the right threshold differs by family (Cas12a ~1300 aa,
            # TnpB ~400 aa), so it cannot be a single global cutoff.
            fam_min = cfg["families"][fam].get("min_aa_len", th["min_aa_len"])
            if len(by_index[idx].aa) < fam_min:
                continue
            fams = sorted({h[0] for h in hs}, key=lambda f: (-prio.get(f, 0), f))
            anchors.append(Anchor(
                orf=by_index[idx], family=fam, model=model, evalue=ievalue, score=score,
                all_models=tuple(dict.fromkeys(h[1] for h in hs)),
                all_families=tuple(fams),
            ))
    return sorted(anchors, key=lambda a: (a.orf.contig, a.orf.start))
