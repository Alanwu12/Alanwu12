#!/usr/bin/env python3
"""
Build and QC the Proteinbase Challenge 1 submission set (anti-EGFR binders).

Scaffold: humanised camelid VHH (nanobody), single chain, ~118-122 aa.
Only the three CDR loops are designed; all four framework regions are
canonical VHH frameworks (hallmark tetrad F42/E49/R50/G52 retained, single
canonical C22-C96 disulfide) so that every submitted chain is near-certain
to fold.

Writes: submission.csv, qc_report.tsv
"""
import csv
import os

# --------------------------------------------------------------------------
# Framework sets. Three canonical VHH frameworks, used to spread the portfolio
# across sequence space (helps the novelty filter) without touching the
# hallmark residues that make a VHH soluble as a single domain.
# --------------------------------------------------------------------------
FRAMEWORKS = {
    "A": dict(
        fr1="QVQLVESGGGLVQAGGSLRLSCAAS",
        fr2="MGWFRQAPGKEREFVA",
        fr3="YYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYC",
        fr4="WGQGTQVTVSS",
    ),
    "B": dict(
        fr1="QVQLQESGGGSVQAGGSLKLSCAAS",
        fr2="LGWFRQAPGKEREGVS",
        fr3="YYADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAVYYC",
        fr4="WGQGTLVTVSS",
    ),
    "C": dict(
        fr1="EVQLVESGGGLVQPGGSLRLSCTAS",
        fr2="MAWFRQAPGQEREGVA",
        fr3="YYVDSVKGRFTISRDNSKNTVYLQMNSLKPEDTAVYYC",
        fr4="WGQGTQVTVSS",
    ),
}

# --------------------------------------------------------------------------
# Designs: (name, objective, framework, cdr1, cdr2, cdr3, note)
#
# Objective 1 - affinity for human EGFR domain III.
#   Paratope is aromatic (Y/W/F) + acidic (D/E): Tyr/Trp to pack the
#   F357/I467/L382 hydrophobic shelf, Asp/Glu to pair with K443/K465, and
#   Ser/Thr to hydrogen bond Q384/Q408/H409 (mature-protein numbering;
#   add 24 for UniProt P00533 numbering).
#
# Objective 2 - human + mouse cross-reactivity.
#   Shorter CDR3, paratope recentred on the ligand-binding core
#   (D355/F357/Q384/Q408/H409), which is conserved because EGF engagement is
#   cross-species. Salt bridges are deliberately depleted so the footprint
#   does not rely on the 465-475 periphery that carries the human-specific
#   determinants (S468; the cetuximab-resistance residue S492 in UniProt
#   numbering) which make cetuximab and panitumumab human-restricted.
#
# Objective 3 - bind at pH 6.5, not at pH 7.4.
#   Histidine switch: 4-5 His placed across the three CDRs and nowhere in the
#   frameworks. At pH 6.5 ~24% of each imidazole is protonated versus ~4% at
#   pH 7.4, so the paratope gains roughly +1 charge on acidification and makes
#   its electrostatic contacts with the acidic domain III subsite. Arg/Lys are
#   kept out of the CDRs so binding energy is His-protonation dependent rather
#   than constitutive; aromatics are retained as a weak pH-independent anchor.
# --------------------------------------------------------------------------
DESIGNS = [
    # ---- Objective 1: human EGFR domain III affinity -------------------
    ("EGFR_D3_AFF_01", 1, "A", "GRTFSSYA", "ISWSGGST", "AADWGYYDSWDYEY",
     "Trp/Tyr shelf + 4 acidic contacts to the K443/K465 face"),
    ("EGFR_D3_AFF_02", 1, "A", "GRTFSDYA", "ISRSGGTT", "AAEYWSGDYYGLDY",
     "CDR2 Arg reaches the D355 subsite; acidic CDR3 tip"),
    ("EGFR_D3_AFF_03", 1, "B", "GFTFSSYA", "ITSGGGDT", "AADRGWYDFDSPYEY",
     "Longer 15-mer CDR3, Pro-kinked tip for a deeper insertion"),
    ("EGFR_D3_AFF_04", 1, "B", "GRTFSTYA", "IYWSGGTT", "AAGYSWDEYYRLDY",
     "Aromatic-dense CDR2; mixed charge tip for shape complementarity"),
    ("EGFR_D3_AFF_05", 1, "C", "GRTFSSWA", "ISWEGGST", "AAWDYGSTYYELDY",
     "CDR1 Trp anchor; Ser/Thr patch for Q384/Q408 hydrogen bonds"),
    ("EGFR_D3_AFF_06", 1, "C", "GFTFSDYA", "ISSGGSTT", "AAEGWYSDYFDLGYDY",
     "16-mer CDR3, extended paratope spanning both lobes of the groove"),
    ("EGFR_D3_AFF_07", 1, "A", "GRTFSSYA", "ISWSGGDT", "AARYWDSGSYYFDLEY",
     "16-mer, Arg/Asp pair pre-organises the loop tip"),
    ("EGFR_D3_AFF_08", 1, "B", "GRSFSTYA", "IDWSGGRT", "AAYDSWGSYYEFDY",
     "Charge-paired CDR2; Phe/Tyr core against the I467 shelf"),
    # ---- Objective 2: human + mouse cross-reactivity --------------------
    ("EGFR_D3_XRC_01", 2, "A", "GRTFSSYA", "ISWSGGST", "AAGYTWYSGSDY",
     "12-mer CDR3, compact footprint on the conserved ligand core"),
    ("EGFR_D3_XRC_02", 2, "A", "GRTFSTYA", "ISSGGSTT", "AASYGTTWYSLDY",
     "13-mer, hydrogen-bond driven, single acidic residue"),
    ("EGFR_D3_XRC_03", 2, "B", "GFTFSSYA", "IYSGGSTT", "AAGSYWSTTYGLDY",
     "Charge-free CDR3; species-invariant aromatic and polar contacts"),
    ("EGFR_D3_XRC_04", 2, "B", "GRTFSSWA", "ITWSGGST", "AATYSGWSSYYLDY",
     "Trp-rich flat paratope, minimal reliance on variable side chains"),
    ("EGFR_D3_XRC_05", 2, "C", "GRTFSTYA", "ISWSGGTT", "AAVQYSGSYWTTDY",
     "Gln mimics the conserved ligand-contact hydrogen bonding"),
    ("EGFR_D3_XRC_06", 2, "C", "GFTFSSWA", "ISSGGSTT", "AAGSTYYWSSPLDY",
     "Pro-rigidified tip; polar/aromatic only, no salt bridges"),
    # ---- Objective 3: pH-selective (binds 6.5, silent at 7.4) -----------
    ("EGFR_D3_PHS_01", 3, "A", "GRTFSHYA", "ISWSGGHT", "AAHDYWSGHYYLDY",
     "4 His spread across all three CDRs; no CDR Arg/Lys"),
    ("EGFR_D3_PHS_02", 3, "A", "GRTFSSYA", "IHWSGGHT", "AAWHSGHYDYSLDY",
     "4 His, CDR2-weighted switch adjacent to the acidic subsite"),
    ("EGFR_D3_PHS_03", 3, "B", "GFTFSHYA", "ISSGGHTT", "AAHGYSHWDHTLDY",
     "5 His for a steeper, more cooperative pH transition"),
    ("EGFR_D3_PHS_04", 3, "B", "GRTFSSHA", "ISWHGGST", "AAHSYWGHDHSLDY",
     "5 His; intramolecular His-Asp pair raises the effective pKa"),
    ("EGFR_D3_PHS_05", 3, "C", "GRTFSHWA", "IHSSGGST", "AAGHYWSHSHDLDY",
     "5 His with a Trp anchor to hold weak pH-independent affinity"),
    ("EGFR_D3_PHS_06", 3, "C", "GRTFSSYA", "IHWSGGHT", "AAHDWSHGHYYLDY",
     "5 His, His stack at the loop tip for maximal switch amplitude"),
]

OBJECTIVE_NAME = {
    1: "binding affinity (human EGFR)",
    2: "mouse cross-reactivity",
    3: "pH-selective binding",
}

AA = set("ACDEFGHIKLMNPQRSTVWY")

# Side chain pKa values (Bjellqvist / EMBOSS style) used for charge and pI.
PKA_POS = {"K": 10.53, "R": 12.48, "H": 6.00}
PKA_NEG = {"D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07}
PKA_NTERM, PKA_CTERM = 9.69, 2.34

# Kyte-Doolittle hydropathy.
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def build(fw, cdr1, cdr2, cdr3):
    f = FRAMEWORKS[fw]
    return f["fr1"] + cdr1 + f["fr2"] + cdr2 + f["fr3"] + cdr3 + f["fr4"]


def net_charge(seq, ph, free_cys=False):
    """Net charge at a given pH. Cys is excluded by default because the two
    canonical cysteines are disulfide bonded and carry no titratable thiol."""
    q = 1.0 / (1.0 + 10 ** (ph - PKA_NTERM))
    q -= 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))
    for aa in seq:
        if aa in PKA_POS:
            q += 1.0 / (1.0 + 10 ** (ph - PKA_POS[aa]))
        elif aa in PKA_NEG:
            if aa == "C" and not free_cys:
                continue
            q -= 1.0 / (1.0 + 10 ** (PKA_NEG[aa] - ph))
    return q


def isoelectric_point(seq):
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def gravy(seq):
    return sum(KD[a] for a in seq) / len(seq)


def liabilities(seq, cdrs):
    """Developability flags. Motifs are only reported when they fall in or
    next to a CDR, where they are actually exposed and mobile."""
    flags = []
    if len(seq) < 10 or len(seq) > 250:
        flags.append(f"LENGTH_OUT_OF_RANGE({len(seq)})")
    bad = sorted(set(seq) - AA)
    if bad:
        flags.append("NON_CANONICAL:" + "".join(bad))
    if seq.count("C") != 2:
        flags.append(f"CYS_COUNT={seq.count('C')}")
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST":
            flags.append(f"NGLYC_SEQUON@{i + 1}:{seq[i:i + 3]}")
    for cdr in cdrs:
        for motif, label in (("NG", "DEAMIDATION"), ("NS", "DEAMIDATION_WEAK"),
                             ("DG", "ISOMERISATION"), ("DP", "ACID_CLEAVAGE")):
            if motif in cdr:
                flags.append(f"{label}:{motif}@CDR")
        if "M" in cdr:
            flags.append("MET_OXIDATION@CDR")
        if "W" in cdr and cdr.count("W") > 3:
            flags.append("TRP_OXIDATION_CLUSTER@CDR")
    return flags or ["none"]


def main():
    out = os.path.dirname(os.path.abspath(__file__))
    rows, qc = [], []
    seen = {}

    for name, obj, fw, c1, c2, c3, note in DESIGNS:
        seq = build(fw, c1, c2, c3)
        if seq in seen:
            raise SystemExit(f"duplicate sequence: {name} == {seen[seq]}")
        seen[seq] = name

        cdrs = [c1, c2, c3]
        cdr_his = sum(c.count("H") for c in cdrs)
        fw_his = seq.count("H") - cdr_his
        # Fraction of an imidazole protonated at each assay pH.
        f65 = 1.0 / (1.0 + 10 ** (6.5 - PKA_POS["H"]))
        f74 = 1.0 / (1.0 + 10 ** (7.4 - PKA_POS["H"]))
        switch = cdr_his * (f65 - f74)

        rows.append(dict(name=name, sequence=seq, molecule_class="nanobody"))
        qc.append(dict(
            name=name,
            objective=f"{obj}:{OBJECTIVE_NAME[obj]}",
            framework=fw,
            length=len(seq),
            cdr1=c1, cdr2=c2, cdr3=c3,
            cdr3_len=len(c3),
            cys=seq.count("C"),
            cdr_his=cdr_his,
            framework_his=fw_his,
            paratope_charge_gain_6p5_vs_7p4=round(switch, 2),
            net_charge_pH6p5=round(net_charge(seq, 6.5), 2),
            net_charge_pH7p4=round(net_charge(seq, 7.4), 2),
            pI=round(isoelectric_point(seq), 2),
            gravy=round(gravy(seq), 3),
            liabilities=";".join(liabilities(seq, cdrs)),
            rationale=note,
        ))

    with open(os.path.join(out, "submission.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(out, "qc_report.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(qc[0]), delimiter="\t")
        w.writeheader()
        w.writerows(qc)

    clean = sum(1 for r in qc if r["liabilities"] == "none")
    print(f"{len(rows)} designs; lengths "
          f"{min(r['length'] for r in qc)}-{max(r['length'] for r in qc)} aa; "
          f"{clean}/{len(qc)} with no developability flags")
    for r in qc:
        if r["liabilities"] != "none":
            print(f"  FLAG {r['name']}: {r['liabilities']}")


if __name__ == "__main__":
    main()
