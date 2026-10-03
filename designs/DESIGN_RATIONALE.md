# Proteinbase Protein Design Competition — Challenge 1

Anti-EGFR binder set. 20 designs, submitted as `designs/submission.csv`.

## Target

Human EGFR extracellular region, UniProt **P00533-1** residues 25–645
(621 aa), structure reference **6ARU chain A**. The recommended epitope is
**domain III**, the lobe that carries the EGF/TGF-α binding surface and the
epitope of cetuximab and panitumumab.

Two numbering schemes appear in the EGFR literature and it is easy to mix
them up. The mature protein starts at P00533 residue 25, so
*mature = UniProt − 24*. The cetuximab contact residues usually quoted as
Q384, Q408, H409, K443, K465, I467, S468 are mature numbering; the famous
cetuximab-resistance mutation S492R is quoted in UniProt numbering and is the
same residue as mature S468. Every epitope residue below is given in mature
numbering.

### Epitope chemistry the designs are complementary to

The domain III ligand-binding face is a shallow groove with three sub-sites:

| Sub-site | Residues (mature) | Paratope answer |
|---|---|---|
| Hydrophobic shelf | F357, L382, I467 | Tyr / Trp / Phe packing |
| Basic rim | K443, K465 | Asp / Glu salt bridges |
| Polar floor | D355, Q384, Q408, H409, S468, N473 | Ser / Thr / Tyr-OH hydrogen bonds |

So the generic paratope recipe is **aromatic + acidic + polar**, which is also
what cetuximab does (its CDR-H3, ALTYYDYEFAY, is exactly aromatic-plus-acidic).
None of the designs reuse cetuximab or panitumumab CDR sequence.

## Scaffold choice

All 20 designs are **humanised camelid VHH single domains (nanobodies)**,
118–122 aa, `molecule_class = nanobody`.

This is a deliberate risk decision. This container has no GPU and no network
access to structure-prediction weights, so nothing here was validated by
RFdiffusion / ProteinMPNN / AlphaFold. A hand-written de novo mini-binder with
no structural validation would be unlikely to fold at all. A VHH localises all
of the uncertainty into the three CDR loops: the framework is a canonical
immunoglobulin fold that will fold and express regardless of what the loops do,
so a design that misses the epitope still produces clean, measurable
no-binding data rather than insoluble aggregate.

Three canonical VHH frameworks (A, B, C) are used to spread the set across
sequence space. All three keep:

- the VHH hallmark tetrad F42 / E49 / R50 / G52 in FR2, which replaces the
  V<sub>H</sub>–V<sub>L</sub> interface and is what makes a single domain soluble;
- the single canonical C22–C96 disulfide, and **no third cysteine** anywhere, so
  there is no scrambling or disulfide-mediated dimerisation risk.

## Objective 1 — binding affinity (8 designs, `EGFR_D3_AFF_*`)

Aromatic-dense, mildly acidic paratopes built on the recipe above. CDR3 lengths
run 14–16 aa, which is the range where VHH CDR3 is long enough to reach into
the domain III groove rather than sitting on the flat surface around it. The
set varies CDR3 length, charge placement and which sub-site the CDR2 reaches,
so the eight designs sample the groove from different approach angles rather
than being eight copies of one idea.

## Objective 2 — mouse cross-reactivity (6 designs, `EGFR_D3_XRC_*`)

Human and mouse EGFR ectodomains are roughly 80–85% identical, and cetuximab
and panitumumab are both human-restricted — the species difference is
concentrated in the peripheral 465–475 loop (mature S468 / UniProt S492 is the
best known determinant), not in the ligand-binding core. The core itself is
conserved, which is why EGF and TGF-α engage both receptors.

The cross-reactive designs therefore:

- use a **shorter CDR3 (12–14 aa)**, giving a smaller footprint that stays over
  the conserved core (D355, F357, Q384, Q408, H409) instead of spreading onto
  the variable periphery;
- are **hydrogen-bond and aromatic driven with salt bridges depleted**, because
  a charged contact is the kind most likely to be lost to a species
  substitution, whereas backbone and aromatic contacts tolerate it;
- avoid contacting the 465–475 region at all.

**Caveat:** this environment's network policy blocked UniProt, NCBI, RCSB and
EBI, so this reasoning rests on published EGFR structural biology rather than a
residue-by-residue human/mouse alignment computed here. If you can run that
alignment (mouse EGFR is UniProt Q01279), it is worth re-checking which
domain III surface positions actually differ before the next submission round.

## Objective 3 — pH-selective binding (6 designs, `EGFR_D3_PHS_*`)

Requirement: bind human EGFR at **pH 6.5** (tumour microenvironment), no
detectable binding at **pH 7.4** (circulation).

This is the inverse of the Chugai-style recycling antibody, and the mechanism
is the same histidine switch run in the other direction. Histidine is the only
side chain with a pKa (~6.0) between the two assay points: it is ~24% protonated
at pH 6.5 but only ~4% at pH 7.4. Placing histidines so that *protonation
creates* a favourable contact with the acidic domain III floor means the
paratope switches on as the pH drops.

Design rules applied:

1. **4–5 His per design, all in CDRs, none in the frameworks.** Several
   histidines make the transition cooperative, so affinity falls off steeply
   across the 0.9 pH unit window instead of drifting gently. The QC table
   reports the paratope charge gained between pH 7.4 and 6.5: **+0.8 to +1.0
   per design.**
2. **No Arg or Lys in the CDRs.** If the paratope carried constitutive positive
   charge it would bind at both pH values and the switch would be swamped.
   Binding energy has to be histidine-dependent to be pH-dependent.
3. **His placed next to Asp in the same loop** in several designs — a
   neighbouring carboxylate raises the effective imidazole pKa, shifting the
   transition midpoint up towards 6.5 and sharpening the on/off contrast.
4. **A retained Trp/Tyr anchor**, so that at pH 7.4 the complex is weak rather
   than the binder being non-functional at both pH values. The aim is a
   measurable affinity window, not a dead design.

## Developability QC

`build_designs.py` regenerates the submission and re-runs the checks; all 20
designs pass clean. Checked: chain length within the 10–250 aa rule, canonical
amino acids only, exactly two cysteines, no N-glycosylation sequons (N-X-S/T),
no Asn deamidation (NG/NS) or Asp isomerisation (DG) or acid-cleavage (DP)
motifs in the CDRs, no methionine oxidation sites in the CDRs, plus pI, net
charge at both assay pH values and GRAVY hydrophobicity. One design
(`AFF_05`) originally carried a DG isomerisation motif in CDR2 and was
corrected to EG, preserving the charge.

See `qc_report.tsv` for the per-design table.

## Known limitations

- **No structural validation.** No docking, no folding prediction, no
  interface scoring. These are knowledge-based designs. A competitive entry
  would run RFdiffusion or BindCraft against the 6ARU domain III surface,
  redesign sequence with ProteinMPNN, and filter on AlphaFold2 pAE_interaction
  and ipTM. That pipeline needs a GPU and model weights neither of which this
  container has.
- **Epitope targeting is unenforced.** Nothing here guarantees a design lands on
  domain III rather than elsewhere on the 621-residue ectodomain.
- **Novelty.** The novelty filter runs automatically on upload. The CDRs are
  original and three different frameworks are used, but nanobody frameworks are
  germline-derived and shared across all VHHs by construction. If any design is
  rejected below 3/4, the fix is to diversify the framework further rather than
  the loops.
