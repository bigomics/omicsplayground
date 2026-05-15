# WGCNA context

## Definitions

- **WGCNA (Weighted Gene Co-expression Network Analysis)**: an unsupervised
  method that groups features (genes, proteins, etc.) into modules based on
  pairwise expression similarity, then summarizes each module by its eigengene
  (first principal component) and correlates that eigengene against sample
  traits.
- **Module**: a cluster of co-expressed genes identified by hierarchical
  clustering of pairwise gene-gene correlations on the normalized expression
  matrix. Each module has a color name (`MEturquoise`, `MEblue`, ...), a member
  gene list, an eigengene (the first principal component of its members), and
  zero or more correlated sample traits.
- **Eigengene**: the first principal component of a module's expression
  matrix; a single per-sample value summarizing the module's expression
  profile.
- **Module-trait correlation**: Pearson correlation between a module's
  eigengene and a sample trait (clinical variable, treatment, time point, etc.).
  In this prompt's data block, the raw `r` is replaced with one of seven
  verbal labels: `strongly correlated` / `correlated` / `moderately correlated`
  / `weakly associated` / `moderately anti-correlated` / `anti-correlated` /
  `strongly anti-correlated`.
- **Hub gene**: a gene with high module membership (`MM`, the correlation
  between the gene's expression and the module eigengene). The most central
  members of a module — useful as anchors for biological interpretation.
- **Module enrichment**: over-representation of a curated gene set
  (Hallmark, Reactome, GO, KEGG, ...) among the module's members. Verbal
  mapping in the data block: `highly significant` (`q < 1e-6`) /
  `significant` / `nominally significant` / `not significant`.
- **Strong / Moderate / Weak signal**: classification provided in the data
  block, derived from significant-enrichment count and trait-correlation
  strength. **Authoritative — do not reclassify.**
- **Top modules**: the subset of modules selected for detailed reporting.
  Out of `n_total_modules` in the network, only `n_top_in_report` are
  described in detail; the remainder are aggregated under "Minor units".
- **Grey module**: the bin of unassigned genes (genes that did not pass
  the inclusion threshold for any module). Always mentioned once in
  Overview, never given its own subsection.

## How to read trait directionality

The "Trait coordination" line in each module's detail block encodes both
directions explicitly: ↑ names the trait where the module is elevated, ↓
names the trait where it is suppressed. State the biological direction in
plain terms — never infer it from trait names alone.

## Hub gene reporting

Always **group functionally and weave** — never present an inventory.
Good: *"replication factors **MCM2**, **MCM4**, **MCM7** alongside the
licensing kinase **CDK1**"*.
Bad: *"Top hub genes: MCM2, MCM4, MCM7, CDK1, RRM1, FEN1, TK1"*.
Each gene mention should travel with its functional role. If the role is not
known from the data block, do not name the gene. **Cap: ≤8 gene names per
paragraph; ≤5 per single sentence.**

## Three rhetorical patterns to use

At least one per Integrated findings section:

1. **Trade-off framing** for opposing units (when two units are anti-correlated):
   *"As the cell-cycle machinery ramps up, the quiescence program is actively
   dismantled."*
2. **Feedback / handoff framing** for sequential or upstream-shared units:
   *"Translation peaks first, proliferation follows, effector maturation closes
   the loop."*
3. **Contrast framing** between data-derived and known biology:
   *"Although MEblue is dominated by replication factors, the absence of E2F-
   target enrichment suggests a non-canonical proliferative axis."*

## Pathway naming

Pathways always in plain English. *"Ribosome biogenesis"*, not `GO_0042254`.
*"Mitotic spindle organization"*, not `R-HSA-68886`. Cluster related pathways
into a 2–3 word **theme** and narrate the theme — not the individual gene-set
identifiers. The data block carries the exact identifiers; the prose carries
the biology.

## Cross-module narrative hierarchy

A WGCNA report's value is in narrating the interplay between
co-expression programs, not in cataloguing them. Follow this hierarchy
when structuring the Integrated findings section:

1. **Opposing programs (anti-correlated module pairs) are the backbone.**
   Every strong anti-correlation (`anti-correlated` or stronger between
   modules) MUST be described as a trade-off, not just listed.
2. **Trait anchoring shows directionality.** Each module's correlation to
   a trait (from the verbal data block) tells the reader which programs
   are activated vs. suppressed. Always state the *biological direction*
   explicitly, never just the verbal label.
3. **Hub gene validation grounds the biology.** When selecting genes to
   weave into a narrative, prefer recognizable canonical markers over
   high-MM but obscure ones. A canonical marker with moderate MM is more
   convincing than a novel gene with the highest MM.
4. **Cross-module pathway convergence (only when present).** If the same
   pathway appears enriched in two modules with opposite trait directions,
   highlight this as an actionable trade-off.

## Writing rules

- Honest about absence: zero-enrichment modules say so up front; low-signal
  datasets say so in Overview.
- Numbers in the data block, not in the prose (with the rare exception of
  the lead module's lead trait correlation).
- See `## Hub gene reporting`, `## Pathway naming`, and `## Three rhetorical
  patterns to use` above for the WGCNA-specific prose conventions.
