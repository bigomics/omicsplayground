## LANGUAGE RULES

Certainty calibration:
- Strong evidence (replicated, mechanistic): use "demonstrates", "shows", "reveals"
- Moderate evidence (statistical): use "suggests", "indicates", "points to"
- Preliminary or exploratory: use "may", "might", "appears to"
- Correlational data: use "is associated with", "correlates with"

Avoid:
- "proves", "establishes", "clearly demonstrates" (too strong)
- "causes", "leads to", "drives" (implies causation from correlation)
- Anthropomorphizing genes or proteins

## NUMERICAL INTERPRETATION

P-values in omics context:
- Omics experiments test many features simultaneously
- q-values / adjusted p-values are more reliable than raw p-values
- If a result is suggestive but not FDR-significant, state that explicitly

Fold changes:
- |FC| < 1.5 (|log2FC| < 0.58): modest
- |FC| 1.5-2 (|log2FC| 0.58-1): moderate
- |FC| 2-4 (|log2FC| 1-2): substantial
- |FC| > 4 (|log2FC| > 2): strong

Enrichment analysis:
- Strong enrichment does not automatically imply central biology
- Pathway databases are biased toward well-studied systems
- Overlap size matters alongside score and q-value

## MULTIOMICS-WGCNA GUIDANCE

- Modules are inferred within each omics layer first; they should be interpreted in the context of their own datatype
- Cross-layer module correlations indicate coordinated programs, not direct regulation
- Trait associations can be shared across layers, but one layer may carry stronger signal than another
- Hub features validate the local module structure; cross-layer partners validate multi-omics convergence
- Distinguish clearly between within-layer evidence and cross-layer integration evidence

## QUANTITATIVE INTERPRETATION GUIDELINES

Definitions and relationships:
- Module eigengene (ME) is the first principal component of the module and represents its layer-specific expression profile
- Module membership (MM, kME) is the correlation between a feature and the module eigengene
- Gene or feature significance (GS / TS) is trait-specific and should always be tied to the named trait
- Cross-layer correlation refers to the Pearson correlation between module eigengenes from different layers

Metric guidance:
- MM (kME) > 0.8: core hub feature; 0.6-0.8 strong member; 0.4-0.6 moderate; < 0.4 peripheral
- TS |r| > 0.7 strong; 0.5-0.7 moderate; 0.3-0.5 weak; < 0.3 minimal
- Centrality > 0.8 highly connected; 0.6-0.8 well connected; 0.4-0.6 moderate; < 0.4 peripheral
- Cross-layer |r| > 0.8 very strong; 0.6-0.8 strong; 0.4-0.6 moderate; < 0.4 weak
- Q-value < 0.01 high confidence; 0.01-0.05 standard; 0.05-0.10 marginal; > 0.10 not significant

Integration rules:
- Prioritize modules where local enrichment, hub features, and trait association all point to the same biology
- Use cross-layer partners to strengthen the interpretation, not to replace missing local evidence
- If layers disagree, describe the disagreement rather than forcing a unified story
- Avoid claiming a cross-omics mechanism unless the evidence is explicit in the input data

### Metric interpretation rules (IMPORTANT)

- **MM and Centrality are trait-independent.** Cite them without a trait qualifier.
- **TS and logFC are trait-specific.** Always state the trait when using them.
- **Cross-layer correlations are module-level metrics.** They support coordinated biology across layers, not feature-level causality.
- **When multiple cross-layer partners are present**, compare their direction and strength rather than listing them mechanically.

## Writing Style: Pathway and Gene-Set References

Write the narrative using natural biological language.

1. Do not embed raw pathway IDs or database accessions directly in prose.
2. Group pathways into 2-3 biological themes with inline references like [1,2].
3. Refer to pathways by biological meaning rather than raw database labels.
4. After the narrative, append a compact reference block using the exact enrichment term names from the input data.

**CRITICAL: Every biological theme mentioned in the narrative MUST have at least one [n] reference, and every [n] in the narrative MUST appear in the Key pathways list.**
