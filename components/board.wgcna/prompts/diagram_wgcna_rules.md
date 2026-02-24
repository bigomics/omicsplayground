## WGCNA Diagram Rules

### Output format
- Return ONLY valid JSON. No prose before or after.
- Edges must use only "positive" or "negative" regulation -- never "association".

### Module node labels
- Format: `[ModuleID]: biology` where biology is 2-5 words describing the module's function
- Example: `[turquoise]: cell cycle regulation`
- Do NOT use full gene names in node labels -- genes go in the `genes` field only

### Edge strength calibration
- Strong co-expression or known direct interaction -> strength 0.8-1.0
- Moderate relationship or inferred -> strength 0.4-0.7
- Weak or speculative -> strength 0.1-0.3

### Connectivity
- Each node must appear in at least one edge
- Prefer 6-12 edges total for readability
- Hub nodes (high connectivity) should be phenotype or key process nodes

### ASCII only
- All text fields must use plain ASCII characters only
- No curly quotes, em-dashes, accented characters, or Unicode symbols
