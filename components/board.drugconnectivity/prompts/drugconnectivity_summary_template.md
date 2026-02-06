## Task

Interpret the drug connectivity analysis results for the contrast **{contrast}** based on the drug enrichment data provided below.

## Analysis Overview

**Contrast:** {contrast}

**Analysis type:** {method}

The drug connectivity analysis correlates the experimental gene expression signature for this contrast with known drug perturbation profiles from the L1000 database. Positive NES indicates drugs that mimic the experimental signature; negative NES indicates drugs that oppose it.

## Top Connected Drugs

The following drugs show the strongest connectivity with the experimental signature:

{top_drugs}

Focus on identifying:
1. The most promising drug candidates (especially those opposing the disease signature)
2. Common pharmacological themes among the top drugs
3. Whether the top drugs suggest specific targetable pathways

## Mechanism of Action Analysis

{moa_summary}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{contrast}** | Drug connectivity summary

Then synthesize the above information to explain:
- The most significant drug-experiment connections and their potential therapeutic implications
- How the top drugs relate to each other through shared mechanisms of action or molecular targets
- Which drug classes or targets are most prominently enriched and what this reveals about the biology
- Potential drug repurposing candidates that oppose the experimental signature
- Use quantitative metrics to support your interpretation (cite specific NES, q-values, and MOA statistics when available)
