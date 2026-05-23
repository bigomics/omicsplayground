### {{contrast}} — {{tier}}

Signal overview (drug counts, direction split, max |NES|):
Drugs tested: {{n_drugs_tested}} | significant (q<0.05): {{n_drugs_sig}} | opposing (neg NES): {{n_neg}} | mimicking (pos NES): {{n_pos}} | max |NES|: {{max_abs_nes}}

Annotation reliability (fraction of drugs carrying MOA / target metadata, plus significance counts at the MOA-class and target levels):
{{frac_annotated_pct}} | sig MOA classes: {{n_sig_moa_classes}} | sig targets: {{n_sig_targets}}

Interpretation evidence summary (primary evidence for report writing):
- Supported opposing MOA terms: {{supported_opposing_moa}}
- Supported mimicking MOA terms: {{supported_mimicking_moa}}
- Corroborating opposing targets: {{corroborating_opposing_targets}}
- Corroborating mimicking targets: {{corroborating_mimicking_targets}}
- Preferred opposing exemplars matching supported MOA terms: {{preferred_opposing_exemplars}}
- Preferred mimicking exemplars matching supported MOA terms: {{preferred_mimicking_exemplars}}
- Annotation confidence: {{annotation_confidence}} ({{frac_annotated_pct}} annotated)

Top opposing drugs (negative NES — candidate signature reversal for L1000 backends; predicted resistance for sensitivity backends):
{{top_opposing_table}}

Top mimicking drugs (positive NES — mechanistic similarity for L1000 backends; predicted vulnerability for sensitivity backends):
{{top_mimicking_table}}

MOA class enrichment (pharmacological mechanisms converging on this signature; fgsea over the drug MOA annotations of the DSEA table):
{{moa_class_table}}

Target enrichment (molecular targets driving the connectivity; fgsea over the drug target annotations of the DSEA table):
{{moa_target_table}}
