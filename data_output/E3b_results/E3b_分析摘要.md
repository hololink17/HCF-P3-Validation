# E3b registered-analysis summary

C1 compares second-wave batch 2 (n=240) with the existing Chinese cohort (n=504).
The margin is |Cliff's delta| < 0.147; it is not a raw mean-difference bound.

| alpha | Cliff delta | 90% CI | lower p | upper p | equivalent |
|---|---:|---|---:|---:|---|
| alpha1 | 0.130721 | [0.061178, 0.200264] | 4.80095e-11 | 0.349968 | False |
| alpha2 | 0.177116 | [0.106026, 0.248207] | 1.11823e-13 | 0.757251 | False |
| alpha3 | 0.047305 | [-0.023405, 0.118015] | 3.56626e-06 | 0.0102589 | True |
| alpha4 | 0.061078 | [-0.010126, 0.132282] | 9.32342e-07 | 0.0236368 | True |

C1 is partially supported: passed=['alpha3', 'alpha4']; not passed=['alpha1', 'alpha2'].
C2 direction-only descriptive replication (both expected negative): direction_consistent=True; no registered model or significance threshold.
C3 comprises 20 items x 5 option positions, 4,800 responses; minimum cell n=2; descriptive only, with no registered sparsity failure threshold.
The batch-1 versus batch-2 mean-scale table is post hoc and is not registered C1.
E3a and order analyses are protocol-specified but unregistered and non-confirmatory.
No item text, mapping, scoring key, parameter, or backend formula is disclosed.
