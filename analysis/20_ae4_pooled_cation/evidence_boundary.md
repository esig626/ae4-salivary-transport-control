# Task 20 evidence boundary

The primary source is Peña-Münzenmayer et al. (2016), *JGP* 147:423–436,
[doi:10.1085/jgp.201611571](https://doi.org/10.1085/jgp.201611571).
The complete [Europe PMC XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4845690/fullTextXML)
was read; its digest and inspected anchors are in `primary_source_record.json`.
The inherited evidence ledger is
`analysis/13_state_resolved_ae4/evidence_freeze.md`, entries E16-01–E16-12,
with additional detail in `analysis/12_ae4_mechanism_reconstruction/experimental_evidence.md`.

Figures 3, 4 and 8 support Na transport, macroscopic electroneutrality,
K transport under imposed gradients, and broad monovalent-cation support.
Cl/HCO3 exchange and reversibility are also supported. They do not identify
independent opposing Na/K cycles in physiological mixed solution.

Figure 10 uses a pooled-cation physiological working model. Its 1:1:2
Cl:cation:HCO3 stoichiometry is assumed. Equal Na/K weights and Task 20's
donor-side partition are modeling choices, not measured microscopic rules.
Figure 9 Hill coefficients do not identify the number of transported ions.
Task 20 tests these declared choices without fitting cation selectivity.

The old ledger's restriction on pooling governed its microscopic state-model
task. The current Task 20 prompt explicitly authorizes a separate coarse-grained
pooled alternative. The ledger and legacy implementation remain unchanged.

The primary paper's illustrative bicarbonate concentration cannot replace the
production bicarbonate obtained from each inherited conserved WT state. Task 20
keeps those states, acid-base chemistry and bath fixed. No genotype data enter
the transport-law tests, feasibility certificate or capacity objective.
