# Task 46E addendum: CBM structural screening before new kinetic model discovery

This addendum is mandatory and supersedes any Task 46 instruction that would begin new kinetic/ODE model construction before the constraint-based screening described here is complete.

## Objective

Before proposing or fitting any new detailed transporter kinetics, build a compact constraint-based model (CBM) of the salivary acinar transport network and use it to determine which flux relationships, capacities and couplings are structurally required to reproduce WT physiology and the AE4-loss phenotype.

The CBM is a fast structural screening layer. It is not the final biological model and must not replace the full dynamic model. Its purpose is to reduce the kinetic model search space and prevent combinatorial model proliferation.

## Starting point

Resume Task 46 from the latest remotely published checkpoint. The Task 44 / checkpoint-01 baseline is authoritative. Do not rerun or re-verify the old dynamic model merely to begin this stage.

Use the preserved model, Tasks 43-44 and the technical report only to define:
- transporter stoichiometries;
- compartment topology;
- experimentally supported flux directions and capacities;
- conservation laws;
- the known failure signature of the old model.

Do not inherit the old kinetic laws as CBM assumptions unless the constraint itself is independently justified.

## Scope of the CBM

Construct the smallest constraint system that preserves the mechanisms relevant to AE4 loss.

At minimum represent conserved balances for:
- Na;
- K;
- Cl;
- total inorganic carbon (TIC);
- total alkalinity (TA);
- charge/electroneutrality where appropriate.

Represent the principal fluxes explicitly, including where relevant:
- NKCC1;
- AE4;
- AE2;
- NHE1;
- NBC-like bicarbonate entry;
- Na/K ATPase;
- apical Cl exit / CaCC;
- K conductances or effective K exit;
- CO2 exchange;
- paracellular transport;
- luminal secretion/outflow bookkeeping.

Use exact transporter stoichiometries. If ATP/energy supply must be represented to prevent impossible pump solutions, include the minimum necessary energetic constraint rather than inventing a metabolic network.

The CBM may use steady-state or quasi-steady-state amount balances. It must state clearly which dynamic storage terms have been set to zero and which conclusions therefore apply only to sustained flux structure.

## Mathematical formulation

Build an explicit stoichiometric/constraint representation of the form

    S v = 0

with clearly documented lower and upper bounds

    l <= v <= u.

Every row and column of S must have a biological interpretation.

Where a fixed physiological state is used to infer feasible transporter direction, document the state and affinity calculation separately from the linear program.

For thermodynamic filtering, use experimentally/physiologically supported concentration ranges or the pinned reference state and impose only sign/direction constraints that are actually justified. Do not pretend that a linear CBM proves global thermodynamic feasibility over unknown concentration ranges.

## Paired WT / AE4-null formulation

Formulate WT and AE4-null conditions jointly rather than as unrelated optimisations:

    S v_WT = 0
    S v_KO = 0
    v_AE4,KO = 0.

Where transporter capacities/abundances are not known to change with genotype, share the same capacity parameter or bound across WT and KO.

In particular, do not give AE4-null NKCC1 an enlarged capacity merely because the intracellular state changed in the old kinetic model.

Treat the isolated salivary NKCC1 result as a principal structural discriminator: the accepted CBM should not require substantial genotype-dependent NKCC1 up-regulation to explain AE4 loss unless protocol-matched evidence supports it.

## Required CBM analyses

Perform the following in this order.

### CBM milestone 1: network construction and verification

Create:
- reaction/transport ledger;
- S matrix in machine-readable form;
- conservation/charge audit;
- bound/provenance table;
- rank and nullspace dimension;
- exact WT/KO shared-capacity definitions.

Do not optimise phenotype yet.

Publish a checkpoint immediately.

### CBM milestone 2: WT feasible-space analysis

Determine whether the WT network admits physiologically plausible sustained flux states under the declared constraints.

Perform flux variability analysis (FVA) for the principal transporters and secretion/output fluxes.

Identify which WT fluxes are tightly constrained and which remain underdetermined.

Do not select a single arbitrary WT optimum unless an experimentally justified objective exists.

Publish a checkpoint immediately.

### CBM milestone 3: AE4-null feasible-space analysis without compensatory assumptions

Set AE4 flux to zero while retaining the same non-genotype-specific transporter capacities.

Determine:
- maximum feasible sustained apical chloride export / secretion proxy;
- minimum and maximum feasible NKCC1 flux;
- feasible pH/alkalinity-support relationships insofar as represented by TIC/TA constraints;
- whether the network can structurally preserve near-WT secretion without increased NKCC1 capacity/activity.

Quantify the best possible KO compensation allowed by the stated biology.

Publish a checkpoint immediately.

### CBM milestone 4: paired WT-KO coupling and NKCC constraint

Add the strongest defensible cross-genotype NKCC constraint supported by the isolated assay/provenance review, expressed as an equality or tolerance only if justified.

For example, if supported:

    v_NKCC,KO <= (1 + epsilon) v_NKCC,WT

with epsilon taken from actual experimental uncertainty/tolerance, not invented.

If no numerical epsilon is justified, use qualitative analyses with no capacity increase and clearly distinguish flux variation caused by altered state from altered transporter capacity.

Determine the maximal KO secretion compatible with the NKCC evidence.

Publish a checkpoint immediately.

### CBM milestone 5: flux coupling and structural redundancy

Analyse only biologically relevant pairwise or small-block flux couplings.

Determine whether AE4 and NKCC1 are genuinely interchangeable once TIC/TA, pump, chloride exit and charge constraints are all included.

Identify exact or interval coupling relationships among:
- AE4 and NKCC1;
- AE4 and NBC/NHE1 alkalinity support;
- total basolateral chloride loading and apical chloride exit;
- pump demand and secretion;
- any other relationship directly implicated by the failure diagnosis.

Do not enumerate all elementary flux modes unless the network is demonstrably tiny and the enumeration is trivial. Exhaustive EFM enumeration is not required and should be avoided if combinatorial.

Publish a checkpoint immediately.

### CBM milestone 6: minimal-relaxation / minimal-repair analysis

Ask the inverse structural question:

What is the smallest set of currently assumed flux bounds/couplings that must change for the CBM to admit both:
- physiological WT secretion; and
- the observed substantial AE4-null secretion deficit without strong NKCC rescue?

Use a sparse/lexicographic relaxation analysis rather than exhaustive combinations.

Prefer:
1. exact feasibility logic;
2. one-at-a-time bound relaxations;
3. L1 or minimal-cardinality relaxation only if necessary.

Do not perform a Cartesian product of transporter variants or arbitrary parameter grids.

Return a ranked list of the smallest mechanistic changes indicated by the CBM.

Publish a checkpoint immediately.

## Anti-combinatorial discipline

This CBM stage exists specifically to reduce model-search complexity.

Therefore:
- build one curated CBM network, not many architecture combinations;
- use FVA and coupling analysis instead of sampling huge parameter grids;
- change one structural constraint at a time unless infeasibility proves that more than one change is jointly necessary;
- do not explore all subsets of candidate mechanisms;
- do not enumerate large flux-mode families;
- do not optimise dozens of arbitrary objectives;
- do not turn the CBM into another fitted model.

The output of the CBM must be a small set of necessary or strongly indicated mechanistic requirements for the later kinetic model.

## Decision gate before kinetic model discovery

Do not proceed to new ODE/kinetic candidate construction until the CBM stage has produced and published a written decision report containing:

1. whether WT physiology is structurally feasible;
2. whether AE4-null near-WT secretion is structurally feasible under fixed/shared transporter capacities;
3. whether strong NKCC compensation is structurally required, optional or excluded;
4. which constraints make AE4 and NKCC1 non-equivalent, if any;
5. which alkalinity/carbon constraints limit AE4 flux;
6. the maximal AE4-null secretion allowed under the independent NKCC evidence;
7. the smallest ranked set of mechanistic changes needed to reconcile the network with the phenotype;
8. exactly which kinetic block(s) should therefore be rebuilt first.

This decision report becomes the search prior for the subsequent new-model discovery stage.

If the CBM shows that the phenotype cannot be obtained within a single-cell steady transport architecture without violating established constraints, stop and publish that result before inventing additional kinetic detail. Identify the missing class of biology instead.

## Required deliverables

Under the Task 46 analysis directory create at minimum:
- `CBM_PROTOCOL.md`;
- reaction/transport ledger;
- machine-readable S matrix;
- bound/provenance table;
- WT FVA output;
- AE4-null FVA output;
- paired WT/KO analysis;
- NKCC-constrained analysis;
- flux-coupling results;
- minimal-relaxation results;
- `CBM_DECISION_REPORT.md`;
- compact scripts sufficient to reproduce every CBM result.

## Checkpoint discipline

Each CBM milestone above is a separate remote checkpoint.

At the end of each milestone:
1. fetch remote state;
2. fast-forward only;
3. verify the milestone;
4. update `CURRENT_STATUS.md`;
5. create immutable `CHECKPOINT_XX.md` or clearly namespaced CBM checkpoint;
6. commit;
7. push immediately;
8. verify the remote SHA;
9. only then continue.

Never force-push.

Do not wait until all CBM analyses are complete before publishing the first results.

## Relationship to the rest of Task 46

Once the CBM decision gate is passed, resume Task 46 new-model discovery.

Use the CBM to narrow the kinetic search to the one or two mechanistic blocks actually implicated by structural necessity.

Do not reopen unrelated transporter families unless the new kinetic model fails a declared holdout and the CBM assumptions are themselves shown to be wrong.

The later full dynamic model must still retain carbon, alkalinity, volumes, lumen, pH, electrical closure, water transport and regulation as required by Task 46.
