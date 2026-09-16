# Task 52 novelty / no-repeat check

Task52 is not a new mechanism-family search. It is a final numerical test of the specific R51H/R51I chloride-reservoir architecture.

## What is genuinely new relative to completed tasks

The candidate differs from prior failed work in three binding ways:

1. **Measured chronic genotype chloride is retained at stimulation onset.** Earlier acute/shared-rest perturbations erased the experimentally measured WT-versus-AE4-KO reservoir difference before stimulation began.
2. **Large state-driven KO NKCC compensation is not allowed to emerge as the central mechanism.** The central paired experiment matches non-AE4 NKCC/AE2 cycle supply between genotypes, motivated by the isolated NKCC assay's lack of detected genotype up-regulation. This is an explicit matched-control idealisation, not an assertion of exact native equality.
3. **The beta/IPR auxiliary anion current is shared demand, not genotype-specific supply or gating.** It can expose an existing chloride-reservoir difference without any AE4-expression term on the channel.

## Distinction from prior tasks

- Tasks16–18: changed AE4 loading share/capacity while retaining common acute rest. Task52 does not repeat this.
- Tasks28–31: NHE/acid-base repair. Task52 does not change NHE or seek chronic-rest closure.
- Task38/40: acute/shared-WT-rest AE4 deletion. Task52 explicitly rejects shared WT chloride initialisation as the central genotype comparison.
- Tasks39–40: Palk-driven state compensation. Task52 excludes Palk.
- Task41/50: AE4-expression-dependent CaCC recruitment/effective coupling. Task52 forbids that multiplier.
- Task42: source/stoichiometric AE4 classes. Task52 leaves AE4 transport law/routing unchanged.
- Task46/47: structural rescue and K/pump limits. Task52 does not revisit those families.
- Task48: chronic genotype rest reconstruction failed under inherited laws. Task52 does not claim to solve chronic resting physiology; it uses a deterministic measurement-constrained stimulation-onset projection.
- Task49: beta*positive-swelling VRAC on beta-blind core. Task52 does not use that gate.
- Task51: beta-NKCC -> swelling -> VRAC. Task52 does not use beta-NKCC or swelling to create genotype specificity; the auxiliary current is a shared source-scale demand term.

## Hard no-repeat gate

Do not reopen:

- Palk/Benjamin NKCC;
- alternative NKCC kinetic families;
- NKCC hard caps selected from phenotype;
- calcium-amplitude changes;
- NHE/CO2 repair;
- NBC family changes;
- AE4 stoichiometry/routing/capacity searches;
- pump/K recycling limits;
- AE4-dependent TMEM16A/VRAC recruitment;
- swelling-gate families;
- transporter subset searches;
- any grid/random/global optimisation;
- exact timing fitting.

The only permitted noncentral perturbations are the explicitly frozen Task52 case sensitivities in `DECISION_LOG.md`.

## Novel scientific claim if successful

If the full model succeeds, the supported claim is not that AE4 molecularly controls an apical channel. It is:

> The measured chronic chloride depletion caused by AE4 loss, together with loss of ordinary beta-activated AE4 chloride supply, is sufficient to create a large stimulated secretory deficit when beta-associated apical chloride demand is shared and large compensatory non-AE4 chloride up-regulation is not imposed.

If it fails, no alternate mechanism is authorised. The final conclusion reverts to the preserved Task50 proof that an additional unidentified beta-conditioned network coupling is required.
