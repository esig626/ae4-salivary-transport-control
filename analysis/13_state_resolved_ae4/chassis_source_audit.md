# Source audit of the whole-cell chassis

## Scope and controlling conclusion

This report audits the fixed Task-12 seven-state chassis against independent
whole-cell transport evidence. It does **not** fit any AE4-null secretion
measurement, choose an AE4 molecular cycle, or claim that a chassis extension
reconstructs the phenotype. The 2015 AE4-null secretion magnitude and time
course remain held out.

The smallest independently supported topology difference is a **coupled
apical cation-current module**:

1. a nonzero apical Ca-activated K conductance; and
2. a nonzero apical Na/K-ATPase, with the same `3 Na out : 2 K in` pump
   orientation as at the basolateral membrane.

The present seven-state chassis has neither. It places all Na/K-pump cycles
and all Ca-activated K efflux at the basolateral membrane. Adding the apical
copies changes the luminal Na and K balances and the apical/basolateral current
closures even if total cellular pump capacity and total K conductance are held
fixed. This difference is directly motivated by Almássy et al. 2012 and 2018,
but those experiments used **mouse parotid** acinar cells, whereas the held-out
AE4 phenotype was measured in **mouse submandibular** glands. The topology is
therefore a source-supported salivary-chassis hypothesis, not yet a directly
validated submandibular fact.

The next smallest source-supported correction is to replace the historical
fitted CO2 quasi-steady formula with a mass-conserving carbon/pH module and to
allow the existing basolateral NHE1 flux to respond to muscarinic stimulation,
cell shrinkage, and acid load using WT evidence. The 2015 AE4 study found no
detectable genotype change in its NHE-dependent alkalinization assay, so an
AE4-null-specific NHE1 multiplier is not admissible. Likewise, it found no
detectable AE4-null increase in isolated NKCC1 activity, so arbitrary NKCC1
upregulation is not an admissible rescue.

## Evidence classes and source boundary

The audit uses the following hierarchy.

- **Direct native evidence:** localization, ion current, ion-sensitive dye,
  gland secretion, or knockout measurements in a stated gland and species.
- **Source-supported mechanism:** transporter stoichiometry, polarity, or
  coupling implied by direct measurements and conservation.
- **Published model result:** a numerical distribution, optimum, or response
  produced by a mathematical model; useful as a prior or sensitivity range,
  not a measurement.
- **Historical implementation evidence:** the internally closed seven-state
  branch reconstructed in Tasks 11-12; useful as a comparison object, not
  biological ground truth.
- **New modelling decision:** a proposed minimal equation or switch below.

Key independent sources are:

1. Peña-Münzenmayer et al. 2015, mouse submandibular gland and native acini,
   [DOI 10.1074/jbc.M114.612895](https://doi.org/10.1074/jbc.M114.612895).
   This is the controlling source for AE4-null chloride, pH, AE2, NKCC1, NHE,
   and held-out secretion observations.
2. Almássy et al. 2018, mouse parotid acinar cells,
   [DOI 10.1007/s00424-018-2109-0](https://doi.org/10.1007/s00424-018-2109-0)
   and the [institutional record](https://real.mtak.hu/85486/). Local apical
   Ca photolysis evoked K current, and Na/K-pump labelling appeared throughout
   the plasma membrane, including the apical pole. Its proposed apical pump
   function and secretion simulations are interpretations/model results.
3. Almássy et al. 2012, mouse parotid acinar cells,
   [DOI 10.1085/jgp.201110718](https://doi.org/10.1085/jgp.201110718).
   Local photolysis and channel blockers directly support functional apical BK
   and IK conductances. An apical Ca flash increased whole-cell K current from
   `1.10 +/- 0.11 nA` to `1.49 +/- 0.16 nA` (`n=12`, `p<0.001`), a reported
   fractional increase of `1.34 +/- 0.02`. This quantifies a functional apical
   contribution, but it does not measure the absolute apical fraction of total
   K conductance.
4. Harmer et al. 2005, mouse submandibular acinar cells,
   [PubMed 15308468](https://pubmed.ncbi.nlm.nih.gov/15308468/). In this
   gland-matched localization experiment, a local apical Ca signal activated
   Cl but not K current; K activation followed broader Ca spread. This is a
   direct warning against treating the parotid apical-K result as a calibrated
   submandibular fraction.
5. Maruyama et al. 1983,
   [DOI 10.1038/302827a0](https://doi.org/10.1038/302827a0), for direct
   basolateral Ca/voltage-activated K-channel evidence in mammalian salivary
   acinar membrane.
6. Park et al. 2001, mouse parotid gland,
   [DOI 10.1074/jbc.M102901200](https://doi.org/10.1074/jbc.M102901200).
   NHE1 was basolateral; its loss suppressed more than 95% of acinar pH
   recovery and reduced 50-minute pilocarpine-stimulated saliva by 34%, with a
   larger late than early effect. These data establish NHE1 importance but do
   not license fitting the AE4-null time course to an NHE1 change.
7. Evans et al. 1999,
   [DOI 10.1074/jbc.274.41.29025](https://doi.org/10.1074/jbc.274.41.29025),
   for muscarinic- and shrinkage-dependent upregulation of NHE1 in mouse
   parotid acini.
8. Nguyen et al. 2000,
   [DOI 10.1111/j.1469-7793.2000.t01-2-00139.x](https://doi.org/10.1111/j.1469-7793.2000.t01-2-00139.x),
   for HCO3-dependent stimulated acidification, carbonic-anhydrase sensitivity,
   and NHE1-dependent pH recovery in mouse sublingual acini. This supports an
   explicit carbon/NHE module while also marking a gland-type transfer.
9. Evans et al. 2000,
   [DOI 10.1074/jbc.M003753200](https://doi.org/10.1074/jbc.M003753200),
   for basolateral acinar NKCC1 and the large salivary defect after NKCC1 loss.
10. Ma et al. 1999,
   [DOI 10.1074/jbc.274.29.20071](https://doi.org/10.1074/jbc.274.29.20071),
   for the requirement for high AQP5-mediated water permeability in saliva
   secretion; AQP5-null mice secreted more than 60% less saliva.
11. Palk et al. 2010,
    [DOI 10.1016/j.jtbi.2010.06.027](https://doi.org/10.1016/j.jtbi.2010.06.027),
    and Vera-Sigüenza et al. 2018,
    [DOI 10.1007/s11538-017-0370-6](https://doi.org/10.1007/s11538-017-0370-6),
    for published model topology and balance conventions, not direct biology.
12. Takano et al. 2021,
    [DOI 10.7554/eLife.66170](https://doi.org/10.7554/eLife.66170), for a later
    model-lineage implementation with KCa channels and Na/K pumps on both
    apical and basolateral membranes. Its distribution sensitivities are model
    results, not new localization measurements.

## What the seven-state chassis currently contains

The Task-12 comparison state is

\[
y=(Na_l,K_l,H,Na_i,K_i,Cl_i,HCO_{3i}),
\]

with intracellular proton concentration eliminated by electroneutrality,
luminal chloride set to `Cl_l=Na_l+K_l`, intracellular CO2 eliminated by a
fitted algebraic expression, and the two membrane potentials solved from a
two-equation quasi-steady current system. The interstitium is an infinite fixed
bath. Lumen volume is not a state. Cell height represents cell volume.

Its non-AE4 topology is:

- basolateral NKCC1, AE2, NHE1, Na/K ATPase, CO2 exchange, water entry, and one
  Ca-activated K conductance;
- apical Ca-activated Cl conductance and water exit;
- paracellular Na, K, and water pathways;
- dynamic luminal Na and K with convective outflow; and
- no apical Na/K pump or apical K conductance.

This is an internally closed historical comparison system. Its absolute flux,
water, and physical-time scales are not certified, and its acid-base reduction
is neither a verified equilibrium chemistry nor a dynamic total-carbon
balance.

## Topology audit

| Module | Seven-state representation | Independent evidence | Audit disposition |
| --- | --- | --- | --- |
| AE2 | Explicit reversible basolateral Cl/HCO3 exchange | Acinus-specific AE2 loss did not detectably change secretion, resting chloride, or initial chloride uptake in the 2015 matched comparison | **Retain explicitly.** A null secretion phenotype is not evidence of absence. Use one shared capacity for WT and AE4-null; set it to zero only in the AE2-null comparison. Require the AE2-null ionic gate as well as the secretion gate. |
| Na/K ATPase | One basolateral pump; its net outward current enters only the basolateral current law | 2018 parotid labelling appeared around the plasma membrane including the apical pole; standard pump orientation implies apical Na delivery and luminal K recovery | **Topology incomplete.** Split into apical and basolateral copies while initially conserving total WT capacity. Do not interpret uniform fluorescence as a 50:50 functional split. |
| Ca-activated K conductance | One basolateral current and one K loss from the cell | Direct apical K currents and apical BK/IK pharmacology in parotid acini; older direct basolateral channel evidence; neither supports an exclusively apical channel population | **Topology incomplete.** Split into apical and basolateral conductances. Add apical K to the lumen balance and apical current closure. |
| NHE1 | Explicit basolateral electroneutral Na-in/H-out flux with a static concentration law | NHE1 is basolateral and dominates pH recovery; WT muscarinic stimulation and shrinkage upregulate it; 2015 found no detectable AE4-genotype difference in stimulated NHE-dependent alkalinization | **Correct protein and polarity, incomplete regulation/chemistry.** Retain one genotype-invariant NHE1 module; estimate its WT stimulus/state response independently. No AE4-null-specific capacity change. |
| Carbon chemistry | Intracellular CO2 is an algebraic fitted function; H is eliminated; the historical buffer rate is scaled and not at chemical equilibrium | HCO3 loss, carbonic-anhydrase sensitivity, and NHE1-dependent recovery establish coupled carbon/pH physiology; the 2015 AE4 study did not directly measure carbonic-anhydrase activity | **Not a validated modern closure.** Use mass-conserving CO2/HCO3/H chemistry or a demonstrably equivalent fast-equilibrium reduction. Add CO3 only if a carbonate AE4 family is tested. |
| NKCC1 | Explicit electroneutral basolateral 1 Na:1 K:2 Cl cotransport | Basolateral localization and large NKCC1-null secretion phenotype; isolated NKCC1 activity was not detectably altered by AE4 loss in 2015 | **Retain and restrict.** Capacity is shared across AE4 genotypes. State-dependent flux may change, but arbitrary knockout upregulation is prohibited. Use the corrected M-to-mM rate convention already established in Task 10. |
| Membrane voltage | Distinct apical and basolateral QSS potentials; currents are apical Cl, basolateral K/pump, and tight-junction Na/K | Polarized current closure is necessary; the apical cation-current evidence changes which currents enter each domain | **Structure useful, current inventory incomplete.** Keep two membrane drops and the tight-junction loop, but insert apical pump and K current before solving the QSS system. Verify signs with charge conservation, not inherited closed forms. |
| Water/cell volume | Apical, basolateral, and paracellular osmotic fluxes; `dV_i/dt=q_b-q_a` | AQP5 and high apical water permeability are directly required for normal saliva; conservation fixes the volume sign | **Retain topology, replace uncertified scaling.** Water flux must use the osmolality of every represented chemical species. Do not fit water gain to AE4-null saliva. |
| Lumen/outflow | Dynamic `Na_l,K_l`; algebraic `Cl_l=Na_l+K_l`; fixed lumen volume; `Q=q_a+q_t`; outflow removes Na/K/Cl | Apical K secretion and apical pump K recovery/Na secretion directly alter lumen balances; carbon-containing saliva invalidates the three-ion electroneutrality reduction when HCO3 is dynamic | **Incomplete after either cation split or carbon extension.** Write complete luminal amount balances and electroneutrality. A fixed-volume control volume is acceptable only with explicit `q_out=q_a+q_t`; otherwise make lumen volume and outflow closure dynamic. |
| Stimulus | One global Ca step on an uncertified time coordinate | Apical K/Cl responses depend on local Ca; AE4 is independently beta/cAMP/PKA regulated; NHE1 responds to muscarinic input and shrinkage | **Protocol incomplete.** At minimum expose separate Ca and beta/PKA inputs and a WT-constrained NHE1 input/state response. Spatial Ca is optional unless the well-mixed reduction fails a source-constrained gate. |

## Modular amount balances and signature vectors

### Sign convention

Use amount states `n_X=V[X]` whenever volume changes. A positive transporter
cycle below has the stated physiological direction. Let

\[
\mathbf b=(Na_i,K_i,Cl_i,HCO_{3i},H_i,CO_{2i};
Na_l,K_l,Cl_l,HCO_{3l};q_a,q_b,q_t).
\]

The first ten entries are mole changes per event. The last three are net
positive-charge equivalents transported from cell to lumen (`q_a`), cell to
bath (`q_b`), and lumen to bath (`q_t`). They are signatures, not currents;
the corresponding conventional current is `F` times event rate times charge
signature. Electroneutral transporters have zero charge entries.

With this convention, the independently supported non-AE4 signatures are:

| Module and positive direction | Signature `sigma` in the coordinates above |
| --- | --- |
| Basolateral NKCC1, bath to cell | `( +1,+1,+2,0,0,0 ; 0,0,0,0 ; 0,0,0 )` |
| Basolateral AE2, Cl into cell/HCO3 to bath | `( 0,0,+1,-1,0,0 ; 0,0,0,0 ; 0,0,0 )` |
| Basolateral NHE1, Na into cell/H out | `( +1,0,0,0,-1,0 ; 0,0,0,0 ; 0,0,0 )` |
| Basolateral Na/K pump | `( -3,+2,0,0,0,0 ; 0,0,0,0 ; 0,+1,0 )` |
| Apical Na/K pump | `( -3,+2,0,0,0,0 ; +3,-2,0,0 ; +1,0,0 )` |
| Basolateral K channel, cell to bath | `( 0,-1,0,0,0,0 ; 0,0,0,0 ; 0,+1,0 )` |
| Apical K channel, cell to lumen | `( 0,-1,0,0,0,0 ; 0,+1,0,0 ; +1,0,0 )` |
| Apical Cl channel, cell to lumen | `( 0,0,-1,0,0,0 ; 0,0,+1,0 ; -1,0,0 )` |
| Effective apical HCO3 loss, cell to lumen | `( 0,0,0,-1,0,0 ; 0,0,0,+1 ; -1,0,0 )` |
| Paracellular Na, bath to lumen | `( 0,0,0,0,0,0 ; +1,0,0,0 ; 0,0,-1 )` |
| Paracellular K, bath to lumen | `( 0,0,0,0,0,0 ; 0,+1,0,0 ; 0,0,-1 )` |
| Carbon hydration `CO2 -> HCO3 + H` | `( 0,0,0,+1,+1,-1 ; 0,0,0,0 ; 0,0,0 )` |

The channel signatures assume positive outward particle flux. Reversal is
represented by a negative event rate, not by a second irreversible module.
The effective HCO3 signature is a boundary balance required when acinar HCO3
loss is represented; the audited evidence does not identify its molecular
carrier, so it must not be relabelled as AE4 or assigned to a specific channel
without a separate source gate.

### Exact reason the membrane split is new

Projection onto the Task-12 intracellular chemical coordinates makes apical
and basolateral pump cycles identical, and likewise makes apical and
basolateral K-channel events identical:

\[
\pi_i\sigma_{P,a}=\pi_i\sigma_{P,b}=(-3,+2,0,0,0,0),
\]

\[
\pi_i\sigma_{K,a}=\pi_i\sigma_{K,b}=(0,-1,0,0,0,0).
\]

The topology difference is visible only in the extended lumen/current space:

\[
\delta\sigma_P:=\sigma_{P,a}-\sigma_{P,b}
=(0_i;+3,-2,0,0;+1,-1,0),
\tag{1}
\]

\[
\delta\sigma_K:=\sigma_{K,a}-\sigma_{K,b}
=(0_i;0,+1,0,0;+1,-1,0).
\tag{2}
\]

These two nonzero directions are linearly independent. They add direct lumen
Na/K loading and redistribute current between the two membrane domains without
adding a new intracellular Na/K stoichiometry. Consequently:

- the Task-12 intracellular source projection could not test this topology;
- redistributing a fixed total pump or fixed total K conductance cannot, by
  itself, repair an intracellular frozen-row source mismatch; but
- it can change luminal composition, tight-junction current, both membrane
  potentials, channel driving forces, osmotic water flux, and therefore the
  relaxed WT/knockout response.

The Task-12 steady cell identity `J_Cl,out=P_total+J_K,total` remains true only
under the same electroneutral cell balances **and no other apical anion loss**,
with

\[
P_{total}=P_a+P_b,
\qquad
J_{K,total}=J_{K,a}+J_{K,b}.
\tag{3}
\]

If an effective apical bicarbonate loss `J_HCO3,a` is admitted, charge and
amount balance instead require

\[
J_{Cl,a}+J_{HCO_3,a}=P_{total}+J_{K,total},
\tag{3a}
\]

for monovalent Cl/HCO3 and no other transcellular current. Thus a carbon/lumen
extension changes the exact anion support identity; it cannot be added only to
osmolality after the chloride projection has been performed.

Equation (3) is location-blind. The new information is the partition of those
totals between membrane and lumen balances. A claimed rescue must therefore
show the full pathway from (1)-(2) through current and luminal closure; a
changed distribution fraction alone is not a mechanism explanation.

### Membrane-split cation equations

Let `P_a,P_b` be apical and basolateral pump-cycle rates, and `J_Ka,J_Kb`
positive K particle effluxes. The cell cation contributions are

\[
\dot n_{Na_i}\big|_{pump}=-3(P_a+P_b),
\qquad
\dot n_{K_i}\big|_{pump,K}=2(P_a+P_b)-J_{K,a}-J_{K,b}.
\tag{4}
\]

The new luminal terms are

\[
\dot n_{Na_l}\big|_{apical\ pump}=+3P_a,
\qquad
\dot n_{K_l}\big|_{apical}=J_{K,a}-2P_a.
\tag{5}
\]

These terms must be added before convective loss. They make the apical pump a
direct transcellular Na source to the lumen and a K-recovery route, precisely
the role proposed in the 2018 parotid study.

A minimal, parameter-economical split conserves the independently calibrated
WT totals:

\[
P_a=f_P P_{tot},\quad P_b=(1-f_P)P_{tot},
\qquad
g_{K,a}=f_K g_{K,tot},\quad g_{K,b}=(1-f_K)g_{K,tot}.
\tag{6}
\]

Equation (6) is a modelling parameterization, not evidence that rates simply
split after gating. A more physical implementation gives each pump its local
external K concentration and each K channel its local voltage and Ca gate,
while splitting only maximum surface capacity. Both `f_P` and `f_K` are
**UNCALIBRATED for mouse submandibular acini**; they are localization
sensitivity parameters until measured in the matched gland and protocol.

### Current closure

To avoid inherited sign ambiguities, define

\[
V_a=\phi_i-\phi_l,\qquad V_b=\phi_i-\phi_e,
\qquad V_l=\phi_l-\phi_e=V_b-V_a.
\]

Let `I_a` be conventional positive current from cell to lumen, `I_b` from cell
to bath, and `I_t` from lumen to bath. With negligible charge accumulation in
the QSS reduction,

\[
I_a-I_t=0,
\qquad
I_b+I_t=0.
\tag{7}
\]

The domain currents must include

\[
I_a=I_{Cl,a}+I_{HCO_3,a}+I_{K,a}+F P_a+I_{other,a},
\]

\[
I_b=I_{K,b}+F P_b+I_{other,b}.
\tag{8}
\]

Here outward Cl particle flux contributes negative conventional current;
outward HCO3 likewise contributes negative current; outward K and each pump
cycle contribute positive current. NKCC1, NHE1, AE2,
and every admitted electroneutral AE4 cycle contribute no direct membrane
current. Equation (7), together with the paracellular current-voltage laws,
should be solved for `V_a,V_b`. The old two-by-two QSS machinery can be reused
only after (8) replaces its incomplete current inventory.

### NHE1 and carbon chemistry

The minimal basolateral NHE1 event is

\[
Na_e+H_i\rightleftharpoons Na_i+H_e,
\qquad J_{NHE1}>0\text{ for Na entry/H exit}.
\tag{9}
\]

It is electroneutral. A genotype-invariant constitutive form may depend on
`Na_i,Na_e,H_i,H_e`, cell volume/shrinkage, and the WT muscarinic input. The
1999/2001 data support regulation by acid load, cell shrinkage, and muscarinic
stimulation. They do not support an AE4-genotype switch. The 2015 native assay
is a gate: the same protocol must not create a large WT-versus-AE4-null NHE
activity difference absent from the data.

For carbon chemistry, use at minimum

\[
CO_2+H_2O\rightleftharpoons H^++HCO_3^-,
\qquad
J_{CA}=k_h[CO_2]_i-k_d[H^+]_i[HCO_3^-]_i,
\tag{10}
\]

and Fickian CO2 entry from the bath and, if retained, exchange with the lumen:

\[
\dot n_{CO_2,i}=J_{CO_2,b}+J_{CO_2,a}-J_{CA},
\]

\[
\dot n_{HCO_3,i}=J_{CA}+\sum J_{HCO_3,transport},
\qquad
\dot n_{H,i}=J_{CA}-J_{NHE1}+\sum J_{H,buffer/transport}.
\tag{11}
\]

A noncarbonate buffer must conserve its total sites, for example
`BH <-> B^- + H+`, or be represented by a derived buffer-capacity reduction.
Fast equilibrium is acceptable only if the reduced equations conserve total
inorganic carbon and total buffer and reproduce (10)-(11) in the appropriate
limit. The historical fitted CO2 denominator does not meet that standard.

If a carbonate-containing AE4 hypothesis is tested, add the distinct reaction

\[
HCO_3^-\rightleftharpoons H^++CO_3^{2-}
\tag{12}
\]

and include `CO3` in charge, osmolality, and every relevant transporter
affinity. Reusing the bicarbonate state for carbonate is inadmissible.

### Lumen, water, and outflow closure

For each luminal solute `X`, write an amount balance

\[
\dot n_{X,l}=J_{X,apical}+J_{X,tight}-q_{out}[X]_l
+J_{X,lumen\ chemistry}.
\tag{13}
\]

If lumen volume is fixed, conservation requires

\[
q_{out}=q_a+q_t.
\tag{14}
\]

If it is not fixed, use

\[
\dot V_l=q_a+q_t-q_{out}
\tag{15}
\]

and give `q_out` an independently stated pressure/duct resistance closure.
Do not use both fixed volume and an independent fitted outflow.

With HCO3/H (and possibly CO3) present, luminal electroneutrality is no longer
`Cl_l=Na_l+K_l`. In the displayed species set it is

\[
[Na]_l+[K]_l+[H]_l
=[Cl]_l+[HCO_3]_l+2[CO_3]_l+\text{other fixed anions}.
\tag{16}
\]

Osmotic water flux must include every represented solute. A generic
conservation-consistent form is

\[
q_a=L_a(\Pi_l-\Pi_i),
\qquad q_b=L_b(\Pi_i-\Pi_e),
\qquad q_t=L_t(\Pi_l-\Pi_e),
\]

\[
\dot V_i=q_b-q_a.
\tag{17}
\]

Pressure terms can be added if independently required. `L_a` represents the
high AQP5-supported apical water permeability; it must be calibrated from WT
water/volume evidence, not from AE4-null secretion. The historical water
coefficients and absolute output scale are not certified physical values.

## Parameter constraints before any knockout secretion reveal

| Quantity | Allowed source constraint | Prohibited interpretation/use |
| --- | --- | --- |
| Apical pump fraction `f_P` | **UNCALIBRATED for submandibular acini.** A strictly nonzero value defines the parotid-informed topology. Qualitative whole-membrane labelling motivates a sensitivity range; local `K_l` must enter apical pump kinetics. Start by conserving total WT pump capacity. | No directly matched submandibular apical-pump capacity fraction was identified in this audit. Do not set `f_P=0.5`: “evenly distributed” does not measure functional capacity or membrane area. Do not fit `f_P` to the AE4-null flow ratio. |
| Apical K fraction `f_K` | **UNCALIBRATED for submandibular acini.** The parotid apical flash increased whole-cell K current by `0.39 nA` (`1.10` to `1.49 nA`; ratio `1.34`), but this perturbation is not a total-conductance partition. Published model values near 20-40% are sensitivity priors only. | The gland-matched submandibular local-Ca experiment did not evoke local apical K current. Do not call 20-40% a measured fraction or force an exclusively apical/basolateral population. |
| Total Na/K-pump capacity | WT Na/K homeostasis, pump kinetics, membrane area, and ATPase evidence. Hold total capacity fixed in the first topology test. | No genotype-specific pump scale chosen to repair AE4-null balances. A later state-dependent change requires independent evidence. |
| Total K conductance and gates | WT currents, Ca dependence, BK/IK pharmacology, and WT membrane potentials. Apical and basal copies use local voltage; separate Ca inputs only if supported. | No conductance chosen from held-out saliva. No inherited historical gate treated as measured biology. |
| NHE1 | WT pH recovery, acid-load, muscarinic, and shrinkage experiments. Same law and parameters in WT and AE4-null. | No AE4-null-specific capacity multiplier; unchanged 2015 NHE readout forbids using one as a rescue. |
| NKCC1 | Basolateral 1:1:2 stoichiometry, corrected rate units, WT chloride uptake, and source-supported stimulation range. Same maximum capacity in WT and AE4-null. | No arbitrary AE4-null upregulation. A changed flux caused by changed concentrations is not the same as changed capacity. |
| AE2 | WT and AE2-null native chloride/pH/secretion comparison; shared capacity in WT and AE4-null. | Do not delete AE2 merely because its knockout flow was negligible. Do not allow a hidden genotype capacity change. |
| Carbon chemistry/buffer | Biochemical equilibrium/rates, WT pH, total inorganic carbon, buffer capacity, and carbonic-anhydrase perturbation where gland-matched. | Do not use the 2015 pH null result as a carbonic-anhydrase measurement; the study explicitly did not measure CA activity. Do not tune buffer capacity to AE4-null saliva. |
| Water/outflow | WT cell volume, osmolarity, AQP5-supported permeability, lumen composition, and conservation. | Do not infer a physical permeability or minute scale from the historical coefficients. Do not fit water gain or outflow resistance to the held-out deficit. |
| Paracellular Na/K | WT transepithelial voltage and lumen composition; apical pump Na delivery may reduce the required paracellular Na flux. | Do not keep the old paracellular current fixed after adding an apical pump; it must be re-solved from voltage and concentrations. |

The historical membrane areas `A_a=36.09 um^2` and `A_b=54.50 um^2` would
give an apical area fraction near `0.398` under equal surface density. This is
historical implementation evidence, not a direct morphometric constraint. It
may be used as a labelled sensitivity initialization, not as the primary
estimate of either `f_P` or `f_K`.

## Minimal ordered extension menu

Each step is switchable and nested. A later step should be attempted only if
the preceding model fails an independently defined gate or if the candidate
AE4 chemistry logically requires it.

### M0. Re-freeze the retained chassis modules

Retain basolateral NKCC1, AE2, NHE1, total Na/K-pump capacity, total
Ca-activated K conductance, apical Cl conductance, the two membrane drops,
paracellular Na/K/water, and cell/lumen water balances. Re-establish WT roots
and AE2-null gates in common physical units. This is not a claim that the old
parameterization is valid; it freezes which existing components may not drift
when a topology switch is tested.

### M1. Add the coupled apical pump/K topology

Implement Eqs. (4)-(8), initially redistributing rather than increasing the
existing total pump and K capacities. The irreducible new parameters are the
two **uncalibrated** localization fractions `f_P,f_K`; local lumen K, local voltages, and local
Ca gates are not optional substitutions. Calibrate only to WT currents,
potentials, lumen Na/K, and source-supported localization/geometry.

This is the first extension because it is the smallest independently supported
topological omission and because it adds the exact extended-space directions
(1)-(2) that were invisible to Task 12. Test apical pump alone, apical K alone,
and the pair as nested controls, but the pair is the biologically coherent
candidate: apical K provides a luminal K source and the apical pump provides a
luminal K recovery route.

### M2. Replace the acid-base reduction and expose NHE1 regulation

Implement Eqs. (9)-(11) with conserved carbon and buffer totals. Fit the shared
WT NHE1 response to independent acid-load, muscarinic, and shrinkage evidence;
then require the 2015 AE4 WT/null NHE readout and resting pH to remain
consistent without genotype retuning. This step supplies acid-equivalent and
Na-loading dynamics, not a free Na source.

### M3. Complete luminal chemistry, osmolality, and outflow

Apply Eqs. (13)-(17) to every represented lumen species. Add luminal HCO3/H
when the carbon module reaches the lumen; add CO3 only for a carbonate AE4
candidate. First test the fixed-volume control-volume closure (14). Promote
`V_l` or pressure-dependent outflow to a state only if fixed-volume closure
fails independent WT volume/composition data.

### M4. Add only source-triggered secondary topology

Two candidates are evidence-supported but not first-line additions here:

- apical NHE2 was detected in mouse parotid acini and NHE2-null parotid glands
  had a late-weighted secretion defect in Park et al. 2001; however, its
  relevance to the 10-minute submandibular AE4 experiment is not established;
- spatially distinct apical/basal Ca signals may be needed to gate the split K
  conductances, but later models show that mean Ca can suffice for total flow
  in some protocols.

Add either only after a gland-matched source or a specific residual demands its
signature. Neither may be selected because it improves held-out secretion.

### M5. Regulatory additions after localization, not before

Only after M1-M3 should one test independently measured pump/K recruitment,
paracellular regulation, AQP5 trafficking, or NKCC1 phosphorylation dynamics.
The first test keeps genotype-independent maximum capacities. An unmeasured
AE4-null capacity change is not a minimal reconstruction.

## Uncertainties and disagreements that must remain visible

1. **Gland mismatch is material.** The AE4 phenotype is submandibular; the
   strongest apical pump/K evidence is parotid. The 2012 parotid study itself
   measured a `0.39 nA` (`34%`) increase in whole-cell K current after local
   apical Ca release (`n=12`, `p<0.001`), whereas Harmer et al. 2005 found that
   localized apical Ca in mouse submandibular acini activated Cl but not K
   until Ca spread more broadly. No directly matched submandibular apical-pump
   capacity fraction was identified in the audited sources. Therefore neither
   parotid result supplies a calibrated submandibular `f_K` or `f_P`; the
   difference is biological heterogeneity, not noise to average away.
2. **Pump labelling is not flux partition.** “Evenly distributed” is a
   qualitative localization result. Apical area, pump surface density,
   substrate access, and functional turnover are all needed to infer `f_P`.
3. **The often-used 40% apical K conductance is model-derived.** Palk/Almássy
   model analyses found efficient secretion with a nonzero partial apical
   fraction and used values near 20-40%; the experiments did not estimate this
   fraction. Later model agreement does not convert it into a measurement.
4. **Basolateral K evidence is not negated by failed local basal photolysis.**
   Almássy et al. 2012 explicitly did not exclude lower-density basolateral
   channels. The defensible topology has both membrane populations.
5. **NHE1 is both important and restricted.** Independent NHE1 loss produces a
   delayed secretion defect and impaired pH recovery, but the AE4-null study
   did not detect altered NHE-mediated alkalinization. This supports shared
   dynamic NHE1 coupling, not genotype-specific retuning.
6. **NHE2 is a live gland-specific omission.** Mouse parotid evidence supports
   apical acinar NHE2, while expression is species/gland dependent. It is not
   silently added to the minimal submandibular chassis.
7. **Carbonic anhydrase remains unmeasured in the AE4 knockout experiment.**
   Normal resting pH does not prove normal carbon flux or buffer capacity.
   Conversely, a carbon module cannot be assigned a knockout change without
   independent evidence.
8. **Unchanged NKCC1 assay activity is a capacity restriction, not a fixed-flux
   rule.** A common NKCC1 law can produce different fluxes at different ionic
   states. What is excluded is arbitrary genotype-dependent activity gain.
9. **Later model lineage is not uniform.** Some descendants included both
   apical and basal cation machinery; another accurate-geometry model avoided
   committing to a K-channel distribution for simplicity. Such choices are
   model reductions, not contradictory experiments.
10. **The old time and water scales remain uncertified.** A topology that
    changes a dimensionless endpoint on the historical coordinate has not yet
    explained minutes of secretion or a physical volume rate.

## Decision handed to the numerical and residual-localization agents

The first nested whole-cell comparison should therefore be:

1. the current total pump and total K conductance, all basolateral;
2. the same totals with pump split only;
3. the same totals with K conductance split only; and
4. the same totals with both split, complete lumen Na/K terms, and the revised
   two-domain current closure.

All four must use identical AE4, NKCC1, AE2, NHE1, carbon, water, and stimulus
parameters. Selection uses WT voltage, lumen Na/K, ionic rest, and AE2-null
evidence only. If the coupled split does not remove the Stage-B ionic residual,
the next test is conserved carbon plus genotype-invariant, WT-constrained NHE1
regulation, followed by full lumen/osmotic closure. This ordering tests the
smallest independently supported topology before adding new total capacities
or regulatory degrees of freedom.
