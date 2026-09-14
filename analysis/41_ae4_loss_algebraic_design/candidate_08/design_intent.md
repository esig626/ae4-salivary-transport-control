Candidate 7 completed all three trajectories and physiological/conservation
checks. Its cumulative deficits were 26.1848% (5%) and 35.6915% (null).
The null is 0.6915 percentage points above the specified range. Its ledger
and parameter remain frozen. Make one algebraic finite-horizon correction
to the SAME CaCC recruitment hypothesis; introduce no new mechanism.

Use the two already computed null data points: Task 40 g=1 and candidate 7
g=b7. Approximate cumulative secretion with a saturating conductance map
Q(g)=Qinf*g/(g+k), analogous to a series conductance. This is a two-point
design surrogate, not an exact identity of the whole-cell ODE. It assumes
zero CaCC-linked secretion at g=0; that untested limit is not asserted as a
biological or model result. For Q1=Q(1) and Qb=Q(b7),

  k=b7*(Q1-Qb)/(Qb-b7*Q1),  Qinf=Q1*(1+k).
  b8=k*Qtarget/(Qinf-Qtarget), Qtarget=0.70*Q_WT_Task40.

Calculate b8 once before rerunning, freeze it, and run WT/5%/null from the
same Task 40 WT rest. There is no exact-fit requirement: both final losses
within 20-35% suffice. These trajectories are target-informed tests, not
held-out validation. No sweep or trajectory optimizer is used.
