# Phase 01. Baseline reproduction

## Goal

Reconstruct the published AE4 secretion model with newly written code and verify its reference physiological state.

## Required checks

- steady state concentrations
- membrane potentials
- cell and luminal osmolarity where applicable
- fluid secretion rate
- mass and charge conventions
- solver independence within numerical tolerance

## Acceptance criterion

The reconstructed baseline agrees with the published reference outputs within prespecified numerical tolerances and passes regression tests.
