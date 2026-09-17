# Dimensional case overlay

These dictionaries are derived from the uploaded `constant.zip` and `system.zip`.
They are a starting configuration for validation, not a complete case and not a
calibration to the paper's normalized example. Keep the actual mesh, regions,
reaction/thermo files, initial fields and boundary conditions from the case.

Use a separate small case first. Review and copy the `constant/` files and include
`system/secondConditioningControls` at the end of that case's `system/controlDict`.
Check that the selected start time is less than endTime and that `application` and
the actual launch command both select the new SPFoam executable.

Changes from the supplied settings:

- Unambiguous `phiMod_m`, `CE` and top-level coupling `condVariable z`.
- Dense mixing uses `(xi,x)`; sparse mixing uses the paper's `(phiModified,x)` option.
- Physical pair limits and kernel support start at the supplied r_i = 0.0026676 m.
  This cap is a numerical starting point, not a measured/calibrated locality scale.
- Closest-group retention is 0.8; see the precise rounding/group definition in
  `docs/secondConditioning.md`.
- Neighbouring-rank pairing avoids globally replicating the particle cloud.
- All mechanism species are coupled so the optional enthalpy kernel can reconstruct
  a temperature consistent with its averaged species and absolute enthalpy.
- Coupling relaxation is fixed at 50 microseconds, equal to the original initial
  value. It no longer changes when the CFD timestep changes.
- Chemistry relTol starts at 1e-4. Repeat with 1e-5 to establish convergence.
- Reference statistics are enabled. All original reaction/transport constants and
  CE=20 remain case inputs requiring physical validation.

Restart migration: `restartMode resume` is the default and requires existing second-
conditioning fields. For a legacy cloud, select `restartMode initialize` once. This
initializes phi from current particle temperature, samples a NEW subset and OU
state, and resets reference age; it preserves species/enthalpy. Restore `resume`
for later runs. A changed R does not silently resample an existing subset.

Rebuild all affected libraries before running: the particle stream now carries
`burnedAge`. Do not mix old and new binaries in one parallel run. Old SC field-file
restarts are supported with burnedAgeRate=0; enabling age on an old checkpoint
requires explicit initialization. Random generator states are not checkpointed,
so restart continuity is statistical, not bitwise.

Before large runs, compare local/global/neighbourPairs on a small case; repeat with
several processor decompositions. Neighbour-pair scheduling is a scalable numerical
approximation, not a claim of decomposition-invariant results. The actual supermesh
must also be checked: all three blockMeshDict files in the supplied archives were
identical, and the region meshes themselves were absent.
