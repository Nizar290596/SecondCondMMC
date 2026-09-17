# Second conditioning: implemented model and validation status

Reference: *A PDF approach to thin premixed flamelets using multiple mapping
conditioning*, DOI 10.1016/j.proci.2016.07.116, section 5 and Appendix II.

The implementation follows the two sequential levels in that paper. It does not
claim to reproduce its figures or to supply closures that the paper leaves open.
The source changes and standalone tests have been exercised; a complete OpenFOAM
v2406 build, parallel case and flame comparison are still required.

## Evolution and conditioning

1. All particles carry phi, evolved by dense `MMCcurl` mixing and
   W(phi)=A_phi*(1-phi)*exp(Z_phi*(phi-1)). Detailed species are not mixed here
   when second conditioning is enabled.
2. Selected particles carry the OU-perturbed reference phiModified and undergo
   sparse `secondCondMMCcurl` mixing followed by detailed chemistry. Selection is
   made at injection/explicit initialization and preserved by clones and transfer.
3. First conditioning localizes in `(xi,x)`. Second conditioning can use
   `(phiModified,x)` with `includeShadowPositions false`, or
   `(phiModified,xi,x)` with that switch true. Shadow entries are resolved by name.
4. Species, absolute enthalpy and transported coupling scalars mix on the sparse
   level. Temperature is reconstructed from the mixed species and enthalpy.
   Neither mixing stage changes shadow position or the selection flag.

`physicalLocalization true` activates physical coordinates in the normalized
k-d-tree metric. `r_i` and `Xim_i` set physical and reference normalizations. A
small normalization means stronger localization. No default can guarantee the
correct hierarchy for every flame; inspect actual pair distances and gradients.
The dimensional overlay chooses the paper's phi/x option to remove unintended
tight shadow conditioning at the sparse level.

`maxPairDistance` is an optional hard physical-distance cap (0 disables it).
`retainedPairFraction` retains the closest candidate groups, ranked by maximum
normalized squared distance within each group. Its default is 0.8 with second
conditioning enabled and 1 otherwise. Whole groups of two or three are retained
within ceil(fraction*N) particles, so rounding can lower the actual fraction. The
paper states closest-80% selection but does not prescribe tie/group rounding;
this is the explicit implementation choice here. It is not a universal calibrated
percentage. Rejection is not compensated by arbitrarily accelerating surviving
pairs; measured variance decay includes the rejected fraction.

## Mixing and reaction times

`mixingTimeScale aISO` uses tau_i=DeltaE_i^2/[CE*vb_i*(D_i+Dt_i)]. With
`meanTimeScale true`, pairs use the harmonic mean; otherwise they use the minimum.
The dense model also retains its previous gradient model as
`mixingTimeScale gradient`. Second conditioning supports aISO and prescribed only.

`mixingTimeScale prescribed; tauMix <seconds>;` gives a constant timescale for
isolated validation. Each level has independent controls. CE=20 in the uploaded
case is preserved in the overlay; it has not been calibrated against the paper.
`C_E` is accepted as a deprecated alias, with conflicts rejected.

The timescale convention is pair-contrast relaxation: an isolated exponentially
mixing pair has variance V(t)=V(0)*exp(-2*t/tauMix). Dense-reference evolution uses
`mixingExtentModel exponential` to avoid finite Poisson jumps as dt tends to zero.
The paper expressly allows jump-free reference evolution.

The sparse default is `mixingExtentModel modifiedCurl`: uniform-extent mixing
events with Poisson rate 3/tauMix. Since E[(1-alpha)^2]=1/3 per event, this has the
same expected variance-decay convention for any finite timestep. All ranks holding
a pair derive identical random events from communicated pair data. The paper does
not specify the finite-step event implementation; this is our explicit choice.

Triplets perform all three edges for dt/2 each. Each particle receives total
exposure dt; all replicas execute the same three exchanges, including remote-
remote exchanges required before a later local exchange. Weighted pair means are
conserved. Triplet variance decay is not identical to that of an isolated pair,
so ensemble calibration still matters.

Progress reaction uses bounded exponential-midpoint substeps with A_phi*h <=
`maxPhiSourceStep` (default 0.05). No post-update clipping conceals overshoot.
Halving this bound provides an integration convergence check. The cheap kinetics
must still be calibrated to the detailed mechanism at the actual inlet state.
The source peaks at phi=1-1/Z_phi (0.95 when Z_phi=20).

The solver reads optional `transportProperties/premixedFlame` values
`laminarFlameSpeed` and `laminarFlameThickness`. Previous hard-coded values remain
the defaults. This makes the assumptions visible; it does not calibrate them.
The original LES flame-speed/vb closure is retained. The paper permits alternative
flame-speed approximations; this is not its normalized 4/7-law benchmark.

## Reference diffusion and optional age

The OU update is exact with unit stationary variance and stable expm1 evaluation.
Its stream is separate from transport and particle selection uses the SC stream.
Shadow-coordinate streams have independent component seeds; clone constructors
preserve state. `randomSeed` controls initial streams. Generator states are not
serialized for bitwise restart equivalence.

The reference is phiModified=(phi+burnedAgeRate*burnedAge)*exp(beta*omega).
With the default burnedAgeRate=0 this is the paper's phi*exp(beta*omega), including
values greater than one. Standard OU is Gaussian with unbounded support; no
arbitrary clipping of omega or phiModified is applied.

Choose either `beta` or `referenceDiffusionRate`, never both. The latter specifies
beta^2/tauOU, so beta=sqrt(referenceDiffusionRate*tauOU). It specifies the coefficient
of multiplicative reference diffusion, not an automatic scalar-dissipation closure.
The mean OU diffusion is beta^2/tauOU times the weighted mean of phiModified^2.
Logs report this and a dense mixing Nphi estimate from weighted variance loss:
Nphi=varianceLoss/(2*dt*totalWeight). These global diagnostics support calibration;
conditional dissipation and consistent mass/PDF weighting must also be assessed.

An optional age extension accumulates
 d(burnedAge)/dt=max(0,(phi-burnedAgeThreshold)/(1-burnedAgeThreshold)).
Age is mixed with the dense progress reference and is initialized to zero at
injection or explicit migration. burnedAgeRate has units 1/s; the threshold is in
[0,1). This is an additional smooth residence-time closure, NOT an equation given
in the paper. Keep it zero for the unextended model and test CO sensitivity before
using it. OU fluctuations above one alone do not represent burned-gas age.

`maxOUTimeStepRatio` limits the solver timestep to that fraction of tauOU when
beta>0 (default 0.1). Exact OU integration does not remove the need to resolve the
changing pairing coordinate. Logs also show the maximum exchange dt/tauMix.

## Parallel pairing and particle management

`pairingMethod local` is a processor-local approximation. `global` is suitable for
small verification cases; it replicates the candidate data and can be prohibitive
on a large cloud. `neighbourPairs` forms disjoint neighbouring processor groups
using a rotating greedy matching of the processor-patch adjacency graph. Each edge
gets priority periodically; a rank exchanges candidates with at most one partner
per call. A common sorted rank order and deterministic ties make pair selection
and ghost updates consistent within each group. This is not a proof of statistical
invariance under mesh decomposition, and periodic coordinate transformations are
not added by this implementation.

The empty/single-particle cases are valid and do not skip MPI collectives. Sparse
exchange reserves based on subset size, not the full cloud. Diagnostics count
unique owned subset particles and report their weight fraction, accepted population
and maximum physical/reference pair separations.

Number control preserves selection flags and reference state through cloning.
Synchronized resampling now uses a common rank/particle order for equal-weight
samples, plus a persistent seeded generator instead of wall-clock reseeding.
Roulette deletion conserves weights in expectation, not exactly in each realization.
R is an injection probability, not an enforced local population fraction. It must
not be forced by promoting particles whose detailed chemistry was not evolved.

## Eulerian feedback

Mixing conditioning and Eulerian coupling conditioning are separate decisions.
The authoritative key is `thermophysicalCoupling/condVariable`. A nested key in
KernelEstimationCoeffs produces a warning identifying the effective setting.
Chemistry always uses the registered mixture fraction z, consistent with the
balanced reaction path.

The kernel measures physical distances over all candidates, enforces rMax,
requires `minKernelParticles` positive-weight samples and
`minEffectiveSamples=(sum w)^2/sum(w^2)`, and logs coverage/support. Failed cells
retain Indicator=0; their transport solution evolves without a fabricated particle
target. Persistent missing coverage is a failed validation, not a solved flame.
The kernel still searches its available local/neighbouring-rank candidates;
coverage must be checked under decomposition changes.

`coupleEnthalpy true` averages absolute enthalpy and all species before reconstructing
temperature. It requires every mechanism species selected and dfMax=0. This preserves
the kernel's mean thermochemical state; it does not make relaxation coupling globally
conservative by itself. The legacy temperature-average mode remains available.

Kernel coupling on phiModified is experimental: a conditional mean evaluated at a
projected mean is generally not an unconditional mean. It requires explicit
`allowExperimentalProgressCoupling true`. The dimensional overlay uses z. A full
conditional-PDF reconstruction is a separate closure, not silently inferred from
the paper.

Use `tauUnits time` for timestep-convergence studies. The supplied 50*timestep
was 50 microseconds initially but varied with adaptive timestep. The overlay uses
50 microseconds as a starting physical value; it is not a calibrated coupling time.

## Configuration, restart, and validation

`mmcVariablesDefinitions` is now read from case constant/, removing the absolute
solver-source path. `phiModified_m` is accepted as an alias for `phiMod_m`; conflicting
aliases fail. Unused nPairSamples/particleFilter/localnessLimited/fullSort controls
are reported. Effective mixing, coupling, and reference settings are printed.

`restartMode resume` preserves valid SC fields. `restartMode initialize` explicitly
initializes missing or existing SC fields from current particle temperature, including
a newly sampled subset. Restore resume afterwards. With SC disabled, legacy clouds
need no SC fields. A new burnedAge field is optional on old SC checkpoints when age
is disabled; it is mandatory to resume an active age model. New stream layout means
all affected libraries/executables must be rebuilt together. Initialization does not
restore omitted detailed-chemistry history and does not constitute calibration.

Run the dependency-free regression checks:

    tests/secondConditioning/run

They compile the production numerical kernels and extract actual pairing/triplet
methods. Checks cover bounded reaction and convergence, exact OU moments/correlation,
Modified Curl variance decay and conservation, N=0..150 pairing, metric dominance,
rejection, disjoint processor scheduling and consistency of triplet replicas, plus actual cloud reference-update paths
(including disabled-noise/source and optional age). They
are NOT a full OpenFOAM template build or an MPI test.

On a configured OpenFOAM v2406 system with existing third-party dependencies:

    tests/secondConditioning/build_openfoam

Then run these gates before production:

1. Fresh/legacy/SC restarts, SC off/on, R=1, balanced and unbalanced chemistry;
   include ranks with 0/1 selected particle and processor-spanning supercells.
2. Homogeneous mixing with chemistry off: weighted species/enthalpy conservation,
   realized variance decay, positivity, and normalization; several timesteps/seeds.
3. Paper benchmark: constant density, unity Le/Sc, normalized D=nu=0.001,
   A=4000, Z=20, dense N=10000 and reacting N=2500; same one-step kinetics for phi
   and Y2. Specify the paper's transport/reference equations and boundary conditions.
   Choose/document tauMix, reference weights, beta and tauOU; their complete numerical
   values are not supplied in the paper. Compare flame location, thickness and PDFs.
   The dimensional overlay is NOT this benchmark; no figure reproduction is claimed.
4. Dimensional LES: calibrate cheap kinetics and flame parameters to the supplied
   multistep mechanism at known pressure, inlet composition and temperature.
5. Compare corrected ordinary SP, SC R=1, then R=0.5/0.25/0.1. R=1 isolates the
   effect of changing conditioning and need not equal ordinary SP. Refine timestep
   with fixed physical coupling time; vary particle number, seed, locality/support,
   chemistry tolerance and MPI decomposition. Evaluate CO separately and profile
   chemistry/mixing/communication, peak memory and total walltime.

The existing particle class still allocates detailed state on dense particles;
subset chemistry reduces chemistry work, not all dense storage. A separate compact
reference-only particle representation is a further performance change, not necessary
for the conditioning algorithm and not included here.
