#set page(margin: 2cm)
#set text(size: 10pt, font: "New Computer Modern")
#set heading(numbering: "1.1.")
#set math.equation(numbering: "(1)")
#set par(justify: true)
#set bibliography(style: "american-physics-society")

#import "@preview/cetz:0.4.0"

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Process-Tensor Variational Quantum Dynamics (PT-VQD):\
    Compressed Kraus Channels for Scalable Multi-Site\
    Open Quantum System Simulation
  ]

  #v(0.5cm)
  #text(size: 11pt, style: "italic")[Ideas Report -- March 2026]
  #v(1cm)
]

= Research Question

Can process tensor matrix product operators (PT-MPOs) provide a _compressed_ Kraus decomposition of non-Markovian bath effects, enabling variational quantum circuits to simulate multi-site open quantum systems (e.g., the 7-chromophore FMO complex) in regimes inaccessible to both classical process tensor methods and existing path-integral-based quantum algorithms?

= Novelty Claim

Existing quantum algorithms for non-Markovian open systems extract Kraus operators from path integrals, requiring construction of the full $d^2 times d^2$ superoperator @Wang2024Kraus. This scales exponentially with system size, limiting demonstrations to 4-site models. Meanwhile, classical process tensor methods @Cygorek2024ACE @Fux2024 efficiently compress bath memory into PT-MPOs but hit a $D^4$ outer bond bottleneck when the system dimension $D$ grows.

*PT-VQD bridges this gap:* by extracting Kraus operators from compressed PT-MPOs, the operator count is bounded _a priori_ by the PT bond dimension $chi$, independent of system Hilbert space dimension. The quantum circuit handles the multi-site system dynamics on $log_2(D)$ qubits, while the PT-MPO handles the non-Markovian bath classically. Neither component alone scales to the full problem.

= Why Now, Why You

*Why now:*
- PT-MPO methods have matured rapidly: efficient construction algorithms reduce scaling from $O(d^8)$ to $O(d^4)$ @Cochin2026, the inner bond structure is now well-understood @Cygorek2025, and production-quality toolkits exist @Cygorek2024ACE @Fux2024.
- The ACS Omega paper @Wang2024Kraus demonstrated Kraus-based non-Markovian dynamics on NISQ hardware, validating the circuit-side approach on spin-boson and 4-site FMO.
- Divide-and-conquer Stinespring dilation @Azevedo2025 addresses the circuit depth problem.
- The gap between these two advances --- compressed PT-MPO extraction and circuit-based Kraus implementation --- is precisely where PT-VQD sits.

*Why you:*
- Background in variational quantum algorithms (p-VQD @Barison2021), density matrix formalism, and manifold optimization.
- Experience with both tensor network (MPSDynamics.jl, ITensors.jl) and quantum circuit (Yao.jl) ecosystems in Julia.
- Positioned at the intersection of quantum information and chemical physics, where multi-site open systems are the motivating application.

= Cross-field Connections

- *Quantum information $arrow.l.r$ chemical physics:* The FMO complex @Adolphs2006FMO is a canonical benchmark in chemical physics; process tensors originate from quantum information theory @Pollock2018. PT-VQD connects these by using information-theoretic compression (PT-MPO) to enable chemical dynamics simulation on quantum circuits.

- *Tensor networks $arrow.l.r$ quantum circuits:* The hybrid architecture mirrors recent trends in TN-assisted VQE and quantum-classical eigensolvers @Fishman2022, but applied to _dynamics_ rather than ground states.

- *Trapped-ion experiments:* Recent work on donor-acceptor exciton transfer with engineered reservoirs @SimChargeTransfer2025 provides experimental benchmarks for the dimer proof of concept.

= Proposed Approach

== Architecture Overview

#figure(
  cetz.canvas({
    import cetz.draw: *

    // Classical side
    rect((-5.5, -1.5), (-0.5, 1.5), name: "classical", stroke: (paint: blue, thickness: 1.2pt), radius: 0.2, fill: blue.lighten(90%))
    content("classical.north", [*Classical (PT-MPO)*], anchor: "south", padding: 0.15)
    content((-3, 0.5), text(size: 8pt)[Per-bath PT-MPO])
    content((-3, 0), text(size: 8pt)[via TEMPO / iTEBD])
    content((-3, -0.5), text(size: 8pt)[Bond dim: $chi$])

    // Quantum side
    rect((1.5, -1.5), (7.5, 1.5), name: "quantum", stroke: (paint: red, thickness: 1.2pt), radius: 0.2, fill: red.lighten(90%))
    content("quantum.north", [*Quantum Circuit*], anchor: "south", padding: 0.15)
    content((4.5, 0.7), text(size: 8pt)[7 system qubits (FMO)])
    content((4.5, 0.2), text(size: 8pt)[$+ ceil(log_2 chi)$ ancilla qubits])
    content((4.5, -0.3), text(size: 8pt)[Variational ansatz $V(theta)$])
    content((4.5, -0.8), text(size: 8pt)[Stinespring dilation])

    // Arrow
    line((-0.5, 0), (1.5, 0), mark: (end: "straight"), stroke: (thickness: 1.5pt))
    content((0.5, 0.35), text(size: 8pt)[Kraus ops], fill: white, frame: "rect", stroke: none, padding: 0.05)
    content((0.5, -0.35), text(size: 8pt)[rank $lt.eq chi$], fill: white, frame: "rect", stroke: none, padding: 0.05)
  }),
  caption: [PT-VQD hybrid architecture. The classical side computes per-bath process tensors as compressed MPOs. Kraus operators extracted from the PT-MPO (rank bounded by $chi$) are implemented on the quantum circuit via Stinespring dilation.],
) <fig:architecture>

== Method

The PT-VQD algorithm proceeds at each time step $t_n arrow t_(n+1)$:

+ *System evolution:* Apply a variational circuit $V(bold(theta)_n)$ encoding the system Hamiltonian (inter-chromophore couplings).

+ *Bath channels (sequential):* For each bath $j = 1, dots, N_"sites"$:
  - Extract the CPTP map $Phi_j$ from bath $j$'s PT-MPO at time step $n$.
  - Decompose via Choi matrix eigendecomposition: $Phi_j(rho) = sum_(k=1)^(r_j) K_j^((k)) rho K_j^((k) dagger)$, with $r_j lt.eq chi_j$.
  - Implement via Stinespring dilation: apply a controlled unitary on ${"site"_j, "ancilla"}$, measure and reset ancilla.

+ *Parameter update:* Optimize $bold(theta)_(n+1)$ via projected variational dynamics @Barison2021 or direct fidelity maximization.

The critical insight: because each bath couples _locally_ to one chromophore, the Kraus operators act as $K_j^((k)) = bb(I)_1 times.circle dots.c times.circle tilde(K)_j^((k)) times.circle dots.c times.circle bb(I)_(N_"sites")$, where $tilde(K)_j^((k))$ is a $2 times 2$ operator on site $j$. This means:
- Baths are applied *sequentially*, not as a tensor product.
- Ancilla qubits are *reused* across baths (with mid-circuit reset).
- Total qubit count: $N_"sites" + ceil(log_2 chi) approx 11$ for 7-site FMO with $chi = 16$.

== Scaling Comparison

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, center),
    stroke: 0.5pt,
    inset: 6pt,
    table.header([], [*Path integral*\ *Kraus* @Wang2024Kraus], [*Classical*\ *PT-MPO* @Cygorek2024ACE], [*PT-VQD*\ *(this work)*]),
    [Bath handling], [$"TNPI"$\ (exp. in memory)], [$"PT-MPO"$\ (linear in $t$)], [$"PT-MPO"$\ (linear in $t$)],
    [System repr.], [$d^2 times d^2$\ superoperator], [$D^4$\ outer bonds], [$log_2(D)$\ qubits],
    [Kraus rank], [Post-hoc\ truncation], [N/A], [Bounded by $chi$\ _a priori_],
    [4-site FMO], [Demonstrated], [Demonstrated], [Target],
    [7-site FMO], [$times$ (Choi too large)], [Borderline\ ($D^4 approx 10^8$)], [*Target*\ ($approx 11$ qubits)],
  ),
  caption: [Scaling comparison of three approaches to non-Markovian multi-site dynamics.],
) <tab:scaling>

== Benchmark Ladder

+ *Spin-boson* (Ohmic, $xi = 0.1$, $omega_c = 7.5$, $beta = 5$): Validate against exact TEMPO and T-TEDOPA. Direct comparison with @Wang2024Kraus parameters.

+ *Donor-acceptor dimer* (2 sites, independent Drude baths, $lambda = 35 "cm"^(-1)$, $gamma = 106.18 "cm"^(-1)$, $T = 300 "K"$): Proof of concept for multi-site PT-VQD. Compare against trapped-ion experimental results @SimChargeTransfer2025.

+ *FMO complex* (7 chromophores, FMO Hamiltonian from @Adolphs2006FMO): The main result. Demonstrate scaling beyond 4-site limit of path-integral Kraus methods.

= Minimum Viable Experiment

*Dimer proof of concept (Step 4 of implementation plan):*

- Build PT-MPO for two independent Drude baths using PT-TEMPO in Julia (ITensors.jl).
- Extract Kraus operators at each time step, verify rank $lt.eq chi$.
- Implement Stinespring dilation in Yao.jl for a 2-qubit system + ancilla.
- Simulate population transfer dynamics $P_(1 arrow 2)(t)$.
- Compare against QuantumDynamics.jl (TEMPO) and MPSDynamics.jl (T-TEDOPA).

*Success criterion:* PT-VQD reproduces exact dynamics with fewer Kraus operators than the path-integral approach, and the Kraus rank is demonstrably bounded by $chi$.

= Success Signal

The problem is solved if:
- PT-VQD simulates 7-site FMO dynamics with accuracy comparable to classical PT methods ($< 1%$ trace distance error).
- The Kraus rank remains bounded by $chi$ (not growing with system size) as confirmed empirically.
- Circuit resource requirements ($approx 11$ qubits, polynomial depth) are achieved.
- A clear computational advantage regime is identified: parameter space where classical PT fails ($D^4$ bottleneck) but PT-VQD succeeds.

= Hope Signal

The approach has hope even if the full vision isn't immediately realized, if:
- The PT-MPO Kraus extraction works correctly and the rank is indeed bounded by $chi$ on the dimer.
- The sequential bath application with ancilla reuse is validated (even if FMO-scale isn't reached yet).
- The Julia PT-TEMPO implementation reproduces OQuPy / ACE results, becoming a useful tool in itself.

= Pivot Signal

Consider abandoning or pivoting if:
- The Kraus rank extracted from PT-MPO is _not_ bounded by $chi$ in practice (e.g., it grows with system size despite compression). This would undermine the core scaling claim.
- The Stinespring dilation circuit depth is too large even for the dimer, making the quantum circuit impractical.
- Classical PT methods (e.g., the $O(d^4)$ algorithm @Cochin2026) continue improving and close the gap before PT-VQD can demonstrate advantage.
- The Trotter splitting between system evolution and bath channels introduces unacceptable errors at feasible time step sizes.

= Open Risks

+ *Mid-circuit measurement vs. deterministic dilation:* Ancilla reuse across baths requires mid-circuit measurement and reset. If deterministic (no-measurement) Stinespring is required, qubit count increases to $N_"sites" times ceil(log_2 chi) approx 35$ for FMO. Still feasible, but less elegant.

+ *Trotter error:* The splitting of system Hamiltonian evolution and bath channel application introduces Trotter error controlled by $Delta t$ and $[H_"sys", H_"sys-bath"]$. This matches the classical PT framework's splitting, but convergence needs empirical verification.

+ *Inter-bath correlations:* The sequential bath application assumes independent baths. If bath-bath correlations are significant (e.g., through common phonon modes), the architecture would need modification.

+ *Classical competition:* The field is moving fast. The March 2026 efficient PT construction @Cochin2026 narrows the quantum advantage window. The target regime (strong coupling, multi-site, non-Markovian) must be carefully identified.

= Target Venue

- *Primary:* Physical Review Research or Quantum (open access, hybrid methods welcome)
- *Alternative:* The Journal of Chemical Physics (computational chemistry audience)
- *Software note:* Julia package could be submitted to Journal of Open Source Software (JOSS)

= Implementation Plan

== Pure Julia Tool Stack

#figure(
  cetz.canvas({
    import cetz.draw: *

    // Boxes
    rect((-5, -0.6), (-2, 0.6), name: "itensors", stroke: 1pt, radius: 0.15, fill: green.lighten(85%))
    content("itensors.center", text(size: 8pt, weight: "bold")[ITensors.jl\ + PT-TEMPO])

    rect((-0.5, -0.6), (2.5, 0.6), name: "linalg", stroke: 1pt, radius: 0.15, fill: yellow.lighten(80%))
    content("linalg.center", text(size: 8pt, weight: "bold")[LinearAlgebra\ Choi $arrow$ Kraus])

    rect((4, -0.6), (7, 0.6), name: "yao", stroke: 1pt, radius: 0.15, fill: red.lighten(85%))
    content("yao.center", text(size: 8pt, weight: "bold")[Yao.jl\ Stinespring])

    // Validation box below
    rect((-5, -2.3), (-2, -1.3), name: "mpsd", stroke: (dash: "dashed"), radius: 0.15, fill: gray.lighten(90%))
    content("mpsd.center", text(size: 8pt)[MPSDynamics.jl\ (validation)])

    rect((-0.5, -2.3), (2.5, -1.3), name: "qd", stroke: (dash: "dashed"), radius: 0.15, fill: gray.lighten(90%))
    content("qd.center", text(size: 8pt)[QuantumDynamics.jl\ (validation)])

    // Arrows
    line("itensors.east", "linalg.west", mark: (end: "straight"), stroke: 0.8pt)
    line("linalg.east", "yao.west", mark: (end: "straight"), stroke: 0.8pt)

    // Labels
    content((-1.25, 0.9), text(size: 7pt)[PT-MPO], fill: white, frame: "rect", stroke: none, padding: 0.03)
    content((3.25, 0.9), text(size: 7pt)[Kraus ops], fill: white, frame: "rect", stroke: none, padding: 0.03)
  }),
  caption: [Pure Julia tool stack for PT-VQD. Solid boxes: core pipeline. Dashed boxes: independent classical validation.],
) <fig:toolstack>

== Timeline

#table(
  columns: (auto, 1fr, auto),
  align: (center, left, left),
  stroke: 0.5pt,
  inset: 6pt,
  table.header([*Weeks*], [*Task*], [*Deliverable*]),
  [1--2], [Classical baselines: spin-boson with QuantumDynamics.jl (TEMPO) and MPSDynamics.jl (T-TEDOPA). Cross-validate. Yao.jl warm-up.], [Validated $angle.l sigma_z (t) angle.r$ dynamics],
  [3--4], [Implement PT-TEMPO in Julia using ITensors.jl. Extend QuantumDynamics.jl TEMPO to output PT-MPO.], [Validated PT-MPO],
  [5--6], [Kraus extraction pipeline: Choi decomposition of PT-MPO slices. Verify rank $lt.eq chi$.], [Kraus extraction module],
  [7--9], [Quantum circuit: Stinespring in Yao.jl. Dimer proof of concept.], [Dimer PT-VQD results],
  [10--13], [Scale to 4-site and 7-site FMO. Scaling analysis.], [FMO results],
  [14--16], [Paper writing and software packaging.], [Manuscript + Julia package],
)

= Key References

#bibliography("references.bib")
