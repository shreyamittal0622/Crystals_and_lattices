# Crystals, Bands and Lattices: Tight-Binding, SSH, Kronig-Penney, DOS and Quantum Confinement

This package covers two notebooks:

| Notebook | Topic |
|---|---|
| `notebooks/02_crystals_and_bands.ipynb` | 1D tight-binding band and effective mass, finite-difference (FD) kinetic energy, particle in a box, SSH chain, Kronig-Penney model, FD with a periodic potential |
| `notebooks/02_crystals_and_lattices_nb2.ipynb` | Square and honeycomb lattice bands, density of states (DOS) in 1D/2D/3D, effective mass versus band gap, quantum confinement (2D/1D/0D DOS), effective-mass curvature fit |

This README explains **every equation used in both notebooks**, shows the notebook figures next to the equation that produced them, tabulates the results, and lists the problems found. Several of them change how the saved figures should be read, so please read Section 7 before reusing any plot.

> Cell numbers below count from 1, top to bottom. Energies in Notebook 1's Kronig-Penney and FD sections use ħ = m = 1, so their "eV" labels are nominal.

---

## 1. Package contents

```
crystals_bands_package/
├── README.md
├── requirements.txt                      numpy, scipy, matplotlib, jupyter
├── notebooks/
│   ├── 02_crystals_and_bands.ipynb
│   └── 02_crystals_and_lattices_nb2.ipynb
└── results/
    ├── notebook1_crystals_and_bands/
    │   ├── figures/    8 plots extracted from the notebook outputs
    │   ├── data/       tb_dispersion, fd_dispersion, box_levels, ssh_bands,
    │   │               ssh_gap_vs_ratio, kronig_penney_band_edges,
    │   │               fd_cosine_potential_bands  (CSV)
    │   └── metrics/    summary_metrics.json, notebook_printed_outputs.txt
    └── notebook2_crystals_and_lattices/
        ├── figures/    6 plots extracted from the notebook outputs
        ├── data/       dos_as_in_notebook, dos_{1D,2D,3D}_band_edge_check,
        │               mass_vs_gap_III-V, confinement_dos, curvature_fit_data  (CSV)
        └── metrics/    summary_metrics.json, notebook_printed_outputs.txt
```

The figures are the notebooks' saved outputs. The CSVs and `summary_metrics.json` were recomputed with clean scripts that follow the notebooks' parameters (the notebooks themselves don't save any data), including the checks that back up the claims in this README.

---

# Part A: Notebook 1, `02_crystals_and_bands.ipynb`

## A1. 1D tight-binding chain and effective mass (Cells 2-5)

### Eq. 1: tight-binding band

$$E(k)=\varepsilon-2t\cos(ka)$$

A chain of atoms with spacing *a*, one orbital per site at energy ε and nearest-neighbour hopping *t*. With plane-wave (Bloch) states ψₙ = e^(ikna), the 1×1 Bloch Hamiltonian is just this number. The band spans ε − 2t to ε + 2t (width 4t) and repeats every 2π/a, so the First Brillouin Zone (FBZ) is −π/a ≤ k ≤ π/a. The notebook sets ε = 0, t = 1, a = 0.5 and uses 500 k-points. (Cell 2 evaluates `E_k = -2*t*cos(k*a)` and does not add ε, which is harmless since ε = 0.)

### Eq. 2: effective mass

$$\frac{1}{m^*}=\frac{1}{\hbar^2}\,\frac{d^2E}{dk^2}$$

The effective mass is the band curvature: a tightly curved band means a light carrier.

### Eq. 3: theoretical effective mass at the band bottom

$$E(k)\approx(\varepsilon-2t)+t a^2k^2\ \Rightarrow\ m^*_{\text{theory}}=\frac{\hbar^2}{2ta^2}=\frac{1}{2(1)(0.25)}=2.0$$

Expanding cos(ka) ≈ 1 − (ka)²/2 gives a parabola E ≈ E₀ + t a²k², which has the free-electron form ħ²k²/2m* with m* = ħ²/(2ta²).

### Eq. 4: numerical extraction

$$E(k)\approx c_2k^2+c_1k+c_0\quad\text{(polynomial fit over }|k|\le0.10\,\pi/a),\qquad m^*_{\text{num}}=\frac{\hbar^2}{2c_2}$$

![Tight-binding fit](results/notebook1_crystals_and_bands/figures/01_tight_binding_effective_mass_fit.png)

| Quantity | Value |
|---|---|
| Fitted c₂ | 0.248238 (ideal t a² = 0.25) |
| m\* numerical | 2.014193 |
| m\* theory | 2.000000 |
| Relative error | **0.7096%** |

**Why the error isn't zero:** the fit window is finite, so the −t a⁴k⁴/12 term in the cosine expansion leaks into c₂. Shrinking the window removes the error, which I checked:

| Fit window (fraction of π/a) | m\* | Error |
|---|---|---|
| 0.10 (notebook) | 2.01419 | 0.710% |
| 0.05 | 2.00324 | 0.162% |
| 0.02 | 2.00054 | 0.027% |
| 0.01 | 2.00007 | 0.003% |

---

## A2. Finite-difference kinetic energy (Cell 6)

### Eq. 5: discretized second derivative

$$-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2}\ \to\ -\frac{\hbar^2}{2ma^2}\left(\psi_{j+1}-2\psi_j+\psi_{j-1}\right),\qquad t_{\text{eff}}=\frac{\hbar^2}{2ma^2}$$

Replacing the continuous kinetic-energy operator with its three-point stencil gives a tridiagonal matrix: **2t_eff on the diagonal and −t_eff on both off-diagonals**. It is the same structure as the tight-binding chain, with hopping t_eff. For m = 1, a = 0.5, ħ = 1: t_eff = 2.

### Eq. 6: Bloch dispersion of the FD Hamiltonian

$$E_{\text{FD}}(k)=2t_{\text{eff}}\bigl(1-\cos ka\bigr)\ \ (\text{the notebook then subtracts }2t_{\text{eff}}),\qquad E\approx\frac{\hbar^2k^2}{2m}\Bigl(1-\frac{(ka)^2}{12}\Bigr)$$

Inserting plane waves into the stencil gives this closed form. The second form shows that the discretization error grows with k²a²/12, so the FD dispersion is accurate only for ka ≪ 1.

### Eq. 7: mass extraction

$$E=c_2k^2+c_0\ \ (\text{fit vs }k^2,\ |k|<0.2),\qquad m^*=\frac{\hbar^2}{2c_2}$$

![FD fit](results/notebook1_crystals_and_bands/figures/02_finite_difference_dispersion_fit.png)

**Result:** input m = 1.0000, extracted m\* = **1.0008** (0.08% high, from the −(ka)²/12 term).

> Note: the cell builds the real-space 100×100 matrix `H`, but never diagonalizes it. E(k) is evaluated from the analytic formula (Eq. 6) instead.

---

## A3. Particle in a box with open boundaries (Cell 7)

### Eq. 8: eigenvalues of the tridiagonal matrix

$$E_n=2t_{\text{eff}}\left[1-\cos\!\frac{n\pi}{N+1}\right]\ \approx\ \frac{\hbar^2}{2m^*}\left(\frac{n\pi}{L}\right)^2,\qquad L=(N+1)a$$

A tridiagonal matrix with zero (Dirichlet) boundaries, with no periodic corner terms, has exactly these eigenvalues. For small n they reduce to the infinite-square-well levels E_n ∝ n².

### Eq. 9: parameters

$$t_{\text{eff}}=\frac{\hbar^2/2m_0}{m^*/m_0\cdot a^2}=\frac{3.810}{0.5\times1.0^2}=7.62\ \text{eV},\qquad N=500,\ a=1\ \text{Å},\ L=501\ \text{Å}$$

ħ²/2m₀ = 3.810 eV·Å² is the standard conversion constant, with m\* = 0.5 m₀.

![Box levels](results/notebook1_crystals_and_bands/figures/03_particle_in_box_levels.png)

| n | E_n (eV) | E_n / E₁ |
|---|---|---|
| 1 | 0.000300 | 1.000 |
| 2 | 0.001198 | 4.000 |
| 3 | 0.002697 | 9.000 |
| 4 | 0.004794 | 15.999 |
| 5 | 0.007490 | 24.998 |
| 10 | 0.029953 | 99.968 |

The numerical eigenvalues match Eq. 8 to 8×10⁻¹⁵ eV, and differ from the continuum n² law by at most 0.033% (n = 10). The printed levels (0.0003, 0.0012, 0.0027, … 0.0300 eV) agree.

> Plot note: the figure adds the *dimensionless* normalized wavefunction (amplitude ≈ ±0.063) to energies of only 0.0003-0.005 eV, so the curves are dominated by the wavefunction amplitude and the vertical axis is not a true energy-level diagram.

---

## A4. SSH (Su-Schrieffer-Heeger) chain (Cells 8-10)

### Eq. 10: Bloch Hamiltonian

$$H(k)=\begin{pmatrix}0&t_1+t_2e^{-ika}\\ t_1+t_2e^{+ika}&0\end{pmatrix}$$

A chain with **two sites per unit cell** and alternating hoppings: t₁ within the cell and t₂ between cells. There is no on-site energy, so the diagonal is zero.

### Eq. 11: energy bands

$$E_\pm(k)=\pm\bigl|t_1+t_2e^{-ika}\bigr|=\pm\sqrt{t_1^2+t_2^2+2t_1t_2\cos ka}$$

A symmetric pair of valence (−) and conduction (+) bands. Maximum |E| = t₁ + t₂ at k = 0 (so 1.4, 2.0 and 2.6 for t₂ = 0.4, 1.0, 1.6, as in the plot); minimum |E| = |t₁ − t₂| at k = π/a.

### Eq. 12: band gap

$$\Delta=2\,|t_1-t_2|$$

The gap is 2×(smallest |E|) and closes only when t₁ = t₂, which is the uniform chain.

![SSH bands](results/notebook1_crystals_and_bands/figures/04_ssh_band_structures.png)

*Solid: conduction bands for three t₂/t₁ values. Dashed: valence bands.*

![SSH gap vs ratio](results/notebook1_crystals_and_bands/figures/05_ssh_gap_vs_hopping_ratio.png)

![SSH gap at three ratios](results/notebook1_crystals_and_bands/figures/06_ssh_gap_at_three_ratios.png)

| t₂/t₁ | Gap (numeric from H(k)) | Gap (2\|t₁−t₂\|) |
|---|---|---|
| 0.4 | 1.2 | 1.2 |
| 1.0 | ≈ 0 (2×10⁻¹⁶) | 0 |
| 1.6 | 1.2 | 1.2 |

Note that t₂/t₁ = 0.4 and 1.6 have **identical gaps** but are different phases (trivial for t₂ < t₁, topological for t₂ > t₁). The gap cannot tell them apart, and the notebook does not compute a winding number or edge states. Also, the plots in Cells 9-10 draw the analytic Eq. 12; only the band plot (Cell 8) diagonalizes H(k).

---

## A5. Kronig-Penney model (Cells 11, 13, 14)

### Eq. 13: the Kronig-Penney dispersion relation

$$\cos(ka)=\cos(qa)+P\,\frac{\sin(qa)}{qa}$$

For a periodic array of delta-function barriers, this is the condition for Bloch states to exist. The dimensionless strength P measures how opaque each barrier is (for a comb V = α Σ δ(x − na), P = mαa/ħ²). A real k exists only where |right-hand side| ≤ 1; those energies form the **allowed bands**, and the rest are **gaps**.

### Eq. 14: energy and wavevector

$$E=\frac{\hbar^2q^2}{2m},\qquad k=\frac{1}{a}\arccos\!\Big[\cos(qa)+P\,\tfrac{\sin(qa)}{qa}\Big]$$

q is the wavevector *inside* the free region. The notebook sweeps q (0.01 to 25), keeps points with |RHS| ≤ 1, and computes k and E. Both ±k are plotted. (Cell 11 writes `E = hbar*q**2/(2*m)` without squaring ħ, which is harmless since ħ = 1.)

### Eq. 15: zone schemes

* **Reduced zone:** all bands folded into −π/a ≤ k ≤ π/a.
* **Repeated (periodic) zone:** the same bands translated by n·2π/a (Cell 11 uses n = −2…2).

Band edges for a = 0.5, ħ = m = 1 (computed on a fine q-grid; the zone-edge free-electron energy is 19.74):

| P | Band 1 (E_low-E_high) | Band 2 | Gap 1 |
|---|---|---|---|
| 0.1 | 0.39-19.74 | 20.53-78.96 | 0.79 |
| 1 | 3.41-19.74 | 26.99-78.96 | 7.25 |
| 5 | 10.44-19.74 | 45.34-78.96 | 25.60 |
| 20 | 16.33-19.74 | 65.56-78.96 | 45.82 |

As P → 0 the bands approach free electrons, and as P grows the bands narrow toward the isolated-well levels, which is the tight-binding limit of Part A1.

![KP zone schemes](results/notebook1_crystals_and_bands/figures/07_kronig_penney_zone_schemes.png)

> **This figure shows only P = 20.** The titles list P = [0.1, 1, 5, 20], but the loop in Cell 11 re-initializes its lists on every P, so only the last value survives. The left panel's two narrow bands (about 16.3-19.7 and 65.6-79.0) match my P = 20 table row exactly. The right panel is a repeated-zone plot made from the +k branch only (the loop reuses the leftover `sign`), not a true extended-zone scheme.

---

## A6. Real-space FD with a periodic potential (Cell 15)

### Eq. 16: Hamiltonian

$$H=-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}+V(x)\ \to\ \text{diag}=2t_{fd}+V(x_i),\quad\text{off-diag}=-t_{fd},\quad V(x)=2V_0\cos\!\frac{2\pi x}{a}$$

One unit cell of length a = 0.5 is split into 100 points (dx = 0.005, t_fd = ħ²/2m dx² = 20 000), with V₀ = 1.5.

### Eq. 17: Bloch boundary condition

$$H_{N,1}=-t_{fd}\,e^{+ika},\qquad H_{1,N}=-t_{fd}\,e^{-ika}$$

The wrap-around corner elements carry a phase e^(±ika) instead of the usual periodic +1. Solving this small Hermitian matrix at each k gives the bands E_n(k). This is the standard way to get a band structure from a real-space grid.

For a pure cosine potential with Fourier component V₀, first-order perturbation theory predicts a first gap of **2V₀ = 3.0** at the zone edge k = π/a, centred on the free-electron energy 19.74.

![KP vs FD](results/notebook1_crystals_and_bands/figures/08_kronig_penney_vs_finite_difference.png)

| Quantity (FD, a = 0.5) | Value |
|---|---|
| E₁ at k = 0 | −0.057 |
| Zone-edge energies (bands 1, 2) | 18.223, 21.223 |
| First gap | **2.9997** (expected 2V₀ = 3.0) |

The FD result is right for its own potential. **However, the Kronig-Penney curve (delta comb, P = 1) and the FD curves (cosine potential) are different physical systems**, so they are not expected to coincide. For example, the KP band bottom is at 3.4, while the cosine potential has zero average and the FD band bottom is at −0.06.

---

# Part B: Notebook 2, `02_crystals_and_lattices_nb2.ipynb`

## B1. Square and honeycomb lattice bands (Cells 2-3)

### Eq. 18: square lattice

$$E(\mathbf k)=-2t\bigl[\cos(k_xa)+\cos(k_ya)\bigr]$$

Nearest-neighbour hopping on a square lattice, with t = 2.7 eV and a = 0.246 nm. One band of total width 8t = 21.6 eV (from −10.8 to +10.8 eV).

### Eq. 19: honeycomb lattice (graphene)

$$E_\pm(\mathbf k)=\pm t\,\bigl|h(\mathbf k)\bigr|=\pm t\sqrt{3+f(\mathbf k)}$$

Two atoms per cell (sublattices A and B) give a conduction band (+) and a valence band (−). Here a_cc = 0.142 nm is the carbon-carbon distance.

### Eq. 20: the structure factor

The standard result for graphene is

$$f_{\text{std}}(\mathbf k)=2\cos\!\left(\sqrt3\,k_ya_{cc}\right)+4\cos\!\left(\tfrac{\sqrt3}{2}k_ya_{cc}\right)\cos\!\left(\tfrac32k_xa_{cc}\right)$$

This gives |h| = 0 at the corner points K, K′ (for example **K = (0, 4π/(3√3 a_cc))**), where the two bands touch in **Dirac cones** and E ≈ ±ħv_F|q| with v_F = 3ta_cc/2ħ. The maximum energy is 3t = 8.1 eV.

The notebook instead uses

$$f_{\text{nb}}=1+4\cos\!\left(\tfrac{\sqrt3}{2}k_ya_{cc}\right)\cos\!\left(\tfrac32k_xa_{cc}\right)+4\cos^2\!\left(\tfrac32k_ya_{cc}\right)$$

![Square vs honeycomb](results/notebook2_crystals_and_lattices/figures/01_square_vs_honeycomb_bands_3d.png)

| | Notebook formula | Standard graphene formula |
|---|---|---|
| Minimum of 3 + f | 1.346 | 0 |
| Minimum \|E\| (band gap/2) | **3.13 eV** (no touching) | 0 (Dirac points) |
| Maximum \|E\| | 9.35 eV | 8.10 eV = 3t |
| Value at K = (0, 4π/(3√3 a_cc)) | 3 + f = 5.13 | 0 |

So the saved honeycomb surfaces do **not** show Dirac cones, and the two bands are separated everywhere by at least 2 × 3.13 ≈ 6.3 eV, despite the plot's title. The square-lattice panel is correct. (The kx, ky grid spans ±4π/(3√3 a_cc) = ±17.03 nm⁻¹, which is larger than the square-lattice FBZ, so that plot shows more than one zone.)

---

## B2. Density of states (Cells 4-8)

### Eq. 21: definition

$$g(E)=\frac1N\sum_{\mathbf k}\delta\bigl(E-E_{\mathbf k}\bigr)$$

The number of states per unit energy, per site. It integrates to 1 over the whole band.

### Eq. 22: Lorentzian broadening

$$\delta(E-E_{\mathbf k})\ \to\ \frac1\pi\,\frac{\eta}{(E-E_{\mathbf k})^2+\eta^2},\qquad\eta=0.05$$

A numerical histogram of discrete k-points is spiky, so each delta function is replaced by a narrow Lorentzian of half-width η. The sum over all k-points then gives a smooth curve.

### Eq. 23: tight-binding bands for the three lattices (t = 1)

$$E_{1D}=-2t\cos k,\quad E_{2D}=-2t(\cos k_x+\cos k_y),\quad E_{3D}=-2t(\cos k_x+\cos k_y+\cos k_z)$$

Sampled on 5000 (1D), 300×300 (2D) and 60×60×60 (3D) k-points. The band bottoms are at −2t, −4t, −6t. The coarse 3D grid combined with small η causes the ripples visible in the 3D curve.

![DOS](results/notebook2_crystals_and_lattices/figures/02_dos_1d_2d_3d_lattices.png)

### Eq. 24: continuum (parabolic-band) DOS near the band bottom

$$g_{1D}(\varepsilon)=\frac1\pi\sqrt{\frac{m^*}{2\hbar^2\varepsilon}},\qquad g_{2D}(\varepsilon)=\frac{m^*}{2\pi\hbar^2},\qquad g_{3D}(\varepsilon)=\frac{1}{4\pi^2}\left(\frac{2m^*}{\hbar^2}\right)^{3/2}\sqrt{\varepsilon}$$

*Spinless, per unit length/area/volume.* ε is the energy measured **from the band bottom**. This gives the key signatures: a 1/√ε divergence in 1D, a constant in 2D and a √ε onset in 3D. For the tight-binding band, m\* = ħ²/(2ta²) = **0.5** (t = a = ħ = 1).

The notebook's "Eq. (4)" curves are:

$$g_{1D}=\frac{1}{\pi\hbar}\sqrt{\frac{2m^*}{E}},\qquad g_{2D}=\frac{m^*}{\pi\hbar^2 t},\qquad g_{3D}=\frac{1}{2\pi^2}\left(\frac{2m^*}{\hbar^2}\right)^{3/2}\sqrt E,\qquad m^*=1$$

![DOS comparison](results/notebook2_crystals_and_lattices/figures/03_dos_numerical_vs_theory.png)

Comparison about 0.5 above the band bottom (numerics from a fine histogram at the ε = 0.475 bin centre, my own check; notebook theory evaluated at E = 0.5):

| Dimension | Numerical TB DOS | Continuum, m\*=0.5 (Eq. 24) | Notebook "theory" |
|---|---|---|---|
| 1D | 0.246 | 0.231 | 0.636 |
| 2D | 0.0848 | 0.0796 | 0.318 |
| 3D | 0.0186 | 0.0175 | 0.101 |

**With correct normalization the lattice DOS agrees with the continuum formulas to about 6-7%** (the remainder is band non-parabolicity). The notebook's red curves disagree by large factors (about 2.6×, 3.8× and 5.5×) for four reasons: they include a spin factor of 2 that the per-site numerics don't have, they use m\* = 1 rather than 0.5, the 2D constant is divided by t, and ε is not measured from the band bottom. See Section 7.

---

## B3. Effective mass versus band gap in III-V semiconductors (Cell 9)

### Eq. 25: two-band k·p (Kane) relation

$$\frac{m_0}{m^*}\approx\frac{E_p}{E_g}\quad\Longrightarrow\quad\frac{m^*}{m_0}=\frac{E_g}{E_p},\qquad E_p=28.8\ \text{eV}$$

In a narrow-gap semiconductor, the conduction and valence bands repel each other through the momentum matrix element (E_p). A smaller gap gives a more strongly curved conduction band and therefore a lighter electron.

![Mass vs gap](results/notebook2_crystals_and_lattices/figures/04_effective_mass_vs_band_gap_III-V.png)

| Material | E_g (eV) | m\*/m₀ measured | m\*/m₀ theory (E_g/E_p) | Measured / theory |
|---|---|---|---|---|
| InSb | 0.235 | 0.014 | 0.0082 | 1.72 |
| InAs | 0.417 | 0.026 | 0.0145 | 1.80 |
| GaSb | 0.812 | 0.039 | 0.0282 | 1.38 |
| InP | 1.424 | 0.079 | 0.0494 | 1.60 |
| GaAs | 1.519 | 0.063 | 0.0527 | 1.19 |

The trend is right (lighter electrons for smaller gaps), but the simplified relation sits **below every data point by 19-80%**. This is expected, since it ignores spin-orbit splitting and the contribution of remote bands. Also, the notebook's GaAs value (0.063) is lower than the widely quoted 0.067, which the same notebook uses in Section B4, so that entry is worth double-checking against its source.

---

## B4. Quantum confinement (Cells 10-13)

### Eq. 26: confinement energy unit

$$E_1=\frac{\hbar^2\pi^2}{2m^*L^2},\qquad m^*=0.067\,m_0\ (\text{GaAs}),\ L=10\ \text{nm}$$

The ground-state energy of an electron in an infinite well of width L. Cell 11 prints **E₁ = 56.12 meV**. I get 56.124 meV.

### Eq. 27: energy levels for each confinement

$$E_n^{\text{well}}=n^2E_1,\qquad E_{n_yn_z}^{\text{wire}}=(n_y^2+n_z^2)E_1,\qquad E_{n_xn_yn_z}^{\text{box}}=(n_x^2+n_y^2+n_z^2)E_1$$

A **quantum well** confines one direction (free in two), a **wire** confines two (free along one), and a **box/dot** confines all three. Each confined direction contributes an n²E₁ term.

### Eq. 28: DOS shape for each case (arbitrary units in the notebook)

$$g_{2D}(E)\propto\sum_n\Theta(E-E_n),\qquad g_{1D}(E)\propto\sum_{n_y,n_z}\frac{1}{\sqrt{E-E_{n_yn_z}}},\qquad g_{0D}(E)\propto\sum\delta(E-E_{n_xn_yn_z})\ \to\ \text{Lorentzian},\ \eta=3\ \text{meV}$$

* **2D:** each subband adds a step (Θ = Heaviside function).
* **1D:** each subband adds a 1/√ spike at its onset (capped at 15 to avoid the infinity).
* **0D:** fully discrete levels, drawn as Lorentzian-broadened peaks.

![Confinement DOS](results/notebook2_crystals_and_lattices/figures/05_confinement_dos_2d_1d_0d.png)

Features that fall inside the plotted 0-400 meV window (computed from Eq. 27):

| Structure | Feature energies (meV) |
|---|---|
| Well (2D) steps | 56.1 (n=1), 224.5 (n=2) |
| Wire (1D) spikes | 112.2 (1,1), 280.6 (1,2 and 2,1) |
| Dot (0D) peaks | 168.4 (111, non-degenerate), 336.7 (211, 121, 112: **3-fold**, so the peak is three times taller) |

The loops run to higher n (the 2D loop to n = 14), but only the features above lie below 400 meV. The infinite-barrier assumption overestimates confinement energies compared with a real finite GaAs/AlGaAs well.

---

## B5. Effective-mass curvature fit (Cell 14)

### Eq. 29: parabolic band and fit

$$E(k)=\frac{\hbar^2}{2m^*}k^2\ \Rightarrow\ E=c_2k^2+c_0,\qquad m^*=\frac{\hbar^2/2m_0}{c_2}\ \ [\text{units of }m_0]$$

Synthetic "band data" E = (3.81/0.067)k² (k = ±0.04 Å⁻¹, 21 points) with Gaussian noise σ = 10⁻⁴ eV (seed 42) is fitted against k².

![Curvature fit](results/notebook2_crystals_and_lattices/figures/06_curvature_fit_effective_mass.png)

| Quantity | Value |
|---|---|
| Fitted c₂ | 56.904 eV·Å² (true 56.866) |
| Extracted m\* | **0.0670 m₀** (0.066955) |
| True m\* | 0.0670 m₀ |

This is a self-consistency test, because the data are generated from the very formula being fitted. In the code, the variable named `hbar_eV_A` holds ħ²/2m₀ = 3.81 eV·Å², not ħ.

---

# 7. Known issues and fixes

### Notebook 1

| # | Where | Issue | Fix |
|---|---|---|---|
| 1 | Cell 11, 13 | Lists are re-created inside the `for p in P` loop, so only **P = 20** is plotted, but titles say P = [0.1, 1, 5, 20] | Move the list initialisation above the loop, plot each P separately |
| 2 | Cell 11, 13 | "Extended zone" uses the leftover `sign` variable and only shifts the +k branch by n·2π/a. That is a repeated-zone picture, not an extended-zone scheme | Build the extended scheme by assigning each band n to \|k\| in [(n−1)π/a, nπ/a] |
| 3 | Cell 15-16 | KP (delta comb) is compared with FD using a cosine potential. These are different systems, so they cannot agree | Use a delta/rectangular-barrier potential in the FD code, or compare with a nearly-free-electron calculation |
| 4 | Cells 6-16 | `k_grid` and `a` are reassigned in several cells. Run top to bottom, `k_grid` comes from Cell 8 (a = 1.0), so the FD curves span only **±0.5 of the BZ** for a = 0.5. The saved figure came from a different run order | Define `k_grid = np.linspace(-np.pi/a, np.pi/a, 300)` right before Cell 15 |
| 5 | Cells 8-9 | The SSH figure is split in two (Cell 8 shows only the left panel, Cell 9 only the right) | Create the figure and call `plt.show()` once after both panels |
| 6 | Cell 7 | Plot adds a dimensionless wavefunction to energies in eV | Scale the wavefunction (for example by the level spacing) or plot probability density separately |
| 7 | Cells 6, 8 | The SSH script overwrites `t2` (a list) with its loop variable; the FD cell builds a matrix `H` that is never used | Rename variables; delete or use `H` |
| 8 | Cells 9-10 | Gap plots draw the analytic Δ = 2\|t₁−t₂\|, not a computed gap | Take the gap from the diagonalization |
| 9 | General | Execution counts are out of order (for example Cell 6 ran as #15, after Cell 7 as #7). Warnings appear for `\AA` and similar escapes | Restart and Run All; use raw strings for LaTeX |

### Notebook 2

| # | Where | Issue | Fix |
|---|---|---|---|
| 1 | Cells 4-8 | **`E_grid` is never defined** (the only definition is commented out). The notebook only worked through a stale kernel variable. Running top to bottom raises `NameError` in Cells 4-8 (I confirmed this) | Uncomment `E_grid = np.linspace(0.005, 3.0, 400)` (or define it properly) |
| 2 | Cell 7 (figure 02) | The DOS curves were computed on E_grid (about 0.005-3) but plotted against `E_grid1` (−6.2 to 6.2), so the **energy axis of figure 02 is wrong**: a 3 eV span is drawn over 12.4 eV, and every feature sits at the wrong energy | Plot against the same grid used to compute |
| 3 | Cells 4-8 | Numerical DOS is computed with E from the **band centre** but labelled "Energy from band bottom". The band bottom is at −2t, −4t, −6t in 1D, 2D, 3D | Evaluate on E = −2dt + ε, or shift the numerics before comparing |
| 4 | Cell 4-6 | The "theory" curves use m\* = 1 (the lattice has m\* = 1/2t = 0.5), include a spin factor 2 that the per-site numerics lack, and 2D is divided by t | Use Eq. 24 with m\* = 0.5; apply the same spin convention to both |
| 5 | Cell 2 | **Honeycomb structure factor is wrong.** It has no Dirac points (min \|E\| = 3.13 eV) and a maximum of 9.35 eV instead of 3t = 8.1 eV | Use f = 2cos(√3 k_y a) + 4cos(√3 k_y a/2)cos(3 k_x a/2) |
| 6 | Cell 9 | Two-band formula ignores spin-orbit and remote bands, so it underestimates m\* by 19-80%. GaAs data point 0.063 differs from the 0.067 used later | Use the fuller k·p expression or state the limitation; check the GaAs value |
| 7 | Cell 12 | DOS in arbitrary units, spike cap 15, and η = 3 meV are arbitrary | State them, or compute with proper prefactors |
| 8 | Cell 14 | Variable `hbar_eV_A` is actually ħ²/2m₀; the test is circular | Rename; for a stronger test fit a band from an actual diagonalization |

**Reproducibility:** Notebook 1 runs top to bottom without errors but gives the half-zone FD plot (Notebook 1, item 4). Notebook 2 stops with a `NameError` until `E_grid` is defined (Notebook 2, item 1).

---

# 8. Results at a glance

| Topic | Result |
|---|---|
| Tight-binding m\* (Notebook 1) | 2.0142 vs theory 2.0000 (0.71% error, window-limited) |
| FD m\* | 1.0008 vs 1.0000 |
| Particle in a box | E₁ = 0.000300 eV; E_n/E₁ → n²; matches exact to 8×10⁻¹⁵ eV |
| SSH gap | 1.2 (t₂/t₁ = 0.4), 0 (1.0), 1.2 (1.6) |
| KP first gap, a = 0.5 | 0.79, 7.25, 25.6, 45.8 for P = 0.1, 1, 5, 20 |
| FD cosine potential first gap | 2.9997 (expected 2V₀ = 3.0) |
| Square lattice bandwidth | 8t = 21.6 eV |
| Honeycomb (as coded) | gap > 6 eV, no Dirac points (formula wrong) |
| DOS band-edge check (corrected) | numerics within ~6-7% of the continuum formulas |
| III-V two-band mass | theory 19-80% below measured |
| Confinement unit E₁ (GaAs, 10 nm) | 56.12 meV |
| Curvature fit m\* | 0.06696 m₀ (true 0.067) |

**Honest summary:** the tight-binding, FD, particle-in-a-box, SSH and confinement results are correct and I could reproduce them independently. The Kronig-Penney figure, the DOS figures and the honeycomb figure contain plotting or formula errors that the notebook's titles hide; Section 7 gives the fixes.
