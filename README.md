# qpc — 2D spin-DFT of a quantum point contact

A Python implementation of two-dimensional Kohn–Sham spin-density-functional theory (LSDA) in a
plane-wave basis, applied to a GaAs quantum point contact (QPC). The goal is to reproduce the
spontaneous **local moment formation** (a net spin of about one electron, S = 1/2, bound at the
QPC) proposed as the origin of the 0.7 conductance anomaly in

> K. Hirose, Y. Meir and N. S. Wingreen, *Local moment formation in quantum point contacts*,
> Phys. Rev. Lett. **90**, 026804 (2003).

This was written as an educational research project. Every function documents the formula it
implements, its units and the reference it comes from.

## Physics

For each spin s = ↑, ↓ the Kohn–Sham equation is solved self-consistently:

```
[ -½∇² + V_ext(x,y) + V_H[n](x,y) + V_xc,s[n↑,n↓](x,y) − σ_s E_Z/2 ] ψ_is = ε_is ψ_is
n_s(x,y) = Σ_i f(ε_is − μ) |ψ_is(x,y)|²
```

| Term | Model | Module |
|---|---|---|
| `V_ext` | parabolic wire `½ω_y² y²` plus the HMW QPC potential with a Pöschl–Teller barrier `V₀/cosh²(x/d)` (HMW Eq. 1) | `potential.py` |
| `V_H` | gate-screened Hartree, kernel `1/ρ − 1/√(ρ² + a²)` with an image charge at a = 100 nm | `hartree.py` |
| `V_xc,s` | 2D LSDA: exact exchange of the 2D electron gas plus Tanatar–Ceperley correlation (PRB 39, 5005, 1989); spin interpolation `quadratic`, `exchange`, `power:p` or `mixed:w` | `xc.py` |
| `E_Z` | Zeeman splitting `g* μ_B B`, used to break spin symmetry | `scf.py` |

The chemical potential μ is fixed by an interacting clean wire (no QPC) at a linear density of
0.028 nm⁻¹. The QPC is then solved at that fixed μ.

## Numerical method

- **Units:** effective GaAs atomic units (`m* = 0.067`, `ε_r = 12.9`): 1 Ha* = 10.956 meV and
  1 a* = 10.189 nm. Physical units are used only at input and output (`units.py`).
- **Basis:** plane waves `exp(iG·r)` with `|G|²/2 ≤ E_cut` on a periodic supercell
  (typically 5000 nm × 320 nm, E_cut = 15 meV, about 3400 plane waves). The kinetic term is
  diagonal in G and the potential is applied with FFTs (`grid.py`, `hamiltonian.py`).
- **Eigensolver:** block LOBPCG (Knyazev 2001, with the basis of Hetmaniuk & Lehoucq 2006) and
  a kinetic preconditioner, or dense diagonalisation for small tests (`solver.py`).
- **Occupations:** Fermi–Dirac at kT = 0.05 meV; μ is found by root-finding at fixed N
  (`occupation.py`).
- **Mixing:** Pulay/DIIS density mixing (Pulay 1980; Kresse & Furthmüller 1996)
  (`mixing.py`).
- **Barrier Fourier analysis:** NumPy FFT, mpmath quadrature and the exact transform of sech²
  for benchmarking (`fourier.py`).

## Package layout

```
qpc/
├── units.py        effective atomic units and conversions
├── grid.py         real-space / plane-wave grid from an energy cutoff
├── fourier.py      Fourier transforms of the barrier (numpy, mpmath, exact sech²)
├── potential.py    HMW QPC potential
├── hamiltonian.py  plane-wave Hamiltonian H = T + V, FFT-based H·X
├── solver.py       block LOBPCG and dense eigensolvers
├── occupation.py   Fermi function, chemical potential, density
├── hartree.py      gate-screened Hartree potential
├── xc.py           2D LSDA exchange-correlation energy and potentials
├── mixing.py       Pulay / DIIS mixing
├── scf.py          Kohn–Sham self-consistency, clean wire, lead models
├── energy.py       Kohn–Sham total energy and grand potential
└── analysis.py     net spin, barrier profile, subbands, local density of states
```

## Requirements

Python ≥ 3.12, `numpy`, `scipy`, `mpmath`.

```bash
pip install numpy scipy mpmath
```

## Minimal example

Run from the folder that contains `qpc/`:

```python
from qpc.units import Units
from qpc.grid import grid_from_cutoff
from qpc.hamiltonian import Hamiltonian
from qpc.hartree import Hartree
from qpc.potential import QPCParams
from qpc.scf import SCFParams, lead_state, run_scf
from qpc.analysis import net_spin

U = Units()
LX_NM, LY_NM, ECUT_MEV, N1D_NM = 5000.0, 320.0, 15.0, 2.8e-2

grid = grid_from_cutoff(U.nm_to_au(LX_NM), U.nm_to_au(LY_NM), U.meV_to_au(ECUT_MEV))
ham = Hamiltonian(grid, U.meV_to_au(ECUT_MEV))
hart = Hartree(grid, SCFParams().a_metal_au(U))
qpc = QPCParams(hbar_wx_meV=1.5)          # barrier curvature: 1.0, 1.5 or 2.0 meV in HMW

p = SCFParams(B_T=1.0, interp="exchange")  # small field to break spin symmetry
lead, mu, V_ext = lead_state(ham, hart, qpc, p, N=N1D_NM * LX_NM)   # clean wire -> mu
res = run_scf(ham, hart, V_ext, p, mu=mu)                           # QPC at fixed mu

print("converged:", res.converged, " net spin M =", net_spin(res, grid))
```

A full run takes on the order of tens of minutes on a laptop. To find the B = 0 local moment, the
field is ramped down (for example 6 T → 0), starting each step from the previous converged
densities through `n_init`. A run that starts exactly spin-symmetric at B = 0 stays symmetric,
because `V_xc,↑ = V_xc,↓` at ξ = 0.

## Notes

- The default spin interpolation in `SCFParams` is `quadratic` (Tanatar–Ceperley's
  prescription). The local moments in the project report were obtained with `exchange`, the
  interpolation of Koskinen, Manninen & Reimann (PRL 79, 1389, 1997). HMW do not state which
  one they used.
- Some docstrings refer to `docs/SPEC.md` and `docs/OPEN_QUESTIONS.md`, which are project notes
  not included in this repository.

## References

1. K. Hirose, Y. Meir, N. S. Wingreen, Phys. Rev. Lett. **90**, 026804 (2003).
2. W. Kohn, L. J. Sham, Phys. Rev. **140**, A1133 (1965).
3. U. von Barth, L. Hedin, J. Phys. C **5**, 1629 (1972).
4. B. Tanatar, D. M. Ceperley, Phys. Rev. B **39**, 5005 (1989).
5. M. Koskinen, M. Manninen, S. M. Reimann, Phys. Rev. Lett. **79**, 1389 (1997).
6. A. V. Knyazev, SIAM J. Sci. Comput. **23**, 517 (2001).
7. P. Pulay, Chem. Phys. Lett. **73**, 393 (1980).
