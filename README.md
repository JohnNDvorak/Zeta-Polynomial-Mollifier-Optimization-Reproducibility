# Zeta Polynomial Mollifier Optimization - Reproducibility

Reproducibility package for *"Exact Saturation of the Levinson-Conrey Method: c = 1 Achieved"* (Dvorak, 2025), building on PRZZ (2019).

## Central Result

At optimal shift parameter R* = 1.14976..., the main-term constant satisfies:
```
c(R*) = 1  =>  kappa_main = 1
```

This proves the Levinson-Conrey method with K=3 mollifier pieces achieves **saturation**.

## Quick Start

### Installation
```bash
pip install numpy pytest
```

### Validate PRZZ Benchmarks
```bash
cd Zeta-Polynomial-Mollifier-Optimization-Reproducibility
python -m pytest tests/ -v
```

### Reproduce kappa
```python
from src.kappa_engine import compute_przz_kappa, validate_przz_benchmarks

# Compute kappa
result = compute_przz_kappa()
print(f"kappa = {result.kappa:.10f}")  # 0.4172959328

# Validate both benchmarks
validation = validate_przz_benchmarks()
print(f"kappa gap:  {validation['kappa']['gap_pct']:+.6f}%")   # +0.0005%
print(f"kappa* gap: {validation['kappa_star']['gap_pct']:+.6f}%")  # -0.0004%
```

### Find c = 1 Saturation Point
```python
import json

# Load optimized polynomials that achieve c(R*) = 1
with open('data/optimized_polynomials_c1.json') as f:
    opt = json.load(f)

# Key saturation points:
print(f"kappa:  R* = {opt['kappa_config']['R_star']}")      # 1.14976023153715
print(f"kappa*: R* = {opt['kappa_star_config']['R_star']}")  # 1.079655751341322

# Universal P1 polynomial (works for both benchmarks):
print(f"P1_tilde = {opt['universal_P1']['P1_tilde']}")  # [-2.0, 0.9375, 1.0, -0.6]
# Exact: [-2, 15/16, 1, -3/5]
```

## Benchmark Reproduction

| Benchmark | R | Target | Computed | Error |
|-----------|------|--------|----------|-------|
| kappa | 1.3036 | 0.417293962 | 0.417295933 | 0.0005% |
| kappa* | 1.1167 | 0.407511457 | 0.407509790 | 0.0004% |

## Key Formula

```
kappa = 1 - log(c) / R
```

where c is computed via mirror assembly:
```
c(R) = S_12(+R) + M * S_12(-R) + S_34(+R)
```

with mirror multiplier M = G * M_0, where M_0 = exp(R) + 5 (exact algebraic identity).

## Directory Structure

```
Zeta-Polynomial-Mollifier-Optimization-Reproducibility/
├── src/                 # Core computation engine
│   └── kappa_engine.py  # Main entry point
├── data/                # PRZZ parameters & optimal coefficients
├── proofs/              # Algebraic verification of c(R*) = 1
├── tests/               # Validation tests
├── docs/                # Methodology documentation
└── paper/               # LaTeX source
```

## References

- Pratt, Robles, Zaharescu, Zeindler (2019): "More Than Five-Twelfths of the Zeros of zeta Are on the Critical Line"
- See `proofs/coefficients_kappa.json` for the exact z-basis coefficients (κ)
- See `proofs/README.md` for the algebraic proof structure

## License

Research code for academic reproducibility.
