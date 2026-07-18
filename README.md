# RSA-4096 Optimization

[![C](https://img.shields.io/badge/C-GCC-A8B9CC?style=flat-square&logo=c&logoColor=white)](https://gcc.gnu.org/)
[![Validation](https://img.shields.io/github/actions/workflow/status/LinkisLethe/rsa4096-optimization/validate.yml?branch=main&style=flat-square&label=validation)](https://github.com/LinkisLethe/rsa4096-optimization/actions/workflows/validate.yml)
[![License](https://img.shields.io/badge/license-MIT-2ea44f?style=flat-square)](LICENSE)

[中文说明](README.zh-CN.md)

A CPU-focused RSA-4096 implementation in C. The project optimizes the supplied big-integer and RSA paths while retaining the fixed round-trip testbench and timing checks.

## Optimization and results

- CIOS Montgomery multiplication reduces repeated modular-reduction cost.
- Sliding-window exponentiation accelerates the general private-exponent path.
- CRT-based private decryption replaces one 4096-bit exponentiation with two smaller operations modulo `p` and `q`.
- A dedicated `e = 65537` path uses 16 squarings and one multiplication for public-key encryption.

The archived server measurements contain three runs per stage. Times are hardware-dependent and are retained as measured results, not universal performance claims.

| Version | Encryption average | Decryption average |
|---|---:|---:|
| Baseline | `0.032837 s` | `9.647460 s` |
| Final stage | `0.004067 s` | `0.686095 s` |

Per-run values and stage-to-stage ratios are available in [`benchmarks/performance_summary.md`](benchmarks/performance_summary.md) and the matching CSV file.

## Build and verify

Requirements: GCC and GNU Make on Linux, macOS, or a compatible environment.

```bash
make
./release/main
```

The executable first verifies a fixed 1,000-byte encrypt/decrypt round trip, then runs ten timed rounds with a 10,019-byte message. `bignum.c` contains large-integer arithmetic, `rsa.c` contains the RSA execution paths, and `main.c` contains the testbench and timing checks. GitHub Actions builds and runs the same program.

## Limits and licensing

This code is derived from a supplied RSA-4096 template and includes a fixed demonstration key in `keys.h`. It has not received a production cryptographic audit and should not protect real data or manage real private keys.

The curated source is available under the [MIT License](LICENSE).
