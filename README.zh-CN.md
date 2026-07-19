# RSA-4096 优化实现

[![C](https://img.shields.io/badge/C-GCC-A8B9CC?style=flat-square&logo=c&logoColor=white)](https://gcc.gnu.org/)
[![Validation](https://img.shields.io/github/actions/workflow/status/LinkisLethe/RSA4096_Optimizer/validate.yml?branch=main&style=flat-square&label=validation)](https://github.com/LinkisLethe/RSA4096_Optimizer/actions/workflows/validate.yml)
[![License](https://img.shields.io/badge/license-MIT-2ea44f?style=flat-square)](LICENSE)

[English](README.md)

这是一个面向 CPU 的 C 语言 RSA-4096 实现。项目在保留固定往返测试和计时检查的前提下，优化了模板中的大整数与 RSA 执行路径。

## 优化与结果

- CIOS Montgomery 乘法减少重复模约减的开销。
- 滑动窗口模幂优化通用私钥指数运算路径。
- 基于 CRT 的私钥解密将一次 4096 位模幂拆成模 `p` 和模 `q` 的两次较小运算。
- `e = 65537` 专用路径通过 16 次平方和 1 次乘法完成公钥加密。

归档的服务器数据中，每个阶段各运行 3 次。耗时会受硬件影响，表中数值只代表已保存的实测结果。

| 版本 | 平均加密时间 | 平均解密时间 |
|---|---:|---:|
| 基线 | `0.032837 s` | `9.647460 s` |
| 最终阶段 | `0.004067 s` | `0.686095 s` |

每轮数据和阶段间倍率见 [`benchmarks/performance_summary.md`](benchmarks/performance_summary.md) 及对应 CSV 文件。

## 构建与验证

环境要求：Linux、macOS 或兼容环境中的 GCC 与 GNU Make。

```bash
make
./release/main
```

程序先验证固定 1,000 字节数据的加密解密往返，再对 10,019 字节消息执行 10 轮计时测试。`bignum.c` 实现大整数运算，`rsa.c` 实现 RSA 执行路径，`main.c` 保存测试与计时逻辑。GitHub Actions 会构建并运行同一程序。

## 局限与许可证

代码基于提供的 RSA-4096 模板，并在 `keys.h` 中包含固定演示密钥。该实现没有经过生产级密码学审计，不能用于保护真实数据或管理真实私钥。

整理后的源码使用 [MIT License](LICENSE)。
