---
layout: page
permalink: /documentation/02_installation/
title: Installation
description: ~
order: 2
---

`fastGMA` is a C++ command-line tool for genomic mixed-model analysis. It is distributed as a prebuilt 64-bit Linux executable.

## 1. Prebuilt Linux executable

A prebuilt 64-bit Linux executable is available here:
[**fastGMA Linux Executable**](https://github.com/chaoning/fastGMA/raw/refs/heads/main/app/linux/fastgma)

After downloading it, make it executable and check that it runs:

```bash
chmod +x fastgma
./fastgma -h
```

This is the quickest option if you only need to run the program. The executable is statically linked against MKL and the Intel OpenMP runtime, so it needs no other libraries on the 64-bit Linux systems it has been tested on. MKL's static libraries still call `dlopen`, so check the binary on your target system. It sets `KMP_AFFINITY=disabled` unless you set `KMP_AFFINITY` yourself, which avoids an affinity error of the Intel OpenMP runtime on some systems.

## 2. Quick validation

Confirm that the executable starts correctly:

```bash
./fastgma -h
```

For an end-to-end test, run the example commands on each module page with the data in `example/`.

## 3. Troubleshooting

- If the executable aborts with an Intel OpenMP affinity error, run it with `KMP_AFFINITY=disabled`.
- If it fails on your system, please open an issue with the error message and the output of `uname -a`.
