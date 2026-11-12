---
title: Developer reference
description: Source layout, command-line usage, and C++ documentation.
---

This page is a starting point for working with the source code. A more detailed
developer reference will be added here later.

## Components

| Component | Role |
| --- | --- |
| OrbitalCore | Orbit determination, dynamics and impact-analysis APIs. |
| TrajectaUI | Qt/QML desktop application. |
| TrajectaCLI | Command-line client for the shared core. |
| StrewnField | Atmospheric entry and strewn-field simulation. |

The repository's [developer guide](https://github.com/dtereshchenko-dev/Trajecta/blob/master/DEVELOPMENT.md)
describes dependencies, runtime data, builds, tests and C++ examples.
The [CLI guide](https://github.com/dtereshchenko-dev/Trajecta/blob/master/TrajectaCLI/README.md)
documents command-line usage, and the
[StrewnField guide](https://github.com/dtereshchenko-dev/Trajecta/blob/master/StrewnField/README.md)
covers atmospheric inputs and outputs.

## API status

Public C++ headers live under `OrbitalCore/api/orbital_core/` and
`StrewnField/api/`. Internal `src/` headers are not a supported third-party
interface. The current API is a same-toolchain C++ interface, not a stable C
ABI across arbitrary compilers.
