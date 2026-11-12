---
title: Command-line interface
description: Status and planned direction of the Trajecta command-line interface.
---

`TrajectaCLI` currently serves primarily as a development and testing tool for
the shared orbit-determination and risk-analysis code. Its command-line options
are not yet designed as a complete, stable workflow for general use.

The planned direction is to configure calculations through JSON files, as
`StrewnCLI` does for atmospheric calculations. The CLI documentation will be
written once that transition is complete.

For the current commands and build instructions, see the
[TrajectaCLI README](https://github.com/dtereshchenko-dev/Trajecta/blob/master/TrajectaCLI/README.md).
The existing JSON-based atmospheric interface is described in the
[StrewnField README](https://github.com/dtereshchenko-dev/Trajecta/blob/master/StrewnField/README.md).
