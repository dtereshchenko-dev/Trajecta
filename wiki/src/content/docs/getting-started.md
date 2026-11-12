---
title: Getting Started
description: Install Trajecta and determine an orbit from observations.
---

This guide takes you from installation to your first orbit solution.

## Install

Download the installer for your platform from
[GitHub Releases](https://github.com/dtereshchenko-dev/Trajecta/releases):

- **Windows:** run `TrajectaSetup-x86_64.exe` and follow the installer.
- **Linux:** make `TrajectaSetup-x86_64.run` executable, then run it and follow
  the installer. For example, in a terminal opened in the download folder:

  ```sh
  chmod +x TrajectaSetup-x86_64.run
  ./TrajectaSetup-x86_64.run
  ```

Installer names and available platforms may differ in older releases.

## Check the data

Open Trajecta and select **Data sources** in the sidebar. Check that the
planetary ephemeris and optical observatory catalogue are connected. These are
used to calculate the example orbit below.

Other sources are needed only for particular observations or settings. For
example, the massive-asteroid ephemeris is needed when those perturbations are
enabled. See [Data Sources](../data-sources/) for details and download options.

## Solve an orbit

Trajecta accepts optical and radar observations in MPC 80-column and ADES
XML/PSV formats. To try a calculation with a supplied example:

1. For your first calculation, use the default fit settings. Before opening the
   file, check in **Settings** that the non-gravitational model is off. You can
   explore the other settings after obtaining your first solution.
2. Select **Open observations** and choose `examples/NEO/2020VV.xml` from the
   installation folder, or drag the file into the Trajecta window. You can
   also use your own observation file.
3. Trajecta starts processing the file when it is opened. If the observation
   review appears, review the observations and continue to the calculation.
   When the calculation finishes, **Solution** displays the orbital elements,
   fit status and observation residuals.
4. Inspect the solution. Open **Covariance** for
   the parameter correlations and **Observations** to review the observations
   included in the fit.

A converged solution is a starting point: check the residuals and the fit
settings before using it for further calculations.

To change the settings after a solution is available, apply the changes and
select **Re-solve** to calculate the orbit again.

## Find the results

After a solution is available, select the folder icon in the top bar to open
its working folder. From the sidebar you can also generate **Ephemerides**,
search for **Close Approaches**, or use **Risk Monitor** for impact analysis.

## Next steps

- [Orbit Determination](../orbit-determination/)
- [Review observations](../guide/observations/)
- [Configure an orbit fit](../guide/fit-settings/)
- [Ephemerides and Close Approaches](../ephemerides-and-close-approaches/)
- [Impact Analysis](../impact-analysis/)
