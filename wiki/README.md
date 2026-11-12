# Trajecta documentation

This directory contains the documentation site, built with Astro Starlight.
The pages are Markdown files in `src/content/docs`; navigation is configured in
`astro.config.mjs`.

With Node.js 24 and npm installed:

```sh
cd wiki
npm ci
npm run dev
```

Open the address printed by Astro, including `/Trajecta/docs/`. To check the
production build and its search index locally:

```sh
npm run build
npm run preview
```

The GitHub Pages publication workflow builds the documentation from the default
branch and adds it under `/Trajecta/docs/` to the same deployment artifact as
the Windows and Linux update repositories. It runs when a stable release is
published or when the update publication workflow is dispatched manually for
the latest stable release. If `wiki` is absent from the default branch, the
workflow publishes only the update repositories; that deployment has no
documentation site. Publishing documentation separately would replace the
Pages deployment and remove the update repositories.
