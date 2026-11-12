import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

export default defineConfig({
  site: 'https://dtereshchenko-dev.github.io',
  base: '/Trajecta/docs',
  redirects: {
    '/verification/dynamics-apophis/': '/Trajecta/docs/numerical-investigations/apophis-2029/',
  },
  integrations: [
    starlight({
      title: 'Trajecta',
      description: 'Documentation for Trajecta orbit determination and impact analysis.',
      logo: { src: '../TrajectaUI/resources/icons/workspace/app-logo.svg' },
      social: [{ icon: 'github', label: 'GitHub', href: 'https://github.com/dtereshchenko-dev/Trajecta' }],
      customCss: ['./src/styles/trajecta.css'],
      components: { Footer: './src/components/Footer.astro' },
      sidebar: [
        { label: 'Start here', collapsed: true, items: [
          { label: 'Getting started', slug: 'getting-started' },
          { label: 'Data sources', slug: 'data-sources' },
        ] },
        { label: 'User guide', collapsed: true, items: [
          { label: 'Orbit determination', slug: 'orbit-determination' },
          { label: 'Observations', slug: 'guide/observations' },
          { label: 'Fit settings', slug: 'guide/fit-settings' },
          { label: 'Ephemerides and close approaches', slug: 'ephemerides-and-close-approaches' },
          { label: 'Impact analysis', slug: 'impact-analysis' },
          { label: 'Orbital Monte Carlo', slug: 'guide/monte-carlo' },
          { label: 'Manifold of Variations', slug: 'guide/mov' },
          { label: 'Atmospheric Monte Carlo', slug: 'guide/atmospheric-mc' },
        ] },
        { label: 'Verification', collapsed: true, items: [
          { label: 'Verification map', slug: 'verification/map' },
          { label: 'Dynamics', collapsed: true, items: [
            { label: 'Integrator round trips', slug: 'verification/integrator-round-trip' },
            { label: 'Holman: dynamics and STM', slug: 'verification/dynamics-holman' },
            { label: 'Parijskij and Ceres', slug: 'verification/dynamics-ceres' },
            { label: '2023 BU', slug: 'verification/dynamics-2023bu' },
          ] },
          { label: 'MOV', slug: 'verification/2019-mo' },
          { label: 'Entry and atmosphere', collapsed: true, items: [
            { label: '2023 CX1', slug: 'verification/2023-cx1' },
            { label: '2008 TC3', slug: 'verification/2008-tc3' },
          ] },
        ] },
        { label: 'Numerical investigations', collapsed: true, items: [
          { label: 'Apophis 2029', slug: 'numerical-investigations/apophis-2029' },
        ] },
        { label: 'Technical reference', collapsed: true, items: [
          { label: 'Dynamics and time', slug: 'reference/dynamics' },
          { label: 'Observation model', slug: 'reference/observation-model' },
          { label: 'Uncertainty', slug: 'reference/uncertainty' },
          { label: 'Atmospheric model', slug: 'reference/atmospheric-model' },
          { label: 'Photometry', slug: 'reference/photometry' },
          { label: 'Files and outputs', slug: 'reference/formats' },
          { label: 'MC archive reference', slug: 'reference/mc-archive' },
          { label: 'CLI', slug: 'reference/cli' },
          { label: 'For developers', slug: 'reference/developers' },
        ] },
      ],
    }),
  ],
});
