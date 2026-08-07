# Interactive tutorial

The public entry point is a JupyterLab tutorial. It provides a guided BHJet and
Gammapy exercise, a widget-based parameter explorer, and links to the advanced
notebooks. The tutorial contains only public, versioned example data.

## Binder launch

[![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rubyduncan/modBHJet/bhjet_tutorial?urlpath=lab/tree/tutorials/00_welcome.ipynb)

The current development launch URL follows the `bhjet_tutorial` branch:

```text
https://mybinder.org/v2/gh/rubyduncan/modBHJet/bhjet_tutorial?urlpath=lab/tree/tutorials/00_welcome.ipynb
```

Before publishing a stable tutorial link, create and push a release tag only
after the `Tutorial image` workflow has passed. Replace `bhjet_tutorial` in the
badge and URL above with that tag (for example, `tutorial-v0.1.0`).

## Classroom use

MyBinder is an excellent public try-it link, but it is an ephemeral volunteer
service and should not be the only platform for a scheduled school. It has a
100-session limit per repository, does not retain user files, and can require
an event quota request. Use the same tested image in a managed JupyterHub for a
classroom event; the exact handoff is in
{doc}`../deployment/tutorial-hosting`.

Participants should download notebooks or results they want to keep before
their session ends. Do not put private data, credentials, or participant data
in this public repository.
