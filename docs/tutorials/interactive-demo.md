# Interactive tutorial

The public entry point is a JupyterLab tutorial. It provides a guided BHJet and
Gammapy exercise, a widget-based parameter explorer, and links to the advanced
notebooks. The tutorial contains only public, versioned example data.

## Binder launch

[![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rubyduncan/modBHJet/tutorial-v0.1.0?urlpath=lab/tree/tutorials/00_welcome.ipynb)

The launch URL deliberately references the immutable `tutorial-v0.1.0` release
tag instead of `main`:

```text
https://mybinder.org/v2/gh/rubyduncan/modBHJet/tutorial-v0.1.0?urlpath=lab/tree/tutorials/00_welcome.ipynb
```

Before publishing this link, create and push the named tag only after the
`Tutorial image` workflow has passed. Each later tutorial release must update
both occurrences of the tag in this page.

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
