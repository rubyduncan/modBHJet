# Tutorial hosting handoff

This page describes the deployment contract for a BHJet school. It is designed
for an institutional JupyterHub, BinderHub, or managed academic JupyterHub;
the organisers choose the host and supported login provider.

## Release artifact

The `Tutorial image` GitHub Actions workflow builds the same environment used
by MyBinder and publishes it to GitHub Container Registry when a
`tutorial-v*` tag is pushed. Deploy the image **by digest**, recorded from the
successful workflow, rather than a mutable tag:

```text
ghcr.io/rubyduncan/modbhjet-tutorial@sha256:<published-digest>
```

The release includes BHJet, a fixed Kariba revision, the public tutorial
notebooks, and the bundled NGC 4261 data. No external download is needed by a
participant server.

## Initial JupyterHub settings

For a 50-person tutorial, configure one isolated single-user server per
participant with the following starting values. The hosting team must validate
these values against the representative tutorial workload before the event.

```yaml
singleuser:
  image: ghcr.io/rubyduncan/modbhjet-tutorial@sha256:<published-digest>
  cpu:
    guarantee: 0.5
    limit: 1
  memory:
    guarantee: 1G
    limit: 2G
  extraEnv:
    OMP_NUM_THREADS: "1"
    OPENBLAS_NUM_THREADS: "1"
    MKL_NUM_THREADS: "1"
    NUMEXPR_NUM_THREADS: "1"
```

Reserve at least 25 CPUs and 50 GiB of allocatable memory for the requested
capacity, plus the hub and system overhead. The limits allow a temporary peak
of up to 50 CPUs and 100 GiB; do not promise that peak unless the host can
provide it.

## Organiser checklist

1. Choose institutional SSO when it is available; otherwise create temporary
   accounts from the event roster. Anonymous access is not the school default.
2. Run a small organiser pilot, then a representative 50-user launch and
   notebook-execution test. Record launch time, memory peak, and any failed
   kernels.
3. Show participants that their server is personal but disposable. Give them a
   download/export step at the end of the school.
4. Name an event administrator who can restart servers and communicate the
   fallback link. Keep the public MyBinder link as a demonstration fallback,
   not as reserved classroom capacity.

## Data and security

The repository data are public and read-only. Do not add credentials, private
observations, or participant uploads to the image. If a later tutorial needs a
large public data set, place it in read-only object storage and provide a
small, tested cache/download step with an in-repository fallback sample.

Further reading: [MyBinder usage guidelines](https://mybinder.readthedocs.io/en/latest/about/user-guidelines.html)
and [JupyterHub capacity planning](https://jupyterhub.readthedocs.io/en/stable/explanation/capacity-planning.html).
