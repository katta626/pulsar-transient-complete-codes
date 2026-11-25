# Linux systems on this Mac

Every Linux system runs from the terminal with `lx`. You never need to open the Docker window. The engine behind it is still Docker Desktop's, so Docker must stay installed; `lx` asks it to start in the background when needed. If a window appears, turn off opening the Dashboard at startup in Docker's settings.

## Set up (once)

From this folder:

```bash
mkdir -p ~/.local/bin && ln -s "$PWD/lx" ~/.local/bin/lx
```

`~/.local/bin` must be on your `PATH`. Without this step, run `docker/lx` from the repository root.

## Use

The folder you are in appears inside Linux as `/work`. Files written there stay on the Mac. Anything written elsewhere inside Linux is discarded on exit.

| Command | What it does |
|---|---|
| `lx list` | the systems and whether each is installed |
| `lx shell presto` | interactive Linux shell |
| `lx run presto accelsearch -zmax 200 x.fft` | one command |
| `lx run pulsarx psrfold_fil --help` | same, for PulsarX |
| `lx run -v ~/data:/data pulsarx ls /data` | also mount another folder (repeatable) |
| `lx run --net-off presto ...` | without network access |
| `lx start`, `lx stop`, `lx status` | engine on, off, state (`stop` frees its memory) |
| `lx old`, `lx resume NAME` | old stopped containers, with the changes made in them |
| `lx build presto` | rebuild from its recipe (downloads packages) |
| `lx where` | images, sizes and the engine's disk file |

Globs such as `*.fft` are expanded on the Mac before Linux sees them, so use relative paths. `lx` never downloads an image by itself; only `lx build` downloads.

## Systems

| Name | Image | CPU type | Contents | Recipe |
|---|---|---|---|---|
| `presto` | `presto-native` | arm64 (native) | PRESTO 6.3.3 | [Dockerfile](images/presto-native/Dockerfile) |
| `pulsarx` | `ypmen/pulsarx` | amd64 (emulated) | PulsarX, TransientX | Docker Hub |
| `fits2fil` | `fits2fil-pulsarx` | amd64 (emulated) | PulsarX + fits2fil | not here; image only |
| `cuda` | `my-linux-env` | amd64 (emulated) | CUDA 12.0, `nvcc`, gcc | not here; image only |
| `singularity` | `quay.io/singularity/singularity:v3.11.5` | amd64 (emulated) | Singularity | quay.io |
| `ubuntu` | `ubuntu` | arm64 (native) | Ubuntu 24.04 | Docker Hub |

- Emulated systems are slower; do not trust their timings. They print a platform warning, which is harmless.
- `fits2fil` and `cuda` are labelled arm64 but hold amd64 programs, so `systems.conf` uses `auto` for them.
- `fits2fil` and `cuda` have no recipe here. If deleted they cannot be rebuilt from this repository. To back one up: `docker save my-linux-env | gzip > my-linux-env.tar.gz` (a few GB).
- There is no GPU on a Mac. `cuda` can compile but not run. To compile our library: `lx run cuda bash -c 'cd /work/filterbank/cuda && make ARCH=-arch=sm_75'` from the repository root (this writes a Linux amd64 library into `cuda/build/`).
- To add a system, add a line to [systems.conf](systems.conf): `name|image|platform|description`. A recipe goes in `images/<image name>/Dockerfile`.

## Where they are

All images and containers are inside one file, `~/Library/Containers/com.docker.docker/Data/vms/0/data/Docker.raw`. Do not move or edit it. Deleting it deletes every system and container.

## Old containers

| Container | From | Kept because |
|---|---|---|
| `my-running-container` | `cuda` | has `/astroflow` sources (75 MB of changes) |
| `gracious_villani` | `ubuntu` | 684 MB of installed packages and an ssh folder |
| `pulsarx_container` | `pulsarx` | mounts `~/sinimg` as `/data`, has shell history |

`fits2fil` contains all of `pulsarx`, but the layers are shared, so it costs little extra disk.

The `presto` recipe was rewritten from the one used to build the installed image (pinned commit, plus a `python` link). It has not been rebuilt since.
