# ROP-Explorer

ROP-Explorer is an interactive tool for **Reaction Order Polyhedra (ROP)** analysis of equilibrium binding networks. It provides a browser UI for constructing binding networks, enumerating structural regimes, visualizing regime graphs, and exploring SISO paths and polyhedral geometry.

![Demo](webapp/demo2.gif)

## Repository Layout

```text
ROP-Explorer/
├── Bnc_julia/    # Local copy of BindingAndCatalysis.jl
├── webapp/       # Julia HTTP backend and frontend assets
├── packaging/    # PackageCompiler build scripts for standalone backend bundles
├── desktop/      # Electron shell for the packaged backend
├── deploy/       # Docker + Nginx deployment files
├── LICENSE
└── README.md
```

## Requirements

- Julia 1.10 or newer
- A modern browser
- Node.js + npm for the Electron desktop shell
- Docker + Docker Compose for server deployment

## Local Development

Clone the repository:

```bash
git clone https://github.com/YanzhangJiang/ROP-Explorer.git
cd ROP-Explorer
```

Install Julia dependencies from the repository root:

```bash
julia --project=webapp -e 'using Pkg; Pkg.develop(path="Bnc_julia"); Pkg.instantiate(); Pkg.precompile()'
```

Start the web app locally:

```bash
cd webapp
./start.sh
```

The server listens on `http://127.0.0.1:8088` by default. To use another port:

```bash
cd webapp
ROP_PORT=8090 julia -t auto --project=. server.jl
```

## Desktop Development

First build the standalone backend bundle:

```bash
julia --project=packaging packaging/build_backend_app.jl
```

This generates:

```text
dist/ROPExplorerBackend/
dist/ROPExplorerBackend/bin/rop-explorer-backend
```

Then start the Electron shell:

```bash
cd desktop
npm install
npm run dev
```

The Electron app launches the bundled backend locally and opens a native window pointed at the local web UI.

## Server Deployment

The `deploy/` directory contains a source-based Docker deployment:

- `deploy/Dockerfile` builds the Julia backend image
- `deploy/docker-compose.yml` runs the Julia app behind Nginx
- `deploy/nginx.conf` serves static assets and proxies API traffic

Build and start the server:

```bash
cd deploy
docker compose build
docker compose up -d
```

The backend runs on port `8088` inside the container, and Nginx exposes the service on ports `80` and `443`.

## Build Release Artifacts

### 1. Standalone Backend Bundle

Build the relocatable backend bundle:

```bash
julia --project=packaging packaging/build_backend_app.jl
```

Output:

```text
dist/ROPExplorerBackend/
```

This bundle contains the Julia runtime, the compiled backend executable, and the frontend assets.

### 2. macOS Desktop App

Build the macOS desktop app from the repository root:

```bash
julia --project=packaging packaging/build_backend_app.jl
cd desktop
npm install
npm run build:mac-app
```

The generated app is placed under:

```text
desktop/release/mac-arm64/ROP-Explorer.app
```

The current Electron packaging script is configured for macOS Apple Silicon (`arm64`).

## What To Upload In Releases

For GitHub Releases, upload only end-user deliverables:

- `desktop/release/mac-arm64/ROP-Explorer.app` as a zipped macOS desktop release
- `dist/ROPExplorerBackend/` as a `.tar.gz` or `.zip` only if you want to publish a backend-only bundle

Do not commit generated build directories such as `desktop/node_modules/`, `desktop/release/`, `desktop/release-packager/`, or `dist/` into the source repository. Build locally, then archive the final `.app` or `ROPExplorerBackend` bundle and upload that archive as the release asset.

## How To Use A Published Release

### Desktop app release

1. Download the packaged macOS release archive.
2. Unzip it to get `ROP-Explorer.app`.
3. Launch the app. It starts its local backend automatically and opens the UI.

### Backend-only release

1. Extract the backend bundle archive.
2. Start the executable:

```bash
cd ROPExplorerBackend
ROP_PORT=8088 ./bin/rop-explorer-backend
```

3. Open `http://127.0.0.1:8088` in a browser.

## API Overview

The backend exposes JSON APIs under `/api/`, including:

- `POST /api/build_model`
- `POST /api/find_vertices`
- `POST /api/build_graph`
- `POST /api/siso_paths`
- `POST /api/siso_polyhedra`
- `POST /api/siso_trajectory`
- `POST /api/rop_cloud`
- `POST /api/vertex_detail`
- `POST /api/fret_heatmap`

Sessions expire after one hour of inactivity.

## Acknowledgment

The core computational engine is [BindingAndCatalysis.jl](https://github.com/Qinguo25/BindingAndCatalysis.jl) by Qinguo Liu. The `Bnc_julia/` directory in this repository is a local copy of that package. Credit for the underlying ROP theory implementation belongs to the original authors.

## License

This repository is released under the MIT License. The vendored `BindingAndCatalysis.jl` copy in `Bnc_julia/` remains separately licensed under MIT by its original author.
