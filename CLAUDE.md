# CLAUDE.md — optimal-control-pages

Interactive essay series on optimal control, published as Pluto notebooks to
GitHub Pages: https://jeremydwong.github.io/optimal-control-pages/
Companion research repo: `jeremydwong/blackmore-opts`.

## How publishing works

- Notebooks (`01-*.jl` … at repo root) are plain Pluto notebooks. On every
  push, `.github/workflows/ExportPluto.yaml` runs them with real solvers
  (PlutoSliderServer v1, Julia 1.12) and deploys static HTML to the
  `gh-pages` branch. GitHub Pages serves that branch.
- `index.html` is **hand-written** (PlutoSliderServer skips generating an
  index when one exists). When adding a page, add its card there, in order.
- **Do not gitignore `*.html`**: the Pages deploy action stages files with
  git, so gitignored exports silently never reach gh-pages (this bit us).
- Interactive demos are self-contained `@htl` JS widgets fed by
  Julia-computed parameter sweeps at export time — **not** `@bind`
  sliders, because PlutoSliderServer's precomputed-bind feature was never
  merged (PR #29); `@bind` would be frozen on the static page.
- Authoring gotcha: Julia's `md""` macro mis-parses multiple `$`
  interpolations mixed with `$...$` math (later ones render literally).
  Any markdown cell that interpolates computed values must use
  `Markdown.parse("...")` with ordinary string interpolation.
- Validate locally before pushing:
  `julia -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add(Pkg.PackageSpec(name="PlutoSliderServer", version="1")); import PlutoSliderServer; PlutoSliderServer.export_directory(".")'`
  then base64-decode `window.pluto_statefile` from the HTML and grep for
  error types and expected rendered text (beware: the notebook *source* is
  also embedded — check for resolved values, not source strings).

## Adding a new page (the common case — no Binder work)

1. Write `NN-page-name.jl` at the repo root, same pattern as existing
   notebooks (frontmatter with `order`/`title`, `Markdown.parse` for
   value-bearing cells, `@htl` widgets for interactivity).
2. Replace its dashed "planned" card in `index.html` with a real link.
3. Push. CI runs the notebook and publishes it.
4. Binder needs touching **only if** the page uses a Julia package outside
   the baked set (InfiniteOpt, Ipopt, JuMP, Clarabel, Plots,
   HypertextLiteral, JSON3) — then follow the runbook below.

## The Binder setup (the "Edit or run" button)

### Plain-English summary

When a reader clicks "Edit or run", Binder rents them a temporary cloud
computer built from a recipe (the `jeremydwong/pluto-on-binder` repo),
which Binder bakes once into a frozen disk snapshot and then copies for
every visitor. Julia translates packages into machine code the first time
they're used ("precompilation", 10+ min for our solvers on Binder's weak
machines) — so the recipe does that translation during the one-time bake
and stores the results in the snapshot. Two gotchas made the first attempt
useless: Julia only reuses stored translations if the Julia *settings* and
the *CPU type* at runtime match the ones used while baking — otherwise it
silently redoes everything live (the "precompiling forever" hang). And
"mybinder.org" is secretly three independent providers behind one address,
each with its own snapshot freezer, so a new recipe must be baked at all
three or some visitors wait ~45 min on a cold one. The details below
exist to keep those three things true.

Every exported page has a built-in button that runs the notebook on the
free mybinder.org service. **It points at a custom environment**:

- Repo: `jeremydwong/pluto-on-binder`, current tag **`opt-control-2`**
- Wired via `Export_binder_url` in `ExportPluto.yaml`:
  `https://mybinder.org/v2/gh/jeremydwong/pluto-on-binder/opt-control-2`

### Why a custom environment

With the stock `fonsp/pluto-on-binder` environment, every session had to
download and precompile InfiniteOpt/Ipopt/Plots/Clarabel from scratch on a
1-CPU Binder pod (10–25 min) — which users experience as an endless
"reconnecting…"/"Precompiling…" hang. The custom image pre-bakes all the
packages the notebooks use, so sessions go click → cells running in ~2–4
minutes.

### How the image is built (and the three traps that broke v1)

The repo is `fonsp/pluto-on-binder` v1.0.1 (Pluto 1.0.1, Julia 1.10.11)
plus a `postBuild` step that `Pkg.add`s + precompiles + **actually loads**
(`using …`) the notebook packages into the image's Julia depot. Three
hard-won rules, each of which independently caused the bake to be useless
in `opt-control-1`:

1. **Compiler flags must match end-to-end.** Julia precompile caches are
   keyed by compiler flags. Upstream's `pluto_server_config.jl` ran
   notebook workers with `pkgimages="no", optimize=1` (a tuning that makes
   *live* precompilation faster — the opposite of our strategy), while Pkg
   operations ran under the server's default flags. Any mismatch → caches
   silently rejected → full re-precompilation at runtime. Fix: the custom
   env removes those options so **bake == server == workers == default
   flags** (which also gives workers native code caches → faster solves).
2. **Pin a portable `JULIA_CPU_TARGET`.** Native caches default to the
   build machine's CPU; mybinder builds and runs on different hardware, so
   native-targeted caches can be rejected at runtime. Both `postBuild`
   (bake time) and the `start` script (runtime) export the portable
   multi-target string official Julia binaries use:
   `generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)`.
3. **End the bake with a real `using` of every package** — verifies
   loadability at build time and bakes any remaining load-time caches.

### mybinder is a federation — warm every member

mybinder.org load-balances across independent members, each with its
**own image registry** (as of 2026-06: `gesis`, `2i2c`, `bids`). An image
built on one member does not exist on the others; users landing on a cold
member sit through a ~45-minute from-scratch build ("Building Pluto.jl…").
After any retag, warm all members:

```sh
for m in gesis 2i2c bids; do
  curl -s -N -m 5400 -H "Accept: text/event-stream" \
    "https://$m.mybinder.org/build/gh/jeremydwong/pluto-on-binder/opt-control-2" \
    | grep -E '"phase": "(ready|failed)"' | head -1 &
done; wait
```

Wait for `"phase": "ready"` from each. (Member list can drift — probe
`https://<member>.mybinder.org/health`.)

### Runbook: notebooks gain a new package dependency

1. Add the package to the `Pkg.add([...])` list **and** the `using` line
   in `pluto-on-binder/postBuild`.
2. Commit, create a new tag (`opt-control-3`, …), push both.
3. Update `Export_binder_url` in `ExportPluto.yaml` to the new tag; push
   (CI re-exports all pages with the new button URL).
4. Warm all federation members (above) — one-time ~45 min each, parallel.
5. Sanity-check a page: `curl -s <page-url> | grep pluto_binder_url`.

### Known limits / fallback

- mybinder's free capacity fluctuates; launches take ~30–60 s on good
  days. The static page + JS widgets never depend on it.
- If Binder ever becomes unworkable: a small VPS running PlutoSliderServer
  (`run_git_directory`) serves these exact notebooks with live `@bind`
  sliders, no Binder involved.
