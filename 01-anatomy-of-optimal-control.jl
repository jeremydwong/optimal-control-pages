### A Pluto.jl notebook ###
# v0.20.0

#> [frontmatter]
#> order = "1"
#> title = "1. The Anatomy of an Optimal Control Problem"
#> description = "Landing a 1-D spacecraft with InfiniteOpt.jl: decision variables, dynamics, boundary conditions, path constraints, and an objective."

using Markdown
using InteractiveUtils

# ╔═╡ a0c0de01-0001-4001-8001-000000000001
md"""
# The Anatomy of an Optimal Control Problem

*Part 1 of a series on thinking in optimal control — from a falling point mass
to a 6-DoF rocket landing.*

Suppose you are falling toward the surface of Mars and you have one engine
pointed at the ground. When should you fire it, and how hard?

This is an **optimal control problem**, and this page dissects one in its
simplest possible form. The take-home is not the rocket — it is the **five
ingredients** that every such problem is built from:

1. **Decision variables** — the things you get to choose, here a thrust profile $a(t)$ *and* the trajectory it produces,
2. **Dynamics** — differential equations linking trajectory to choices,
3. **Boundary conditions** — where you start and where you must end,
4. **Path constraints** — what must stay true the whole way (don't be underground),
5. **An objective** — the single number that says which feasible choice is *best*.

We will write these five ingredients almost verbatim in
[InfiniteOpt.jl](https://github.com/infiniteopt/InfiniteOpt.jl), a Julia
package that lets you state the *continuous-time* problem and handles the
discretization for you.

> 🏃 **You can run all of this.** Press **Edit or run** (top right) and choose
> **Binder** to get a live Julia session for this exact notebook — edit the
> code, move the numbers, break things. (It is free, but takes a few minutes
> to boot.) The sliders further down work instantly, right on this page.
"""

# ╔═╡ a0c0de01-0002-4002-8002-000000000002
md"""
## The problem

A point mass at altitude $h_0$ is descending at speed $v_0$. Gravity pulls
down at $g$; the engine pushes up with acceleration $a(t)$, limited to
$0 \le a(t) \le a_{\max}$. We must reach the ground at rest at a fixed final
time $t_f$:

```math
\begin{aligned}
\dot h &= v, \qquad \dot v = a - g \\
h(0) &= h_0, \quad v(0) = v_0, \qquad h(t_f) = 0, \quad v(t_f) = 0 \\
h(t) &\ge 0, \qquad 0 \le a(t) \le a_{\max}
\end{aligned}
```

For the objective we minimize the **control effort** $\int_0^{t_f} a^2\,dt$
(a smooth stand-in for "don't work the engine hard"). Why not fuel,
$\int a\,dt$? Hold that thought — there is a small surprise below.
"""

# ╔═╡ a0c0de01-0003-4003-8003-000000000003
begin
    using InfiniteOpt, Ipopt, Plots, HypertextLiteral, JSON3
    const g_mars = 3.71     # m/s²
    const h₀ = 1000.0       # m
    const v₀ = -50.0        # m/s (descending)
    const a_max = 10.0      # m/s²
end

# ╔═╡ a0c0de01-0004-4004-8004-000000000004
md"""
## The five ingredients, in code

Below is the whole problem. Read it next to the math above — the
correspondence is nearly one-to-one. `t` is an *infinite parameter* (a
continuous variable the states are functions of); `∂(h, t)` is the
derivative $\dot h$; InfiniteOpt turns the continuous problem into a finite
one by enforcing everything on a grid of `num_supports` points
(*transcription* — we'll have much more to say about this later in the
series).
"""

# ╔═╡ a0c0de01-0005-4005-8005-000000000005
function solve_landing(tf; h0 = h₀, v0 = v₀, g = g_mars, amax = a_max,
                       n = 101)
    model = InfiniteModel(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)

    # time: the "infinite parameter" everything else is a function of
    @infinite_parameter(model, t ∈ [0, tf], num_supports = n)

    # (1) decision variables: the trajectory AND the control
    @variable(model, h, Infinite(t))                 # altitude [m]
    @variable(model, v, Infinite(t))                 # velocity [m/s]
    @variable(model, 0 <= a <= amax, Infinite(t))    # thrust accel [m/s²]

    # (2) dynamics
    @constraint(model, ∂(h, t) == v)
    @constraint(model, ∂(v, t) == a - g)

    # (3) boundary conditions
    @constraint(model, h(0) == h0)
    @constraint(model, v(0) == v0)
    @constraint(model, h(tf) == 0)
    @constraint(model, v(tf) == 0)

    # (4) path constraint: stay above ground the whole way
    @constraint(model, h >= 0)

    # (5) objective: minimize control effort
    @objective(model, Min, ∫(a^2, t))

    optimize!(model)

    ok = termination_status(model) in (OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED)
    ts = supports(t)
    av = value(a)
    (ok = ok, status = termination_status(model), t = ts,
     h = value(h), v = value(v), a = av,
     fuel = sum(av[2:end]) * tf / (n - 1))  # ∫a dt under this transcription
end

# ╔═╡ a0c0de01-0006-4006-8006-000000000006
md"""
That's it. No Riccati equations, no calculus of variations by hand — you
*declare* the problem and a solver (here [Ipopt](https://github.com/coin-or/Ipopt))
searches for the best feasible trajectory. Let's solve a comfortable case,
$t_f = 40$ s:
"""

# ╔═╡ a0c0de01-0007-4007-8007-000000000007
begin
    sol40 = solve_landing(40.0)
    sol40_status = string(sol40.status)
    sol40_fuel = round(sol40.fuel, digits = 1)
    md"Solver says: **$sol40_status** — fuel (Δv) delivered: **$sol40_fuel m/s**"
end

# ╔═╡ a0c0de01-0008-4008-8008-000000000008
begin
    p1 = plot(sol40.t, sol40.h, lw = 2.5, color = :goldenrod, legend = false,
        ylabel = "h [m]", title = "Effort-optimal landing, tf = 40 s")
    p2 = plot(sol40.t, sol40.v, lw = 2.5, color = :deepskyblue, legend = false,
        ylabel = "v [m/s]")
    p3 = plot(sol40.t, sol40.a, lw = 2.5, color = :orangered, legend = false,
        ylabel = "a [m/s²]", xlabel = "t [s]", ylim = (-0.5, a_max * 1.1))
    hline!(p3, [a_max], ls = :dash, color = :red, alpha = 0.6)
    plot(p1, p2, p3, layout = (3, 1), size = (680, 540))
end

# ╔═╡ a0c0de01-0009-4009-8009-000000000009
md"""
With plenty of time, the optimizer spreads the work out: a gentle,
smooth thrust profile, comfortably below the $a_{\max}$ line.

**Now squeeze the schedule.** Drag the slider: as $t_f$ shrinks, the same
five ingredients produce a qualitatively different answer — the thrust
saturates against its ceiling and the profile sharpens into *coast, then
brake hard*: the famous **suicide burn**. Nobody programmed that strategy;
it fell out of the constraints.
"""

# ╔═╡ a0c0de01-0010-4010-8010-000000000010
tf_sweep = let
    tfs = round.(collect(range(23.0, 60.0, length = 20)), digits = 1)
    runs = []
    for tf in tfs
        s = solve_landing(tf; n = 81)
        s.ok && push!(runs, (tf = tf, t = round.(s.t, digits = 2),
                           h = round.(s.h, digits = 1),
                           v = round.(s.v, digits = 2),
                           a = round.(s.a, digits = 3),
                           fuel = round(s.fuel, digits = 1),
                           apeak = round(maximum(s.a), digits = 2)))
    end
    runs
end

# ╔═╡ a0c0de01-0011-4011-8011-000000000011
@htl("""
<div style="background:#0e1117; border-radius:10px; padding:14px; color:#ddd; font-family:system-ui;">
  <div style="margin-bottom:8px;">
    <b style="color:goldenrod">final time tf</b>:
    <input type="range" min="0" max="$(length(tf_sweep)-1)" value="$(length(tf_sweep)-1)" style="width:55%; vertical-align:middle;">
    <span class="tfval" style="color:goldenrod; font-weight:bold;"></span>
    &nbsp;&nbsp; Δv delivered: <span class="fuelval" style="color:#7fdfff"></span>
    &nbsp;&nbsp; peak a: <span class="apeak" style="color:#ff9966"></span>
  </div>
  <canvas class="c_a" width="640" height="170" style="width:100%"></canvas>
  <canvas class="c_v" width="640" height="150" style="width:100%"></canvas>
  <canvas class="c_h" width="640" height="150" style="width:100%"></canvas>
  <script>
    const root = currentScript.parentElement;
    const DATA = $(HypertextLiteral.JavaScript(JSON3.write(tf_sweep)));
    const AMAX = $(a_max);
    const slider = root.querySelector("input");
    function chart(cv, xs, ys, opts) {
      const ctx = cv.getContext("2d");
      const W = cv.width, H = cv.height, L = 46, R = 10, T = 18, B = 22;
      ctx.clearRect(0, 0, W, H);
      const xmin = 0, xmax = opts.xmax, ymin = opts.ymin, ymax = opts.ymax;
      const px = x => L + (x - xmin) / (xmax - xmin) * (W - L - R);
      const py = y => H - B - (y - ymin) / (ymax - ymin) * (H - T - B);
      ctx.strokeStyle = "#444"; ctx.lineWidth = 1;
      ctx.strokeRect(L, T, W - L - R, H - T - B);
      ctx.fillStyle = "#999"; ctx.font = "11px system-ui";
      ctx.fillText(opts.label, L + 6, T - 5);
      for (let i = 0; i <= 4; i++) {
        const yv = ymin + (ymax - ymin) * i / 4;
        ctx.fillText(yv.toFixed(0), 4, py(yv) + 4);
        ctx.strokeStyle = "#2a2a2a";
        ctx.beginPath(); ctx.moveTo(L, py(yv)); ctx.lineTo(W - R, py(yv)); ctx.stroke();
      }
      if (opts.hline !== undefined) {
        ctx.strokeStyle = "#cc3333"; ctx.setLineDash([5, 4]);
        ctx.beginPath(); ctx.moveTo(L, py(opts.hline)); ctx.lineTo(W - R, py(opts.hline)); ctx.stroke();
        ctx.setLineDash([]);
      }
      ctx.strokeStyle = opts.color; ctx.lineWidth = 2.2;
      ctx.beginPath();
      for (let i = 0; i < xs.length; i++) {
        const X = px(xs[i]), Y = py(ys[i]);
        i === 0 ? ctx.moveTo(X, Y) : ctx.lineTo(X, Y);
      }
      ctx.stroke();
    }
    function draw() {
      const d = DATA[+slider.value];
      const xmax = DATA[DATA.length - 1].tf;
      root.querySelector(".tfval").textContent = d.tf.toFixed(1) + " s";
      root.querySelector(".fuelval").textContent = d.fuel.toFixed(0) + " m/s";
      root.querySelector(".apeak").textContent = d.apeak.toFixed(1) + " m/s2";
      chart(root.querySelector(".c_a"), d.t, d.a,
        {xmax: xmax, ymin: -0.5, ymax: AMAX * 1.15, color: "#ff7043", label: "thrust accel a(t)  [m/s2]", hline: AMAX});
      chart(root.querySelector(".c_v"), d.t, d.v,
        {xmax: xmax, ymin: -110, ymax: 10, color: "#4fc3f7", label: "velocity v(t)  [m/s]"});
      chart(root.querySelector(".c_h"), d.t, d.h,
        {xmax: xmax, ymin: 0, ymax: 1050, color: "#d4af37", label: "altitude h(t)  [m]"});
    }
    slider.addEventListener("input", draw);
    draw();
  </script>
</div>
""")

# ╔═╡ a0c0de01-0012-4012-8012-000000000012
md"""
## The surprise: an objective that can't tell solutions apart

Watch the **Δv delivered** readout as you move the slider — it is *exactly*
$|v_0| + g\,t_f$, no matter what the thrust profile looks like. In 1-D with
upward-only thrust, integrating $\dot v = a - g$ between the fixed endpoints
pins the fuel completely:

```math
\int_0^{t_f} a\,dt \;=\; v(t_f) - v(0) + g\,t_f \;=\; |v_0| + g\,t_f.
```

So if we had minimized *fuel* at fixed $t_f$, **every feasible trajectory
would have cost the same** — the optimizer would have shrugged and returned
an arbitrary one. This is the first habit of thinking in optimal control:

> **Ask whether your objective can actually distinguish your decisions.**
> If it can't, the solver's answer is an accident of the algorithm, not a
> design.

(The fix here was to minimize effort instead. In 3-D — next pages — thrust
direction matters, gravity losses depend on the path, and fuel becomes a
real, honest objective.)

Notice also what the formula *does* say: fuel grows linearly with $t_f$, so
the fuel-optimal landing is the **minimum-time** landing. Let's find it.
"""

# ╔═╡ a0c0de01-0013-4013-8013-000000000013
t_min = let
    # bisect on feasibility: the optimizer itself tells us what's possible
    lo, hi = 5.0, 60.0
    for _ in 1:12
        mid = (lo + hi) / 2
        solve_landing(mid; n = 61).ok ? (hi = mid) : (lo = mid)
    end
    round(hi, digits = 2)
end

# ╔═╡ a0c0de01-0014-4014-8014-000000000014
md"""
Minimum feasible time: **$(t_min) s**. Drag the slider down toward it and
watch the thrust pin itself to the ceiling — at $t_f = t_{\min}$ the
*constraints alone* dictate the trajectory (maximum braking the whole way
down, no choice left), which is exactly what "minimum time" means.

Notice *how* we found it: we asked the solver "is this $t_f$ feasible?"
twelve times and bisected. Wrapping an optimizer in an outer search is a
standard and honorable move — we will use it again on page 4 for
free-final-time landings.

## Take-home

Every optimal control problem in this series — up to and including a 6-DoF
rocket — is these same five ingredients:

| Ingredient | Here | In InfiniteOpt |
|---|---|---|
| Decision variables | $h, v, a$ as functions of $t$ | `@variable(model, h, Infinite(t))` |
| Dynamics | $\dot h = v,\ \dot v = a - g$ | `@constraint(model, ∂(h,t) == v)` |
| Boundary conditions | start falling, end at rest | `@constraint(model, h(0) == h0)` |
| Path constraints | $h \ge 0$, $0 \le a \le a_{\max}$ | bounds & `@constraint` |
| Objective | $\min \int a^2 dt$ | `@objective(model, Min, ∫(a^2, t))` |

The art is rarely in the solver. It is in *posing* the problem: an objective
that distinguishes solutions, constraints that capture reality, and — as
we'll see next — a formulation the solver can actually digest.

**Next: [Part 2 — Nonconvexity and the lossless trick](02-lossless-convexification.html)**,
where a real rocket engine (which cannot throttle below ~30%) punches a hole
in our feasible set, and a beautiful piece of 2007-era theory patches it.
"""

# ╔═╡ Cell order:
# ╟─a0c0de01-0001-4001-8001-000000000001
# ╟─a0c0de01-0002-4002-8002-000000000002
# ╠═a0c0de01-0003-4003-8003-000000000003
# ╟─a0c0de01-0004-4004-8004-000000000004
# ╠═a0c0de01-0005-4005-8005-000000000005
# ╟─a0c0de01-0006-4006-8006-000000000006
# ╠═a0c0de01-0007-4007-8007-000000000007
# ╠═a0c0de01-0008-4008-8008-000000000008
# ╟─a0c0de01-0009-4009-8009-000000000009
# ╟─a0c0de01-0010-4010-8010-000000000010
# ╟─a0c0de01-0011-4011-8011-000000000011
# ╟─a0c0de01-0012-4012-8012-000000000012
# ╠═a0c0de01-0013-4013-8013-000000000013
# ╟─a0c0de01-0014-4014-8014-000000000014
