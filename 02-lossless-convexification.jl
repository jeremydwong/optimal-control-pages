### A Pluto.jl notebook ###
# v0.20.0

#> [frontmatter]
#> order = "2"
#> title = "2. Nonconvexity and the Lossless Trick"
#> description = "A real engine cannot throttle to zero. The feasible set becomes a ring with a hole — and a beautiful relaxation repairs it, losslessly."

using Markdown
using InteractiveUtils

# ╔═╡ b0c0de02-0001-4001-8001-000000000001
md"""
# Nonconvexity and the Lossless Trick

*Part 2 — [back to Part 1](01-anatomy-of-optimal-control.html)*

On page 1 our engine could produce any acceleration from $0$ to $a_{\max}$.
Real rocket engines cannot. Once lit, a typical descent engine can throttle
down to perhaps 30% of full thrust — below that the combustion goes unstable.
And you do not relight a landing engine on the way down.

So the honest constraint on the thrust *vector* $\mathbf{u}$ is

```math
\rho_1 \;\le\; \|\mathbf{u}\|\;\le\; \rho_2,
```

with $\rho_1 > 0$. Geometrically the feasible set of thrust vectors is an
**annulus** — a disk with a hole punched in the middle. And that hole breaks
the single property that makes optimization reliable: **convexity**.

> 🏃 As before, everything here is runnable — **Edit or run** (top right) →
> **Binder**. Sliders work instantly on this page.
"""

# ╔═╡ b0c0de02-0002-4002-8002-000000000002
md"""
## Feel the hole

A set is **convex** if the straight line between any two of its points stays
inside the set. Drag the slider to move the second thrust choice around the
ring, and watch their average:
"""

# ╔═╡ b0c0de02-0003-4003-8003-000000000003
begin
    using JuMP, Clarabel, LinearAlgebra, Plots, HypertextLiteral, JSON3
    const ρ2_demo = 10.0
end

# ╔═╡ b0c0de02-0004-4004-8004-000000000004
@htl("""
<div style="background:#0e1117; border-radius:10px; padding:14px; color:#ddd; font-family:system-ui;">
  <div style="margin-bottom:6px;">
    <b style="color:goldenrod">throttle floor ρ₁</b>:
    <input class="s_rho" type="range" min="0" max="8" step="0.5" value="5" style="width:34%; vertical-align:middle;">
    <span class="rhoval" style="color:goldenrod; font-weight:bold;"></span>
    &nbsp;&nbsp;<b style="color:#7fdfff">second point</b>:
    <input class="s_th" type="range" min="0" max="360" step="2" value="120" style="width:34%; vertical-align:middle;">
    <span class="verdict" style="font-weight:bold; margin-left:10px;"></span>
  </div>
  <canvas width="640" height="360" style="width:100%"></canvas>
  <script>
    const root = currentScript.parentElement;
    const cv = root.querySelector("canvas"), ctx = cv.getContext("2d");
    const sR = root.querySelector(".s_rho"), sT = root.querySelector(".s_th");
    const R2 = 10.0;
    const cx = 320, cy = 185, scale = 15.5;
    function draw() {
      const r1 = +sR.value, th = +sT.value * Math.PI / 180;
      ctx.clearRect(0, 0, cv.width, cv.height);
      // annulus
      ctx.beginPath(); ctx.arc(cx, cy, R2 * scale, 0, 7);
      ctx.arc(cx, cy, r1 * scale, 0, 7, true);
      ctx.fillStyle = "rgba(212,175,55,0.25)"; ctx.fill("evenodd");
      ctx.strokeStyle = "#d4af37"; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(cx, cy, R2 * scale, 0, 7); ctx.stroke();
      if (r1 > 0) { ctx.beginPath(); ctx.arc(cx, cy, r1 * scale, 0, 7); ctx.stroke(); }
      // two feasible choices + their average
      const rA = (r1 + R2) / 2;
      const A = [cx + rA * scale, cy];
      const B = [cx + rA * scale * Math.cos(th), cy - rA * scale * Math.sin(th)];
      const M = [(A[0] + B[0]) / 2, (A[1] + B[1]) / 2];
      const mr = Math.hypot(M[0] - cx, M[1] - cy) / scale;
      const inside = (mr >= r1 - 1e-9) && (mr <= R2 + 1e-9);
      ctx.strokeStyle = "#888"; ctx.setLineDash([5, 4]);
      ctx.beginPath(); ctx.moveTo(A[0], A[1]); ctx.lineTo(B[0], B[1]); ctx.stroke();
      ctx.setLineDash([]);
      function dot(P, color) {
        ctx.beginPath(); ctx.arc(P[0], P[1], 6, 0, 7);
        ctx.fillStyle = color; ctx.fill();
      }
      dot(A, "#d4af37"); dot(B, "#4fc3f7");
      dot(M, inside ? "#66bb6a" : "#ef5350");
      ctx.fillStyle = "#999"; ctx.font = "12px system-ui";
      ctx.fillText("two feasible thrust choices (gold, blue) and their average", 14, 18);
      root.querySelector(".rhoval").textContent = r1.toFixed(1);
      const v = root.querySelector(".verdict");
      if (r1 === 0) { v.textContent = "no hole: any average is feasible (convex)"; v.style.color = "#66bb6a"; }
      else if (inside) { v.textContent = "this average is feasible"; v.style.color = "#66bb6a"; }
      else { v.textContent = "average INSIDE THE HOLE — infeasible!"; v.style.color = "#ef5350"; }
    }
    sR.addEventListener("input", draw); sT.addEventListener("input", draw);
    draw();
  </script>
</div>
""")

# ╔═╡ b0c0de02-0005-4005-8005-000000000005
md"""
With $\rho_1 = 0$ the set is a disk — convex, and optimizing over it is a
solved problem in the strongest sense: convex solvers find the *global*
optimum, fast, with a certificate. With $\rho_1 > 0$, two perfectly legal
thrust choices can average to something illegal. Algorithms that work by
blending and interpolating candidate solutions — which is what optimizers
*do* — keep getting pulled into the hole. You lose the guarantee, and in
practice you lose convergence too.

## The lift

Here is the trick, from Açıkmeşe & Ploen (2007), later flight-proven as the
heart of Mars pinpoint-landing guidance (Blackmore et al., 2013 — the paper
this whole series reimplements). Introduce one extra scalar variable
$\sigma$, and replace the annulus with:

```math
\|\mathbf{u}\| \le \sigma, \qquad \rho_1 \le \sigma \le \rho_2 .
```

Look at what happened. In the **lifted** space $(\mathbf{u}, \sigma)$, both
constraints are convex (a cone and a box — no holes anywhere). The hole has
not been removed; it has been *hidden in a projection*. The original set is
the slice $\|\mathbf{u}\| = \sigma$; the lifted set also contains points with
$\|\mathbf{u}\| < \sigma$ — thrust magnitudes the engine cannot actually fly.
We have **relaxed** the problem: every honest solution is still feasible,
but so are some impostors.
"""

# ╔═╡ b0c0de02-0006-4006-8006-000000000006
begin
    plt = plot(size = (640, 360), xlabel = "‖u‖", ylabel = "σ",
        xlim = (0, 11.5), ylim = (0, 11.5), legend = :topleft,
        title = "The lifted feasible set (cross-section)")
    ρ1_fig = 5.0
    xs = range(0, ρ2_demo, length = 100)
    plot!(plt, Shape(vcat([0.0, ρ1_fig], collect(range(ρ1_fig, ρ2_demo, length = 40)), [ρ2_demo, 0.0]),
                     vcat([ρ1_fig, ρ1_fig], collect(range(ρ1_fig, ρ2_demo, length = 40)), [ρ2_demo, ρ2_demo])),
        color = :goldenrod, alpha = 0.3, lw = 0, label = "relaxed set (convex)")
    plot!(plt, [ρ1_fig, ρ2_demo], [ρ1_fig, ρ2_demo], lw = 4, color = :orangered,
        label = "original set  ‖u‖ = σ ∈ [ρ₁, ρ₂]")
    plot!(plt, [0, 11.5], [0, 11.5], ls = :dash, color = :gray, label = "‖u‖ = σ")
    hline!(plt, [ρ1_fig, ρ2_demo], ls = :dot, color = :red, label = "ρ₁, ρ₂")
    plt
end

# ╔═╡ b0c0de02-0007-4007-8007-000000000007
md"""
So when is this not cheating? The famous result is the **losslessness
theorem** (Açıkmeşe & Ploen 2007; Blackmore et al. 2013): for the rocket's
minimum-fuel problem — *under hypotheses*, notably on the final time — the
optimal solution of the relaxed convex problem satisfies
$\|\mathbf{u}(t)\| = \sigma(t)$. The impostors are feasible but never
optimal, so solving the convex problem solves the nonconvex one.

A theorem is a contract, and its hypotheses are the fine print. The way to
respect fine print in numerical work is to **check the conclusion on every
solve**: compute $\max_t \big|\sigma - \|\mathbf{u}\|\big|$ and look at
it. It costs one line. Let's do exactly that on a 2-D lander with proper
fuel bookkeeping ($z = \ln m$, $\dot z = -\alpha\sigma$, thrust-per-mass
bounds $\rho\,e^{-z}$ — the page-3 formulation in miniature), discretizing
time directly into `N` steps — *transcription by hand* this time — and
handing the convex problem to
[Clarabel](https://github.com/oxfordcontrol/Clarabel.rs), a conic solver:
"""

# ╔═╡ b0c0de02-0008-4008-8008-000000000008
function solve_planar(ρ1; tf = 20.0, ρ2 = ρ2_demo, N = 60,
                      r0 = [-300.0, 400.0], v0 = [30.0, -20.0], g = 3.71,
                      α = 4.5e-4, scvx_iters = 5)
    # z = ln(m/m₀) with ż = -ασ; thrust-per-mass bounds ρ·e^{-z} are
    # linearized about z_ref and re-solved (successive convexification)
    dt = tf / N
    z_ref = zeros(N)
    out = nothing
    for _ in 1:scvx_iters
        model = Model(Clarabel.Optimizer); set_silent(model)

        @variable(model, r[1:N+1, 1:2])          # position (x, z)
        @variable(model, v[1:N+1, 1:2])          # velocity
        @variable(model, u[1:N, 1:2])            # thrust acceleration
        @variable(model, σ[1:N])                 # the lifted slack
        @variable(model, z[1:N+1])               # log-mass ratio

        for k in 1:N, i in 1:2
            gi = i == 2 ? -g : 0.0
            @constraint(model, r[k+1, i] == r[k, i] + dt * v[k, i] + 0.5dt^2 * (gi + u[k, i]))
            @constraint(model, v[k+1, i] == v[k, i] + dt * (gi + u[k, i]))
        end
        @constraint(model, [k = 1:N], z[k+1] == z[k] - α * dt * σ[k])
        @constraint(model, z[1] == 0)
        @constraint(model, [i = 1:2], r[1, i] == r0[i])
        @constraint(model, [i = 1:2], v[1, i] == v0[i])
        @constraint(model, [i = 1:2], r[N+1, i] == 0)
        @constraint(model, [i = 1:2], v[N+1, i] == 0)
        @constraint(model, [k = 1:N+1], r[k, 2] >= 0)     # stay above ground

        for k in 1:N
            @constraint(model, [σ[k]; u[k, :]] in SecondOrderCone())  # ‖u‖ ≤ σ
            z0 = z_ref[k]; ez = exp(-z0)
            @constraint(model, σ[k] >= ρ1 * ez * (1 - (z[k] - z0) + 0.5(z[k] - z0)^2))
            @constraint(model, σ[k] <= ρ2 * ez * (1 - (z[k] - z0)))
        end

        @objective(model, Min, dt * sum(σ))      # fuel ∝ ∫σ dt
        optimize!(model)
        termination_status(model) in (OPTIMAL, ALMOST_OPTIMAL) || return nothing

        zv = value.(z); z_ref = zv[1:N]
        un = [norm(value.(u[k, :])) for k in 1:N]
        σv = value.(σ)
        out = (t = collect(range(0, tf - dt, length = N)), r = value.(r),
               unorm = un, σ = σv, gap = maximum(abs.(σv .- un)),
               fuel = dt * sum(σv))
    end
    out
end

# ╔═╡ b0c0de02-0009-4009-8009-000000000009
begin
    snug = solve_planar(5.0; tf = 20.0)
    roomy = solve_planar(5.0; tf = 30.0)
    snug_gap = round(snug.gap, sigdigits = 2)
    roomy_gap = round(roomy.gap, sigdigits = 3)
    Markdown.parse("""
    Same lander, same throttle floor ρ₁ = 5, two flight times:

    | flight time | max gap, σ vs ‖u‖ | verdict |
    |---|---|---|
    | tf = 20 s (snug — ~25% above minimum) | **$snug_gap** | ✅ lossless |
    | tf = 30 s (roomy — ~90% above minimum) | **$roomy_gap** | ❌ **impostor** |

    The check passes when time is snug and **fires when time is roomy**.
    Why? With lots of spare time, the optimizer wants long stretches of
    *less-than-minimum* thrust. On those arcs σ is pinned at its floor
    while the actual ‖**u**‖ the dynamics need sits strictly below it —
    and since the objective only sees σ, every u beneath the floor costs
    exactly the same. Many optima, most of them impostors; the
    interior-point solver hands you a centered one. (Page 1's lesson in a
    new costume: an objective that cannot tell the candidates apart.)

    This is precisely the regime the theorem's fine print excludes — its
    hypotheses bound the final time. The structure of the failure even
    tells you the remedy: don't hand the problem a wasteful tf. Choosing
    tf *well* is page 4's whole story.
    """)
end

# ╔═╡ b0c0de02-0010-4010-8010-000000000010
rho_sweep = let
    tfs = [18.0, 20.0, 24.0, 28.0, 32.0]
    rhos = collect(0.0:1.0:8.0)
    grid = []
    for tf in tfs
        row = []
        for ρ1 in rhos
            s = solve_planar(ρ1; tf = tf)
            push!(row, s === nothing ? nothing :
                (t = round.(s.t, digits = 2),
                 x = round.(s.r[:, 1], digits = 1), z = round.(s.r[:, 2], digits = 1),
                 unorm = round.(s.unorm, digits = 3), sigma = round.(s.σ, digits = 3),
                 gap = round(s.gap, sigdigits = 2), fuel = round(s.fuel, digits = 1)))
        end
        push!(grid, row)
    end
    (tfs = tfs, rhos = rhos, grid = grid)
end

# ╔═╡ b0c0de02-0011-4011-8011-000000000011
@htl("""
<div style="background:#0e1117; border-radius:10px; padding:14px; color:#ddd; font-family:system-ui;">
  <div style="margin-bottom:8px;">
    <b style="color:goldenrod">throttle floor ρ₁</b>:
    <input class="s_r" type="range" min="0" max="$(length(rho_sweep.rhos)-1)" value="5" style="width:30%; vertical-align:middle;">
    <span class="rv" style="color:goldenrod; font-weight:bold;"></span>
    &nbsp;&nbsp;<b style="color:#7fdfff">flight time tf</b>:
    <input class="s_t" type="range" min="0" max="$(length(rho_sweep.tfs)-1)" value="1" style="width:26%; vertical-align:middle;">
    <span class="tv" style="color:#7fdfff; font-weight:bold;"></span>
    &nbsp;&nbsp; Δv: <span class="fv" style="color:#aaa"></span>
    &nbsp;&nbsp; max gap: <span class="gv" style="font-weight:bold;"></span>
  </div>
  <canvas class="c_tr" width="640" height="220" style="width:100%"></canvas>
  <canvas class="c_u"  width="640" height="190" style="width:100%"></canvas>
  <script>
    const root = currentScript.parentElement;
    const SW = $(HypertextLiteral.JavaScript(JSON3.write(rho_sweep)));
    const R2 = $(ρ2_demo);
    const sr = root.querySelector(".s_r"), st = root.querySelector(".s_t");
    function axes(ctx, W, H, L, R, T, B, xmin, xmax, ymin, ymax, label) {
      ctx.clearRect(0, 0, W, H);
      ctx.strokeStyle = "#444"; ctx.lineWidth = 1; ctx.strokeRect(L, T, W - L - R, H - T - B);
      ctx.fillStyle = "#999"; ctx.font = "11px system-ui"; ctx.fillText(label, L + 6, T - 5);
      const px = x => L + (x - xmin) / (xmax - xmin) * (W - L - R);
      const py = y => H - B - (y - ymin) / (ymax - ymin) * (H - T - B);
      for (let i = 0; i <= 4; i++) {
        const yv = ymin + (ymax - ymin) * i / 4;
        ctx.fillText(yv.toFixed(0), 4, py(yv) + 4);
        ctx.strokeStyle = "#2a2a2a";
        ctx.beginPath(); ctx.moveTo(L, py(yv)); ctx.lineTo(W - R, py(yv)); ctx.stroke();
      }
      return [px, py];
    }
    function polyline(ctx, px, py, xs, ys, color, width) {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.beginPath();
      for (let i = 0; i < xs.length; i++) {
        i === 0 ? ctx.moveTo(px(xs[i]), py(ys[i])) : ctx.lineTo(px(xs[i]), py(ys[i]));
      }
      ctx.stroke();
    }
    function draw() {
      const tfv = SW.tfs[+st.value], rho = SW.rhos[+sr.value];
      const d = SW.grid[+st.value][+sr.value];
      root.querySelector(".rv").textContent = rho.toFixed(0);
      root.querySelector(".tv").textContent = tfv.toFixed(0) + " s";
      const gv = root.querySelector(".gv"), fv = root.querySelector(".fv");
      const c1 = root.querySelector(".c_tr"), c2 = root.querySelector(".c_u");
      if (!d) {
        fv.textContent = "—";
        gv.textContent = "INFEASIBLE (not enough time)"; gv.style.color = "#ef9a3c";
        c1.getContext("2d").clearRect(0, 0, c1.width, c1.height);
        c2.getContext("2d").clearRect(0, 0, c2.width, c2.height);
        return;
      }
      fv.textContent = d.fuel.toFixed(1) + " m/s";
      gv.textContent = d.gap.toExponential(1) + (d.gap < 1e-3 ? "  ✓ lossless" : "  ✗ impostor!");
      gv.style.color = d.gap < 1e-3 ? "#9ccc65" : "#ef5350";
      const x1 = c1.getContext("2d");
      const [pxa, pya] = axes(x1, c1.width, c1.height, 46, 10, 18, 22, -320, 40, 0, 420, "trajectory  z vs x  [m]");
      polyline(x1, pxa, pya, d.x, d.z, "#d4af37", 2.2);
      x1.fillStyle = "#ef5350"; x1.beginPath();
      x1.arc(pxa(0), pya(0), 5, 0, 7); x1.fill();
      const x2 = c2.getContext("2d");
      const [pxb, pyb] = axes(x2, c2.width, c2.height, 46, 10, 18, 22, 0, SW.tfs[SW.tfs.length-1], 0, R2 * 1.3,
        "thrust:  ‖u(t)‖ (gold)  σ(t) (blue dashed)  floor/ceiling (red)");
      x2.strokeStyle = "#cc3333"; x2.setLineDash([5, 4]);
      [rho, R2].forEach(function(b) {
        x2.beginPath(); x2.moveTo(pxb(0), pyb(b)); x2.lineTo(pxb(tfv), pyb(b)); x2.stroke();
      });
      x2.setLineDash([]);
      polyline(x2, pxb, pyb, d.t, d.unorm, "#d4af37", 2.6);
      x2.setLineDash([6, 5]);
      polyline(x2, pxb, pyb, d.t, d.sigma, "#4fc3f7", 1.8);
      x2.setLineDash([]);
    }
    sr.addEventListener("input", draw);
    st.addEventListener("input", draw);
    draw();
  </script>
</div>
""")

# ╔═╡ b0c0de02-0012-4012-8012-000000000012
md"""
Things to try:

- **Find the failure boundary.** At $t_f = 20$ s the check passes for every
  feasible floor; push $t_f$ to 28–32 s with $\rho_1 \ge 4$ and watch the
  gold $\|\mathbf{u}\|$ sag below the blue dashed $\sigma$ on the
  min-thrust arcs — the readout flips to **✗ impostor**. The theorem's
  final-time fine print, made visible.
- Where the check passes, $\|\mathbf{u}(t)\|$ rides the **bounds** — full
  thrust or minimum thrust with brief transitions. This *bang-bang*
  structure is the maximum principle showing through, just as the
  saturation on page 1 was.
- Drop $t_f$ to 18 s with a high floor and you hit **infeasible** — the
  three-way frontier between impossible, optimal, and wasteful is all
  visible in one slider.
- Fuel **rises** with $\rho_1$: a pickier engine is a costlier engine.

## When the check "fails": a diagnostics war story

This entire series exists because of a measurement like that gap check going
wrong. An earlier study of this exact problem class concluded that a whole
transcription method (orthogonal collocation, 17 configurations swept)
"fundamentally cannot" achieve a tight gap — best case $2.61$, which is
$10^9$ times worse than what you just saw.

The number $2.61$ should have rung an alarm. It equals $\rho_1/m_{\text{wet}}$
*exactly* — the lower bound of $\sigma$ at $t = 0$. Printing the gap **per
time point** instead of as a max revealed everything: at every interior
point the gap was $\sim 10^{-9}$; the entire failure was *one* point,
$t = 0$, where (under that discretization) the control variable appears in
no dynamics equation at all. The optimizer had parked a meaningless free
variable at zero, and a max-norm metric reported it as a method failure.
One extra constraint tying that endpoint to its neighbor, and the "broken"
method beat the working one by three orders of magnitude.

> **When an optimality check fails, look at *where* it fails before
> concluding *why*.** A failure metric that exactly equals a problem
> constant is a confession, not a coincidence.

## Take-home

- Convexity is the property that buys you global optima and certificates;
  the throttle floor destroys it.
- A **lift** ($\sigma$) can hide nonconvexity in a projection; a
  **relaxation** is honest only if the impostors it admits are never
  *optimal*. That is a theorem **with hypotheses** — and the hypotheses
  are about *your problem data* (here, the final time), not just the
  formulas.
- **Always verify** — the tightness check costs one line, it genuinely
  fires when a hypothesis is violated, and *where* it fires tells you the
  story (impostors live on the min-thrust arcs that a wasteful $t_f$
  creates).

**Next: [Part 3 — Landing on Mars in 3-DoF](03-mars-landing.html)**, where
we add the third dimension, glide-slope safety cones, and — the real
complication — the rocket gets lighter as it burns.
"""

# ╔═╡ Cell order:
# ╟─b0c0de02-0001-4001-8001-000000000001
# ╟─b0c0de02-0002-4002-8002-000000000002
# ╠═b0c0de02-0003-4003-8003-000000000003
# ╟─b0c0de02-0004-4004-8004-000000000004
# ╟─b0c0de02-0005-4005-8005-000000000005
# ╟─b0c0de02-0006-4006-8006-000000000006
# ╟─b0c0de02-0007-4007-8007-000000000007
# ╠═b0c0de02-0008-4008-8008-000000000008
# ╟─b0c0de02-0009-4009-8009-000000000009
# ╟─b0c0de02-0010-4010-8010-000000000010
# ╟─b0c0de02-0011-4011-8011-000000000011
# ╟─b0c0de02-0012-4012-8012-000000000012
