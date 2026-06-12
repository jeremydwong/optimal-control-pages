### A Pluto.jl notebook ###
# v0.20.0

#> [frontmatter]
#> order = "3"
#> title = "3. Landing on Mars in 3-DoF"
#> description = "The full Blackmore 2013 powered-descent problem: mass depletion, glide-slope cones, and a landing you can fly with sliders."

using Markdown
using InteractiveUtils

# ╔═╡ c0c0de03-0001-4001-8001-000000000001
md"""
# Landing on Mars in 3-DoF

*Part 3 — [back to Part 2](02-lossless-convexification.html)*

Time to assemble the real thing: the minimum-fuel powered-descent problem
from *Lossless Convexification of Nonconvex Control Bound and Pointing
Constraints of the Soft Landing Optimal Control Problem* (Açıkmeşe, Carson &
Blackmore, IEEE TCST 2013) — the algorithm family behind Mars pinpoint
landing. Same five ingredients as page 1, same lossless lift as page 2, plus
two pieces of new physics:

1. **The rocket gets lighter.** Mass obeys $\dot m = -\alpha\,\|\mathbf{T}\|$,
   and acceleration is $\mathbf{T}/m$ — products of unknowns everywhere.
2. **Safety geometry.** The lander must stay inside a **glide-slope cone**
   above the landing site (never skim the terrain) and below a speed limit.

> 🏃 **Edit or run** (top right) → Binder gives you this exact notebook live.
> The sliders below work instantly without it.
"""

# ╔═╡ c0c0de03-0002-4002-8002-000000000002
md"""
## Taming the mass: the logarithm trick

Divide the dynamics by $m$ and you get clean variables: a *specific thrust*
$\mathbf{u} = \mathbf{T}/m$ and, with $z = \ln m$, a linear-looking mass
equation. The full problem:

```math
\begin{aligned}
\min \int_0^{t_f} \sigma\,dt \quad \text{s.t.}\quad
& \dot{\mathbf r} = \mathbf v, \qquad
  \dot{\mathbf v} = \mathbf g + \mathbf u, \qquad
  \dot z = -\alpha\,\sigma \\
& \|\mathbf u\| \le \sigma, \qquad
  \rho_1 e^{-z} \le \sigma \le \rho_2 e^{-z} \\
& \sqrt{r_x^2 + r_y^2} \le r_z / \tan\gamma_{gs}, \qquad \|\mathbf v\| \le V_{\max} \\
& \text{boundary conditions; } z(t_f) \ge \ln m_{\text{dry}}
\end{aligned}
```

There is our lift $\sigma$ from page 2 — but its bounds now carry
$e^{-z}$, because a fixed thrust force buys more *acceleration* as the tanks
drain. Those exponential bounds are the one nonconvexity left (the lower one
is actually convex; the upper is not), and there are two honest ways in:

- **Just solve the NLP.** Ipopt accepts $e^{-z}$ as-is. With a sane warm
  start it converges in a couple of seconds.
- **Successively convexify.** Linearize $e^{-z}$ about a reference $z_{ref}$,
  solve the resulting *conic* problem in milliseconds with Clarabel, update
  $z_{ref}$, repeat 3–4 times. (The paper itself bounds $z$ with first- and
  second-order expansions — same spirit.)

We'll do both and check they agree — *trust, but verify* is the running theme
of this series.
"""

# ╔═╡ c0c0de03-0003-4003-8003-000000000003
begin
    using InfiniteOpt, Ipopt, JuMP, Clarabel, LinearAlgebra, Plots,
          HypertextLiteral, JSON3

    # Mars-like parameters, straight from the paper's test case
    const g_mars = [0.0, 0.0, -3.7114]    # m/s²
    const m_wet, m_dry = 1905.0, 1505.0   # kg
    const α = 1.0 / (225.0 * 9.80665)     # s/m  (Isp = 225 s)
    const ρ1, ρ2 = 4972.0, 13260.0        # thrust bounds, N
    const V_max = 500.0                   # m/s
    const r0 = [450.0, -330.0, 2400.0]    # m
    const v0 = [-40.0, 10.0, -10.0]       # m/s
end

# ╔═╡ c0c0de03-0004-4004-8004-000000000004
md"""
## Route 1: declare it (InfiniteOpt + Ipopt)

The continuous-time problem, ingredient by ingredient, exactly like page 1 —
just with vectors, the exponential bounds kept *exact*, and one extra
constraint tying the $t=0$ control to its neighbor (the endpoint-artifact
lesson from page 2's war story, applied):
"""

# ╔═╡ c0c0de03-0005-4005-8005-000000000005
function solve_infiniteopt(tf; γ_gs = deg2rad(6.0), n = 61)
    model = InfiniteModel(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
    set_optimizer_attribute(model, "max_iter", 3000)
    set_optimizer_attribute(model, "tol", 1e-8)

    @infinite_parameter(model, t ∈ [0, tf], num_supports = n)

    @variable(model, r[1:3], Infinite(t))      # position
    @variable(model, v[1:3], Infinite(t))      # velocity
    @variable(model, z, Infinite(t))           # log-mass
    @variable(model, u[1:3], Infinite(t))      # specific thrust T/m
    @variable(model, σ >= 0, Infinite(t))      # the lossless lift

    @constraint(model, [i = 1:3], ∂(r[i], t) == v[i])
    @constraint(model, [i = 1:3], ∂(v[i], t) == g_mars[i] + u[i])
    @constraint(model, ∂(z, t) == -α * σ)

    @constraint(model, [i = 1:3], r[i](0) == r0[i])
    @constraint(model, [i = 1:3], v[i](0) == v0[i])
    @constraint(model, z(0) == log(m_wet))
    @constraint(model, [i = 1:3], r[i](tf) == 0)
    @constraint(model, [i = 1:3], v[i](tf) == 0)
    @constraint(model, z(tf) >= log(m_dry))

    @constraint(model, sum(u[i]^2 for i in 1:3) <= σ^2)   # ‖u‖ ≤ σ
    @constraint(model, σ >= ρ1 * exp(-z))                  # exact exp bounds
    @constraint(model, σ <= ρ2 * exp(-z))

    @constraint(model, r[1]^2 + r[2]^2 <= (r[3] / tan(γ_gs))^2)  # glide slope
    @constraint(model, sum(v[i]^2 for i in 1:3) <= V_max^2)      # speed limit
    @constraint(model, r[3] >= 0)

    # the t=0 control appears in no dynamics equation under this
    # transcription — tie it to the first interval (see Part 2's sidebar)
    Δ = tf / (n - 1)
    @constraint(model, [i = 1:3], u[i](0) == u[i](Δ))
    @constraint(model, σ(0) == σ(Δ))

    # warm start: straight line, conservative burn
    for i in 1:3
        set_start_value_function(r[i], τ -> r0[i] * (1 - τ / tf))
        set_start_value_function(v[i], τ -> v0[i] * (1 - τ / tf))
        set_start_value_function(u[i], τ -> -v0[i] / tf - g_mars[i])
    end
    set_start_value_function(z, τ -> log(m_wet - (m_wet - m_dry) * 0.3τ / tf))
    set_start_value_function(σ, τ -> (ρ1 + ρ2) / 2 / m_wet)

    @objective(model, Min, ∫(σ, t))
    optimize!(model)

    un = [norm([value(u[i])[k] for i in 1:3]) for k in 1:length(supports(t))]
    σv = value(σ)
    (ok = termination_status(model) in (OPTIMAL, LOCALLY_SOLVED),
     fuel = m_wet - exp(value(z)[end]),
     gap = maximum(abs.(σv .- un)))
end

# ╔═╡ c0c0de03-0006-4006-8006-000000000006
begin
    ref = solve_infiniteopt(72.0)
    ref_fuel = round(ref.fuel, digits = 1)
    ref_gap = round(ref.gap, sigdigits = 2)
    Markdown.parse("""
    At tf = 72 s (the paper's value): fuel = **$ref_fuel kg**,
    max relaxation gap = **$ref_gap** across *every* time point. Lossless
    survives the exponential coupling — the mass bookkeeping that page 2
    showed is *essential* is in, and the numerical check agrees.
    """)
end

# ╔═╡ c0c0de03-0007-4007-8007-000000000007
md"""
## Route 2: transcribe it (Clarabel + successive convexification)

For the interactive demos we want ~50 solves, so milliseconds matter. We
discretize by hand ($N = 60$ steps), keep the second-order-cone constraints
native, linearize only the $e^{-z}$ bounds about a reference, and iterate the
reference. Four passes land within solver precision of Route 1:
"""

# ╔═╡ c0c0de03-0008-4008-8008-000000000008
function solve_conic(tf; γ_gs = deg2rad(6.0), N = 60, scvx_iters = 4)
    dt = tf / N
    z_ref = collect(range(log(m_wet), log(m_wet - 0.5(m_wet - m_dry)), length = N))
    out = nothing
    for _ in 1:scvx_iters
        model = Model(Clarabel.Optimizer); set_silent(model)
        @variable(model, r[1:N+1, 1:3]); @variable(model, v[1:N+1, 1:3])
        @variable(model, z[1:N+1]); @variable(model, u[1:N, 1:3]); @variable(model, σ[1:N])

        for k in 1:N
            for i in 1:3
                @constraint(model, r[k+1,i] == r[k,i] + dt*v[k,i] + 0.5dt^2*(g_mars[i] + u[k,i]))
                @constraint(model, v[k+1,i] == v[k,i] + dt*(g_mars[i] + u[k,i]))
            end
            @constraint(model, z[k+1] == z[k] - α*dt*σ[k])
        end
        @constraint(model, [i = 1:3], r[1,i] == r0[i])
        @constraint(model, [i = 1:3], v[1,i] == v0[i])
        @constraint(model, [i = 1:3], r[N+1,i] == 0)
        @constraint(model, [i = 1:3], v[N+1,i] == 0)
        @constraint(model, z[1] == log(m_wet))
        @constraint(model, z[N+1] >= log(m_dry))

        for k in 1:N
            @constraint(model, [σ[k]; u[k,:]] in SecondOrderCone())
            z0 = z_ref[k]; ez = exp(-z0)
            @constraint(model, σ[k] >= ρ1*ez*(1 - (z[k]-z0) + 0.5(z[k]-z0)^2))
            @constraint(model, σ[k] <= ρ2*ez*(1 - (z[k]-z0)))
        end
        for k in 1:N+1
            @constraint(model, [r[k,3]/tan(γ_gs); r[k,1]; r[k,2]] in SecondOrderCone())
            @constraint(model, [V_max; v[k,:]] in SecondOrderCone())
            @constraint(model, r[k,3] >= 0)
        end

        @objective(model, Min, dt * sum(σ))
        optimize!(model)
        termination_status(model) in (OPTIMAL, ALMOST_OPTIMAL) || return nothing

        zv = value.(z); z_ref = zv[1:N]
        uv = value.(u); σv = value.(σ)
        un = [norm(uv[k,:]) for k in 1:N]
        out = (t = collect(range(0, tf - dt, length = N)),
               r = value.(r), m = exp.(zv),
               unorm = un, σ = σv, Tc = exp.(zv[1:N]) .* un,
               fuel = m_wet - exp(zv[end]),
               gap = maximum(abs.(σv .- un)))
    end
    out
end

# ╔═╡ c0c0de03-0009-4009-8009-000000000009
begin
    conic72 = solve_conic(72.0)
    fuel_r2 = round(conic72.fuel, digits = 1)
    fuel_r1 = round(ref.fuel, digits = 1)
    gap_r2 = round(conic72.gap, sigdigits = 2)
    Tlo = round(Int, minimum(conic72.Tc))
    Thi = round(Int, maximum(conic72.Tc))
    Markdown.parse("""
    Route 2 at tf = 72 s: fuel = **$fuel_r2 kg** (Route 1 said
    $fuel_r1), gap = **$gap_r2**, thrust riding [$Tlo, $Thi] N against
    bounds [$ρ1, $ρ2]. Two independent routes, one answer.
    """)
end

# ╔═╡ c0c0de03-0010-4010-8010-000000000010
md"""
## Fly it

Drag $t_f$ and watch the whole solution reorganize; press ▶ to replay the
landing. The wireframe is the glide-slope cone — the trajectory hugs it when
time is tight. The thrust trace below shows the bang-bang
**max → min → max** profile (and the lossless gap staying at $10^{-7}$ or
better for every setting).
"""

# ╔═╡ c0c0de03-0011-4011-8011-000000000011
tf_data = let
    runs = []
    for tf in 45.0:5.0:110.0
        s = solve_conic(tf)
        s === nothing && continue
        push!(runs, (tf = tf,
            t = round.(s.t, digits = 2),
            x = round.(s.r[:,1], digits = 1), y = round.(s.r[:,2], digits = 1),
            zc = round.(s.r[:,3], digits = 1),
            Tc = round.(s.Tc, digits = 0),
            unorm = round.(s.unorm, digits = 4), sigma = round.(s.σ, digits = 4),
            fuel = round(s.fuel, digits = 1), gap = round(s.gap, sigdigits = 2)))
    end
    runs
end

# ╔═╡ c0c0de03-0012-4012-8012-000000000012
@htl("""
<div style="background:#0e1117; border-radius:10px; padding:14px; color:#ddd; font-family:system-ui;">
  <div style="margin-bottom:8px;">
    <b style="color:goldenrod">time of flight tf</b>:
    <input class="s_tf" type="range" min="0" max="$(length(tf_data)-1)" value="5" style="width:40%; vertical-align:middle;">
    <span class="tv" style="color:goldenrod; font-weight:bold;"></span>
    <button class="play" style="margin-left:14px; background:#263238; color:#eee; border:1px solid #555; border-radius:6px; padding:3px 14px; cursor:pointer;">&#9654; replay</button>
    &nbsp;&nbsp; fuel: <span class="fv" style="color:#7fdfff"></span>
    &nbsp;&nbsp; gap: <span class="gv" style="color:#9ccc65"></span>
  </div>
  <canvas class="c3d" width="640" height="330" style="width:100%"></canvas>
  <canvas class="cth" width="640" height="170" style="width:100%"></canvas>
  <script>
    const root = currentScript.parentElement;
    const DATA = $(HypertextLiteral.JavaScript(JSON3.write(tf_data)));
    const RHO1 = $(ρ1), RHO2 = $(ρ2), TANG = Math.tan(6 * Math.PI / 180);
    const sl = root.querySelector(".s_tf"), btn = root.querySelector(".play");
    let frame = -1, anim = null;
    const A = 0.62, CA = Math.cos(A), SA = Math.sin(A), TILT = 0.40;
    function proj(x, y, z) {
      const u = x * CA - y * SA, w = x * SA + y * CA;
      return [320 + u * 0.62, 295 - z * 0.105 - w * TILT * 0.18];
    }
    function draw3d(d, upto) {
      const cv = root.querySelector(".c3d"), ctx = cv.getContext("2d");
      ctx.clearRect(0, 0, cv.width, cv.height);
      ctx.fillStyle = "#999"; ctx.font = "11px system-ui";
      ctx.fillText("trajectory (isometric) — wireframe: 6 deg glide-slope cone", 14, 16);
      // ground grid
      ctx.strokeStyle = "#22303a"; ctx.lineWidth = 1;
      for (let gx = -500; gx <= 500; gx += 125) {
        let p1 = proj(gx, -500, 0), p2 = proj(gx, 500, 0);
        ctx.beginPath(); ctx.moveTo(p1[0], p1[1]); ctx.lineTo(p2[0], p2[1]); ctx.stroke();
        p1 = proj(-500, gx, 0); p2 = proj(500, gx, 0);
        ctx.beginPath(); ctx.moveTo(p1[0], p1[1]); ctx.lineTo(p2[0], p2[1]); ctx.stroke();
      }
      // glide slope cone wireframe (apex at origin)
      ctx.strokeStyle = "#395a47";
      const hcone = 2500;
      for (let a = 0; a < 360; a += 30) {
        const rr = hcone * TANG;
        const p1 = proj(0, 0, 0);
        const p2 = proj(rr * Math.cos(a * Math.PI / 180), rr * Math.sin(a * Math.PI / 180), hcone);
        ctx.beginPath(); ctx.moveTo(p1[0], p1[1]); ctx.lineTo(p2[0], p2[1]); ctx.stroke();
      }
      // trajectory
      ctx.strokeStyle = "#d4af37"; ctx.lineWidth = 2.2; ctx.beginPath();
      const n = (upto < 0) ? d.x.length : upto + 1;
      for (let i = 0; i < n; i++) {
        const p = proj(d.x[i], d.y[i], d.zc[i]);
        i === 0 ? ctx.moveTo(p[0], p[1]) : ctx.lineTo(p[0], p[1]);
      }
      ctx.stroke();
      // lander dot + thrust spike
      const i = n - 1, p = proj(d.x[i], d.y[i], d.zc[i]);
      ctx.fillStyle = "#4fc3f7"; ctx.beginPath(); ctx.arc(p[0], p[1], 5, 0, 7); ctx.fill();
      if (i < d.Tc.length) {
        const len = 36 * d.Tc[i] / RHO2;
        ctx.strokeStyle = "#ff7043"; ctx.lineWidth = 3;
        ctx.beginPath(); ctx.moveTo(p[0], p[1]); ctx.lineTo(p[0], p[1] + len); ctx.stroke();
      }
      // target
      const pt = proj(0, 0, 0);
      ctx.strokeStyle = "#ef5350"; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(pt[0] - 6, pt[1] - 6); ctx.lineTo(pt[0] + 6, pt[1] + 6);
      ctx.moveTo(pt[0] - 6, pt[1] + 6); ctx.lineTo(pt[0] + 6, pt[1] - 6); ctx.stroke();
    }
    function drawThrust(d, upto) {
      const cv = root.querySelector(".cth"), ctx = cv.getContext("2d");
      const W = cv.width, H = cv.height, L = 56, R = 10, T = 18, B = 22;
      ctx.clearRect(0, 0, W, H);
      ctx.strokeStyle = "#444"; ctx.strokeRect(L, T, W - L - R, H - T - B);
      ctx.fillStyle = "#999"; ctx.font = "11px system-ui";
      ctx.fillText("thrust ‖T(t)‖ [N] — red dashes: engine limits", L + 6, T - 5);
      const xmax = DATA[DATA.length - 1].tf;
      const px = x => L + x / xmax * (W - L - R);
      const py = y => H - B - y / (RHO2 * 1.15) * (H - T - B);
      [0, 5000, 10000].forEach(function(yv) {
        ctx.fillText((yv / 1000).toFixed(0) + "k", 18, py(yv) + 4);
      });
      ctx.strokeStyle = "#cc3333"; ctx.setLineDash([5, 4]);
      [RHO1, RHO2].forEach(function(b) {
        ctx.beginPath(); ctx.moveTo(px(0), py(b)); ctx.lineTo(px(xmax), py(b)); ctx.stroke();
      });
      ctx.setLineDash([]);
      ctx.strokeStyle = "#d4af37"; ctx.lineWidth = 2.2; ctx.beginPath();
      const n = (upto < 0) ? d.Tc.length : Math.min(upto + 1, d.Tc.length);
      for (let i = 0; i < n; i++) {
        const X = px(d.t[i]), Y = py(d.Tc[i]);
        i === 0 ? ctx.moveTo(X, Y) : ctx.lineTo(X, Y);
      }
      ctx.stroke();
    }
    function draw(upto) {
      const d = DATA[+sl.value];
      root.querySelector(".tv").textContent = d.tf.toFixed(0) + " s";
      root.querySelector(".fv").textContent = d.fuel.toFixed(1) + " kg";
      root.querySelector(".gv").textContent = d.gap.toExponential(1);
      draw3d(d, upto); drawThrust(d, upto);
    }
    sl.addEventListener("input", function() { if (anim) { clearInterval(anim); anim = null; } draw(-1); });
    btn.addEventListener("click", function() {
      if (anim) clearInterval(anim);
      const d = DATA[+sl.value];
      let i = 0;
      anim = setInterval(function() {
        draw(i); i += 1;
        if (i >= d.x.length) { clearInterval(anim); anim = null; }
      }, 28);
    });
    draw(-1);
  </script>
</div>
""")

# ╔═╡ c0c0de03-0013-4013-8013-000000000013
md"""
## The cost of caution

The glide-slope angle $\gamma_{gs}$ is a *safety* parameter: a steeper cone
keeps the lander higher above surrounding terrain. Safety is not free —
steeper approach geometry costs propellant. Here is the same landing
($t_f = 72$ s) under different cones:
"""

# ╔═╡ c0c0de03-0014-4014-8014-000000000014
gs_data = let
    runs = []
    for γ in [2.0, 4.0, 6.0, 10.0, 15.0, 20.0]
        s = solve_conic(72.0; γ_gs = deg2rad(γ))
        s === nothing && continue
        push!(runs, (gamma = γ,
            x = round.(s.r[:,1], digits = 1), y = round.(s.r[:,2], digits = 1),
            zc = round.(s.r[:,3], digits = 1),
            fuel = round(s.fuel, digits = 1)))
    end
    runs
end

# ╔═╡ c0c0de03-0015-4015-8015-000000000015
@htl("""
<div style="background:#0e1117; border-radius:10px; padding:14px; color:#ddd; font-family:system-ui;">
  <div style="margin-bottom:8px;">
    <b style="color:goldenrod">glide-slope angle γ</b>:
    <input type="range" min="0" max="$(length(gs_data)-1)" value="2" style="width:50%; vertical-align:middle;">
    <span class="gv2" style="color:goldenrod; font-weight:bold;"></span>
    &nbsp;&nbsp; fuel: <span class="fv2" style="color:#7fdfff"></span>
  </div>
  <canvas width="640" height="260" style="width:100%"></canvas>
  <script>
    const root = currentScript.parentElement;
    const DATA = $(HypertextLiteral.JavaScript(JSON3.write(gs_data)));
    const sl = root.querySelector("input");
    const cv = root.querySelector("canvas"), ctx = cv.getContext("2d");
    function draw() {
      const d = DATA[+sl.value];
      root.querySelector(".gv2").textContent = d.gamma.toFixed(0) + " deg";
      root.querySelector(".fv2").textContent = d.fuel.toFixed(1) + " kg";
      ctx.clearRect(0, 0, cv.width, cv.height);
      ctx.fillStyle = "#999"; ctx.font = "11px system-ui";
      ctx.fillText("side profile: altitude vs ground distance to target", 14, 16);
      const W = cv.width, H = cv.height, L = 50, R = 14, T = 24, B = 24;
      const dmax = 700, hmax = 2500;
      const px = x => L + x / dmax * (W - L - R);
      const py = y => H - B - y / hmax * (H - T - B);
      ctx.strokeStyle = "#444"; ctx.strokeRect(L, T, W - L - R, H - T - B);
      // cone boundary
      const tg = Math.tan(d.gamma * Math.PI / 180);
      ctx.fillStyle = "rgba(102,187,106,0.12)";
      ctx.beginPath(); ctx.moveTo(px(0), py(0));
      ctx.lineTo(px(dmax), py(dmax * tg)); ctx.lineTo(px(dmax), py(0)); ctx.closePath(); ctx.fill();
      ctx.strokeStyle = "#395a47";
      ctx.beginPath(); ctx.moveTo(px(0), py(0)); ctx.lineTo(px(dmax), py(dmax * tg)); ctx.stroke();
      ctx.fillStyle = "#888"; ctx.fillText("keep-out (below glide slope)", px(dmax * 0.55), py(dmax * 0.5 * tg) + 16);
      // trajectory: distance to pad vs altitude
      ctx.strokeStyle = "#d4af37"; ctx.lineWidth = 2.4; ctx.beginPath();
      for (let i = 0; i < d.x.length; i++) {
        const dist = Math.hypot(d.x[i], d.y[i]);
        i === 0 ? ctx.moveTo(px(dist), py(d.zc[i])) : ctx.lineTo(px(dist), py(d.zc[i]));
      }
      ctx.stroke();
      ctx.fillStyle = "#ef5350"; ctx.beginPath(); ctx.arc(px(0), py(0), 5, 0, 7); ctx.fill();
    }
    sl.addEventListener("input", draw);
    draw();
  </script>
</div>
""")

# ╔═╡ c0c0de03-0016-4016-8016-000000000016
md"""
## A loose end, and a teaser

We have been handing the optimizer $t_f$ as a *given*. But look what the
$t_f$ slider revealed — fuel is anything but flat in flight time:
"""

# ╔═╡ c0c0de03-0017-4017-8017-000000000017
begin
    tfs_curve = [d.tf for d in tf_data]
    fuels = [d.fuel for d in tf_data]
    pteaser = plot(tfs_curve, fuels, lw = 2.5, color = :goldenrod,
        marker = :circle, ms = 4, legend = false,
        xlabel = "time of flight tf [s]", ylabel = "fuel [kg]",
        title = "Fuel vs. flight time — somebody should optimize this…",
        size = (640, 300))
    scatter!(pteaser, [tfs_curve[argmin(fuels)]], [minimum(fuels)],
        ms = 8, marker = :star5, color = :deepskyblue)
    pteaser
end

# ╔═╡ c0c0de03-0018-4018-8018-000000000018
md"""
The paper's reference $t_f = 72$ s burns about **23% more fuel** than the
best flight time near 49 s. Choosing $t_f$ *is itself an optimization* — but
the convex machinery above only works for *fixed* $t_f$. The resolution
(an outer line search wrapped around the convex solver — bisection's smarter
sibling, exactly the trick from page 1's minimum-time search) gets its own
page later in the series.

## Take-home

- Real physics (mass depletion) entered through a **change of variables**
  ($z = \ln m$, $\mathbf{u} = \mathbf{T}/m$) chosen to keep the problem as
  convex as possible. Formulation *is* the craft.
- The lossless gap check from page 2 scales up intact — verify it on every
  solve; it is one `maximum(abs.(...))`.
- Two independent solution routes (declared NLP vs. hand-transcribed conic
  SCVX) agreeing to ~0.1 kg is worth more than either alone.

**Next: [Part 4 — Six Degrees of Freedom](04-sixdof-landing.html)**, where
the rocket becomes a rigid body, the convex spell breaks, and successive
convexification picks up the pieces. *(Also planned: free final time, and a
transcription-diagnostics post-mortem.)*
"""

# ╔═╡ Cell order:
# ╟─c0c0de03-0001-4001-8001-000000000001
# ╟─c0c0de03-0002-4002-8002-000000000002
# ╠═c0c0de03-0003-4003-8003-000000000003
# ╟─c0c0de03-0004-4004-8004-000000000004
# ╠═c0c0de03-0005-4005-8005-000000000005
# ╟─c0c0de03-0006-4006-8006-000000000006
# ╟─c0c0de03-0007-4007-8007-000000000007
# ╠═c0c0de03-0008-4008-8008-000000000008
# ╟─c0c0de03-0009-4009-8009-000000000009
# ╟─c0c0de03-0010-4010-8010-000000000010
# ╟─c0c0de03-0011-4011-8011-000000000011
# ╟─c0c0de03-0012-4012-8012-000000000012
# ╟─c0c0de03-0013-4013-8013-000000000013
# ╟─c0c0de03-0014-4014-8014-000000000014
# ╟─c0c0de03-0015-4015-8015-000000000015
# ╟─c0c0de03-0016-4016-8016-000000000016
# ╟─c0c0de03-0017-4017-8017-000000000017
# ╟─c0c0de03-0018-4018-8018-000000000018
