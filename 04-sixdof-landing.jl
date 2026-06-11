### A Pluto.jl notebook ###
# v0.20.0

#> [frontmatter]
#> order = "4"
#> title = "4. Six Degrees of Freedom: When Convexity Runs Out"
#> description = "Attitude dynamics break the convex spell. Successive convexification — and an honest account of every way it failed before it worked."

using Markdown
using InteractiveUtils

# ╔═╡ d0c0de04-0001-4001-8001-000000000001
md"""
# Six Degrees of Freedom: When Convexity Runs Out

*Part 4 — [back to Part 3](03-mars-landing.html)*

The 3-DoF lander of page 3 was a point with a magic thrust vector that could
point anywhere instantly. A real rocket is a rigid body: the engine is
bolted to the bottom, it can only gimbal a few degrees, and to redirect
thrust the whole vehicle must *rotate* — subject to inertia, torque limits,
and a tilt envelope. This is the SpaceX-booster-landing problem, and it is
where the beautiful convex story of pages 2–3 runs out.

New state: a quaternion $q$ (attitude) and body rates $\omega$. New physics:

```math
\dot q = \tfrac12\,\Omega(\omega)\,q, \qquad
J\dot\omega = \tau - \omega \times J\omega, \qquad
\tau = r_{cp} \times m\,\mathbf{u}_B,
```

and the acceleration the trajectory feels is $C(q)\,\mathbf{u}_B$ — thrust
expressed in the *body* frame, rotated into the world by the attitude. New
constraints: a **gimbal cone** (thrust within 20° of the body axis), a
**tilt cone** (body axis within 60° of vertical), and a **rate limit**
$\|\omega\| \le \omega_{\max}$.

> 🏃 As always: **Edit or run** → Binder for a live session. The widgets
> below are instant.
"""

# ╔═╡ d0c0de04-0002-4002-8002-000000000002
md"""
## Why there is no lossless trick this time

Look at the terms the new physics introduced: $C(q)\,\mathbf{u}_B$ multiplies
*unknown attitude* by *unknown thrust*. $\Omega(\omega)\,q$ multiplies
*unknown rate* by *unknown attitude*. $\omega \times J\omega$ is quadratic
in the unknown rate. These are **bilinear** terms — products of decision
variables — and no change of variables is known that hides them losslessly
inside a cone, the way $z = \ln m$ and the $\sigma$-lift handled mass and
the throttle annulus.

So we change strategy. If the problem won't become convex *globally*, make
it convex *locally* and iterate — **successive convexification (SCvx)**:

1. Take a reference trajectory $\bar x(t), \bar u(t)$ (an honest guess).
2. **Linearize** the nonconvex dynamics about it. Keep everything that is
   already convex — the $\sigma$-lift from page 2 stays, untouched.
3. Solve the resulting conic problem, but only trust it **near the
   reference** (a *trust region* $\|x - \bar x\| \le \eta$), and let the
   linearized dynamics be violated at a steep price (*virtual controls*
   $\nu$, penalized in the objective) so the subproblem is never infeasible.
4. Score the step: did the *true nonlinear* dynamics improve as much as the
   linear model predicted? Accept and grow $\eta$, or reject and shrink.
5. Repeat until the virtual controls vanish and the defect against true
   dynamics is negligible.

Convex optimization is demoted from *the answer* to *the inner loop* — but
each subproblem still solves in milliseconds with a certificate, and that
is what makes the outer loop practical.
"""

# ╔═╡ d0c0de04-0003-4003-8003-000000000003
begin
    using JuMP, Clarabel, LinearAlgebra, HypertextLiteral, JSON3

    # Mars lander, same trajectory problem as page 3, now with a body
    const g_I = [0.0, 0.0, -3.7114]
    const m_wet, m_dry = 1905.0, 1505.0
    const α_fuel = 1.0 / (225.0 * 9.80665)
    const ρ₁, ρ₂ = 4972.0, 13260.0
    const J_B = Diagonal([4400.0, 4400.0, 2000.0])   # inertia [kg m²]
    const J_inv = inv(Matrix(J_B))
    const r_cp = [0.0, 0.0, -1.5]                    # engine lever arm [m]
    const δ_max = deg2rad(20.0)                      # gimbal cone
    const θ_max = deg2rad(60.0)                      # tilt cone
    const ω_max = deg2rad(60.0)                      # rate limit
    const γ_gs = deg2rad(6.0)
    const r₀ = [450.0, -330.0, 2400.0]; const v₀ = [-40.0, 10.0, -10.0]
    const q₀ = [1.0, 0.0, 0.0, 0.0];    const ω₀ = [0.0, 0.0, 0.0]
    const t_f = 72.0; const N_steps = 40; const Δt = t_f / N_steps

    # Scaling: every variable O(1). This is not cosmetic — the trust region
    # and the virtual-control penalty compare these numbers across blocks.
    const S_r = 2400.0; const S_v = 50.0; const S_ω = ω_max
    const S_u = ρ₂ / m_wet
end

# ╔═╡ d0c0de04-0004-4004-8004-000000000004
begin
    skew(w) = [0 -w[3] w[2]; w[3] 0 -w[1]; -w[2] w[1] 0]

    function quat_to_dcm(q)   # body → inertial
        w, x, y, z = q
        [1-2(y^2+z^2)  2(x*y-z*w)    2(x*z+y*w);
         2(x*y+z*w)    1-2(x^2+z^2)  2(y*z-x*w);
         2(x*z-y*w)    2(y*z+x*w)    1-2(x^2+y^2)]
    end

    # ∂(C(q)u)/∂q at (q̄, ū): how attitude changes the felt acceleration
    function dcbi_dq_times_u(q, u)
        w, x, y, z = q
        e = [x, y, z]
        exu = cross(e, u)
        d_dw = 2exu
        d_dx = 2w*cross([1,0,0],u) + 2cross([1,0,0],exu) + 2cross(e,cross([1,0,0],u))
        d_dy = 2w*cross([0,1,0],u) + 2cross([0,1,0],exu) + 2cross(e,cross([0,1,0],u))
        d_dz = 2w*cross([0,0,1],u) + 2cross([0,0,1],exu) + 2cross(e,cross([0,0,1],u))
        hcat(d_dw, d_dx, d_dy, d_dz)
    end

    function omega_mat(w)
        x, y, z = w
        [0 -x -y -z; x 0 z -y; y -z 0 x; z y -x 0]
    end

    qnorm(q) = (n = norm(q); n > 1e-15 ? q / n : [1.0, 0, 0, 0])

    function slerp(qa, qb, w)
        d = clamp(dot(qa, qb), -1.0, 1.0)
        if d < 0; qb = -qb; d = -d; end
        d > 0.9995 && return qnorm(qa + w * (qb - qa))
        θ = acos(d)
        (sin((1-w)θ) * qa + sin(w*θ) * qb) / sin(θ)
    end

    # Defect of a candidate against the TRUE nonlinear discrete dynamics,
    # accumulated in the same scaled units as the virtual controls ν —
    # so the trust-region ratio compares like with like
    function compute_defects(r, v, q, ω, z, u_B, σ)
        total = 0.0
        for k in 1:N_steps
            C = quat_to_dcm(q[k,:]); mk = exp(z[k])
            a = g_I + C * u_B[k,:]
            r_p = r[k,:] + Δt*v[k,:] + 0.5Δt^2*a
            v_p = v[k,:] + Δt*a
            q_p = qnorm(q[k,:] + Δt*0.5*omega_mat(ω[k,:])*q[k,:])
            τ = cross(r_cp, mk*u_B[k,:]); gyro = cross(ω[k,:], J_B*ω[k,:])
            ω_p = ω[k,:] + Δt*J_inv*(τ - gyro)
            z_p = z[k] - α_fuel*Δt*σ[k]
            total += norm(r[k+1,:]-r_p)/S_r + norm(v[k+1,:]-v_p)/S_v
            total += norm(q[k+1,:]-q_p) + norm(ω[k+1,:]-ω_p)/S_ω + abs(z[k+1]-z_p)
        end
        total
    end
end

# ╔═╡ d0c0de04-0005-4005-8005-000000000005
md"""
## The reference: lift the 3-DoF answer into 6-DoF

SCvx needs a starting reference, and the smartest cheap one is the solution
of the **simpler model we already trust**: solve the page-3 problem (3-DoF,
milliseconds), then *attach an attitude to it* — point the body axis along
each thrust vector.

One subtlety here cost a day of debugging, so it gets a code comment and a
rule: **the reference must satisfy the hard constraints.** The thrust at
$t=0$ is ~29° off vertical, but the boundary condition pins the initial
attitude upright — and with the body upright, that thrust direction violates
the 20° gimbal cone. A trust region is a ball *around the reference*; center
it on an infeasible point and every small-η subproblem is infeasible by
construction, which reads as mysterious solver failures. So: blend the
attitude from the required $q(0)$ into the thrust-aligned profile, and clamp
the body-frame thrust onto the gimbal cone.
"""

# ╔═╡ d0c0de04-0006-4006-8006-000000000006
function solve_3dof()
    N = N_steps
    z_ref = collect(range(log(m_wet), log(m_wet - 0.5(m_wet - m_dry)), length = N))
    local rs, vs, zs, us, ss
    for _ in 1:4
        model = Model(Clarabel.Optimizer); set_silent(model)
        @variable(model, r[1:N+1,1:3]); @variable(model, v[1:N+1,1:3])
        @variable(model, z[1:N+1]); @variable(model, u[1:N,1:3]); @variable(model, s[1:N])
        for k in 1:N, i in 1:3
            @constraint(model, r[k+1,i] == r[k,i] + Δt*v[k,i] + 0.5Δt^2*(g_I[i]+u[k,i]))
            @constraint(model, v[k+1,i] == v[k,i] + Δt*(g_I[i]+u[k,i]))
        end
        for k in 1:N; @constraint(model, z[k+1] == z[k] - α_fuel*Δt*s[k]); end
        for i in 1:3
            @constraint(model, r[1,i] == r₀[i]); @constraint(model, v[1,i] == v₀[i])
            @constraint(model, r[N+1,i] == 0);   @constraint(model, v[N+1,i] == 0)
        end
        @constraint(model, z[1] == log(m_wet)); @constraint(model, z[N+1] >= log(m_dry))
        for k in 1:N
            @constraint(model, [s[k]; u[k,:]] in SecondOrderCone())
            z0 = z_ref[k]; ez = exp(-z0)
            @constraint(model, s[k] >= ρ₁*ez*(1-(z[k]-z0)+0.5(z[k]-z0)^2))
            @constraint(model, s[k] <= ρ₂*ez*(1-(z[k]-z0)))
        end
        for k in 1:N+1
            @constraint(model, [r[k,3]/tan(γ_gs); r[k,1]; r[k,2]] in SecondOrderCone())
            @constraint(model, r[k,3] >= 0)
        end
        @objective(model, Min, Δt*sum(s))
        optimize!(model)
        rs = value.(r); vs = value.(v); zs = value.(z); us = value.(u); ss = value.(s)
        z_ref = zs[1:N]
    end
    rs, vs, zs, us, ss
end

# ╔═╡ d0c0de04-0007-4007-8007-000000000007
function init_6dof(r3, v3, z3, u3, s3)
    N = N_steps
    q = zeros(N+1, 4); ω = zeros(N+1, 3); uB = zeros(N, 3)
    for k in 1:N           # attitude that points the body axis along thrust
        T = u3[k,:]; Tn = norm(T)
        T̂ = Tn > 1e-8 ? T/Tn : [0.0,0,1]
        cr = cross([0.0,0,1], T̂); d = dot([0.0,0,1], T̂)
        q[k,:] = d > 0.9999 ? [1,0,0,0] : qnorm([1+d; cr])
    end
    q[N+1,:] = q₀
    n_blend = 4            # ...but honor the q(0) boundary condition
    q_tgt = q[n_blend+1,:]
    for k in 1:n_blend
        q[k,:] = slerp(q₀, q_tgt, (k-1)/n_blend)
    end
    for k in 1:N           # body thrust, clamped into the gimbal cone
        u = quat_to_dcm(q[k,:])' * u3[k,:]
        n = norm(u)
        if n > 1e-9 && u[3]/n < cos(δ_max)
            perp = u .- [0.0,0,u[3]]; pn = norm(perp)
            dir = pn > 1e-12 ? perp/pn : [1.0,0,0]
            u = n * (sin(δ_max)*dir .+ [0.0,0,cos(δ_max)])
        end
        uB[k,:] = u
    end
    for k in 1:N           # rates from finite-differenced attitude
        dq = q[k+1,:] - q[k,:]
        qc = [q[k,1], -q[k,2], -q[k,3], -q[k,4]]
        qv = qc[1]*dq[2:4] + dq[1]*qc[2:4] + cross(qc[2:4], dq[2:4])
        west = 2qv/Δt; wn = norm(west)
        ω[k,:] = wn > 0.8ω_max ? west*0.8ω_max/wn : west
    end
    ω[1,:] = ω₀
    (r3/S_r, v3/S_v, copy(z3), q, ω/S_ω, uB/S_u, s3/S_u)
end

# ╔═╡ d0c0de04-0008-4008-8008-000000000008
md"""
## The convex subproblem

Everything below is one conic program. Three details carry the whole method,
and each one is a scar:

- **Anchor the linearization at the *propagated* reference** $f(\bar x_k,
  \bar u_k)$ — *not* at the reference's own next state $\bar x_{k+1}$. With
  the wrong anchor, setting the variables equal to the reference satisfies
  the constraint with $\nu \equiv 0$ *identically*, so the subproblem is
  structurally blind to how dynamically inconsistent the reference is. The
  optimizer then "fixes" everything by exploiting linearization error (it
  swung the attitude 70° on iteration 1), the predicted-vs-actual ratio is
  garbage, and no step is ever accepted.
- **Virtual controls** $\nu$ on every dynamics row, with an exact-penalty
  weight: the subproblem can always fall back to "reference + $\nu$", so it
  is never infeasible — and $\|\nu\|$ *is* the linear model's estimate of
  the dynamics defect.
- **One trust region over all scaled states** — which only means anything
  because everything was scaled to O(1) first. A radius of 1.0 is modest
  for position (2.4 km) and absurd for a quaternion (~120°); the shared η
  must be sized for the most fragile block (the attitude linearization).
"""

# ╔═╡ d0c0de04-0009-4009-8009-000000000009
function solve_subproblem(r_ref, v_ref, z_ref, q_ref, ω_ref, u_ref, s_ref, η;
                          λ_vc = 1e5, λ_att = 5e2)
    N = N_steps
    model = Model(Clarabel.Optimizer); set_silent(model)

    @variable(model, r[1:N+1,1:3]); @variable(model, v[1:N+1,1:3])
    @variable(model, q[1:N+1,1:4]); @variable(model, ω[1:N+1,1:3])
    @variable(model, z[1:N+1]); @variable(model, u_B[1:N,1:3]); @variable(model, σ[1:N])
    @variable(model, ν_r[1:N,1:3]); @variable(model, ν_v[1:N,1:3])
    @variable(model, ν_q[1:N,1:4]); @variable(model, ν_w[1:N,1:3])
    @variable(model, ν_rn[1:N] >= 0); @variable(model, ν_vn[1:N] >= 0)
    @variable(model, ν_qn[1:N] >= 0); @variable(model, ν_wn[1:N] >= 0)

    R_vr = Δt*S_v/S_r; R_ar = 0.5Δt^2/S_r; R_av = Δt/S_v
    g_r = 0.5Δt^2*g_I/S_r; g_v = Δt*g_I/S_v

    for k in 1:N
        qr = q_ref[k,:]; ωr = ω_ref[k,:]*S_ω; ur = u_ref[k,:]*S_u
        mk = exp(z_ref[k]); C_k = quat_to_dcm(qr)
        D_k = dcbi_dq_times_u(qr, ur)              # attitude sensitivity

        for i in 1:3                                # position & velocity
            acc_u = sum(C_k[i,j]*u_B[k,j] for j in 1:3)
            acc_q = sum(D_k[i,j]*(q[k,j]-qr[j]) for j in 1:4)
            @constraint(model, r[k+1,i] == r[k,i] + R_vr*v[k,i] + g_r[i]
                        + R_ar*S_u*acc_u + R_ar*acc_q + ν_r[k,i])
            @constraint(model, v[k+1,i] == v[k,i] + g_v[i]
                        + R_av*S_u*acc_u + R_av*acc_q + ν_v[k,i])
        end

        # quaternion: anchored at the PROPAGATED reference (see above)
        Ω_k = omega_mat(ωr)
        dqdω = 0.5*[-qr[2] -qr[3] -qr[4]; qr[1] -qr[4] qr[3];
                     qr[4]  qr[1] -qr[2]; -qr[3] qr[2]  qr[1]]
        A_q = Matrix{Float64}(I,4,4) + Δt*0.5Ω_k
        B_q = Δt*dqdω*S_ω
        q_prop = qnorm(qr + Δt*0.5Ω_k*qr)
        for i in 1:4
            rhs = q_prop[i]
            for j in 1:4; rhs += A_q[i,j]*(q[k,j]-qr[j]); end
            for j in 1:3; rhs += B_q[i,j]*(ω[k,j]-ω_ref[k,j]); end
            @constraint(model, q[k+1,i] == rhs + ν_q[k,i])
        end

        # body rates: same anchoring rule
        Jωr = J_B*ωr; τ_ref = cross(r_cp, mk*ur); gyro_ref = cross(ωr, Jωr)
        A_w = Matrix{Float64}(I,3,3) - Δt*J_inv*(skew(ωr)*J_B - skew(Jωr))
        B_wu = Δt*J_inv*(mk*skew(r_cp))
        ω_prop = (ωr + Δt*J_inv*(τ_ref - gyro_ref))/S_ω
        for i in 1:3
            val = ω_prop[i]
            for j in 1:3
                val += A_w[i,j]*(ω[k,j]-ω_ref[k,j])
                val += (S_u/S_ω)*B_wu[i,j]*(u_B[k,j]-u_ref[k,j])
            end
            @constraint(model, ω[k+1,i] == val + ν_w[k,i])
        end

        @constraint(model, z[k+1] == z[k] - α_fuel*Δt*S_u*σ[k])
        @constraint(model, [ν_rn[k]; ν_r[k,:]] in SecondOrderCone())
        @constraint(model, [ν_vn[k]; ν_v[k,:]] in SecondOrderCone())
        @constraint(model, [ν_qn[k]; ν_q[k,:]] in SecondOrderCone())
        @constraint(model, [ν_wn[k]; ν_w[k,:]] in SecondOrderCone())
    end

    for i in 1:3
        @constraint(model, r[1,i] == r₀[i]/S_r); @constraint(model, v[1,i] == v₀[i]/S_v)
        @constraint(model, ω[1,i] == ω₀[i]/S_ω)
        @constraint(model, r[N+1,i] == 0); @constraint(model, v[N+1,i] == 0)
        @constraint(model, ω[N+1,i] == 0)
    end
    for i in 1:4; @constraint(model, q[1,i] == q₀[i]); end
    @constraint(model, z[1] == log(m_wet)); @constraint(model, z[N+1] >= log(m_dry))

    for k in 1:N                       # LCvx lift — pages 2-3, unchanged
        @constraint(model, [σ[k]; u_B[k,:]] in SecondOrderCone())
        z0 = z_ref[k]; ez = exp(-z0)
        @constraint(model, S_u*σ[k] >= ρ₁*ez*(1-(z[k]-z0)+0.5(z[k]-z0)^2))
        @constraint(model, S_u*σ[k] <= ρ₂*ez*(1-(z[k]-z0)))
        @constraint(model, cos(δ_max)*σ[k] <= u_B[k,3])      # gimbal cone
    end
    for k in 1:N+1
        @constraint(model, [r[k,3]/tan(γ_gs); r[k,1]; r[k,2]] in SecondOrderCone())
        @constraint(model, r[k,3] >= 0)
        @constraint(model, [1.0; ω[k,:]] in SecondOrderCone())   # ‖ω‖ ≤ ω_max
        qr = q_ref[k,:]                                          # tilt (linearized)
        tv = 1 - 2(qr[2]^2 + qr[3]^2)
        @constraint(model, tv - 4qr[2]*(q[k,2]-qr[2]) - 4qr[3]*(q[k,3]-qr[3]) >= cos(θ_max))
    end

    for k in 1:N+1                     # the trust region
        δ = vcat(r[k,:]-r_ref[k,:], v[k,:]-v_ref[k,:],
                 q[k,:]-q_ref[k,:], ω[k,:]-ω_ref[k,:])
        @constraint(model, [η; δ] in SecondOrderCone())
    end
    for k in 1:N
        @constraint(model, [η; vcat(u_B[k,:]-u_ref[k,:], [σ[k]-s_ref[k]])] in SecondOrderCone())
    end

    @variable(model, att_err >= 0)
    @constraint(model, [att_err; q[N+1,:]-q₀] in SecondOrderCone())
    vc = sum(ν_rn[k]+ν_vn[k]+ν_qn[k]+ν_wn[k] for k in 1:N)
    @objective(model, Min, Δt*S_u*sum(σ) + λ_vc*vc + λ_att*att_err)

    optimize!(model)
    termination_status(model) in (OPTIMAL, ALMOST_OPTIMAL) || return nothing
    qv = value.(q); for k in 1:N+1; qv[k,:] = qnorm(qv[k,:]); end
    (r = value.(r), v = value.(v), q = qv, ω = value.(ω), z = value.(z),
     u = value.(u_B), σ = value.(σ),
     vc = sum(value(ν_rn[k])+value(ν_vn[k])+value(ν_qn[k])+value(ν_wn[k]) for k in 1:N))
end

# ╔═╡ d0c0de04-0010-4010-8010-000000000010
md"""
## The outer loop: trust, but score every step

The acceptance test is the conscience of the method: after each subproblem,
re-propagate the candidate through the **true** nonlinear dynamics and ask
whether reality improved as much as the linear model promised
($\rho = \text{actual}/\text{predicted}$). Honest steps grow the trust
region; dishonest ones shrink it and are rejected. Crucially, both sides of
that ratio must be in the **same units** — the third scar below.
"""

# ╔═╡ d0c0de04-0011-4011-8011-000000000011
function run_scvx(; max_iters = 30, tol_vc = 1e-5, tol_defect = 1e-2)
    r3, v3, z3, u3, s3 = solve_3dof()
    fuel3 = m_wet - exp(z3[end])
    r_ref, v_ref, z_ref, q_ref, ω_ref, u_ref, s_ref = init_6dof(r3, v3, z3, u3, s3)

    snap(rr, q, z, u, s, defect, η, ρr, act) = (
        x = round.(rr[:,1]*S_r, digits = 1), y = round.(rr[:,2]*S_r, digits = 1),
        zc = round.(rr[:,3]*S_r, digits = 1),
        tilt = [round(rad2deg(acos(clamp(1-2(q[k,2]^2+q[k,3]^2), -1, 1))), digits = 1)
                for k in 1:N_steps+1],
        Tc = round.([exp(z[k])*norm(u[k,:])*S_u for k in 1:N_steps], digits = 0),
        fuel = round(m_wet - exp(z[end]), digits = 1),
        defect = round(defect, sigdigits = 3), eta = round(η, digits = 3),
        rho = ρr === nothing ? nothing : round(ρr, digits = 2), action = act)

    defect_prev = compute_defects(r_ref*S_r, v_ref*S_v, q_ref, ω_ref*S_ω,
                                  z_ref, u_ref*S_u, s_ref*S_u)
    J_prev = Δt*sum(s_ref*S_u) + 1e5*defect_prev
    η = 0.25
    history = Any[snap(r_ref, q_ref, z_ref, u_ref*S_u, s_ref*S_u, defect_prev, η, nothing, "init")]

    final = nothing
    for it in 1:max_iters
        sol = solve_subproblem(r_ref, v_ref, z_ref, q_ref, ω_ref, u_ref, s_ref, η)
        if sol === nothing
            η = max(1e-3, 0.5η)
            push!(history, snap(r_ref, q_ref, z_ref, u_ref*S_u, s_ref*S_u,
                                defect_prev, η, nothing, "infeasible→shrink"))
            η <= 1e-3 && break
            continue
        end
        defect = compute_defects(sol.r*S_r, sol.v*S_v, sol.q, sol.ω*S_ω,
                                 sol.z, sol.u*S_u, sol.σ*S_u)
        fuel_Δv = Δt*S_u*sum(sol.σ)
        J_cand = fuel_Δv + 1e5*defect
        J_model = fuel_Δv + 1e5*sol.vc
        predicted = J_prev - J_model; actual = J_prev - J_cand
        ρr = abs(predicted) > 1e-6 ? actual/predicted : (actual >= 0 ? 1.0 : 0.0)

        if ρr >= 0.0     # accept
            r_ref, v_ref, q_ref, ω_ref, z_ref = sol.r, sol.v, sol.q, sol.ω, sol.z
            u_ref, s_ref = sol.u, sol.σ
            J_prev = J_cand; defect_prev = defect
            act = ρr >= 0.7 ? "accept↑" : ρr >= 0.1 ? "accept" : "accept↓"
            η = ρr >= 0.7 ? min(3.0, 1.5η) : ρr >= 0.1 ? η : max(1e-3, 0.5η)
        else
            act = "REJECT"
            η = max(1e-3, 0.5η)
        end
        push!(history, snap(sol.r, sol.q, sol.z, sol.u*S_u, sol.σ*S_u, defect, η, ρr, act))

        if sol.vc < tol_vc && defect < tol_defect && startswith(act, "accept")
            final = (r = r_ref*S_r, q = q_ref, ω = ω_ref*S_ω, z = z_ref,
                     u = u_ref*S_u, σ = s_ref*S_u)
            break
        end
        η <= 1e-3 && break
    end
    (final = final, history = history, fuel3 = fuel3)
end

# ╔═╡ d0c0de04-0012-4012-8012-000000000012
scvx = run_scvx()

# ╔═╡ d0c0de04-0013-4013-8013-000000000013
begin
    fin = scvx.final
    un6 = [norm(fin.u[k,:]) for k in 1:N_steps]
    gap6 = round(maximum(abs.(fin.σ .- un6)), sigdigits = 2)
    fuel6 = round(m_wet - exp(fin.z[end]), digits = 1)
    fuel3r = round(scvx.fuel3, digits = 1)
    tilt6 = round(maximum(h -> h, [rad2deg(acos(clamp(1-2(fin.q[k,2]^2+fin.q[k,3]^2),-1,1)))
                                   for k in 1:N_steps+1]), digits = 1)
    gim6 = round(maximum([rad2deg(acos(clamp(fin.u[k,3]/(norm(fin.u[k,:])+1e-12),-1,1)))
                          for k in 1:N_steps]), digits = 1)
    om6 = round(maximum([rad2deg(norm(fin.ω[k,:])) for k in 1:N_steps+1]), digits = 1)
    nit = length(scvx.history) - 1
    Markdown.parse("""
    Converged in **$nit iterations**. The scorecard:

    | quantity | value | limit / reference |
    |---|---|---|
    | fuel | **$fuel6 kg** | 3-DoF optimum: $fuel3r kg |
    | lossless gap, σ vs ‖u‖ | **$gap6** | the page-2 check, alive inside SCvx |
    | max tilt | $tilt6 ° | 60° |
    | max gimbal | $gim6 ° | 20° |
    | max body rate | $om6 °/s | 60°/s |

    Full rigid-body attitude costs about **3 kg** over the point-mass ideal —
    the price of having to *rotate* to redirect thrust. And the lossless
    relaxation stays tight inside every subproblem: the convex core of
    pages 2–3 survives intact at the heart of a nonconvex method.
    """)
end

# ╔═╡ d0c0de04-0014-4014-8014-000000000014
md"""
## Watch it converge

This is the part no equation conveys: drag through the SCvx iterations and
watch a guessed trajectory *negotiate* with physics. Early iterations move
boldly (large trust region), rejected steps pull η down, and the defect —
the gap between the linear story and true dynamics — collapses by orders of
magnitude.
"""

# ╔═╡ d0c0de04-0015-4015-8015-000000000015
@htl("""
<div style="background:#0e1117; border-radius:10px; padding:14px; color:#ddd; font-family:system-ui;">
  <div style="margin-bottom:8px;">
    <b style="color:goldenrod">SCvx iteration</b>:
    <input type="range" min="0" max="$(length(scvx.history)-1)" value="0" style="width:46%; vertical-align:middle;">
    <span class="iv" style="color:goldenrod; font-weight:bold;"></span>
    &nbsp; <span class="av" style="font-weight:bold;"></span>
    &nbsp;&nbsp; defect: <span class="dv" style="color:#7fdfff"></span>
    &nbsp;&nbsp; η: <span class="ev" style="color:#aaa"></span>
    &nbsp;&nbsp; fuel: <span class="fv" style="color:#9ccc65"></span>
  </div>
  <canvas class="c3d" width="640" height="300" style="width:100%"></canvas>
  <canvas class="cti" width="640" height="160" style="width:100%"></canvas>
  <script>
    const root = currentScript.parentElement;
    const H = $(HypertextLiteral.JavaScript(JSON3.write(scvx.history)));
    const sl = root.querySelector("input");
    const A = 0.62, CA = Math.cos(A), SA = Math.sin(A);
    function proj(x, y, z) {
      const u = x*CA - y*SA, w = x*SA + y*CA;
      return [320 + u*0.62, 268 - z*0.095 - w*0.072];
    }
    function draw() {
      const d = H[+sl.value];
      root.querySelector(".iv").textContent = sl.value;
      const av = root.querySelector(".av");
      av.textContent = d.action;
      av.style.color = d.action.startsWith("accept") ? "#9ccc65" :
                       d.action === "init" ? "#aaa" : "#ef5350";
      root.querySelector(".dv").textContent = Number(d.defect).toExponential(1);
      root.querySelector(".ev").textContent = d.eta;
      root.querySelector(".fv").textContent = d.fuel.toFixed(1) + " kg";
      const cv = root.querySelector(".c3d"), ctx = cv.getContext("2d");
      ctx.clearRect(0, 0, cv.width, cv.height);
      ctx.fillStyle = "#999"; ctx.font = "11px system-ui";
      ctx.fillText("candidate trajectory at this iteration (isometric)", 14, 16);
      ctx.strokeStyle = "#22303a";
      for (let g = -500; g <= 500; g += 125) {
        let p1 = proj(g,-500,0), p2 = proj(g,500,0);
        ctx.beginPath(); ctx.moveTo(p1[0],p1[1]); ctx.lineTo(p2[0],p2[1]); ctx.stroke();
        p1 = proj(-500,g,0); p2 = proj(500,g,0);
        ctx.beginPath(); ctx.moveTo(p1[0],p1[1]); ctx.lineTo(p2[0],p2[1]); ctx.stroke();
      }
      ctx.strokeStyle = "#d4af37"; ctx.lineWidth = 2.2; ctx.beginPath();
      for (let i = 0; i < d.x.length; i++) {
        const p = proj(d.x[i], d.y[i], d.zc[i]);
        i === 0 ? ctx.moveTo(p[0], p[1]) : ctx.lineTo(p[0], p[1]);
      }
      ctx.stroke();
      const pt = proj(0,0,0);
      ctx.strokeStyle = "#ef5350"; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(pt[0]-6,pt[1]-6); ctx.lineTo(pt[0]+6,pt[1]+6);
      ctx.moveTo(pt[0]-6,pt[1]+6); ctx.lineTo(pt[0]+6,pt[1]-6); ctx.stroke();
      const c2 = root.querySelector(".cti"), x2 = c2.getContext("2d");
      const W = c2.width, Hh = c2.height, L = 46, R = 10, T = 18, B = 22;
      x2.clearRect(0, 0, W, Hh);
      x2.strokeStyle = "#444"; x2.strokeRect(L, T, W-L-R, Hh-T-B);
      x2.fillStyle = "#999"; x2.font = "11px system-ui";
      x2.fillText("tilt angle over flight [deg] — red dash: 60 deg cone", L+6, T-5);
      const px = i => L + i/(d.tilt.length-1)*(W-L-R);
      const py = yv => Hh - B - yv/70*(Hh-T-B);
      x2.strokeStyle = "#cc3333"; x2.setLineDash([5,4]);
      x2.beginPath(); x2.moveTo(L, py(60)); x2.lineTo(W-R, py(60)); x2.stroke();
      x2.setLineDash([]);
      x2.strokeStyle = "#4fc3f7"; x2.lineWidth = 2.2; x2.beginPath();
      for (let i = 0; i < d.tilt.length; i++) {
        i === 0 ? x2.moveTo(px(i), py(d.tilt[i])) : x2.lineTo(px(i), py(d.tilt[i]));
      }
      x2.stroke();
    }
    sl.addEventListener("input", draw);
    draw();
  </script>
</div>
""")

# ╔═╡ d0c0de04-0016-4016-8016-000000000016
md"""
## The landed solution, with attitude

Press ▶: the white segment is the body axis, the orange spike is the
gimbaled thrust (note it is *not* always along the body — that small wiggle
is the gimbal doing its 15-degree best while the body swings through its
much larger tilt).
"""

# ╔═╡ d0c0de04-0017-4017-8017-000000000017
replay = let
    f = scvx.final
    bz = zeros(N_steps+1, 3); Td = zeros(N_steps, 3)
    for k in 1:N_steps+1
        bz[k,:] = quat_to_dcm(f.q[k,:]) * [0.0, 0, 1]
    end
    for k in 1:N_steps
        TI = quat_to_dcm(f.q[k,:]) * f.u[k,:]
        n = norm(TI); Td[k,:] = n > 1e-9 ? TI/n : [0.0,0,1]
    end
    (x = round.(f.r[:,1], digits = 1), y = round.(f.r[:,2], digits = 1),
     zc = round.(f.r[:,3], digits = 1),
     bx = round.(bz[:,1], digits = 3), by = round.(bz[:,2], digits = 3),
     bz = round.(bz[:,3], digits = 3),
     tx = round.(Td[:,1], digits = 3), ty = round.(Td[:,2], digits = 3),
     tz = round.(Td[:,3], digits = 3),
     Tc = round.([exp(f.z[k])*norm(f.u[k,:]) for k in 1:N_steps], digits = 0))
end

# ╔═╡ d0c0de04-0018-4018-8018-000000000018
@htl("""
<div style="background:#0e1117; border-radius:10px; padding:14px; color:#ddd; font-family:system-ui;">
  <button class="play" style="background:#263238; color:#eee; border:1px solid #555; border-radius:6px; padding:4px 16px; cursor:pointer;">&#9654; replay landing</button>
  <span class="tv" style="margin-left:12px; color:goldenrod;"></span>
  <canvas width="640" height="360" style="width:100%; margin-top:8px;"></canvas>
  <script>
    const root = currentScript.parentElement;
    const D = $(HypertextLiteral.JavaScript(JSON3.write(replay)));
    const RHO2 = $(ρ₂); const DT = $(Δt);
    const cv = root.querySelector("canvas"), ctx = cv.getContext("2d");
    const A = 0.62, CA = Math.cos(A), SA = Math.sin(A);
    function proj(x, y, z) {
      const u = x*CA - y*SA, w = x*SA + y*CA;
      return [320 + u*0.62, 330 - z*0.115 - w*0.075];
    }
    let anim = null;
    function draw(i) {
      ctx.clearRect(0, 0, cv.width, cv.height);
      root.querySelector(".tv").textContent = "t = " + (i*DT).toFixed(1) + " s";
      ctx.strokeStyle = "#22303a";
      for (let g = -500; g <= 500; g += 125) {
        let p1 = proj(g,-500,0), p2 = proj(g,500,0);
        ctx.beginPath(); ctx.moveTo(p1[0],p1[1]); ctx.lineTo(p2[0],p2[1]); ctx.stroke();
        p1 = proj(-500,g,0); p2 = proj(500,g,0);
        ctx.beginPath(); ctx.moveTo(p1[0],p1[1]); ctx.lineTo(p2[0],p2[1]); ctx.stroke();
      }
      ctx.strokeStyle = "#d4af37"; ctx.lineWidth = 1.6; ctx.globalAlpha = 0.55;
      ctx.beginPath();
      for (let k = 0; k <= i; k++) {
        const p = proj(D.x[k], D.y[k], D.zc[k]);
        k === 0 ? ctx.moveTo(p[0], p[1]) : ctx.lineTo(p[0], p[1]);
      }
      ctx.stroke(); ctx.globalAlpha = 1.0;
      // body axis segment
      const L = 130;
      const c = [D.x[i], D.y[i], D.zc[i]];
      const b = [D.bx[i], D.by[i], D.bz[i]];
      const top = proj(c[0]+b[0]*L/2, c[1]+b[1]*L/2, c[2]+b[2]*L/2);
      const bot = proj(c[0]-b[0]*L/2, c[1]-b[1]*L/2, c[2]-b[2]*L/2);
      // gimbaled flame from the base, opposite thrust
      const ki = Math.min(i, D.Tc.length - 1);
      const fl = 200 * D.Tc[ki] / RHO2;
      const fp = proj(c[0]-b[0]*L/2 - D.tx[ki]*fl, c[1]-b[1]*L/2 - D.ty[ki]*fl,
                      c[2]-b[2]*L/2 - D.tz[ki]*fl);
      ctx.strokeStyle = "#ff7043"; ctx.lineWidth = 3.5;
      ctx.beginPath(); ctx.moveTo(bot[0], bot[1]); ctx.lineTo(fp[0], fp[1]); ctx.stroke();
      ctx.strokeStyle = "#fff"; ctx.lineWidth = 4;
      ctx.beginPath(); ctx.moveTo(bot[0], bot[1]); ctx.lineTo(top[0], top[1]); ctx.stroke();
      ctx.fillStyle = "#4fc3f7";
      ctx.beginPath(); ctx.arc(top[0], top[1], 3.5, 0, 7); ctx.fill();
      const pt = proj(0,0,0);
      ctx.strokeStyle = "#ef5350"; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(pt[0]-6,pt[1]-6); ctx.lineTo(pt[0]+6,pt[1]+6);
      ctx.moveTo(pt[0]-6,pt[1]+6); ctx.lineTo(pt[0]+6,pt[1]-6); ctx.stroke();
    }
    root.querySelector(".play").addEventListener("click", function() {
      if (anim) clearInterval(anim);
      let i = 0;
      anim = setInterval(function() {
        draw(i); i += 1;
        if (i >= D.x.length) { clearInterval(anim); anim = null; }
      }, 70);
    });
    draw(D.x.length - 1);
  </script>
</div>
""")

# ╔═╡ d0c0de04-0019-4019-8019-000000000019
md"""
## What was actually hard — a debugging confession

This page's algorithm did not work on the first try, or the fifth. The first
running version made *zero* progress: every step rejected, trust region
shrinking to nothing, subproblems eventually "infeasible". For the record —
because the failure modes are more instructive than the successes — the four
bugs, in the order they were found:

1. **The linearization was anchored at the reference's next state**
   ($\bar x_{k+1}$) instead of the propagated reference $f(\bar x_k, \bar
   u_k)$. Consequence: virtual controls were identically zero *at the
   reference*, so the subproblem literally could not see the dynamics
   defect. The tell was a flat contradiction in the logs: $\|\nu\| \approx
   10^{-9}$ ("dynamics perfect!") while the true defect *grew* from 0.1
   to 196. **A model that cannot represent its own error will lie to you
   with perfect confidence.**
2. **The initial reference violated the hard constraints** (thrust-aligned
   attitude at $t=0$ vs. the upright boundary condition and gimbal cone).
   Every small trust region around it was infeasible *by geometry*, which
   surfaced as inscrutable solver failures only after the trust region had
   shrunk — far from the actual cause.
3. **The trust ratio compared different units**: the model's $\|\nu\|$ in
   scaled coordinates vs. the true defect in meters. The acceptance test
   was dimensionally meaningless; whether a step was "good" depended on
   the units of the position vector.
4. **The conic solver itself was fragile**: the original returned
   NUMERICAL_ERROR / INFEASIBLE on small-trust-region subproblems that are
   provably feasible (reference + ν is always available). Swapping
   solvers (ECOS → Clarabel) made every subproblem solve cleanly —
   sometimes the bug is in the part you assumed was solid.

None of these are exotic. All four produced the *same* outward symptom —
"SCvx doesn't converge" — which is why the diagnosis discipline from page 2
matters more here than anywhere: **instrument the model and the truth
separately, and chase contradictions between them.**

## Principles for taking convexification to harder problems

Distilled from this series, in roughly the order you should reach for them:

1. **Spend formulation effort before algorithm effort.** The log-mass trick
   and the σ-lift bought pages 2–3 a *globally* solvable problem. Hunt for
   changes of variables and lossless lifts first; they are worth more than
   any outer loop.
2. **Keep the convex core convex.** When you do iterate, don't linearize
   what is already a cone. SCvx here linearizes *only* the attitude
   coupling; the LCvx machinery rides inside unchanged — and its tightness
   check keeps passing, which is free verification of the inner loop.
3. **Linearize about a propagated reference, and give the model an escape
   valve** (virtual controls with exact penalty). The subproblem must be
   able to *express* the defect, or it will spend linearization error to
   hide it.
4. **Trust regions are promises about linearization validity** — size them
   for the most nonlinear block (here, attitude: a quaternion perturbation
   of 1 is a ~120° swing), which is only possible if you have **scaled
   every variable to O(1)** and measure every penalty in the same units.
5. **Initialize from the model one notch simpler.** The 3-DoF solution
   lifted into 6-DoF (with constraints *enforced on the lift*) starts the
   loop a few honest iterations from the answer.
6. **Score steps against the truth, not the model.** Re-propagate the real
   dynamics every iteration; accept on actual-vs-predicted improvement.
   This single habit is what turns "linearize and hope" into an algorithm.
7. **Verify at every level, and read *where* failures occur.** Tightness
   gaps, defects, constraint margins — each is one line to compute, and
   each failure's *location* (which arc, which timestep, which block)
   points at its cause. Every wrong conclusion this series corrected was
   caught by exactly such a check.

*Planned next: free final time (wrapping a line search around the convex
solver), and a transcription-diagnostics post-mortem (the full story of the
"2.61 gap" from page 2).*
"""

# ╔═╡ Cell order:
# ╟─d0c0de04-0001-4001-8001-000000000001
# ╟─d0c0de04-0002-4002-8002-000000000002
# ╠═d0c0de04-0003-4003-8003-000000000003
# ╟─d0c0de04-0004-4004-8004-000000000004
# ╟─d0c0de04-0005-4005-8005-000000000005
# ╠═d0c0de04-0006-4006-8006-000000000006
# ╠═d0c0de04-0007-4007-8007-000000000007
# ╟─d0c0de04-0008-4008-8008-000000000008
# ╠═d0c0de04-0009-4009-8009-000000000009
# ╟─d0c0de04-0010-4010-8010-000000000010
# ╠═d0c0de04-0011-4011-8011-000000000011
# ╠═d0c0de04-0012-4012-8012-000000000012
# ╟─d0c0de04-0013-4013-8013-000000000013
# ╟─d0c0de04-0014-4014-8014-000000000014
# ╟─d0c0de04-0015-4015-8015-000000000015
# ╟─d0c0de04-0016-4016-8016-000000000016
# ╟─d0c0de04-0017-4017-8017-000000000017
# ╟─d0c0de04-0018-4018-8018-000000000018
# ╟─d0c0de04-0019-4019-8019-000000000019
