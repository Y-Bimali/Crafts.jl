"""
    bondcolors(poly::Polyform)

The binding site colors `(c1, c2)` of every bond of `poly`, in the order `Roly.bonds` gives them.
"""
bondcolors(poly::Polyform) = ((_sitecolor(poly, p1, s1), _sitecolor(poly, p2, s2))
                              for ((p1, s1), (p2, s2)) in bonds(poly))

_sitecolor(poly::Polyform, particle, site) =
    Roly.color(bindingrules(poly), Roly.SpeciesSiteLoc(poly.particles[particle].speciesindex, site))

"""
    twistcorrection(poly::Polyform, bond)

Return the rotation that carries the partner site of `bond` into the twist the potential is written
for.

A site with a rotational symmetry about its normal is only determined up to that symmetry, so a bond
may sit in any of `Roly.twistfreedom` twists, all of them physically the same contact. The potentials
here are written for the reference twist alone, and evaluating one at a bond that sits in another
branch puts the bonded configuration at a saddle rather than a minimum.

Throws if the sites are in no admissible twist at all, which means the contact is not a contact.
"""
function twistcorrection(poly::Polyform, bond)
    (p1, s1), (p2, s2) = bond
    D = dimension(poly)
    F = Roly.numtype(poly)
    sys = bindingrules(poly)

    b1 = bindingsite(poly.particles[p1], sys, s1)
    b2 = bindingsite(poly.particles[p2], sys, s2)
    ntwists = Roly.twistfreedom(b1, b2)
    t = Roly.twist(b1, b2)
    isnothing(t) && throw(ArgumentError(
        "sites ($p1, $s1) and ($p2, $s2) are bonded but their frames match no admissible twist " *
        "of the $ntwists this bond allows, so the bonded configuration is not a minimum"))

    return Roly.standard_rotation(F, Val(D), t, ntwists) * inv(Roly.standard_rotation(F, Val(D)))
end

"""
    map_potential(bond_potential, poly::Polyform; embed3d=false)

Build the function that gives the total bond energy of `poly` as a function of a `dtot x nparticles`
matrix of displacements from the bonded configuration.
"""
function map_potential(bond_potential, poly::Polyform; embed3d=false)
    sys = bindingrules(poly)
    bs = collect(bonds(poly))
    twists = [twistcorrection(poly, b) for b in bs]

    d = embed3d ? 3 : dimension(poly)
    dtot = d * (d+1) ÷ 2

    function energy_fn(ξs::AbstractMatrix)
        E = 0.0
        rottype = d == 2 ? Angle2d{eltype(ξs)} : RotXYZ{eltype(ξs)}
        for (bi, ((p1, s1), (p2, s2))) in enumerate(bs)
            particle1 = poly.particles[p1]
            particle2 = poly.particles[p2]
            spcs1 = particle1.speciesindex
            spcs2 = particle2.speciesindex

            particle_pose1 = embed3d ? Pose{3}(particle1.pose) : particle1.pose
            particle_pose2 = embed3d ? Pose{3}(particle2.pose) : particle2.pose
            site_pose1 = let sp=bindingsite(species(sys, spcs1), s1).pose
                embed3d ? Pose{3}(sp) : sp
            end
            site_pose2 = let sp=bindingsite(species(sys, spcs2), s2).pose * twists[bi]
                embed3d ? Pose{3}(sp) : sp
            end

            # transform coordinates such that minimum is at zero.
            pose1 = Pose(ξs[1:d, p1] + particle_pose1.x, rottype(ξs[d+1:dtot, p1]...) * particle_pose1.psi)
            pose2 = Pose(ξs[1:d, p2] + particle_pose2.x, rottype(ξs[d+1:dtot, p2]...) * particle_pose2.psi)

            x1 = pose1 * site_pose1
            x2 = pose2 * site_pose2

            E += bondenergy(bond_potential, x1, x2, _sitecolor(poly, p1, s1),
                            _sitecolor(poly, p2, s2))
        end
        return E
    end

    return energy_fn
end

"""
    hessian(bond_potential, poly; embed3d=false)

Hessian of [`map_potential`](@ref) at the bonded configuration, of shape
`(dtot*nparticles, dtot*nparticles)`.
"""
function hessian(bond_potential, poly; embed3d=false)
    d = embed3d ? 3 : dimension(poly)
    n = nparticles(poly)
    dtot = d * (d+1) ÷ 2

    energy_fn = map_potential(bond_potential, poly; embed3d)
    x0 = zeros(Roly.numtype(poly), dtot, n)
    return ForwardDiff.hessian(energy_fn, x0 .+ sqrt(eps(eltype(x0))))
end

function orientational_volume(d::Integer)
    d == 2 && return 2π
    d == 3 && return 8π^2
    throw(ArgumentError("orientational volume requires d=2 or d=3"))
end

function inertia_tensor(poly::Polyform; embed3d=false)
    d = dimension(poly)
    com = sum(let x=p.pose.x; x end for p in poly.particles) / nparticles(poly)
    inertia = sum(let x=p.pose.x - com; norm(x)^2 * I(d) - x * x' end for p in poly.particles)
    embed3d || return inertia
    return [inertia zeros(2); zeros(2)' tr(inertia)]
end

"""
    EntropySolver

Supertype of the methods for turning a structure and a bond potential into a partition function.
See [`TreeLike`](@ref), [`MeanField`](@ref), [`TetheredLaplace`](@ref) and [`COMLaplace`](@ref).
"""
abstract type EntropySolver end

"""
    TetheredLaplace()

Entropy solver that expands the bond potential to second order about the bonded configuration,
holding one particle fixed to remove the global translations and rotations.
"""
struct TetheredLaplace <: EntropySolver end

"""
    COMLaplace()

Entropy solver that expands the bond potential to second order about the bonded configuration,
working in the frame of the structure's centre of mass.
"""
struct COMLaplace <: EntropySolver end

"""
    MeanField()

Entropy solver that applies mean field theory and averages over all particle tetherings.
"""
struct MeanField <: EntropySolver end

"""
    MeanFieldIsotropic()
    MeanFieldTopological()

Entropy solver that is [`MeanField`](@ref) with every single-particle Hessian replaced by `zᵢ` copies
of the single-bond one, reducing the answer to the tree-like result times a penalty fixed by the
degree sequence alone.

`MeanFieldTopological` substitutes the structure mean coordination number `z̄ = 2 n_b / n` for every `zᵢ`,
resulting in

    Ω_s = |SO(d)| ω^(n-1) (n / 2 n_b)^(dtot (n-1) / 2) / σ_s .

This underestimates `Ω_s`, by an amount growing with the variance of the degrees.
"""
struct MeanFieldIsotropic{uniform} <: EntropySolver end
MeanFieldIsotropic() = MeanFieldIsotropic{false}()

const MeanFieldTopological = MeanFieldIsotropic{true}

"""
    TreeLike(omega=1)

Entropy solver that treats every structure as a tree. A structure that really is a tree has exactly `n - 1` bonds
and the answer is exact. For non-trees, the bond volumes are averaged over the bonds. 
Given a bond potential, the bondvolumes are computed from [`bondvolume`](@ref), with no potential each 
bond contributes `omega`.
"""
struct TreeLike <: EntropySolver
    omega::Float64
end
TreeLike() = TreeLike(1.0)

"""
    EntropyModel(potential, solver=COMLaplace(); embed3d=false)
    EntropyModel(solver=TreeLike(); embed3d=false)

An EntropyModel is a bond potential together with a solver that computes bond volumes.

`embed3d=true` embeds a two-dimensional structure in 3d.
"""
struct EntropyModel{P,S<:EntropySolver}
    potential::P
    solver::S
    embed3d::Bool
end
EntropyModel(potential, solver::EntropySolver=COMLaplace(); embed3d=false) =
    EntropyModel(potential, solver, embed3d)

function EntropyModel(solver::EntropySolver=TreeLike(); embed3d=false)
    solver isa TreeLike || throw(ArgumentError(
        "$(nameof(typeof(solver))) needs a bond potential; use `EntropyModel(potential, solver)`."))
    return EntropyModel(nothing, solver, embed3d)
end

"""
    entropy(model::EntropyModel, poly::Polyform)

Partition function of the structure `poly` under `model`.
"""
entropy(m::EntropyModel, poly::Polyform) = _entropy(m.potential, poly, m.solver; m.embed3d)

"""
    entropydimension(m::EntropyModel, x)

The dimension in which `m` scores the structures of `x`, a `BindingRules` or a `Polyform`: 3 whenever
`m` embeds, and the dimension of `x` itself otherwise.
"""
entropydimension(m::EntropyModel, x) = m.embed3d ? 3 : dimension(x)

# A Laplace solver reads off the curvature at the pose Roly bonds the sites in, so it is only an
# expansion of anything if the potential is actually stationary there.
# Only a rigid spring reports how far the bonded pose sits above its own minimum, so that is the
# only potential whose expansion point can be checked here rather than taken on trust.
_checkbondsrelaxed(bond_potential, poly::Polyform) = nothing

function _checkbondsrelaxed(rsp::RigidSpringPotential, poly::Polyform)
    for (c1, c2) in bondcolors(poly)
        excess = contactexcess(rsp, c1, c2)
        excess > sqrt(eps(typeof(excess))) && throw(ArgumentError(
            "invalid minimizer for Laplace's method: binding site colors ($c1, $c2) bond energy is $excess above the potential's minimum/"))
    end
    return nothing
end

function _entropy(::Nothing, poly::Polyform, solver::TreeLike; embed3d)
    d = embed3d ? 3 : dimension(poly)
    return solver.omega^(nparticles(poly) - 1) * orientational_volume(d) / symmetrynumber(poly)
end

function _entropy(bond_potential, poly::Polyform, ::TreeLike; embed3d)
    d = embed3d ? 3 : dimension(poly)
    Ω = orientational_volume(d) / symmetrynumber(poly)
    nparticles(poly) == 1 && return Ω

    # A tree of `n` particles has exactly `n - 1` bonds, so the geometric mean is just their
    # product; anything with cycles is sees `n - 1` "typical" bonds instead.
    logω = mean(logbondvolume(bond_potential, c1, c2; d) for (c1, c2) in bondcolors(poly))
    return Ω * exp((nparticles(poly) - 1) * logω)
end

# The coordination number of every particle of `poly`, in one pass over its bonds.
function _degrees(poly::Polyform)
    degrees = zeros(Int, nparticles(poly))
    for bond in bonds(poly), loc in bond
        degrees[loc.particle] += 1
    end
    return degrees
end

function _entropy(bond_potential, poly::Polyform, ::MeanFieldIsotropic{uniform};
                  embed3d) where {uniform}
    treelike = _entropy(bond_potential, poly, TreeLike(1); embed3d)
    n = nparticles(poly)
    n == 1 && return treelike

    d = embed3d ? 3 : dimension(poly)
    dtot = d * (d + 1) ÷ 2
    uniform && return treelike * (2 * nbonds(poly) / n)^(-dtot * (n - 1) / 2)
    return treelike * exp(-dtot * (n - 1) / (2n) * sum(log, _degrees(poly)))
end

function _entropy(bond_potential, poly::Polyform, ::MeanField; embed3d)
    d = embed3d ? 3 : dimension(poly)
    n = nparticles(poly)
    Ω = orientational_volume(d) / symmetrynumber(poly)
    n == 1 && return Ω

    _checkbondsrelaxed(bond_potential, poly)
    dtot = d * (d + 1) ÷ 2
    H = hessian(bond_potential, poly; embed3d)
    logω = 0.0
    for i in 1:n
        block = ((dtot * (i - 1)) + 1):(dtot * i)
        logω += dtot / 2 * log(2π) - logdet(Symmetric(H[block, block])) / 2
    end
    return Ω * exp((n - 1) / n * logω)
end

function _entropy(rsp::RigidSpringPotential, poly::Polyform, ::MeanField; embed3d)
    d = embed3d ? 3 : dimension(poly)
    n = nparticles(poly)
    Ω = orientational_volume(d) / symmetrynumber(poly)
    n == 1 && return Ω

    _checkbondsunstressed(rsp, poly)
    logω = sum(_logomega_compatible(k, S) for (k, S) in _patchcovariances(rsp, poly, embed3d))
    return Ω * exp((n - 1) / n * logω)
end

_sitepose(poly::Polyform, particle, site, embed3d) =
    let sp = bindingsite(species(bindingrules(poly), poly.particles[particle].speciesindex),
                         site).pose
        embed3d ? Pose{3}(sp) : sp
    end

# The total stiffness kᵢ carried by each particle, and Σᵢ, the covariance of every patch on it in
# its own frame. In an unstressed structure a neighbour's averaged patches land exactly on the ones
# facing them, so ⟨Sᵢ⟩ = Σᵢ Rᵢᵀ and Φ_d sees only the eigenvalues of Σᵢ.
function _patchcovariances(rsp::RigidSpringPotential{D,F}, poly::Polyform, embed3d;
                           bondlist=bonds(poly)) where {D,F}
    n = nparticles(poly)
    ks = zeros(F, n)
    means = fill(zero(SVector{D,F}), n)
    seconds = fill(zero(SMatrix{D,D,F}), n)

    for ((p1, s1), (p2, s2)) in bondlist
        slots = _slots(rsp, _sitecolor(poly, p1, s1), _sitecolor(poly, p2, s2))
        k = rsp.k[slots...]
        for (p, s, slot) in ((p1, s1, slots[1]), (p2, s2, slots[2]))
            pose = _sitepose(poly, p, s, embed3d)
            ks[p] += k
            means[p] += k * pose.x
            # the cloud is centred on its own site, so the site sits at its mean patch position and
            # carrying it over to the particle's frame is the parallel axis theorem
            seconds[p] += k * (pose.psi * rsp.Scc[slot] * pose.psi' + pose.x * pose.x')
        end
    end
    return ((ks[i], seconds[i] / ks[i] - (means[i] / ks[i]) * (means[i] / ks[i])') for i in 1:n)
end

# Mean field takes the neighbours' patches at their averaged positions, which are the glued ones
# only if every bond really is relaxed where it sits.
function _checkbondsunstressed(rsp::RigidSpringPotential, poly::Polyform)
    _checkbondsrelaxed(rsp, poly)
    for (c1, c2) in bondcolors(poly)
        strain = strainenergy(rsp, c1, c2)
        strain > sqrt(eps(typeof(strain))) && throw(ArgumentError(
            "binding site colors ($c1, $c2) carry a strain energy of $strain, but mean field " *
            "theory is only worked out here for unstressed structures; see `strainenergy`"))
    end
    return nothing
end

function _entropy(bond_potential, poly::Polyform, ::TetheredLaplace; embed3d)
    d = embed3d ? 3 : dimension(poly)
    dtot = d * (d+1) ÷ 2

    _checkbondsrelaxed(bond_potential, poly)
    H = hessian(bond_potential, poly; embed3d)
    λs = eigvals(H[dtot+1:end, dtot+1:end])
    S_vib = -0.5 * sum(log, λs / (2π); init=0)
    return orientational_volume(d) * exp(S_vib) / symmetrynumber(poly)
end

function _entropy(bond_potential, poly::Polyform, ::COMLaplace; embed3d)
    d = embed3d ? 3 : dimension(poly)
    N = nparticles(poly)
    dtot = d * (d+1) ÷ 2

    _checkbondsrelaxed(bond_potential, poly)
    H = hessian(bond_potential, poly; embed3d)
    λs = eigvals(H)
    S_vib = -0.5 * sum(log, λ / (2π) for λ in @view λs[dtot+1:end]; init=0)

    inertia = inertia_tensor(poly; embed3d)
    if d == 2
        detG = N^d * (tr(inertia) + N)
    elseif d == 3
        detG = N^d * det(inertia + N * I(dtot - d))
    else
        throw(ArgumentError("dimension must be 2 or 3"))
    end
    return orientational_volume(d) * sqrt(detG) * exp(S_vib) / symmetrynumber(poly)
end
