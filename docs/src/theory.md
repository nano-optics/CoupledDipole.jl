# Coupled dipole theory

```@contents
```

From a microscopic viewpoint, as taken in the discrete dipole approximation, a collection of small inclusions such as molecules dispersed in a solvent may be described as a set of point dipoles in a vaccuum, whereby not only the molecules of interest but also those comprising the medium itself are described as discrete dipoles. This approach is numerically prohibitive, due to the large number of dipoles required to model the environment, which is in fact infinite in most situations of interest. 

With a more pragmatic approach, we are led to consider the response of dipoles within the framework of the macroscopic Maxwell equations, where the optical response of the solvent is encompassed in its refractive index ``n=\sqrt{\varepsilon}``. Waves propagating in this medium have a phase velocity ``c/n``. Light scattering by a collection of particles (or molecules) will thus be described as the interaction between dipole moments induced by the incident light and in mutual interaction via their own scattered field in the surrounding homogeneous medium, as illustrated in Fig.~\ref{fig:00schematic}.


We define the dipole moments ``\mathbf{P}`` to be proportional to the local (macroscopic) electric field,

```math
\mathbf{P} = \mathbb{\alpha} \mathbf{E}
```

where ``\mathbb{\alpha}`` is a ``3\times 3`` polarisability tensor, which will be discussed in section~\ref{sec:XXX}. The standard CD equations take the form,

```math
\mathbf{P}^{i}= \mathbb{\alpha}_i \left( \mathbf{E_\text{inc}} + \sum_{j\neq i} \mathbb{G}_{ij} \mathbf{P}^j\right),
```

yielding a linear system of 3N equations to be solved for the unknown polarisation vector ``\mathbf{P}``.

## Green's tensor in a medium

We derive in the Appendix the macroscopic field associated with a dipole in a medium,

```math
\mathbf{E}(\mathbf{r}) =\kappa^{-1} \frac{\mathrm{exp}(ikr)}{r}\left\{k_1^2\left[\mathbf{P} - (\mathbf{\hat r}\cdot\mathbf{P})\mathbf{\hat r} \right] + \left(\frac{1}{r^2} - \frac{ik_1}{r}\right) \left[3\mathbf{\hat r}(\mathbf{\hat r}\cdot\mathbf{P}) - \mathbf{P}\right] \right\},\qquad \kappa = 4\pi\varepsilon_0\varepsilon
```

which we may write as

```math
\mathbf{E} = \mathbb{G} \mathbf{P}
```

with, 

```math
\mathbb{G} = \kappa^{-1} \frac{\mathrm{exp}(ikr)}{r}\left\{k_1^2\left[\mathbb{I} - \mathbf{\hat r}\otimes\mathbf{\hat r} \right] - \left(\frac{1}{r^2} - \frac{ik_1}{r}\right) \left[\mathbb{I} - 3\mathbf{\hat r}\otimes\mathbf{\hat r}\right] \right\}.
```

Note the constant prefactor ``\kappa`` which will also appear in the definition of the polarisability, allowing us to diregard it throughout the code implementation, as discussed in section~\ref{sec:equivalence}.

## Coupled dipole equations 

Since our code aims to model not only particles but also molecules, which can have a non-invertible tensor (for example a uniaxial response), this formulation is not suitable\footnote{note that for uniaxial or isotropic dipoles, a simplified version of the CDA could be derived cf Govorov, Berova, but we are aiming for a fully general coupled dipole method}: solving for the polarisation ``\mathbf{P}``, one needs to form a matrix with block diagonal elements ``\mathbb{\alpha}^{-1}``, where two of the principal components would become infinite. Instead, we recast the coupled dipole equations in a different form, solving for the self-consistent macroscopic fields at the dipole positions. This field ``\mathbf{E}^{i}`` is the sum of the incident field and the field scattered by all other dipoles,

```math
\mathbf{E}^i=\mathbf{E_\text{inc}}(\mathbf{r}_i) + \sum_{j\neq i} \mathbb{G}_{ij} \mathbb{\alpha}_j \mathbf{E}^j,
```

where ``\mathbb{G}_{ij} = \mathbb{G}(\mathbf{r}_i,\mathbf{r}_j)`` is the Green's tensor for the field created by dipole ``i`` at the location of dipole ``j`` in the infinite surrounding medium, 

```math
\mathbb{G}_{ij} = \kappa^{-1}\frac{\mathrm{exp}(ikr)_{ij}}{r_{ij}}\left\{k_1^2\left[\mathbb{I} - \mathbf{\hat r}_{ij}\otimes\mathbf{\hat r}_{ij} \right] - \left(\frac{1}{r_{ij}^2} - \frac{ik_1}{r_{ij}}\right) \left[\mathbb{I} - 3\mathbf{\hat r}_{ij}\otimes\mathbf{\hat r}_{ij}\right] \right\}.
```

Note that optical reciprocity imposes the symmetry ``\mathbb{G}_{ji} = \mathbb{G}^\intercal_{ij}``. Grouping the fields ``\mathbf{E}^i`` and ``\mathbf{E_\text{inc}}(\mathbf{r}_i)`` into vectors ``\mathbf{E}`` and ``\mathbf{E_\text{inc}}`` of length ``3N``, where ``N`` is the number of dipoles, the coupled system takes the matrix form, 

```math
\mathbb{F} \mathbf{E}  = \mathbf{E_\text{inc}}
```

defining the free-space interaction matrix ``\mathbb F  = (\mathbb I - \mathbb G\mathbb{\alpha} )``.  The ``3\times 3`` diagonal blocks of ``\mathbb{F}`` are identity matrices, and each off-diagonal block is of the form ``\mathbb F_{ij} = -\mathbb G_{ij}\mathbb{\alpha}_{j}``. Note that the product ``\mathbb G\mathbb{\alpha}`` implies that ``\mathbb F`` is not longer (block-)symmetric, for an arbitrary set of polarisabilies (either unlike particles, or simply rotated along arbitrary orientations). For a three-dipoles system, the matrix representation can be written,

```math
\left[
\begin{pmatrix}
\mathbb{I}&\mathbb{O}&\mathbb{O}\\
\mathbb{O}&\mathbb{I}&\mathbb{O}\\
\mathbb{O}&\mathbb{O}&\mathbb{I}
\end{pmatrix} -
\begin{pmatrix}
\mathbb{O}&\mathbb{G}_{12}\mathbb{\alpha}_2&\mathbb{G}_{13}\mathbb{\alpha}_3\\
\mathbb{G}_{21}\mathbb{\alpha}_1&\mathbb{O}&\mathbb{G}_{23}\mathbb{\alpha}_3\\
\mathbb{G}_{31}\mathbb{\alpha}_1&\mathbb{G}_{32}\mathbb{\alpha}_2&\mathbb{O}
\end{pmatrix}
\right]
\begin{pmatrix}
\mathbf{E}^1\\
\mathbf{E}^2\\
\mathbf{E}^3
\end{pmatrix} = 
\begin{pmatrix}
\mathbf{E_\text{inc}}^1\\
\mathbf{E_\text{inc}}^2\\
\mathbf{E_\text{inc}}^3
\end{pmatrix} .
```

The linear system is solved by standard linear algebra routines (LU factorisation), yielding the self-consistent macroscopic electric field ``\mathbf{E}`` at each dipole location. From these local fields we can calculate the polarisation ``\mathbf{P} = \mathbb{\alpha} \mathbf{E}``.

## Order-of-scattering solution

An alternative approach to solving the system of coupled dipole equations is to follow an order-of-scattering approximation scheme. The linear system of equations is formally re-written as

```math
\begin{aligned}
\mathbb{F} \mathbf{E} &= \mathbf{E_\text{inc}}\\
 \mathbf{E} &= (\mathbb{I} - \mathbb{\alpha}\mathbb{G})^{-1} \mathbf{E_\text{inc}}\\
 \mathbf{E} &= (\mathbb{I} + \mathbb{\alpha}\mathbb{G} + (\mathbb{\alpha}\mathbb{G})^2 + (\mathbb{\alpha}\mathbb{G})^3 +\dots) \mathbf{E_\text{inc}}\\
 \mathbf{E} &= \underbrace{\mathbf{E_\text{inc}}}_{\mathbf{E}^0} + \underbrace{\mathbb{\alpha}\mathbb{G}\mathbf{E}^0}_{\mathbf{E}^1} + \underbrace{\mathbb{\alpha}\mathbb{G}\mathbf{E}^1}_{\mathbf{E}^2} + \underbrace{\mathbb{\alpha}\mathbb{G}\mathbf{E}^2}_{\mathbf{E}^3} +\dots
\end{aligned}.
```

The system ``(\mathbb{I} - \mathbb{\alpha}\mathbb{G}) \mathbf{E} = \mathbf{E_\text{inc}}`` is thus solved iteratively by defining a temporary variable ``\mathbf{E}^n`` and taking the following steps,

```math
\begin{aligned}
\mathbf{E}^0 &= \mathbf{E_\text{inc}}, & \mathbf{E} &= \mathbf{E}^0\\
\mathbf{E}^1 &= \mathbb{\alpha}\mathbb{G}\mathbf{E}^0,  & \mathbf{E} &= \mathbf{E}^0 + \mathbf{E}^1\\
\mathbf{E}^2 &= \mathbb{\alpha}\mathbb{G}\mathbf{E}^1,         & \mathbf{E} &= \mathbf{E}^0 + \mathbf{E}^1 + \mathbf{E}^2\\
\dots &
\end{aligned}
```

At each step we update the net field by adding a component corresponding to an higher order of multiple scattering. The procedure is followed iteratively until convergence of ``\mathbf{E}`` is achieved with sufficient accuracy. Typically, the convergence criterion may be chosen as the relative change in the value of a far-field cross-section from one iteration to the next. 

If the interaction between dipoles ``\mathbb{\alpha}\mathbb{G}`` is no stronger than the contribution due to the incident field alone (``\mathbb{I}``), convergence may be obtained and after a few iterations the polarisation and electric fields have reached their self-consistent value.

This procedure was implemented, and can be useful for relatively large systems of dipoles provided they do not interact too strongly. 


## Cross-sections

Physically-relevant quantities accessible to experiments that may be derived from the solution of the CD equations include the far-field cross-sections for scattering, absorption, and extinction.

### Scattering

The scattering cross-section can be computed by integrating the flux of the Poynting vector over all scattering directions in the far-field. For a single dipole, the result is,

```math
\sigma_\text{sca} = \frac{k_1^4} {6\pi\varepsilon} |\alpha|^2.
```

And for a collection of dipoles,

```math
\sigma_\text{sca} = \kappa^{-2} k_1^4 \iint_\Omega \left|\sum_i \left(\mathbb{I} - \mathbf{\hat n}\otimes\mathbf{\hat n}\right) \mathbf{P}_i \mathrm{exp}(-ik_1 \mathbf{r}_i\cdot\mathbf{\hat n})\right|^2 \text{d}\Omega.
```

As an alternative, the scattering cross-section may be obtained from the difference between extinction and absorption cross-sections\footnote{Markel derived this alternative expression by direct integration of the scattering cross-section over all scattering directions.}.

### Extinction

The extinction cross-section may be obtained from the optical theorem as the interference between incident and scattered fields\cite{mishchenko1},

```math
\sigma_\text{ext} = 4\pi\kappa^{-1} k_1 \Im\left(\mathbf{P} \cdot \mathbf{E_\text{inc}}^* \right).
```

### Absorption 

The absorption cross-section is obtained by evaluating the work done by the total field, (but excluding self-reaction), on the dipoles:

```math
\sigma_\text{abs} =  4\pi \kappa^{-1} k_1 \left[\Im\left(\mathbf{P} \cdot \mathbf{E}^* \right) - \frac 2 3 \kappa^{-1} k_1^3 |\mathbf{P}|^2\right].
```

## Polarisability prescriptions

We defined the polarisability by the following relation between dipole moment and electric field,

```math
\mathbf{P}=\mathbb{\alpha}\mathbf{E}
```

noting that different authors include various prefactors such as ``4\pi``, ``\varepsilon``, ``\varepsilon_0`` or combinations thereof. Polarisability prescriptions that are compatible with our chosen conventions are reviewed in this section. First, we consider the case of a subwavelength nanoparticle, with the standard static polarisability of a small sphere ("Clausius-Mossotti" polarisability), and a long-wavelength approximation for elongated particles that includes corrections for radiative damping and dynamic depolarisation.

We also discuss the polarisability of dye molecules, and its link to a microscopic (intrinsic) polarisability that may be obtained from first principle calculations.

### Effective scaling prefactor

The prefactor ``\kappa`` appears in several equations and simplifies in the final expression for the cross-sections. We can therefore simplify the formalism throughout by defining reduced variables ``\bar\alpha, \bar \mathbf{P}, \bar\mathbb{G}``. The correspondence is summarised in Table~\ref{tab:equivalence}. Only the scaled quantities are used in the code, but for simplicity we refer to them with conventional variable names.

Equivalence between theory and rescaled equations used in the code. The prefactor ``\kappa=4\pi\varepsilon_0\varepsilon`` can be simplified throughout by defining a suitably scaled polarisability.


|  | Theory (with ``\kappa = 4\pi\varepsilon_0\varepsilon``) | Equivalent formulation used in the code  |
|-----------------|:-------------|:--------------------------------|
| polarisability | ``\alpha``  | ``\bar\alpha = \kappa^{-1}\alpha``    |
| sphere (CM)     | ``\alpha_\text{cm} = \kappa a^3\frac{\varepsilon - \varepsilon}{\varepsilon +2 \varepsilon}``          | ``\bar\alpha_\text{cm} = a^3\frac{\varepsilon - \varepsilon}{\varepsilon +2 \varepsilon}``            |
| dipole moment      | ``\mathbf{P} = \alpha\mathbf{E}``         | ``\bar\mathbf{P} = \kappa^{-1}\mathbf{P}``            |
| Green's function      | ``\mathbb{G}= \kappa^{-1} \frac{\mathrm{exp}(ikr)}{r}\left\{\dots\right\}``         | ``\bar\mathbb{G} = \kappa\mathbb{G} = \frac{\mathrm{exp}(ikr)}{r}\left\{\dots\right\}``            |
| ``\sigma_\text{ext}``      | ``4\pi \kappa^{-1} k_1   \Im\left(\mathbf{P} \cdot \mathbf{E_\text{inc}}^* \right)``         | ``4\pi k_1  \Im\left(\bar\mathbf{P} \cdot \mathbf{E_\text{inc}}^* \right)``            |
| ``\sigma_\text{abs}``      | ``4\pi \kappa^{-1} k_1 \left[ \Im\left(\mathbf{P} \cdot \mathbf{E}^* \right) - \frac{2}{3}\kappa^{-1}k^3\|\mathbf{P}\|^2 \right]``         | ``4\pi k_1 \left[ \Im\left(\bar\mathbf{P} \cdot \mathbf{E}^* \right) - \frac{2}{3}k^3\|\bar\mathbf{P}\|^2 \right]``            |
| ``\sigma_\text{sca}``      | ``\kappa^{-2}k_1^4 \iint_\Omega \left\|\sum_i \left(\mathbb{I} - \mathbf{\hat n}\otimes\mathbf{\hat n}\right) \mathbf{P}_i e^{-ik_1 \mathbf{r}_i\cdot\mathbf{\hat n}}\right\|^2 \text{d}\Omega``         | ``k_1^4 \iint_\Omega \left\|\sum_i \left(\mathbb{I} - \mathbf{\hat n}\otimes\mathbf{\hat n}\right) \bar\mathbf{P}_i e^{-ik_1 \mathbf{r}_i\cdot\mathbf{\hat n}}\right\|^2 \text{d}\Omega``            |




### Particle polarisability

The static polarisability of a sphere of radius ``a`` in a medium is given by (ref. Jackson, Griffith)\footnote{note that the dielectric constant ``\varepsilon`` is sometimes factored in the dipole moment instead.},

```math
\mathbb{\alpha} = \kappa a^3 \frac{\varepsilon - \varepsilon'_p}{\varepsilon +2 \varepsilon'_p}\mathbb{I}, \qquad \text{(with ``\mathbf{P}=\mathbb{\alpha}\mathbf{E}``)}.
```

This equation describes the response of a spherical particle to a static electric field. In nano-optics, the nanoparticles may be much smaller than the wavelength (and therefore are excited by an essentially constant field at any given time), but the time-variation of the fields implies that the induced dipole *radiates*. This requires a correction to the polarisability to satisfy energy conservation, namely an imaginary component accounting for radiative reaction. This subject will be further discussed in section~\ref{sec:rr}. Meier and Wokaun proposed\cite{} that larger nanoparticles may also benefit from a dynamic-depolarisation correction, which takes into account the dephasing of the internal field across finite nanoparticles and provides a more accurate description of their optical response. 

For homogeneous nanospheres and spherical nanoshells, the Mie theory readily provides an exact polarisability, but in practice it is often sufficient to use an approximate (but closed-form) formula\cite{dmitri}. For other geometries, no analytical solution is known. The static polarisability of ellipsoids may be derived analytically, and long-wavelength corrections have been proposed following the strategy of Meier and Wokaun\ref{schatz}. The accuracy is however rather poor in practice\cite{moroz,luke}, and we use instead the semi-empirical formula proposed by Kuwata et al. based on fits of full-wave simulations\cite{kuwata}. 

```math
\alpha_i =  \frac{V}{4\pi}\frac{1}{L_i + \frac{\varepsilon}{\varepsilon'_p-\varepsilon} + 
              A \varepsilon x^2 +  B \varepsilon^2 x^4 - i\frac{4\pi^2\varepsilon^{3/2}}{3}  V /\lambda^3}
```
with ``V`` the volume of the particle, and ``x = \frac{2\pi a_i}{\lambda}`` the size parameter along the semi-axis ``a_i`` (``a``, ``b``, or ``c`` in the case of spheroids).

A numerical test of the validity of these polarisability prescriptions against exact results is presented in Appendix\ref{sec:kuwata}.

### Dye polarisability

In our derivation of the CD equations we have taken a fully macroscopic viewpoint. In the case of molecular dipoles, the embedded dipoles (external inclusions), however, are still microscopic entities responding to the microscopic (local) field\footnote{Note that this strategy relies on relatively low concentrations of inclusions, with the solvent concentration being so large in comparison that it largely dictates the optical response (refractive index) of the medium (?). In the DDA such a requirement is not applicable, since dense particles are modelled with discrete dipoles. The strategy used to consider an embedding medium is different: one considers an equivalent problem in vaccuum, rescaling the permittivity of each consituent material and the incident wavevector. Perhaps this strategy could also be used for molecular dipoles, but it is not clear how to scale the polarisability. One idea might be to consider the three orthogonal dipoles as describing a static ellipsoid, for which depolarisation factors are known, reverse-engineer the link between alpha and epsilon for this effective particle, and apply the equivalence strategy to this effective epsilon before performing the inverse transformation to an alpha. Does that make sense? Are we losing generality by considering an ellipsoid? (I wouldn't think so -- it's the only shape that has homogeneous internal field, and we want an homogeneous field to define a macroscopic polarisability over that effective volume).}.

At a microscopic level, a dipole reacts to the applied field with an intrinsic polarisability ``\alpha^\mu``. The applied field is enhanced by a local-field correction arising from the polarisation of the embedding medium, thus we may define an effective polarisability for the dipole responding to a macroscopic field,

```math
 p = L \alpha^\mu \mathbf{E}
```

where ``L = \dfrac{\varepsilon +2}{3}``. In turn, such a dipole moment produces a field which is enhanced by the same local field factor, from a reciprocity argument. For our purposes, the effective macroscopic polarisability is therefore defined as

```math
\alpha = L^2 \alpha^\mu.
```

where ``\alpha^\mu`` may have been obtained from first-principle calculations (e.g. DFT). With this link, the dipole moments all respond to macroscopic fields, and the coupling between dipoles follows the same equations as above, with macroscopic fields propagating in the ambient medium. We note that if the molecular response is derived from an experimental UV-vis absorption measurement, these considerations are not important: we can directly link the extinction (or absorption) cross-section to a macroscopic polarisability. This strategy is further discussed in the Appendix.

## Rotations and orientation averaging

The simulations are set up with respect to a fixed frame of reference ("lab frame""). The particle positions and orientations, as well as the incident field, are all specified in this common frame. 


Rotation matrices are defined via Euler angles. Such matrices are used to specify particle orientations by rotation of their polarisability tensor. 

Similarly, the incident wavevector and electric field are rotated to an arbitrary direction of incidence.  

For the simulation of optical activity, or of disordered systems, we are often required to perform an average over the directions of incident light (incidence and polarisation). A given incident field is characterised by a wave-vector ``\mathbf{k}`` and an electric vector ``{\mathbf{E}_\text{inc}}`` describing the light polarisation, and both of these vectors can be rotated using a rotation matrix ``R``. The spectra are averaged over two orthogonal polarisations, with incident wave-vectors that span the full range of ``\varphi\in[0,2\pi], \theta\in[0,\pi]``. This averaging is performed by numerical integration\footnote{note: a few papers propose analytical angular-averaged cross-sections, e.g. Khlebstov, but there doesn't appear to be any gain in efficiency as the formulas are relatively cumbersome and costly to evaluate numerically.}

```math
\left<\sigma\right>_\Omega=\frac1 {4\pi} \int_0^{2\pi}\int_{0}^{\pi} \sigma(\varphi,\psi) \sin \theta d\theta d\varphi.
```

using a quadrature scheme (Gauss-Legendre, Monte-Carlo, etc.). A lighter alternative is often used in such calculations, averaging only three particular directions of the incident field. Note that by choosing ``\theta`` points such that ``\cos\theta`` is uniformly distributed in ``[-1,1]``, the associated change of variable automatically embeds the ``\sin\theta`` factor.

