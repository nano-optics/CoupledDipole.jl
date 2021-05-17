# Far field cross-sections


### Extinction

The extinction cross-section may be obtained from the optical theorem as the interference between incident and scattered fields,

```math
\sigma_\text{ext} = 4\pi\kappa^{-1} k_1 \Im\left(\mathbf{P} \cdot \mathbf{E_\text{inc}}^* \right).
```

```@docs
extinction!
```

### Absorption 

The absorption cross-section is obtained by evaluating the work done by the total field, (but excluding self-reaction), on the dipoles:

```math
\sigma_\text{abs} =  4\pi \kappa^{-1} k_1 \left[\Im\left(\mathbf{P} \cdot \mathbf{E}^* \right) - \frac 2 3 \kappa^{-1} k_1^3 |\mathbf{P}|^2\right].
```


```@docs
absorption!
```

### Scattering

The scattering cross-section can be computed in two ways:

- as the difference between extinction and absorption
- by integrating the flux of the Poynting vector over all scattering directions in the far-field. 

Comparing both results might be useful as a consistency check.

```math
\sigma_\text{sca} = \kappa^{-2} k_1^4 \iint_\Omega \left|\sum_i \left(\mathbb{I} - \mathbf{\hat n}\otimes\mathbf{\hat n}\right) \mathbf{P}_i \mathrm{exp}(-ik_1 \mathbf{r}_i\cdot\mathbf{\hat n})\right|^2 \text{d}\Omega.
```

```@docs
scattering!
```
