# Background and Motivation

The coupled dipole method is a convenient approximation to describe multiple scattering in relatively sparse clusters of particles (or molecules) smaller than the wavelength. 

The central approximation is that the response of a given particle to the incident field is described in the electric dipole approximation as 
```math
\mathbf{P} = \mathbb{\alpha} \mathbf{E}
```
where ``\mathbf{E}`` is the net electric field at that location. More details on the theory are summarised [here](./theory).

## Code design and conventions

From a user's perspective the code provides 2 high-level functions to perform calculations of:

1. Fixed orientation far-field cross-sections (extinction, absorption and scattering), for multiple wavelengths and incidence directions
2. Orientation-averaged cross-sections (extinction, absorption and scattering) and associated circular dichroism spectra, using numerical cubature over the full solid angle

The high-level functions require at least 2 inputs: 

- A `Cluster` object, describing the geometry of the particle cluster
- A `Material` object, describing the wavelength-dependent optical properties of the various media

### Geometry description

`Cluster` is a structure comprising 4 fields,

- `positions`, a vector of N particle positions stored as 3-vectors storing cartesian coordinates x,y,z

### Material description

### High level functions

#### Fixed orientation cross-sections

#### Angular averaging and circular dichroism

### Rotations

### Angular averaging

