# Signal models

## Moving lens

The moving-lens painter evaluates a truncated-NFW lensing deflection profile for each halo and projects the transverse velocity onto the local angular direction. Individual halo templates are summed linearly on the flat-sky canvas.

Required halo quantities are redshift, concentration, NFW scale radius, characteristic density, comoving distance, two transverse velocity components, and angular position.

## Thermal SZ

The tSZ painter uses a circular Gaussian profile. The thermodynamic-temperature conversion is frequency dependent through

```math
g(x)=x\coth(x/2)-4,
\qquad
x=h\nu/(k_B T_{\rm CMB}).
```

The default settings are 150 GHz, intrinsic Gaussian FWHM 5.83 arcmin, and amplitude multiplier 1.328.

## Kinetic SZ

The kSZ painter uses a Gaussian optical-depth profile and the relation

\[
\Delta T_{\rm kSZ}/T_{\rm CMB} = -\tau v_{\rm los}/c.
\]

Positive line-of-sight velocity is defined as motion away from the observer. The profile amplitude is supplied as a dimensionless peak optical depth.

## Point and impulse sources

The point-source layer assigns each catalog object to the nearest pixel. Unit amplitudes produce a source-location map; optional amplitudes allow weighted catalogs.

The dust/CIB test layer uses the same operation in temperature units. Its default source amplitude is 4.25 µK. It should be interpreted as a compact test-source layer, not as a complete physical CIB model.

## Localization disks

A small unit disk can be painted around every source location. This is useful for checking catalog-to-map coordinate placement.

## Circular masks

The mask painter fills a disk around every catalog position using a per-object angular radius. Additive and binary-union modes are both available.
