# Spectral mass construction

Scripts for the hadron-spectrum chapter of *The Cosserat Supersolid* monograph.
Run any script from inside this directory; the modules import their siblings by name.

https://doi.org/10.5281/zenodo.18636501

- `meson_baryon_docking.py`: kaon and rho clusters docked on the nucleon and Delta hosts at every coherent translate and orientation; shows that a compound with no void pair sits at or above its two-body threshold (the Lambda(1405), Sigma(1670) and Delta(1920) compound readings).
- `channel_blocks.py`: the one-dimensional channels of the shell and Born as exact 2x2 blocks [[d+c^2, c],[c, mu+1]] with mu the graph-Laplacian eigenvalue and c the curl norm, equal to one on every channel.
- `static_distortion_map.py`: eigenstrain relaxation of Shockley offsets (three partials at 120 degrees give zero) and each host's mechanisms re-carried axially; the shell and void cluster return the integer rungs 8 and 9, the cap hosts do not.
