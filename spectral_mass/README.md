# Spectral mass construction

Scripts for the hadron-spectrum chapter of *The Cosserat Supersolid* monograph.
Run any script from inside this directory; the modules import their siblings by name.

https://doi.org/10.5281/zenodo.18636501

- `meson_baryon_docking.py`: kaon and rho clusters docked on the nucleon and Delta hosts at every coherent translate and orientation; shows that a compound with no void pair sits at or above its two-body threshold (the Lambda(1405), Sigma(1670) and Delta(1920) compound readings).
- `channel_blocks.py`: the one-dimensional channels of the shell and Born as exact 2x2 blocks [[d+c^2, c],[c, mu+1]] with mu the graph-Laplacian eigenvalue and c the curl norm, equal to one on every channel.
- `static_distortion_map.py`: eigenstrain relaxation of Shockley offsets (three partials at 120 degrees give zero) and each host's mechanisms re-carried axially; the shell and void cluster return the integer rungs 8 and 9, the cap hosts do not.
- `kk_cluster_matrix.py`: the Kaluza-Klein reduced 4D Cosserat cluster matrix (u, u_4, phi, phi_i4) at compact wavenumber n = 0, +/-1; factorisation check at n = 0 and the compact rung lambda = 3 + boundary count that carries the accommodated hyperons at N(m_0 - m_e).
- `d4_shell_cluster.py`: the coordination shell as a D_4 cluster with all 24 bonds; severed inter-layer bonds reproduce the slice A_2u pair exactly, clamped ones would put the nucleon at 1050 MeV.
