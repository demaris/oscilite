# cmlMinimal

This code implements a hybrid diffusively coupled map lattice (CML) and globally coupled map (GCM).
It is intended that you have a full scipy distribution, but really relies only on scipy.signal and scipy.stats
This is a macro-scale simulation of balanced networks of excitatory and inhibitory neurons, with a small world network hopping between
(small world dynamics are shown to be similar to GCM).
This code is currently under test and unsupported. Stay tuned!
It is roughly working, but some values of slider settings will cause a numerical problem, numpy warnings, and graphics cease to function.

