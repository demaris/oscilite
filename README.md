# cml_mini

This package contains two types of dynamical simulation kernels which can be combined.  They are abstract models studied within the
computational physics and nonlinear dynamics community, and may have applicability to neuroscience and ecological modeling.

The diffusiveCML class implements a hybrid diffusively coupled map lattice (CML) and globally coupled map (GCM) and some statistical utilities,
including an ising model with a momentum based conception of spin of each lattice site.  The spin variable can be displayed and
measures of frustration used to control the dynamics to reach an equilibrium.

The competitiveCML class implements a competitive CML where sites with high values.

By setting the local coupling to zero, the hybrid map reduces to globally coupled map (GCM) formulations.
Conversely, setting global coupling to zero reduces the dynamics to the CML formulation.
While most of the literature has a fixed setting for the nonlinearity parameter alpha, the simulator and GUI supports
alpha being a field variable, with a range of values from low (current alpha setting)

In addition you can observe the effects of cycling the alpha, local and global coupling parameters. In earlier work
on image and pattern recognition using these systems I varied these parameters.

It is intended that you have a full scipy distribution, but really relies only on scipy.signal, scipy.stats,
and matplotlib for display and user interface controls.

The diffusive cml can be considered as simulation of columns of balanced networks of excitatory and inhibitory neurons.
The statistics of globally coupled maps resembles a small world connectity matrix between sites.

The combination of recurrent dynamics and convolution may create mapping manifolds similar to deep feedforward networks,
with subspaces in the dynamics acting as the nodes in hidden layers and as output units. In my image recognition work,
I searched for dynamics which would create similar distributions across subspaces (partitions) for different views of an object,
with a distance function comparing distribution values similar to comparision nodes in Siamese networks.
The dynamics were constrained to a few iterations, with time varying alpha and gamma parameters.


