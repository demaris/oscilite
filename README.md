# oscilite
David DeMaris
demaris@well.com

davidl.demaris@gmail.com

This package contains two types of discrete dynamical simulation classes which can be made to interact, and which serve as models of
cooperative and competitive spatial dynamical  systems, and requires only a scientific python distribution. The 'lite'
designator indicates that many contemplated components may be available in a later 'heavy weight' distribution.

There are two testing and visualization demo environments provided to give some intutition for the behavior of coupled map lattice
(CML) dynamics.

python demo_diffusive_cml.py

python demo_competitive_cml.py

CML are abstract complex systems models originally studied within the computational physics and nonlinear dynamics community,
and may have applicability to neuroscience, social community formation, and ecological modeling.  In my masters and Ph.D. work
at University of Texas Austin EE department, a matlab version of the cooperative dynamics code was used for both unsupervised and
supervised computer vision tasks, and modeling of psychophysics  (Necker Cube, Muller Lyer Illusions).
This work was informed by emerging theories of oscillations, balanced networks, and cross frequency coupling, and presented at conferences including INNS and Computational Neuroscience, and in particular the then small community interested in oscillations and neural field theory. This work was unusual in introducing nonstationarity (time varying, dynamic during inference) networks, spatio-temporal interactions, and a view of perception as recurrent activity in short frames where the representation is statistical and coded over the full range of ensemble activity rates, rather than high frequency (activity) units.  
Several papers are available here. 

https://daviddemaris.com/pubs/

In supervised and semi-supervised learning, coarsed grained representations of lattice activity are used to adjust a small number of
parameters applied everywhere in a field; this procedure recently been described as downward causation.

https://royalsocietypublishing.org/doi/full/10.1098/rsta.2016.0338

But more straightforwardly representations are points on a manifold of the ensemble activty. 

With respect to neuroscience, the individual sites or cells can be considered as abstractions of a local E-I network
or pool, then coupled to its neighbors. It can be considered an efficient neural field or mass action model, in the sense that spikes
are not treated, but only ensemble average frequencies about a default rate. The state variables model a recurrent pool in a column
High coupling increases the synchronization over the network, where synchronization is used as in the statistical physics literature - 
high synchronization means the distribution is focused (possibly multi-modally), low synchronization means that it is dispersed. 
There is no notion of phase synchronization per se, though there may be phase synchronization at the level of spatiotemporal pattersn 

The diffusive class implements a hybrid  of diffusively coupled map lattice (CML) and globally coupled map (GCM)
and some statistical utilities, including an ising model with a momentum based conception of spin of each lattice site.
The spin field can be displayed, and measures of frustration used to control the dynamics to reach an equilibrium.
Some work at IBM in the 80s used this self organizing annealing-like method to discriminate defect images. 

In addition the diffusive CML implementation supports two forms of control via masks.  There is 'pinning', where some states are
fixed while others vary according to local dynamics, allowing the pinned states to act as 'lattice defects' which may
control the dynamics.  There is also an activation mask and an example given of lateral spreading dynamics.

By setting the local coupling to zero, the hybrid map reduces to globally coupled map (GCM) formulations.
Conversely, setting global coupling to zero reduces the dynamics to the CML formulation.
While most of the literature has a fixed setting for the nonlinearity parameter alpha, the simulator and GUI supports
alpha being a field variable, with a range of values from low (current alpha setting)

In addition you can observe the effects of cycling the alpha, local and global coupling parameters. In earlier work
on image and pattern recognition using these systems, these parameters are varied to perform unsupervised clustering,
or via semi-supervised learning using a genetic algorithm to search the paramter space. 

The competitive_cml class implements a competitive CML where sites with high state values supress neighbors.

It is intended that you have a full scipy distribution, but this code relies only on scipy.signal, scipy.stats,
and matplotlib for display and user interface controls.

The diffusive cml can be considered as simulation of columns of balanced networks of excitatory and inhibitory neurons.
The statistics of globally coupled maps resembles a small world connectivity matrix between sites.

The combination of recurrent dynamics and convolution may create mapping manifolds similar to deep feedforward networks,
with subspaces visited by the dynamics acting as the nodes in hidden layers and as output units. In my image recognition work,
I searched for dynamics which would create similar distributions across subspaces (partitions) for different views of an object,
with a distance function comparing distribution values similar to comparision nodes in Siamese networks.
The dynamics were constrained to a few iterations, with time varying nonlinearity and coupling parameters.
The working hypothesis was that slow oscillations (alpha, delta) serve as nonlinear control parameters for fast dynamics in the gamma range,
with the convolution and coupling occuring periodically in low gamma. 
Learning was performed with genetic algorithms to search parameters for each class, with a family of recognizers
creating the overall manifold.  It was shown that even though the GA objective (loss) functions only handled invariance over views and
maximizing distance from other representations via KL distance term, projections into a low dimensional space exhibited clustering
of the recognized objects according to the number of jointed segments.

A Brief History of Excitable Map-Based Neurons and Neural Networks
https://arxiv.org/pdf/1303.0256.pdf
