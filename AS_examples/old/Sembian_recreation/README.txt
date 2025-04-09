Example to recreate results, experiments, and simulations from Sembian Phys. of Fluids 2016 (https://pubs.aip.org/aip/pof/article-abstract/28/5/056102/1079453/Plane-shock-wave-interaction-with-a-cylindrical?redirectedFrom=fulltext)

Note:
Sembian uses the Tait EOS for liquid phase, this simulation currently uses Stiffened gas EOS as Tait is not yet implemented. 11/19/2024
The droplet location in relation to the shock from the Sembian experiments and simulations are yet to be confirmed. 


Simulation setup:
High pressure region to model exploding wire facility.
Mesh stretching in x at the beginning of the domain to avoid shock reflections. This can be removed by setting nx_stretch to zero in input file.
Wall boundary conditions on the left, top, and right boundaries (can be changed in input file).
Symmetry boundary condition applied at bottom of domain by applying a slip wall (inviscid simulation with wall at bottom, set droplet center on bottom of domain).
