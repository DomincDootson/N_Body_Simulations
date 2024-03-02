# N-Body Simulations 
This file contains all the code that I wrote to study the evolution of 2D disc using N-body methods. This N body takes a particle perturbation approach: where a perturbation and unperturbed simulation are run, and the different is outputted. This approach allows for direct comparison to my Linear response work. 

## File Overviews
### Bodies
This file contains two important files:
- `Bodies.h`: This class contains all the information about a particle and the methods to calculate its evolution
- `particleSampling.cpp`: This file samples the particles and writes them to file. 

### Box
This file contains the methods to calculate the gravitational field sourced by the particles.

### N Body Classes
The main N body class, `NBody.h`, holds all the particles and perturbation fields and is in charge of running the experiment. The rest of the files are devoted to specific implimentations with different perturbations, e.g. making sure that the bar is correctly updated at each time step.

### Orbit Sections 
Code to calculate the surface of a section after a simulation

### Perturbation Grids 
To study the effects of a perturbation to the disc, we have a virtual class, `PerturbationGrid`, that calculates the extra forces felt by the particles. The other files in this folder calculate implement the methods to calculate the force that is felt with given perturbations. 