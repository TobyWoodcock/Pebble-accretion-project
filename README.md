# Pebble-accretion-project
## Project description
Project supervisor: [Sebastiaan Krijt](https://www.skrijt.com/)

MPhys project to investigate the formation of planets in protoplanetary discs as a result of pebble and gas accretion. Pebble accretion is a mechanism by which a protoplanets gains solid material through the accretion of centimetre-scale 'pebbles'. Gas accretion, on the other hand, is the development of a gaseous envelope around a sufficiently massive protoplanet. 

To model the formation of a planet through pebble and gas accretion, we use [Python](https://www.python.org). We make use of [`pebble-predictor`](https://github.com/astrojoanna/pebble-predictor.git) by [Dr. Joanna Drążkowska](https://www2.mps.mpg.de/homes/drazkowska/Home.html) to predict the populations of pebbles and their sizes in the protoplanetary disc and [`epsilon`](https://staff.fnwi.uva.nl/c.w.ormel/software.html) by [Chris Ormel](https://staff.fnwi.uva.nl/c.w.ormel/index.html) and [Beibei Liu](https://sites.google.com/view/beibei-liu/) to estimate the efficiency of Pebble accretion onto a protoplanet. 

Using these Python modules, we incorporate three main processes
 1. Pebble accretion driven by the drift of small solids from the outer disc which we forecast using `pebble-predictor` as a 'pebble flux'. `pebble-predictor` also predicts the size of the pebble in the disc which is used by `epsilon` to calculate the efficiency of pebble accretion onto a planet
 2. Planetary migration as a result of gravitational forces exerted on the protoplanet by the disc, following the approach of [Nielsen et al. 2023](https://arxiv.org/abs/2308.15504)
 3. Gas accretion through the collapse of a gas cloud surrounding a protoplanet, also following Nielsen et al. 2023

We implement three principle physical phenomena that affect the rate of accretion
- The water snowline, the radius of the disc at which water transitions from being found in liquid form to as solid ice
- The pebble isolation mass, the mass of a protoplanet at which gas cloud in its vicinity is dense enough to deflect incoming pebbles
- The corotation radius, the inner radius of the disc at which gas no longer lies in a central 'midplane', and protoplanets stop migrating

Combining these processes and phenomena, we are able to create a basic picture of a protoplanet's evolution in a protoplanetary disc, beginning from an intial embryo. We use numerical integration to simulate how the mass, orbital radius and density of a protoplanet would change over time as a result of these processes. 

![image](https://github.com/user-attachments/assets/ea947af8-6ff9-49f4-abcc-9eb06061f5e5)

^ an example plot showing the growth and migration of a protoplanet in two different types of discs

## Code
To capture the complexity of astronomical objects, we use the Python class system. Two basic classes are implemented, which could be transferrable to other projects focused on simulating planetary physics, these are the star and planet classes. 
- The star class has basic functionaliy relevant to the project. A star object is created with a mass and period of rotation. The former is used in a method to calculate the Keplerian orbital frequency of any other body in its vicinity, the latter to calculate the corotation radius in a later stage. 
- The planet class has features that allow for a more detailed analysis. A planet object is created with an initial mass, chemical composition, orbital radius and solid fraction. Given the planet's mass, the chemical composition and solid fraction enable the calculation of the average density and radius of its solid core. Given a star, a planet's orbital radius reveals its orbital period. All of these quantities are encoded as properties of the class. Crucially, this means they are calculated upon retrieval from an object.

In the project files, the basic classes are defined in systems.py

(Sidenote: using the properties feature has advantages and drawbacks when it comes to efficiency. It is best used for quantities that requre computationally untensive calculations, and that change frequently, as each time a property is retrieved, its value is calculated again in the 'getter' function. For more involved quantities, it is better to do the calculation before beginning a loop, or at the start of each stage of the loop, and assigning the constant calculated value to the object. This way, the program does not need to repeatedly calculate the same quantity)

Two more specialised classes are now added, the disc class, and the protoplanet class, which inherits functionality from the planet class.
- the disc class has a broad range of attributes, including a central star object, a total mass and a dust-to-gas ratio. These quantities enable the calculation of density and temperature gradients in the disc, which in turn may be used to locate the snowline. Most notably, the disc class has a radial grid and a temporal grid, which `pebble predictor` utilises in order to forecast the pebble flux in the disc. `pebble predictor` returns a two-dimensional array of fluxes over the combined grid
- the protoplanet class adds attributes to the planet class, enabling the simulation of protoplanetary evolution in a disc. Fundamentally, it has a disc object: its host protoplanetary disc. It uses 



![image](https://github.com/user-attachments/assets/a7d4d866-4ae9-41ab-b9f2-543c1d6a5643)


