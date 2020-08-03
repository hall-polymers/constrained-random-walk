# Constrained-random-walk
Python script to generate unconfined or confined random walks (RW) in the presence of walls. The fundamental assumption is that when the RW would pass through a wall, it instead reflects off

## Motivation
<p align="center">
	<img src="demo/motivation.png"  width="580" height="311"/>
</p>

For now, it is also assumed that the box is a unit cube with appropriate periodic boundary conditions
(or equivalently, we only start in the unit cube, but the system expands into all space with the appropriate symmetry)

input paramters are 
 ```M``` - number of trials 
 ```N``` - number of steps in each trial
 ```L``` - boxsize (set relative to the unit step size)
 ```constraint``` - 
   ```none```: no constraint 
   ```lamellae```: assumes that the z=0 and z=L planes are walls 
   ```cylinder```: assumes that the wall is a cylinder surface along z-dir, and its center is (L/2, L/2) and radius is 2/L
   ```gyroid```: assumes that the wall is the surface of a double gyroid unit cell

optional parameters
  ```x0, y0, z0``` - the inital location of the RW, a negative number will choose random location between 0 and L FIXME: does not check that these choices are reasonable, so an infinite loop can happen if you force the system to start inside a particle
  ```outputchains``` - boolean to determine if printing to lammpstrj format

<p align="center">
	<img src="demo/crw_demo.gif" width="300" height="317"/>
</p>