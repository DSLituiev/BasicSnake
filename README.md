Snake : Active Contour
======================
by Dirk-Jan Kroon
08 Jul 2010 (Updated 08 Nov 2011)

Snake Segmentation (Kass et al), 2D/3D including GVF and balloon force (Easy to read code)

These functions implements the basic snake segmentation contour, as introduced by Kass et al.

## Algorithm:
A snake is an active (moving) contour, in which the points are attracted by edges and other image boundaries. To keep the contour smooth, a membrane and thin plate energy is used as contour regularization.

## Implementation:
All code is well commented, and is probably easy readable. It includes,
- A balloon force, which is not in the original paper
- Binary image segmentation result,
- The "Gradient Vector Flow" (GVF) method.
- 3D implementation, including parameter to keep the mesh from self-intersecting

Try the example in the help of function Snake.m !

## Literature:
- Michael Kass, Andrew Witkin and Demetri TerzoPoulos "Snakes Active Contour Models", 1987
- Jim Invins and John Porril, "Everything you always wanted to know about snakes (but were afraid to ask)
- Chenyang Xu and Jerry L. Prince, "Gradient Vector Flow: A New external force for Snakes
- Christoph Lurig, Leif Kobbelt, Thomas Ertl, "Hierachical solutions for the Deformable Surface Problem in Visualization"

Please report bugs, successes and other comments 
