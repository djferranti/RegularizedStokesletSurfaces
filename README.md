# Regularized Stokeslet Surfaces
This is a library that implements the regularized Stokeslet surfaces method as described in the Journal of Computational Physics paper https://www.sciencedirect.com/science/article/pii/S0021999124002535.
A preprint of the paper is available https://arxiv.org/abs/2310.14470 with minor expository changes at the end. 

# Code basics
To understand how to use the code, read through the examples in the **tests** folder. There are four example scripts:  

1. *triangleStructTest.m* : evaluates the basic geometric data and stores in a Triangle struct for a particular triangle. The triangle struct data is a basic component of how the code works so it is important to know what is needed. 
2. *testTbaseCases.m*: Creates a random set of points and evaluates the integrals used in the method analytically and compares to the quadrature. For most points, these will agree well. However, for points on the triangle or in the interior, there will be a large discrepancy. This occurs because of the regularized Stokeslet kernel is nearly singular, i.e. it is highly peaked near the force points. To demonstrate this, we add to the list of evaluation points the vertices of the triangle.
3. *testTranslatingSphere.m* : Used the method to evaluate the velocity of a triangulated sphere in Stokes flow given forces that should approximately impose a unit velocity in the x-direction.
4. *testDragSphere.m* : Does the inverse of the previous script. Calculates the drag on a triangulated sphere in Stokes flow given an imposed velocity field (unit in x-direction) at the triangle vertices.

5. The last two scripts use a Delaunay triangulation of the sphere from John Burkardt's code (see **sphere_delaunay**) distributed under the LGPL license https://www.gnu.org/licenses/lgpl-3.0.en.html. See John Burkardt's website at
6. https://people.math.sc.edu/Burkardt/m_src/sphere_delaunay/sphere_delaunay.html for more information.

# Future development 
More test scripts may be added, but this will likely not change much in the near future. Feel free to reach out to me (dferranti--at--wpi--dot--edu) if you have questions.
