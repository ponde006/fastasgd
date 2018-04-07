# fastasgd
This code uses a variant of the fastsg bijection to calculate Clenshaw-Curtis type adaptive sparse grids.

In order to use this program, include it in your Julia session using
include("[your filepath]/fastasgd.jl")
using fastasgd

Basic usage is as follows

```julia
d = 2 # Set the dimensions of the function
states = 1 # Set the number of discrete states, a new grid will be generated for each state
ord = 11 # Set the maximum depth of approximation
lb = [-1.0, -1.0] # Set the lower bound for each dimension
ub = [1.0, 1.0] # Set the upper bound for each dimension

myFunc = ASG(d, states, ord, lb, ub); # Initialize the sparse grid object
# Generate the adaptive sparse grid hierarchal surpluses
# Note that the function must be vector valued. The third argument controls local adaptivity
# The fourth argument controls spatial adaptivity (see Jakeman et. al. "Local and Dimension Adaptive Sparse Grid Interpolation and Quadrature"
# The final argument controls which discrete state is being updated. In this example there is only 
# one discrete state so it is always 1
@time adaptivegrid( myFunc, x -> (x[1] + x[2]/2.0)^2, 0.01, 0.01, 1);

```

Once the interpolation object is created we can evaluate it by calling
```julia
@show myFunc(0.5, 0.66, 1)
```
where again, the last argument controls which discrete state is evaluated. We can also evaluate the derivative of each argument using the following expression
```julia
@show myFunc(1, 0.5, 0.66, 1) # Derivative with respect to the first argument
@show myFunc(2, 0.5, 0.66, 1) # Derivative with respect to the second argument
```

Finally, I have included the function iinterpolate to evaluate the integral of our function over the relevant domain
```julia
u1 = [1.0, 1.0]
@show myFunc(u1, myFunc, myFunc.l, 1) # The third arugment controls the depth of the approximation. Using myFunc.l makes it
# the maximum depth, i.e. myFunc.l = ord above.
```
