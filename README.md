# Project 1: Brush

## Design
For the design of this project, I placed all my fields and methods
inside `canvas2d.cpp`. 

In terms of the general implementation,
everything pivots around the `paintBrush` function. This takes a pair of coordinates
and paints there by iteratively editing the RGBA values in the radius around there.

`changeMask` regenerates the mask vector everytime the `brushType`, `brushRadius`, or `brushDensity` (for spray brush)
is changed. It takes a `(float -> float)` function, which gives it the intensity values
at each pixel on the mask vector. For efficiency, purposes, the mask vector is not
regenerated if only the `brushColour` is changed, as that would be unnecessary.

The smudge brush works by "picking up" the canvas pixels at each point, and storing them in a vector
`smudge_data`. Then, whenever the mouse is dragged, it puts these pixels down at the next location, and
picks up new pixels, etc. etc.

## Extra Credit
I implemented the spray brush, which works as follows: before painting the pixel, a random number
between 1 and 400 (inclusive) is generated. If the random number is less than `brushDensity`, then the 
intensity remains at 1. If not, the intensity becomes 0 (i.e. that pixel is not painted).


# Project 2: Filter
## Design
For the design of the project, I placed all my fields and methods in `canvas2d.cpp`.

For the design of the project, I created methods for each filter type, and had `filterImage()` 
call one depending on the user selection.

`filterEdgeDetect()` implements the Sobel edge detection. It first converts the image to grayscale,
and then convolves the separated Sobel kernels in the X and Y directions, first horizontally and then
vertically. These gradients are stored in separate vectors, and finally `setGradient()`
computes the magnitude of the gradients and modifies the canvas to show this.

`filterBlur()` implements the blur, using a triangle filter. First it generates a separated 1D blur kernel based on the
user's selected radius, then it calls the general-purpose functions `convolveX()` and `convolveY()` to convolve by the kernel in each direction.

`filterScale()` implements the scaling. It calls `scaleXDirection()` and `scaleYDirection()` sequentially, and scales by the desired factors in each
direction. Each of these depends on the `triangleWeight()` function, which provides the weight of a pixel according to the triangle filter; `backMap()`, which
takes an output coordinate and computes its equivalent in the input image; and the `h_prime_x()` and `h_prime_y()` functions, respectively, which computes the colour
of an output pixel. Unfortunately there is some redundancy between `h_prime_x()` and `h_prime_y()`, as each one accesses pixels differently. I tried to have it take a function
as an argument, which would help it decide how to get pixels, but unfortunately Qt kept complaining and I couldn't do it.

## Extra Credit

I implemented the median filter in `filterMedian()`, which iterates through the image pixels and calculates the median for each colour channel with `calculateMedian()`.
This, in turn, calls `insertSorted()` to insert the new pixels in a sorted manner, so that the median can be easily indexed.
