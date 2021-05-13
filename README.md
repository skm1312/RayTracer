## Parallel RayTracer

The code presented performs ray tracing and rendering spheres. This was a project created for the CS 802 - Parallel Algorithm class. 
We can work with different shapes as well. The files presented in the work are Parallel.cpp, Parallel_RT.cpp, and a header file to support the workflow named geometry.h. 
(The references for the code is mentioned in the respective files)

## How to run the code?

I have used a Cygwin compiler to run the code so that we have a system-based performance idea for the working of the code. We have used the OpenMP API to implement parallelism
For compilation we need to run:
Parallel: g++ -fopenmp Parallel.cpp –o parallel
Parallel Run time: g++ -fopenmp Parallel_RT.cpp –o parallelrt

For execution, we need to run: (we can also prepend the execution statements with the time keyword so that we obtain the actual runtime of processors)
Parallel: time ./parallel
Parallel Run time: time ./parallelrt

We obtain a '.ppm' file as output for this program which can be converted to a '.png' or '.jpg' file. Or simply use a PPM viewer to verify the rendering.

## Understanding and Conclusion

As simple serial ray tracer runs in 14868.5 milliseconds. The time taken for it to run in parallel is equal to 3240.6 milliseconds. 
By observing this we can see that the parallel is running in 11628 milliseconds less than the serial code. This code satisfies Amdahl’s law that states even if we provide too many threads or processes the speedup will increase only slowly and saturates at around 20X speedup. (In this case, it happens around the thread value of 7).
