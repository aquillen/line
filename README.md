# soft spheres in a line

This does a bunch of soft spheres with downward gravity

The mass-spring model code is built upon the N-body code <i> rebound </i>, 
see https://github.com/hannorein/rebound

Slightly edited and older (not up to date) rebound code is in src/

Spring routines (not part of rebound and used for our mass/spring model spin stuff) are in src_spring/

To find something to run go into myexamples/line
>  make 

>  ./rebound_spring


There is a parameter file b1.  To run the code with the parameter file
>  ./rebound_spring b1


I found that you need to install the glfw or glfw3 library for the opengl display to work

