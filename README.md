# soft spheres in a line

The mass-spring model code is built upon the N-body code rebound, 
see https://github.com/hannorein/rebound

Slightly edited rebound code is in src/

Spring routines (not part of rebound) are in src_spring/

To find something to run go into myexamples/line
>  make 
>  ./rebound_spring


There is a parameter file b1.  To run the code with the parameter file
>  ./rebound_spring b1
