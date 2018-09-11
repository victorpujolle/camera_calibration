**This code is used to calibrate a camera using 4 markers placed on the corners
of a 10cm 10cm square**


You need to provide a csv file (../Calib_list/1.csv here) which contains some data
used for the calibration.

- 1:8 	: detected marker corner pixels (clock wise) (x,y)x4
- 9:11	: expected marker center position
- 12:13	: focal length x,y in pixels
- 14:16	: camera focal point
- 17:19	: camera angle, degrees (z1,x2,z3) --> R=Rz1.Rx2.Rz3, ("-z" direction is looking focal point)
- 20	: camera distance to focal point

First the line where no markers were detected are removed.

Then a parameters struct is created, it is used to store the variables used in different functions.

Then the algo tries to minimize the cost function which returns the standard deviation of pixel offsets.
The minimization algo is nmsimplex. It is a non-derivative optimization algorithm. It is really sensitive to local
minimums, this is why we launch it many times.

At the end the function optimizer calculate two vectors u,v and some other parameters linked with u and v
u and v are the deviation between calculated positions and detected position
