# fetch-limited-model
Enables user to compute fetch at cell centers of an Unstructured Grid. Utilizes stompy to create and UnstructuredGrid object. Currently, set up to only interact with UnTRIM08 subclass of UnstructuredGrid, but can be easily generalized.

Takes thetaw (direction of wind source CCW and relative to the positve x-axis. This is not typical, so be cautious. Will fix in the future.

Generally, the cartesian convention is CCW relative to the positive x-axis but requires the direction of the wind vector (ie 180 deg clockwise of what is herein requested). Another convention is the nautical convention which requests the direction of the wind source; however, it is CW relative to the positive y-axis (more specifically, relative to True North). 
