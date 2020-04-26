import numpy as np
import pytropicana as trop

A = np.array ( 
    [
        [ -1.0, -1.0,  0.0 ],
        [  1.0, -1.0,  0.0 ],
        [ -1.0,  1.0,  0.0 ],
        [  1.0,  1.0,  0.0 ],
        [  0.0,  0.0,  1.0 ]
    ], 
    dtype=float )

print (trop.compute_volume(A))