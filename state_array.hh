#include <Eigen/Core>

class ArrayState
{
protected:

    double* data;
    int n;

public:

    ArrayState (int dim) 
    { 
	assert (dim > 0);

	n = dim;
	data = new double [dim * dim];

	for (int j = 0; j < n; ++j)
	    for (int i = 0; i < n; ++i)
		data [i + j*n] = (i == j) ? 1.0 : 0.0;
    }

    void pivot (int k, const double* v)
    {
	assert (k >= 0 && k < n);

	double r = 0.0;
	
	for (int i = 0; i < n; ++i)
	    r += v[i] * data [i + k*n]; 
	for (int i = 0; i < n; ++i)
	    data [i + k*n] /= r;
	for (int j = 0; j < n; ++j) {
	    if (j != k) {
		r = 0.0;
		for (int i = 0; i < n; ++i)
		    r += v[i] * data [i + j*n];
		for (int i = 0; i < n; ++i)
		    data [i + j*n] -= r * data [i + k*n];
	    }
	}
    }
};

