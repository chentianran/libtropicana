#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "trop.hh"
#include "debug.hh"

namespace py = boost::python;
namespace np = boost::python::numpy;

long long compute_volume (const np::ndarray& S)
{
    np::dtype ftype = np::dtype::get_builtin<double>();
    if (ftype != S.get_dtype()) {
        return -1;
    }

    int m = S.shape(0);         // n.o. rows
    int n = S.shape(1);         // n.o. columns
    int row_s = S.strides(0);    // row stride
    int col_s = S.strides(1);       // column stride
    char* data = S.get_data();      // raw data pointer

    RowMatrix A (m, n);

    for (int i = 0; i < m; ++i) {
        char* row = data + (i * row_s);       
        for (int j = 0; j < n; ++j) {
            char*   rp = row + (j * col_s);
            double* fp = reinterpret_cast<double*>(rp);
            A (i, j) = *fp;
        }
    }

    Trop trop;
    trop.compute (A);

    return trop.volume;
}

BOOST_PYTHON_MODULE(pytropicana)
{
    Py_Initialize();
    np::initialize();

    py::def("compute_volume", compute_volume);
}
