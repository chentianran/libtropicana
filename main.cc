#include "trop.hh"

#include <iostream>

using namespace std;

int main ()
{
    try {
        const int W = 2, H = 1;
        RowMatrix A ((W+1)*(H+1),2);
        for (int i = 0; i <= H; ++i) {
            for (int j = 0; j <= W; ++j) {
                A.row(i*(W+1)+j)(0) = i;
                A.row(i*(W+1)+j)(1) = j;
            }
        }

        cout << A << endl;

        Trop trop;
        trop.compute (A);

        cout << "Volume: " << trop.volume << endl;

        //cout << "B*D = \n" << (B*D) << endl << endl;
    } catch (const char* err) {
        cerr << err << endl;
        return 1;
    }
    return 0;
}
