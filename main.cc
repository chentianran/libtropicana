#include <iostream>
#include <algorithm>
#include <utility>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include "trop.hh"
#include "debug.hh"


using namespace std;

bool is_this_empty (const std::string& str) {return str.empty();}

vector<double> get_vector (const std::string& line)
{
    vector<double> v;
    vector<string> tokens;
    boost::split (tokens, line, boost::is_any_of("\t "), boost::token_compress_on);
    vector<string>::iterator new_end = remove_if (tokens.begin(), tokens.end(), is_this_empty);
    std::transform (tokens.begin(), new_end, back_inserter(v), boost::bind (boost::lexical_cast<double,string>,_1));
    return v;
}

int main ()
{
    try {

        vector< vector<double> > rows;
        int n = 0;                                                       // keep track of maximum size of each line
        for (string line; getline (cin, line);) {                            // for each line in the stream
            if (line.empty())
                continue;
            vector<double> row = get_vector (line);                                //   turn each line into a Term object
            n = std::max (n, static_cast<int>(row.size()));                                  //   update max line size
            rows.push_back (row);                                         //   add the term to the list
        }

        int m = rows.size();
        RowMatrix A (m, n);

        LOGVAR(1,m);
        LOGVAR(1,n);

        for (int i = 0; i < m; ++i)                                  // for each term we created
            for (int j = 0; j < rows[i].size(); ++j)
                A(i,j) = rows[i][j];

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
