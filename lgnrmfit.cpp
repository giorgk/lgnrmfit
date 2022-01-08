// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*

    This is an example illustrating the use the general purpose non-linear 
    least squares optimization routines from the dlib C++ Library.

    This example program will demonstrate how these routines can be used for data fitting.
    In particular, we will generate a set of data and then use the least squares  
    routines to infer the parameters of the model which generated the data.
*/


#include <dlib/optimization.h>
#include <iostream>
#include <vector>


using namespace std;
using namespace dlib;

// ----------------------------------------------------------------------------------------

typedef matrix<double,2,1> input_vector;
typedef matrix<double,3,1> parameter_vector;

const double sqrt2pi = std::sqrt(2*pi);
const double sqrtpi = std::sqrt(pi);
// ----------------------------------------------------------------------------------------

// We will use this function to generate data.  It represents a function of 2 variables
// and 3 parameters.   The least squares procedure will be used to infer the values of 
// the 3 parameters based on a set of input/output pairs.
double model (
    const input_vector& input,
    const parameter_vector& params
)
{
    const double a = params(0);
    const double b = params(1);

    const double x = input(0);

    /*
               /               2 \
               |   (a - log(x))  |
    sqrt(2) exp| - ------------- |
               |           2     |
               \        2 b      /
    ------------------------------
            b x sqrt(pi) 2
     */

    const double logx = log(x);
    const double nominator = sqrt_2 * exp(- ((a - logx)*(a - logx))/(2*b*b));
    const double denominator = b * sqrtpi * 2;
    return nominator/denominator;

}

// ----------------------------------------------------------------------------------------

// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double residual (
    const std::pair<input_vector, double>& data,
    const parameter_vector& params
)
{
    return model(data.first, params) - data.second;
}

// ----------------------------------------------------------------------------------------

// This function is the derivative of the residual() function with respect to the parameters.
parameter_vector residual_derivative (
    const std::pair<input_vector, double>& data,
    const parameter_vector& params
)
{
    /*
    wrt a (mean)
             /               2 \
             |   (a - log(x))  |
  sqrt(2) exp| - ------------- | (2 a - 2 log(x))
             |           2     |
             \        2 b      /
- -----------------------------------------------
                   3
                  b  x sqrt(pi) 4


    wrt b (standard deviation)
               /               2 \
               |   (a - log(x))  |        2                 2    2
    sqrt(2) exp| - ------------- | (log(x)  - 2 a log(x) + a  - b )
               |           2     |
               \        2 b      /
    ---------------------------------------------------------------
                             4
                            b  x sqrt(pi) 2
    */
    parameter_vector der;

    const double a = params(0);
    const double b = params(1);

    const double x = data.first(0);
    const double logx = log(x);
    const double sqb = b*b;
    const double exp_fraction = sqrt_2 * exp(- ((a - logx)*(a - logx))/(2*sqb));

    const double numeratorA = exp_fraction * (2 * a - 2 * logx);
    const double denominatorA = b*sqb * sqrtpi * 4;

    const double numeratorB = exp_fraction * (logx*logx - 2 * a * logx + a*a - sqb);
    const double denominatorB = sqb*sqb * x * sqrtpi * 2;

    der(0) = numeratorA/denominatorA;
    der(1) = numeratorB/denominatorB;

    return der;
}

// ----------------------------------------------------------------------------------------

bool read_inputfile(std::string filename, std::vector<double>& data, std::string& errmsg){
    bool out = false;
    std::ifstream datafile(filename.c_str());
    if (!datafile.good()) {
        errmsg = "Can't open the file" + filename;
    }
    else{
        std::string line;
        double d;
        while (getline(datafile, line)){
            if (line.size() > 0){
                std::istringstream inp(line.c_str());
                inp >> d;
                data.push_back(d);
            }
        }
        out = true;
    }
    return out;
}

int main(int argc, char* argv[])
{
    std::string infile = argv[1];
    std::vector<double> data;
    std::string msg;
    read_inputfile(infile, data, msg);
    try
    {
        // randomly pick a set of parameters to use in this example
        const parameter_vector params = 10*randm(3,1);
        cout << "params: " << trans(params) << endl;


        // Now let's generate a bunch of input/output pairs according to our model.
        std::vector<std::pair<input_vector, double> > data_samples;
        input_vector input;
        for (int i = 0; i < 1000; ++i)
        {
            input = 10*randm(2,1);
            const double output = model(input, params);

            // save the pair
            data_samples.push_back(make_pair(input, output));
        }

        // Before we do anything, let's make sure that our derivative function defined above matches
        // the approximate derivative computed using central differences (via derivative()).  
        // If this value is big then it means we probably typed the derivative function incorrectly.
        cout << "derivative error: " << length(residual_derivative(data_samples[0], params) - 
                                               derivative(residual)(data_samples[0], params) ) << endl;





        // Now let's use the solve_least_squares_lm() routine to figure out what the
        // parameters are based on just the data_samples.
        parameter_vector x;
        x = 1;

        cout << "Use Levenberg-Marquardt" << endl;
        // Use the Levenberg-Marquardt method to determine the parameters which
        // minimize the sum of all squared residuals.
        solve_least_squares_lm(objective_delta_stop_strategy(1e-7).be_verbose(), 
                               residual,
                               residual_derivative,
                               data_samples,
                               x);

        // Now x contains the solution.  If everything worked it will be equal to params.
        cout << "inferred parameters: "<< trans(x) << endl;
        cout << "solution error:      "<< length(x - params) << endl;
        cout << endl;




        x = 1;
        cout << "Use Levenberg-Marquardt, approximate derivatives" << endl;
        // If we didn't create the residual_derivative function then we could
        // have used this method which numerically approximates the derivatives for you.
        solve_least_squares_lm(objective_delta_stop_strategy(1e-7).be_verbose(), 
                               residual,
                               derivative(residual),
                               data_samples,
                               x);

        // Now x contains the solution.  If everything worked it will be equal to params.
        cout << "inferred parameters: "<< trans(x) << endl;
        cout << "solution error:      "<< length(x - params) << endl;
        cout << endl;




        x = 1;
        cout << "Use Levenberg-Marquardt/quasi-newton hybrid" << endl;
        // This version of the solver uses a method which is appropriate for problems
        // where the residuals don't go to zero at the solution.  So in these cases
        // it may provide a better answer.
        solve_least_squares(objective_delta_stop_strategy(1e-7).be_verbose(), 
                            residual,
                            residual_derivative,
                            data_samples,
                            x);

        // Now x contains the solution.  If everything worked it will be equal to params.
        cout << "inferred parameters: "<< trans(x) << endl;
        cout << "solution error:      "<< length(x - params) << endl;

    }
    catch (std::exception& e)
    {
        cout << e.what() << endl;
    }
}
