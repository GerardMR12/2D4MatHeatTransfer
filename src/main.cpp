#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

// Variables -------------------------------------------------------------
// Define all the variables, including the j object to read the json file
json j;
int i, t, N, M, max_iter, wrong_time_iter;
double T0, mytime, V, VP, dx, dy, delta, rf, dt, ti, tf, nt, L, H, W, TwBottom, TwRight_0, TwRight_slope, qwTop, Text, alpha_ext, Se, Sw, Sn, Ss, SwallH, SwallL;
vector<double> p1, p2, p3, rho, cp, lambda, xP, yP;
vector<vector<int>> mat;
vector<vector<double>> aP, aE, aW, aN, aS, bP;
vector<vector<vector<double>>> T;

// Functions -------------------------------------------------------------
// Overload the - operator for matrices
vector<vector<double>> operator-(const vector<vector<double>>& a, const vector<vector<double>>& b) {
    // Create the result matrix
    vector<vector<double>> result = a;

    // Subtract the matrices
    for(size_t i = 0; i < a.size(); i++) {
        for(size_t j = 0; j < a[i].size(); j++) {
            result[i][j] -= b[i][j];
        }
    }

    return result;
}

// Overload the + operator for matrices
vector<vector<double>> operator+(const vector<vector<double>>& a, const vector<vector<double>>& b) {
    // Create the result matrix
    vector<vector<double>> result = a;

    // Add the matrices
    for(size_t i = 0; i < a.size(); i++) {
        for(size_t j = 0; j < a[i].size(); j++) {
            result[i][j] += b[i][j];
        }
    }

    return result;
}

// Overload the * operator for matrices times scalars
vector<vector<double>> operator*(const vector<vector<double>>& a, const double& b) {
    // Create the result matrix
    vector<vector<double>> result = a;

    // Multiply the matrix by the scalar
    for(size_t i = 0; i < a.size(); i++) {
        for(size_t j = 0; j < a[i].size(); j++) {
            result[i][j] *= b;
        }
    }

    return result;
}

// Absolute value of matrix elements stored in a vector
vector<double> abs_mat_to_vec(const vector<vector<double>>& a) {
    // Create the result matrix
    vector<double> result;

    // Calculate the absolute value of the matrix
    for(size_t i = 0; i < a.size(); i++) {
        for(size_t j = 0; j < a[i].size(); j++) {
            result.push_back(abs(a[i][j]));
        }
    }

    return result;
}

// Get the correct lambda
double harmonic_lambda(int a, int b, string position) {
    // Check the position and return the correct lambda
    if(position == "E") {
        return 2*lambda[mat[a][b]]*lambda[mat[a][b+1]]/(lambda[mat[a][b]] + lambda[mat[a][b+1]]);
    } else if(position == "W") {
        return 2*lambda[mat[a][b-1]]*lambda[mat[a][b]]/(lambda[mat[a][b-1]] + lambda[mat[a][b]]);
    } else if(position == "N") {
        return 2*lambda[mat[a][b]]*lambda[mat[a-1][b]]/(lambda[mat[a][b]] + lambda[mat[a-1][b]]);
    } else if(position == "S") {
        return 2*lambda[mat[a+1][b]]*lambda[mat[a][b]]/(lambda[mat[a+1][b]] + lambda[mat[a][b]]);
    } else {
        cout << "Invalid position!" << endl;
        exit(0);
        return 0.0;
    }
}

// Actions ---------------------------------------------------------------
// Assign the constants from the json file to the variables
void set_constants() {
    // Read the file
    ifstream json_file("src/data.json");
    json_file >> j;

    // Set the Geometry
    p1.resize(2, 0.0);
    p2.resize(2, 0.0);
    p3.resize(2, 0.0);
    p1 = j["Geometry"]["p1"].get<vector<double>>();
    p2 = j["Geometry"]["p2"].get<vector<double>>();
    p3 = j["Geometry"]["p3"].get<vector<double>>();

    // General dimensions
    L = p3[0];
    H = p3[1];

    // Set the Material Properties
    rho.resize(4, 0.0);
    cp.resize(4, 0.0);
    lambda.resize(4, 0.0);
    for(i = 0; i < 4; i++) {
        rho[i] = j["Material Properties"][i]["rho"];
        cp[i] = j["Material Properties"][i]["cp"];
        lambda[i] = j["Material Properties"][i]["lambda"];
    }

    // Set the Boundary Conditions
    TwBottom = j["Boundary Conditions"]["TwBottom"];
    TwRight_0 = j["Boundary Conditions"]["TwRight"][0];
    TwRight_slope = j["Boundary Conditions"]["TwRight"][1];
    qwTop = j["Boundary Conditions"]["qwTop"];
    Text = j["Boundary Conditions"]["Text"];
    alpha_ext = j["Boundary Conditions"]["alpha_ext"];

    // Set the Control Volumes
    N = j["Control Volumes"]["N"];
    M = j["Control Volumes"]["M"];
    dx = L/M;
    dy = H/N;

    // Set the Time Parameters
    dt = j["Time Parameters"]["dt"];
    ti = j["Time Parameters"]["ti"];
    tf = j["Time Parameters"]["tf"];
    nt = (tf-ti)/dt;

    // Set the Solver Parameters
    T0 = j["Solver Parameters"]["T0"];
    rf = j["Solver Parameters"]["Relaxation Factor"];
    max_iter = j["Solver Parameters"]["Maximum Iterations"];
    delta = j["Solver Parameters"]["delta"];
}

// Find the material properties for each node
void find_materials() {
    // Resize the matrix
    mat.resize(N+2, vector<int>(M+2, 0));

    // Find the materials
    for(int i = 0; i < N+2; i++) {
        for(int j = 0; j < M+2; j++) {
            if(xP[j] <= p1[0] && yP[i] <= p1[1]) {
                mat[i][j] = 0;
            } else if(xP[j] > p1[0] && yP[i] <= p2[1]) {
                mat[i][j] = 1;
            } else if(xP[j] <= p2[0]) {
                mat[i][j] = 2;
            } else {
                mat[i][j] = 3;
            }
        }
    }

    // Create the materials file
    ofstream mater;
    mater.open("output/materials.txt");

    // Print the materials
    for(int i = 0; i < N+2; i++) {
        for(int j = 0; j < M+2; j++) {
            mater << mat[i][j] << " ";
        }
        mater << endl;
    }
    mater.close();
}

// Set the discretization constants that do not change over time
void discretization_constants() {
    // Initialize values
    aP.resize(N+2, vector<double>(M+2, 0.0));
    aE.resize(N+2, vector<double>(M+2, 0.0));
    aW.resize(N+2, vector<double>(M+2, 0.0));
    aN.resize(N+2, vector<double>(M+2, 0.0));
    aS.resize(N+2, vector<double>(M+2, 0.0));
    bP.resize(N+2, vector<double>(M+2, 0.0));

    // i = 1; j = 2, ..., M+1 (TOP)
    for(int j = 1; j < M+1; j++) {
        aS[0][j] = harmonic_lambda(0, j, "S")/(abs(yP[0]-yP[1]))*Ss;
        aP[0][j] = aS[0][j];
        bP[0][j] = qwTop*Sn;
    }

    // i = 2, ..., N+1; j = 2, ..., M+1 (MIDDLE)
    for(int i = 1; i < N+1; i++) {
        for(int j = 1; j < M+1; j++) {
            aW[i][j] = harmonic_lambda(i, j, "W")/(abs(xP[j]-xP[j-1]))*Sw;
            aE[i][j] = harmonic_lambda(i, j, "E")/(abs(xP[j+1]-xP[j]))*Se;
            aN[i][j] = harmonic_lambda(i, j, "N")/(abs(yP[i]-yP[i-1]))*Sn;
            aS[i][j] = harmonic_lambda(i, j, "S")/(abs(yP[i+1]-yP[i]))*Ss;
            aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] + rho[mat[i][j]]*cp[mat[i][j]]*VP/dt;
        }
    }

    // i = 2, ..., N+1; j = 1 (LEFT)
    for(int i = 1; i < N+1; i++) {
        aE[i][0] = harmonic_lambda(i, 0, "E")/(abs(xP[1]-xP[0]))*Se;
        aP[i][0] = aE[i][0] + alpha_ext*Sw;
        bP[i][0] = alpha_ext*Sw*Text;
    }
}

// Perform initial calculations
void previous_calculations() {
    // Geometry properties
    V = L*H;
    SwallH = H;
    SwallL = L;
    xP.resize(M+2, 0.0);
    yP.resize(N+2, 0.0);
    xP[0] = 0.0;
    yP[0] = H;
    xP[1] = dx/2;
    yP[1] = H - dy/2;
    xP[M+1] = L;
    yP[N+1] = 0.0;

    // Control volume x positions
    for(i = 2; i < M+1; i++) {
        xP[i] = xP[i-1] + dx;
    }

    // Control volume y positions
    for(i = 2; i < N+1; i++) {
        yP[i] = yP[i-1] - dy;
    }

    // Control volume properties (structured mesh)
    VP = dx*dy;
    Sn = dx;
    Ss = dx;
    Sw = dy;
    Se = dy;

    // Find material properties for each node
    find_materials();

    // Set the discretization constants that do not change over time
    discretization_constants();
}

// Solve the problem with Gauss-Seidel
void solve_with_GS() {
    // Create the T_unsolved matrix
    vector<vector<double>> T_unsolved;
    T_unsolved.resize(N+2, vector<double>(M+2, T0));

    // Create the error variables
    double error;
    vector<double> error_vec;

    // Current map of the temperature
    T[t] = T[t-1];

    // i = N+2; j = 1, ..., M+2 (BOTTOM)
    for(int j = 0; j < M+2; j++) {
        T[t][N+1][j] = TwBottom;
    }

    // i = 1, ..., N+1; j = M+2 (RIGHT)
    for(int i = 0; i < N+1; i++) {
        T[t][i][M+1] = TwRight_0 + TwRight_slope*mytime;
    }

    // Perform the iterations
    for(int iter = 0; iter < max_iter; iter++) {
        // Store the previous temperature
        T_unsolved = T[t];

        // i = 1; j = 2, ..., M+1 (TOP)
        for(int j = 1; j < M+1; j++) {
            T[t][0][j] = (aS[0][j]*T[t][1][j] + bP[0][j])/aP[0][j];
        }

        // i = 2, ..., N+1; j = 2, ..., M+1 (MIDDLE)
        for(int i = 1; i < N+1; i++) {
            for(int j = 1; j < M+1; j++) {
                bP[i][j] = rho[mat[i][j]]*cp[mat[i][j]]*VP/dt*T[t-1][i][j];
                T[t][i][j] = (aE[i][j]*T[t][i][j+1] + aW[i][j]*T[t][i][j-1] + aN[i][j]*T[t][i-1][j] + aS[i][j]*T[t][i+1][j] + bP[i][j])/aP[i][j];
            }
        }

        // i = 2, ..., N+1; j = 1 (LEFT)
        for(int i = 1; i < N+1; i++) {
            T[t][i][0] = (aE[i][0]*T[t][i][1] + bP[i][0])/aP[i][0];
        }

        // Compute the error
        error_vec = abs_mat_to_vec(T[t] - T_unsolved);
        error = *max_element(error_vec.begin(), error_vec.end());

        // Update the temperature
        T[t] = T_unsolved + (T[t] - T_unsolved)*rf;

        // Check the convergence
        if(error < delta) {
            cout << "Converged in " << iter << " iterations! ";
            break;
        }
    }

    // i = 1; j = 1 (TOP LEFT)
    T[t][0][0] = (T[t][0][1] + T[t][1][0])/2.0;

    // Check if the solution converged
    if(error >= delta) {
        cout << "Solution did not converge!" << endl;
    }
}

// Check the energy balance
void energy_balance() {
    // Calculate the energy balance
    double Q = 0.0;
    double U = 0.0;

    // i = 1; j = 2, ..., M+1 (TOP)
    Q += qwTop*SwallL;

    // i = 2, ..., N+1; j = 1 (LEFT)
    for(int i = 1; i < N+1; i++) {
         Q += alpha_ext*Sw*(Text - T[t][i][0]);
    }

    // i = 2, ..., N+1; j = M+2 (RIGHT)
    for(int i = 1; i < N+1; i++) {
        Q += harmonic_lambda(i, M+1, "W")*(T[t][i][M+1] - T[t][i][M])/(abs(xP[M+1]-xP[M]))*Sw;
    }

    // i = N+2; j = 2, ..., M+2 (BOTTOM)
    for(int j = 1; j < M+1; j++) {
        Q += harmonic_lambda(N+1, j, "N")*(T[t][N+1][j] - T[t][N][j])/(abs(yP[N+1]-yP[N]))*Sn;
    }

    // Internal energy of every control volume
    for(int i = 1; i < N+1; i++) {
        for(int j = 1; j < M+1; j++) {
            U += rho[mat[i][j]]*cp[mat[i][j]]*VP*(T[t][i][j] - T[t-1][i][j])/dt;
        }
    }

    cout << "(U - Q)/U: " << abs((U - Q)/U) << endl;

    // Check the energy balance
    if(abs((U-Q)/U) > 1e-5) {
        wrong_time_iter += 1;
    }
}

// Solve the transitory state
void solve_transitory() {
    // Set the initial temperature
    T.resize(nt+1, vector<vector<double>>(N+2, vector<double>(M+2, T0)));

    // Energy balance
    wrong_time_iter = 0;

    // i = 1, ..., N+1; j = M+2 (RIGHT)
    for(int i = 0; i < N+1; i++) {
        T[0][i][M+1] = TwRight_0;
    }

    // i = N+2; j = 1, ..., M+2 (BOTTOM)
    for(int j = 0; j < M+2; j++) {
        T[0][N+1][j] = TwBottom;
    }

    // Perform the iterations
    for(t = 1; t <= nt; t++) {
        // Real time variable
        mytime = ti + t*dt;

        // Current time iteration
        cout << "Time step " << t << " of " << nt << ". ";

        // Solve the problem with Gauss-Seidel
        solve_with_GS();

        // Check the energy balance
        energy_balance();
    }

    // Print the energy balance
    cout << "There have been " << wrong_time_iter << " iterations where the energy balance was not satisfied!" << endl;
}

// Display results
void display_results(string method) {
    // Create output file
    ofstream file;
    file.open("output/" + method + " N=" + to_string(N) + " M=" + to_string(M) + " nt=" + to_string(int(nt)) + " results.txt");

    // Print the results
    for(int i = 0; i < N+2; i++) {
        for(int j = 0; j < M+2; j++) {
            file << T[nt][i][j] - 273 << " ";
        }
        file << endl;
    }
    file.close();
}

// Main ------------------------------------------------------------------
int main() { // Routine to run the program
    // Read the json file and assign values to the variables
    set_constants();

    // Perform initial calculations such as control volumes, discretization constants...
    previous_calculations();

    // Start time
    auto start = chrono::high_resolution_clock::now();

    // Solve the transitory state
    solve_transitory();

    // End time and duration
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    // Print the duration
    int minutes = int(duration.count()/60);
    double seconds = duration.count() - minutes*60;
    cout << "Elapsed time to solve transitory state: " << minutes << " minutes and " << seconds << " seconds." << endl;

    // Display the results
    display_results("GS");

    return 0;
}