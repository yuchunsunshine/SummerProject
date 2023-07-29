#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

// Define the functions f(x, y) and g(x, y), and the exact solution u(x,y)
double f(double x, double y) {
    return x * exp(y);
}

double g(double x, double y) {
    if (x == 0.0) {
        return 0.0;  // u(0, y) = 0
    } else if (y == 0.0) {
        return x;  // u(x, 0) = x
    } else if (y == 1.0) {
        return M_E * x;  // u(x, 1) = ex
    } else if (x == 2.0) {
        return 2.0 * exp(y);  // u(2, y) = 2e^{y}
    } else {
        return 0.0;  
    }
}

double exactSolution(double x, double y) {
    return x * exp(y);
}


int main() {
    // Define the input values
    double a,b,c,d;
    a = 0.0;
    b = 2.0;
    c = 0.0; 
    d = 1.0;  
    // endpoints
    int m,n; //5,6
    m = 10;
    n = 12;          
    // grid dimensions
    double TOL = 1e-10;        
    // tolerance
    int N = 10000;             
    // maximum number of iterations

    // Step 1: Set h and k
    double h = (b - a) / n; // 2/6=0.3333
    double k = (d - c) / m; // 1/5=0.2

    // Step 2: Generate mesh points
    double* x = new double[n + 1];
    for (int i = 0; i < n + 1; i++) {
        x[i] = a + i * h;
    }

    // Step 3: Generate mesh points
    double* y = new double[m + 1];
    for (int j = 0; j < m + 1; j++) {
        y[j] = c + j * k;
    }

    // Step 4: Initialize w matrix
    double** w = new double*[n + 1];
    for (int i = 0; i < n + 1; i++) {
        w[i] = new double[m];
        for (int j = 0; j < m + 1; j++) {
            w[i][j] = g(x[i], y[j]);
        }
    }

    // Step 5: Set lambda, mu, and l
    double lambda = pow(h, 2) / pow(k, 2);
    double mu = 2 * (1 + lambda);
    int l = 1;

    // Step 6: Gauss-Seidel iterations
    while (l <= N) {
        // Step 7: Update w[1, m-1]
        double z = (-pow(h, 2) * f(x[1], y[m - 1]) + g(a, y[m - 1]) + lambda * g(x[1], d) + lambda * w[1][m - 2] + w[2][m - 1]) / mu;
        double NORM = std::abs(z - w[1][m - 1]);
        w[1][m - 1] = z;

        // Step 8: Update w[2, m-1] to w[n-2, m-1]
        for (int i = 2; i < n - 1; i++) {
            z = (-pow(h, 2) * f(x[i], y[m - 1]) + lambda * g(x[i], d) + w[i - 1][m - 1] + w[i + 1][m - 1] + lambda * w[i][m - 2]) / mu;
            if (std::abs(z - w[i][m - 1]) > NORM) {
                NORM = std::abs(z - w[i][m - 1]);
            }
            w[i][m - 1] = z;
        }

        // Step 9: Update w[n-1, m-1]
        z = (-pow(h, 2) * f(x[n - 1], y[m - 1]) + g(b, y[m - 1]) + lambda * g(x[n - 1], d) + lambda * w[n - 1][m - 2] + w[n - 2][m - 1]) / mu;
        if (std::abs(z - w[n - 1][m - 1]) > NORM) {
            NORM = std::abs(z - w[n - 1][m - 1]);
        }
        w[n - 1][m - 1] = z;

        // Step 10: Update w[1, j] to w[n-1, j] for j = m-2 to 2
        for (int j = m - 2; j >= 2; j--) {
            // Step 11: Update w[1, j]
            z = (-pow(h, 2) * f(x[1], y[j]) + g(a, y[j]) + lambda * w[1][j + 1] + lambda * w[1][j - 1] + w[2][j]) / mu;
            if (std::abs(z - w[1][j]) > NORM) {
                NORM = std::abs(z - w[1][j]);
            }
            w[1][j] = z;

            // Step 12: Update w[2, j] to w[n-2, j]
            for (int i = 2; i < n - 1; i++) {
                z = (-pow(h, 2) * f(x[i], y[j]) + lambda * w[i][j + 1] + w[i - 1][j] + w[i + 1][j] + lambda * w[i][j - 1]) / mu;
                if (std::abs(z - w[i][j]) > NORM) {
                    NORM = std::abs(z - w[i][j]);
                }
                w[i][j] = z;
            }

            // Step 13: Update w[n-1, j]
            z = (-pow(h, 2) * f(x[n - 1], y[j]) + g(b, y[j]) + lambda * w[n - 1][j + 1] + lambda * w[n - 1][j - 1] + w[n - 2][j]) / mu;
            if (std::abs(z - w[n - 1][j]) > NORM) {
                NORM = std::abs(z - w[n - 1][j]);
            }
            w[n - 1][j] = z;
        }

        // Step 14: Update w[1, 1]
        z = (-pow(h, 2) * f(x[1], y[1]) + g(a, y[1]) + lambda * g(x[1], c) + lambda * w[1][2] + w[2][1]) / mu;
        if (std::abs(z - w[1][1]) > NORM) {
            NORM = std::abs(z - w[1][1]);
        }
        w[1][1] = z;

        // Step 15: Update w[2, 1] to w[n-2, 1]
        for (int i = 2; i < n - 1; i++) {
            z = (-pow(h, 2) * f(x[i], y[1]) + lambda * g(x[i], c) + w[i - 1][1] + lambda * w[i][2] + w[i + 1][1]) / mu;
            if (std::abs(z - w[i][1]) > NORM) {
                NORM = std::abs(z - w[i][1]);
            }
            w[i][1] = z;
        }

        // Step 16: Update w[n-1, 1]
        z = (-pow(h, 2) * f(x[n - 1], y[1]) + g(b, y[1]) + lambda * g(x[n - 1], c) + lambda * w[n - 1][2] + w[n - 2][1]) / mu;
        if (std::abs(z - w[n - 1][1]) > NORM) {
            NORM = std::abs(z - w[n - 1][1]);
        }
        w[n - 1][1] = z;

        // Step 17: Check convergence
        if (NORM <= TOL) {
            // Step 18: Output results
            for (int i = 0; i < n + 1; i++) {
                for (int j = 1; j < m + 1; j++) {
                    double difference = std::abs(w[i][j] - exactSolution(x[i], y[j]));
                    // std::cout << i << '&' << j << '&' << x[i] << '&' << y[j] << '&' << w[i][j] << '&' << exactSolution(x[i], y[j]) << '&' << difference << '\\' << std::endl;
                    // std::cout << "x = " << x[i] << ", y = " << y[j] << ", w[i][j] = " << w[i][j] << ", Absolute Difference = " << difference << std::endl;
                     std::cout <<  " Absolute Difference = " << difference << std::endl;

                }
            }
            // Print the final value of l
            std::cout << "Final value of l: " << l << std::endl;


           // compare with the exact solution u(x,y) = x*e^{y}
           // std::cout << "Value of Tol: " << TOL << std::endl;

            // Step 19: Stop
            break;
        }

        // Step 20: Increment l
        l++;
    }

    // Step: Save the data to a file
   /* std::ofstream dataFile("data.txt");
    if (dataFile.is_open()) {
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m + 1; j++) {
                double xval = x[i];
                double yval = y[j];
                double wval = w[i][j];
                dataFile << xval << " " << yval << " " << wval << "\n";
            }
            dataFile << "\n";
        }
        dataFile.close();
    } else {
        std::cout << "Unable to open data file." << std::endl;
    }  */
    

    // Step 21: Output failure message
    if (l > N) {
        std::cout << "Maximum number of iterations exceeded." << std::endl;
    }

    // Cleanup: Deallocate memory
    delete[] x;
    delete[] y;
    for (int i = 0; i < n; i++) {
        delete[] w[i];
    }
    delete[] w;

    return 0;
}





