#include <iostream>
#include <cmath>

// Define the functions f(x, y) and g(x, y) according to your specific problem
double f(double x, double y) {
    return 0.0;
}

double g(double x, double y) {
    if (x == 0.0) {
        return 0.0;  // u(0, y) = 0
    } else if (y == 0.0) {
        return 0.0;  // u(x, 0) = 0
    } else if (y == 0.5) {
        return 200.0 * x;  // u(x, 0.5) = 200x
    } else if (x == 0.5) {
        return 200.0 * y;  // u(0.5, y) = 200y
    } else {
        return 0.0;  // Default value, can be modified based on your problem
    }
}

int main() {
    // Define the input values
    double a,b,c,d;
    a = 0.0;
    b = 0.5;
    c = 0.0; 
    d = 0.5;  
    // endpoints
    int m,n;
    m = 4;
    n = 4;          
    // grid dimensions
    double TOL = 0.00001;        
    // tolerance
    int N = 1000;             
    // maximum number of iterations

    // Step 1: Set h and k
    double h = (b - a) / n;
    double k = (d - c) / m;

    // Step 2: Generate mesh points
    double* x = new double[n];
    for (int i = 0; i < n; i++) {
        x[i] = a + i * h;
    }

    // Step 3: Generate mesh points
    double* y = new double[m];
    for (int j = 0; j < m; j++) {
        y[j] = c + j * k;
    }

    // Step 4: Initialize w matrix
    double** w = new double*[n];
    for (int i = 0; i < n; i++) {
        w[i] = new double[m];
        for (int j = 0; j < m; j++) {
            w[i][j] = 0.0;
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
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    std::cout << "x = " << x[i] << ", y = " << y[j] << ", w = " << w[i][j] << std::endl;
                }
            }
            // Print the final value of l
            std::cout << "Final value of l: " << l << std::endl;
            // Print the Tol
            std::cout << "Value of Tol: " << TOL << std::endl;

            // Step 19: Stop
            break;
        }

        // Step 20: Increment l
        l++;
    }

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







