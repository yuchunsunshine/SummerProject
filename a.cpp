#include <iostream>
#include <fstream>
#include <cstdlib>

int main() {
    // Step 1: Generate or load your data in suitable arrays or vectors
    // Assuming you have the following arrays or vectors representing your data
    int n = 10;  // Number of data points along x-axis
    int m = 10;  // Number of data points along y-axis
    double* x = new double[n];
    double* y = new double[m];
    double** w = new double*[n];
    
    // Step 2: Populate x, y, and w arrays or vectors with your data

    // Step 3: Save the data to a text file
    std::ofstream dataFile("data.txt");
    if (!dataFile) {
        std::cerr << "Error opening data file." << std::endl;
        return 1;
    }

    for (int i = 0; i < n + 1; i++) {
        dataFile << x[i] << " ";
    }
    dataFile << std::endl;

    for (int j = 0; j < m + 1; j++) {
        dataFile << y[j] << " ";
        for (int i = 0; i < n + 1; i++) {
            dataFile << w[i][j] << " ";
        }
        dataFile << std::endl;
    }

    dataFile.close();

    // Step 4: Generate the Gnuplot script
    std::ofstream scriptFile("plot.gp");
    if (!scriptFile) {
        std::cerr << "Error opening script file." << std::endl;
        return 1;
    }

    scriptFile << "set pm3d interpolate 20,20" << std::endl;
    scriptFile << "set xlabel 'x'" << std::endl;
    scriptFile << "set ylabel 'y'" << std::endl;
    scriptFile << "set zlabel 'w'" << std::endl;
    scriptFile << "set title '3D Surface Plot'" << std::endl;
    scriptFile << "set view 60,30" << std::endl;
    scriptFile << "splot 'data.txt' using 1:2:3 with pm3d" << std::endl;

    scriptFile.close();

    // Step 5: Execute the Gnuplot script using system command
    std::string command = "gnuplot plot.gp";
    int status = system(command.c_str());

    if (status == -1) {
        std::cerr << "Error executing Gnuplot script." << std::endl;
        return 1;
    }

    // Clean up dynamically allocated memory
    delete[] x;
    delete[] y;
    for (int i = 0; i < n + 1; i++) {
        delete[] w[i];
    }
    delete[] w;

    return 0;
}
