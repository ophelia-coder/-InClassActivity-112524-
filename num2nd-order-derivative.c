#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Function to compute the second-order derivative
void compute_second_derivative(double (*f)(double), double a, double b, int N, double* approx, double* exact) {
    double h = (b - a) / (N - 1);
    for (int i = 0; i < N; i++) {
        double x = a + i * h;
        if (i == 0 || i == N - 1) {
            approx[i] = 0.0; // Boundary points not computed
        } else {
            approx[i] = (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
        }
        exact[i] = -sin(x); // Exact second derivative
    }
}

// Function to compute the L2-norm of the error
double compute_L2_error(double* exact, double* approx, int N, double h) {
    double error_sum = 0.0;
    for (int i = 1; i < N - 1; i++) { // Ignore boundary points
        double diff = exact[i] - approx[i];
        error_sum += diff * diff;
    }
    return sqrt(error_sum * h);
}

// Main function
int main() {
    double a = 0.0, b = M_PI;
    int N = 100; // Number of grid points
    double h = (b - a) / (N - 1);

    double* approx = (double*)malloc(N * sizeof(double));
    double* exact = (double*)malloc(N * sizeof(double));

    // Compute the derivatives
    compute_second_derivative(sin, a, b, N, approx, exact);

    // Compute the L2 error
    double L2_error = compute_L2_error(exact, approx, N, h);

    printf("L2 error: %e\n", L2_error);



    return 0;
}
