#include <random>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>  // For std::sort
using namespace std;

// Function to compute adaptive threshold
double compute_adaptive_threshold(complex fft_result[], int N) {
    vector<double> amplitudes(N);
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        amplitudes[i] = ComplexMul(fft_result[i], ComplexConjugate(fft_result[i])).re;
    }

    sort(amplitudes.begin(), amplitudes.end());


    return amplitudes[int(0.75 * N)];
}

complex* denoise_signal(complex input_seq[], int N) {
    
    complex* fft_result = DIT_FFT_reordered(input_seq, N);

  
    double adaptive_threshold = compute_adaptive_threshold(fft_result, N);
    cout << "Adaptive Threshold: " << adaptive_threshold << endl;

    // Apply low-pass filter by zeroing out components with small amplitudes
    for (int i = 0; i < N; i++) {
        double PSD = ComplexMul(fft_result[i], ComplexConjugate(fft_result[i])).re;
        if (PSD < adaptive_threshold) {
            fft_result[i].re = 0;
            fft_result[i].im = 0;
        }
    }

  
    complex* denoised_signal = DIT_IFFT_reordered(fft_result, N);

    delete[] fft_result;  
    return denoised_signal;
}
