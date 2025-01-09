#ifndef COMPLEX_H
#define COMPLEX_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <ctime>
#define PI 3.1415926
using namespace std;

typedef struct Complex {
    double re;
    double im;
    Complex() {
        re = 0;
        im = 0;
    };
    Complex(double a,double b) {
        re = a;
        im = b;
    };
} complex;


complex* append_seq(complex seq_1[], complex seq_2[], int N);
complex* reorder_seq(complex input_seq[], int N);
complex* Calc_WN(int N);
int reverse_bit(int value, int N);

// Multiplier
inline complex ComplexMul(complex c1, complex c2)
{
    complex r;
    r.re = c1.re*c2.re - c1.im*c2.im;
    r.im = c1.re*c2.im + c1.im*c2.re;
    return r;
}

// Adder
inline complex ComplexAdd(complex c1, complex c2)
{
    complex r;
    r.re = c1.re + c2.re;
    r.im = c1.im + c2.im;
    return r;
}

inline complex ComplexSub(complex a, complex b) {
    complex result;
    result.re = a.re - b.re;  
    result.im = a.im - b.im;  
    return result;
}

// -c
inline complex ReverseComplex(complex c)
{
    c.re = -c.re;
    c.im = -c.im;
    return c;
}

inline complex ComplexConjugate(complex a){
       complex result;
       result.re= a.re;
       result.im= -a.im;
       return result;
}

// Append [seq_1] & [seq_2]
inline complex* append_seq(complex seq_1[], complex seq_2[], int N) {
    complex* total_seq = new complex[N*2];
    #pragma omp parallel for 
    for (int i = 0; i < N; i++) {
        total_seq[i] = seq_1[i];
        total_seq[i+N] = seq_2[i];
    }
    return total_seq;
}

// Reorder the input_seq
inline complex* reorder_seq(complex input_seq[], int N) {
    complex* reordered_seq = new complex[N];
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int k = reverse_bit(i, static_cast<int>(log2(N)));
        reordered_seq[k] = input_seq[i];
    }
    return reordered_seq;
}

// Reverse Bit
inline int reverse_bit(int value, int N) {
    int ret = 0;
    int i = 0;
    while (i < N) {
        ret <<= 1;
        ret |= (value>>i) & 1;
        i++;
    }
    return ret;
}

// Calc WN[]
inline complex* Calc_WN(int N) {
    complex* WN = new complex[N];
    complex WN_unit; 
    WN_unit.re = cos(2*PI/N); 
    WN_unit.im = -sin(2*PI/N);
    WN[0].re=1; WN[0].im=0;
    for (int i = 1; i < N; ++i) {
        WN[i] = ComplexMul(WN[i-1], WN_unit);
    }
    return WN;
}

inline complex* Calc_Inverse_WN(int N) {
    complex* WN = new complex[N];
    complex WN_unit; 
    WN_unit.re = cos(2*PI/N); 
    WN_unit.im = sin(2*PI/N);
    WN[0].re=1; WN[0].im=0;
    for (int i = 1; i < N; ++i) {
        WN[i] = ComplexMul(WN[i-1], WN_unit);
    }
    return WN;
}

#endif

