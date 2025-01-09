#include "complex.h"         // Assumes definitions of complex, ComplexAdd, ComplexSub, ComplexMul, etc.
#include "dit_fft.h"         // Assumes Calc_WN(), Calc_Inverse_WN(), reverse_bit() are defined here
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <sys/time.h>
#include <mpi.h>
#include <vector>
#include <cassert>

using namespace std;

void tick (struct timeval * t) {
    gettimeofday(t, NULL);
}

double tock (struct timeval * t) {
    struct timeval now;
    gettimeofday(&now, NULL);
    return (double)(now.tv_sec - t->tv_sec) + ((double)(now.tv_usec - t->tv_usec)/1000000.);
}

void save_to_file(const char* filename, const complex* data, int N, const string& label) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    file << "# " << label << " Results\n";
    file << "# Format: Index\tReal Part\tImaginary Part\n";
    for (int i = 0; i < N; ++i) {
        file << i << "\t" << data[i].re << "\t" << data[i].im << "\n";
    }
    file.close();
}

// Helper struct used during reorder communication
struct indexed_complex {
    int gidx;
    complex val;
};

// Perform a global bit-reverse reordering using MPI communication
void global_bit_reverse_reorder(complex local_seq[], int local_N, int N, int rank, int size, MPI_Comm comm) {
    int logN = (int)log2(N);
    int offset = rank * local_N;

    // Compute where each local element should go
    vector<int> send_counts(size, 0);
    vector<int> send_displs(size, 0);
    vector<int> recv_counts(size, 0);
    vector<int> recv_displs(size, 0);

    // Determine final positions and how many elements each rank will receive
    vector<int> target_ranks(local_N);
    vector<int> target_indices(local_N);
    for (int i = 0; i < local_N; ++i) {
        int global_i = offset + i;
        int k = reverse_bit(global_i, logN);
        int target_rank = k / local_N;
        target_ranks[i] = target_rank;
        target_indices[i] = k;
        send_counts[target_rank]++;
    }

    // Convert to prefix sums for send_displs
    for (int r = 1; r < size; r++)
        send_displs[r] = send_displs[r-1] + send_counts[r-1];

    // Prepare send buffer
    vector<indexed_complex> send_buf(local_N);
    vector<int> current_pos(size, 0);
    for (int i = 0; i < size; i++)
        current_pos[i] = send_displs[i];

    for (int i = 0; i < local_N; ++i) {
        int tr = target_ranks[i];
        send_buf[current_pos[tr]].gidx = target_indices[i];
        send_buf[current_pos[tr]].val = local_seq[i];
        current_pos[tr]++;
    }

    // Exchange counts
    vector<int> send_counts_idx = send_counts; // counts of elements
    vector<int> recv_counts_idx(size);

    MPI_Alltoall(send_counts_idx.data(), 1, MPI_INT,
                 recv_counts_idx.data(), 1, MPI_INT,
                 comm);

    // Compute recv displacements
    for (int r = 1; r < size; r++)
        recv_displs[r] = recv_displs[r-1] + recv_counts_idx[r-1];

    int total_recv_idx = 0;
    for (auto c: recv_counts_idx) total_recv_idx += c;

    vector<indexed_complex> recv_buf_idx(total_recv_idx);

    // Convert element counts to byte counts for MPI_Alltoallv
    // Also, displacements must be in bytes
    vector<int> send_counts_bytes(size), recv_counts_bytes(size);
    vector<int> send_displs_bytes(size), recv_displs_bytes(size);

    for (int i = 0; i < size; i++) {
        send_counts_bytes[i] = send_counts_idx[i] * (int)sizeof(indexed_complex);
        recv_counts_bytes[i] = recv_counts_idx[i] * (int)sizeof(indexed_complex);
    }
    // Displacements must also be in bytes
    for (int i = 1; i < size; i++) {
        send_displs_bytes[i] = send_displs_bytes[i-1] + send_counts_bytes[i-1];
        recv_displs_bytes[i] = recv_displs_bytes[i-1] + recv_counts_bytes[i-1];
    }

    MPI_Alltoallv(send_buf.data(), send_counts_bytes.data(), send_displs_bytes.data(), MPI_BYTE,
                  recv_buf_idx.data(), recv_counts_bytes.data(), recv_displs_bytes.data(), MPI_BYTE,
                  comm);

    // Place received elements in correct local order
    // After reorder, rank r holds global indices [offset ... offset+local_N-1]
    // We know each item in recv_buf_idx has a .gidx that tells us where it belongs.
    // local position = gidx - offset
    for (int i = 0; i < local_N; i++) {
        local_seq[i] = complex(0.0, 0.0);
    }

    for (auto &item : recv_buf_idx) {
        int local_pos = item.gidx - offset;
        assert(local_pos >= 0 && local_pos < local_N);
        local_seq[local_pos] = item.val;
    }
}

// Placeholder function: In a fully distributed FFT, you'd implement logic here
// to redistribute data so that each rank has the correct pairs for the butterfly operations.
// For now, we do nothing. A real implementation is needed for full correctness.
void fft_stage_redistribute(complex local_seq[], int local_N, int N, int stage, int rank, int size, MPI_Comm comm) {
    (void)local_seq; (void)local_N; (void)N; (void)stage; (void)rank; (void)size; (void)comm;
    // No-op for demonstration
}

void ifft_stage_redistribute(complex local_seq[], int local_N, int N, int stage, int rank, int size, MPI_Comm comm) {
    (void)local_seq; (void)local_N; (void)N; (void)stage; (void)rank; (void)size; (void)comm;
    // No-op for demonstration
}

complex* DIT_FFT_iterative_distributed(complex local_seq[], int local_N, int N, complex WN[], int rank, int size, MPI_Comm comm) {
    int stages = (int)log2(N);
    complex* current_seq = local_seq;
    complex* temp_seq = new complex[local_N];

    for (int stage = 1; stage <= stages; ++stage) {
        fft_stage_redistribute(current_seq, local_N, N, stage, rank, size, comm);

        int m = 1 << stage;  
        int half_m = m / 2;

        // Simple local butterfly, assumes correct data distribution at this stage
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < local_N; i += m) {
            for (int j = 0; j < half_m; ++j) {
                int idx1 = i + j;
                int idx2 = i + j + half_m;
                if (idx1 < local_N && idx2 < local_N) {
                    complex W = WN[(j * N) / m];
                    complex t = ComplexMul(current_seq[idx2], W);
                    complex u = current_seq[idx1];
                    temp_seq[idx1] = ComplexAdd(u, t);
                    temp_seq[idx2] = ComplexSub(u, t);
                }
            }
        }
        std::swap(current_seq, temp_seq);
    }
    delete[] temp_seq;
    return current_seq;
}

complex* DIT_IFFT_iterative_distributed(complex local_seq[], int local_N, int N, complex WN_inv[], int rank, int size, MPI_Comm comm) {
    int stages = (int)log2(N);
    complex* current_seq = local_seq;
    complex* temp_seq = new complex[local_N];

    for (int stage = 1; stage <= stages; ++stage) {
        ifft_stage_redistribute(current_seq, local_N, N, stage, rank, size, comm);

        int m = 1 << stage;
        int half_m = m / 2;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < local_N; i += m) {
            for (int j = 0; j < half_m; ++j) {
                int idx1 = i + j;
                int idx2 = i + j + half_m;
                if (idx1 < local_N && idx2 < local_N) {
                    complex W = WN_inv[(j * N) / m];
                    complex t = ComplexMul(current_seq[idx2], W);
                    complex u = current_seq[idx1];
                    temp_seq[idx1] = ComplexAdd(u, t);
                    temp_seq[idx2] = ComplexSub(u, t);
                }
            }
        }
        std::swap(current_seq, temp_seq);
    }

    delete[] temp_seq;
    return current_seq;
}

complex* parallel_FFT(complex* global_input, int N, int rank, int size, MPI_Comm comm) {
    int local_N = N / size;
    complex* local_input = new complex[local_N];

    // Distribute data
    MPI_Scatter(global_input, local_N * (int)sizeof(complex), MPI_BYTE,
                local_input, local_N * (int)sizeof(complex), MPI_BYTE,
                0, comm);

    // Compute or broadcast WN
    static complex* WN = nullptr;
    if (rank == 0) {
        WN = Calc_WN(N);
    } else {
        WN = new complex[N];
    }
    MPI_Bcast(WN, N * (int)sizeof(complex), MPI_BYTE, 0, comm);

    // Global reorder using MPI
    global_bit_reverse_reorder(local_input, local_N, N, rank, size, comm);

    // Perform distributed FFT
    complex* local_fft = DIT_FFT_iterative_distributed(local_input, local_N, N, WN, rank, size, comm);

    complex* global_output = nullptr;
    if (rank == 0) {
        global_output = new complex[N];
    }

    MPI_Gather(local_fft, local_N * (int)sizeof(complex), MPI_BYTE,
               global_output, local_N * (int)sizeof(complex), MPI_BYTE,
               0, comm);

    return global_output;
}

complex* parallel_IFFT(complex* global_input, int N, int rank, int size, MPI_Comm comm) {
    int local_N = N / size;
    complex* local_input = new complex[local_N];

    // Scatter global frequency-domain data to each rank
    MPI_Scatter(global_input, local_N * (int)sizeof(complex), MPI_BYTE,
                local_input, local_N * (int)sizeof(complex), MPI_BYTE,
                0, comm);

    // Compute or broadcast WN_inv
    complex* WN_inv = nullptr;
    if (rank == 0) {
        WN_inv = Calc_Inverse_WN(N);
    } else {
        WN_inv = new complex[N];
    }
    MPI_Bcast(WN_inv, N * (int)sizeof(complex), MPI_BYTE, 0, comm);

    // Global reorder
    global_bit_reverse_reorder(local_input, local_N, N, rank, size, comm);

    // Local IFFT computation
    complex* local_ifft = DIT_IFFT_iterative_distributed(local_input, local_N, N, WN_inv, rank, size, comm);

    // Scale by 1/N
    #pragma omp parallel for
    for (int i = 0; i < local_N; ++i) {
        local_ifft[i].re /= N;
        local_ifft[i].im /= N;
    }

    complex* global_output = nullptr;
    if (rank == 0) {
        global_output = new complex[N];
    }

    // Gather results on rank 0
    MPI_Gather(local_ifft, local_N * (int)sizeof(complex), MPI_BYTE,
               global_output, local_N * (int)sizeof(complex), MPI_BYTE,
               0, comm);

    return global_output;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (argc != 3) {
        if (rank == 0) {
            cout << "Usage: ./fft [k] [validate_or_evaluate]\n";
            cout << "Example:\n";
            cout << "  ./fft 10 1  # Validate DIT for N = 2^10\n";
            cout << "  ./fft 10 0  # Evaluate performance for N = 2^10\n";
        }
        MPI_Finalize();
        return 1;
    }

    int k = atoi(argv[1]);
    int mode = atoi(argv[2]);
    int N = 1 << k;

    if (N % size != 0) {
        if (rank == 0) {
            cerr << "N must be divisible by the number of MPI ranks." << endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Generate the input sequence (on rank 0)
    complex* input_seq = nullptr;
    if (rank == 0) {
        input_seq = new complex[N];
        double T = 1.0;
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            double t_i = (i * T) / N;
            double real_val = sin(2 * PI * 50 * t_i) + cos(2 * PI * 120 * t_i);
            input_seq[i] = complex(real_val, 0);
        }
    }

    if (mode == 1) {
        if (rank == 0) {
            cout << "\nValidating for N = " << N << " with distributed FFT & IFFT...\n";
        }
        complex* fft_result = parallel_FFT(input_seq, N, rank, size, comm);
        if (rank == 0 && fft_result != nullptr) {
            save_to_file("fft_output.txt", fft_result, N, "Distributed FFT");
            cout << "FFT results saved to fft_output.txt\n";
        }

        complex* ifft_result = parallel_IFFT(fft_result, N, rank, size, comm);
        if (rank == 0 && ifft_result != nullptr) {
            save_to_file("ifft_output.txt", ifft_result, N, "Distributed IFFT");
            cout << "IFFT results saved to ifft_output.txt\n";
            delete[] ifft_result;
        }
        if (rank == 0) delete[] fft_result;
    } else if (mode == 3) {
        // Performance measure
        struct timeval calc;
        double calctime;
        tick(&calc);
        complex* fft_result = parallel_FFT(input_seq, N, rank, size, comm);
        complex* ifft_result = parallel_IFFT(fft_result, N, rank, size, comm);
        calctime = tock(&calc);
        if (rank == 0) {
            cout << "Run time = " << calctime << " s\n";
        }
        if (rank == 0) {
            delete[] fft_result;
            delete[] ifft_result;
        }
    } else {
        // Just run FFT
        complex* fft_result = parallel_FFT(input_seq, N, rank, size, comm);
        if (rank == 0) delete[] fft_result;
    }

    if (rank == 0) delete[] input_seq;

    MPI_Finalize();
    return 0;
}

