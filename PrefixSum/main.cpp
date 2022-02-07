#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <omp.h>

#define DTYPE int


const int NUM_THREADS = 1;
const size_t N = 1<<28;
timeval start, end;

// Sequential prefix-sum algorithm implementation
void get_groundtruth(std::vector<DTYPE>& gt, const std::vector<DTYPE>& in);

// Parallel prefix-sum algorithm implementations
void parallel_prefix_sum_1(const std::vector<DTYPE>& in, std::vector<DTYPE>& out);


// Helper functions
inline void print_vector(std::vector<DTYPE>& in);
inline float get_sec();
void check_result(std::vector<DTYPE>& gt, std::vector<DTYPE>& out);

int main(void) {

    /*******************************************
     * Data initialization
     * in : the input sequence.
     * out : the result of parallel algorithm
     * gt : the ground truth, calculated by the get_groundtruth() function
     ********************************************/
    std::cout << "************************************************************************\n";
    std::cout << "Parallel Prefix Sum Algorithm Start\n";
    std::cout << "--- The length of sequence is " << N <<"\n";
    std::cout << "--- The total memory size of sequence is " << N*sizeof(DTYPE)/1024.0/1024/1024 <<" GB \n";
    std::cout << "************************************************************************\n";
    srand(time(NULL));
    std::vector<DTYPE> in(N); 
    std::vector<DTYPE> out(N, 0); 
    std::generate(in.begin(), in.end(), [](){return std::rand()%200-100;});
    std::vector<DTYPE> gt(N, 0);


    /*******************************************
     * Getting the ground truth with a naive sequential algorithm
     ********************************************/
    std::cout << "\nSequential Prefix Sum Algorithm Start\n";
    gettimeofday(&start, NULL);
    get_groundtruth(gt, in);
    gettimeofday(&end, NULL);
    std::cout << "--- Total elapsed time : " << get_sec() << " s\n\n";


    /*******************************************
     * Getting the result of prefix-sum with parallel algorithm
     ********************************************/
    std::cout << "\nParallel Prefix Sum Algorithm Start\n";
    gettimeofday(&start, NULL);
    parallel_prefix_sum_1(in, out);
    gettimeofday(&end, NULL);
    std::cout << "--- Total elapsed time : " << get_sec() << " s\n\n";

    

    /*******************************************
     * Checking the correctness
     ********************************************/
    check_result(gt, out);
    return 0;

}

DTYPE prefix_sum_1_first_path(const std::vector<DTYPE>& in, std::vector<DTYPE>& out, size_t start, size_t end) {

    if (end-start <= (N>>4)) {
        out[start] = in[start];
        for (size_t i=start+1; i<end; i++)
            out[i] = out[i-1] + in[i];
        return out[end-1];
    }

    size_t mid = (start+end)/2;
    DTYPE left=0, right=0;

    #pragma omp task shared(out, left)
    left = prefix_sum_1_first_path(in, out, start, mid);

    right = prefix_sum_1_first_path(in, out, mid, end);

    #pragma omp taskwait 
    out[end-1] = left + right;
    return out[end-1];
}

void prefix_sum_1_second_path(const std::vector<DTYPE>& in, std::vector<DTYPE>& out, size_t start, size_t end, DTYPE left) {


    if (end-start <= (N>>4)) {
        for (size_t i=start; i<end-1; i++)
            out[i] += left;
        out[end-1] = out[end-2] + in[end-1];
        return;
    }

    size_t mid = (start+end)/2;
    DTYPE temp = out[mid-1];

    #pragma omp task shared(out)
    prefix_sum_1_second_path(in, out, start, mid, left);
    prefix_sum_1_second_path(in, out, mid, end, left+temp);

    #pragma omp taskwait

}

void parallel_prefix_sum_1(const std::vector<DTYPE>& in, std::vector<DTYPE>& out) {


    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp single
        {
            prefix_sum_1_first_path(in, out, 0, in.size());
            prefix_sum_1_second_path(in, out, 0, in.size(), 0);
        }
    }
}




void get_groundtruth(std::vector<DTYPE>& gt, const std::vector<DTYPE>& in) {
    /*
     * The simple sequenTal algorithm: accumulate the sum from left to right
     * Work : O(N), Span : O(N)
     */
    gt[0] = in[0];
    for (size_t i=1; i<gt.size(); i++) {
        gt[i] = gt[i-1] + in[i];
    }
    
}

void print_vector(std::vector<DTYPE>& in) {
    for(auto i=in.begin(); i!=in.end(); i++)
        std::cout << *i << " ";
    std::cout << std::endl;
}

float get_sec() {
    return (end.tv_sec-start.tv_sec) + (end.tv_usec-start.tv_usec)*1e-6;
}

void check_result(std::vector<DTYPE>& gt, std::vector<DTYPE>& out) {

    std::cout << "Start Checking result...\n";
    bool result = true;
    for (size_t i=0; i<gt.size(); i++) {
        if (gt[i] != out[i]) {
            result = false;
            std::cout << "--- [[ERR]] Incorrect at out["<<i<<"] with result="<<out[i]<<" but, gt="<<gt[i]<<"\n";
            break;
        }
    }

    if (result) {
        std::cout << "--- The result is correct!!!\n";
    } else {
        std::cout << "--- The result is incorect!!!\n";
    }

}








/*******************************************
 * Legacy
 ********************************************/


/*
void parallel_prefix_sum_1(std::vector<DTYPE> in, std::vector<DTYPE>& out) {

    out[0] = in[0];
    std::vector<DTYPE>::iterator p[2] = {in.begin(), out.begin()};
    bool target = true;
    for (int d=0; d<log2(N); d++) {
        //#pragma omp parallel for num_threads(8) schedule(static)
        for (int i=0; i<N; i++) {
            if (i >= pow(2, d)) {
                p[(int)target][i] = p[(int)(!target)][i-pow(2,d)] + p[(int)(!target)][i];
            } else {
                p[target][i] = p[!target][i];
            }
        }
        target = !target;
    }

    if (target) {
        std::copy(in.begin(), in.end(), out.begin());
    }

}
*/