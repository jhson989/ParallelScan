#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#define DTYPE float
const size_t N = 1<<20;
timeval start, end;

// Sequential prefix-sum algorithm implementation
void get_groundtruth(std::vector<DTYPE>& gt, std::vector<DTYPE> in);

// Parallel prefix-sum algorithm implementations
void parallel_prefix_sum_1(std::vector<DTYPE> in, std::vector<DTYPE>& out);


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
    std::cout << "Parallel Prefix Sum algorithm Start\n";
    std::cout << "--- The length of sequence is " << N <<"\n";
    std::cout << "--- The total memory size of sequence is " << N*sizeof(DTYPE)/1024.0/1024/1024 <<" GB \n";
    srand(time(NULL));
    std::vector<DTYPE> in(N); 
    std::vector<DTYPE> out(N); 
    std::generate(in.begin(), in.end(), [](){return std::rand()%200-100;});
    std::vector<DTYPE> gt(N, 0);


    /*******************************************
     * Getting the ground truth with a naive sequential algorithm
     ********************************************/
    std::cout << "Sequential Prefix Sum algorithm Start\n";
    gettimeofday(&start, NULL);
    get_groundtruth(gt, in);
    gettimeofday(&end, NULL);
    std::cout << "--- Total elapsed time : " << get_sec() << " s\n\n";


    /*******************************************
     * Getting the result of prefix-sum with parallel algorithm
     ********************************************/
    std::cout << "Parallel Prefix Sum algorithm Start\n";
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

void parallel_prefix_sum_1(std::vector<DTYPE> in, std::vector<DTYPE>& out) {

    std::copy(in.begin(), in.end(), out.begin());
    std::vector<DTYPE>::iterator p[2] = {in.begin(), out.begin()};
    for (int d=0; d<log2(N); d++) {
        #pragma omp parallel for num_threads(8) schedule(dynamic)
        for (int i=0; i<in.size(); i++) {
            if (i >= pow(2, d)) {
                out[i] = in[i-pow(2,d)] + out[i];
            }
        }
        std::copy(out.begin(), out.end(), in.begin());
    }
}



void get_groundtruth(std::vector<DTYPE>& gt, std::vector<DTYPE> in) {
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

void check_result(std::vector<DTYPE>& gt, std::vector<DTYPE>& in) {

    std::cout << "Start Checking result...\n";
    bool result = true;
    for (size_t i=0; i<gt.size(); i++) {
        if (gt[i] != in[i]) {
            result = false;
            break;
        }
    }

    if (result) {
        std::cout << "--- The result is correct!!!\n";
    } else {
        std::cout << "--- The result is incorect!!!\n";
    }

}