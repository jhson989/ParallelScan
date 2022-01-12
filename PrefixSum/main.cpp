#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <sys/time.h>
#define DTYPE float
const size_t N = 1<<18;
timeval start, end;

// Sequential prefix-sum algorithm implementation
void get_groundtruth(std::vector<DTYPE>& gt, std::vector<DTYPE>& in);

// Helper functions
inline void print_vector(std::vector<DTYPE>& in);
inline float get_sec();
void check_result(std::vector<DTYPE>& gt, std::vector<DTYPE>& in);

int main(void) {

    /*******************************************
     * Data initialization
     * in : input sequence. Prefix-sum is a in-place operation
     * gt : ground truth, calculated by the get_groundtruth() function
     ********************************************/
    std::cout << "Parallel Prefix Sum algorithm Start\n";
    std::cout << "--- The length of sequence is " << N <<"\n";
    std::cout << "--- The total memory size of sequence is " << N*sizeof(DTYPE)/1024.0/1024/1024 <<" GB \n";
    srand(time(NULL));
    std::vector<DTYPE> in(N); 
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
     * Checking the correctness
     ********************************************/
    check_result(gt, in);
    return 0;

}

void get_groundtruth(std::vector<DTYPE>& gt, std::vector<DTYPE>& in) {

    for (size_t i=0; i<gt.size(); i++) {
        for (size_t j=0; j<=i; j++) {
            gt[i] += in[j];
        }
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