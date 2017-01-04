#ifndef INCLUDED_BASIC_STATS
#define INCLUDED_BASIC_STATS

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace std;

template<typename T>
class Basic_Stats {
    vector<T> v_;
    vector<T> vs_;
    int N;
    T sum_;
    T mean_;
    T var_; //sample variance
    T std_; //sample standard deviation
    T min_;
    T max_;
    T median_;
    void compute_basic_stats();
public:
    Basic_Stats<T>(vector<T> vec);
    int size() {return N;}
    T sum() {return sum_;}
    T mean() {return mean_;}
    T var() {return var_;}
    T std() {return std_;}
    T median() {return median_;}
    T ptile(int percent); //return percentile value
};

template<typename T>
Basic_Stats<T>::Basic_Stats(vector<T> vec){
    N= vec.size();
    v_= vec;
    vs_= v_;
    sort(vs_.begin(),vs_.end());
    min_= vs_[0];
    max_= vs_[N-1];
    compute_basic_stats();
}

template<typename T>
void Basic_Stats<T>::compute_basic_stats(){
    sum_= accumulate(v_.begin(),v_.end(),0.0);
    mean_= sum_/N;
    T sq_sum= inner_product(v_.begin(),v_.end(),v_.begin(),0.0); // squared_sum
    T sq_mean= sq_sum/N;
    var_= (sq_mean-mean_*mean_)*N*1.0/(N-1.0);
    std_= sqrt(var_);
    median_= N%2 ? vs_[N/2] : (vs_[N/2-1]+vs_[N/2])/2.0;
}

template<typename T>
T Basic_Stats<T>::ptile(int percent){
    int index= N*percent/100;
    if(index > N) index= N-1;
    if(index < 0) index= 0;
    return vs_[index];
}


#endif
