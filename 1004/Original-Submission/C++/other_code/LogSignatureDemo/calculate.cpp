#include<array>
#include<iostream>
#include<random>
#include<vector>

#include "bch.h"

using namespace std;

const size_t D = DIM; //dimension of path
const size_t M = LEVEL; //level of log signature to calculate

//generates a random D-dimensional path with a new seed every time
vector<array<float,D>> generateRandomPath(){
  mt19937 randomNumbers{random_device()()};
  vector<array<float,D>> path;
  uniform_real_distribution<float> urd(-1,1);
  for(int i=0;i!=100; ++i){
    array<float, D> point;
    for(int d=0; d<D; ++d){
      point[d]=urd(randomNumbers);
    }
    path.push_back(point);
  }
  return path;
}

//returns the logsignature of a path
vector<double> calculateLogSignature(
		      const vector<array<float,D>>& path){
  LogSignature<D,M> logsig{};
  if (path.size()>1)
    for(int j=0; j<D; ++j)
      logsig[j] = path[1][j]-path[0][j];
  for(size_t i=2; i<path.size(); ++i){
    Segment<D> displacement;
    for(int j=0; j<D; ++j)
      displacement[j] = path[i][j]-path[i-1][j];
    joinSegmentToSignatureInPlace<D,M>(logsig,displacement);
  }
  return {logsig.begin(), logsig.end()}; 
}

int main(){
  auto path = generateRandomPath();
  vector<double> logsignature = calculateLogSignature(path);
  for(double p : logsignature)
    cout<<p<<" ";
  cout<<"\n";
}
