/* EDGE: Ensemble Dimensionality Reduction */
/* calculating the affinity matrix */
/* 08/26/2019 */

//[[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include <cfloat>
#include <boost/random.hpp>
#include <algorithm>
#include <numeric>
#include <utility>
using namespace Rcpp;


// random number generator
typedef boost::mt19937 Rng;
typedef boost::uniform_real<float> uniform;
typedef boost::uniform_int<int> uniformInt;
typedef boost::uniform_int<unsigned> uniformUnsigned;

bool sortdesc(const std::pair<double, int>& a,
              const std::pair<double, int>& b) {
  return (std::get<0>(a) > std::get<0>(b));
}

// [[Rcpp::export]]
NumericMatrix affinity(NumericMatrix& dat, int n_wl, int n_dm, int nk, int H=1017881, int seed=5489){
  int i, j, k, l, m, n, count;
  int data_cell;
  int data_gene;
  int zero_ct;
  float min_, max_;
  double data_min, data_max;
  uint32_t index, it, dt, wt, bi;
  float tt;
  std::vector<int> idx;
  std::vector<int> eles;
  std::vector<std::vector<int> > dims;
  std::vector<std::vector<double> > thresholds;
  std::vector<std::vector<unsigned int> > weights;
  std::vector<std::vector<std::vector<int> > > bins;

  dims.resize(n_wl);
  thresholds.resize(n_wl);
  weights.resize(n_wl);
  bins.resize(n_wl);
  Rng rng(seed);

  min_ = FLT_MAX;
  max_ = -1*FLT_MAX;
  data_cell = dat.rows();
  data_gene = dat.cols();
  std::vector<std::vector<int> > amat(data_cell, std::vector<int>(data_cell, 0));
  std::vector<std::vector<int> > kmat(data_cell, std::vector<int>(data_cell, 0));
  // calculating maximum and minimum of whole data
  for(i=0; i<data_cell; i++){
    for(j=0; j<data_gene; j++){
      if(dat(i,j) > max_) max_ = dat(i,j);
      if(dat(i,j) < min_) min_ = dat(i,j);
    }
  }
  data_max = max_;
  data_min = min_;

  // generating weights and thresholds
  for(i=0; i<n_wl; i++){
    dims[i].resize(n_dm);
    thresholds[i].resize(n_dm);
    weights[i].resize(n_dm);
    for(j=0; j<n_dm; j++){
      dims[i][j] = boost::variate_generator<Rng &, uniformUnsigned>(rng, uniformUnsigned(0, data_gene-1))();
      if(data_min != data_max)
        thresholds[i][j] = boost::variate_generator<Rng &, uniform>(rng, uniform(data_min, data_max))();
      else
        thresholds[i][j] = data_min;
    }
    for(j=0; j<n_dm; j++)
      weights[i][j] = rng();
  }

  // calculating bins
  for(i=0; i<n_wl; i++){
    bins[i].resize(H);
    for(j=0; j<data_cell; j++){
      index = 0;
      for(k=0; k<n_dm; k++){
        dt = dims[i][k];
        tt = thresholds[i][k];
        wt = weights[i][k];
        it = (dat(j, dt) > tt)?1:0;
        index += (wt * it);
      }
      index = index % H;
      bins[i][index].push_back(j);
    }
  }

  // calculating the affinity matrix
  for(i=0; i<n_wl; i++){
    std::vector<std::vector<int> > tmat(data_cell, std::vector<int>(data_cell, 0));
    for(j=0; j<data_cell; j++){
      index = 0;
      for(k=0; k<n_dm; k++){
        dt = dims[i][k];
        tt = thresholds[i][k];
        wt = weights[i][k];
        it = (dat(j, dt) > tt)?1:0;
        index += (wt * it);
      }
      index = index % H;
      for(l=0; l<bins[i][index].size(); l++){
        bi = bins[i][index][l];
        tmat[j][bi] = 1;
      }
    }
    for(m=0; m<data_cell; m++){
      for(n=0; n<data_cell; n++){
        amat[m][n] = amat[m][n] + tmat[m][n];
      }
    }
  }

  // symmetric matrix
  std::vector<std::vector<int> > smat(data_cell, std::vector<int>(data_cell, 0));
  for(m=0; m<data_cell; m++){
    for(n=0; n<data_cell; n++){
      smat[m][n] = amat[m][n] + amat[n][m];
      smat[n][m] = amat[m][n] + amat[n][m];
    }
  }

  // nearest neighbors
  for(i=0; i<data_cell; i++){
    std::vector<std::pair<int, int> > ixp;
    for(int gg=0; gg < data_cell; gg++){
      ixp.push_back(std::make_pair(smat[i][gg], gg));
    }
    std::stable_sort(ixp.begin(), ixp.end(), sortdesc);
    idx.resize(data_cell);
    for(int gg=0; gg<data_cell; gg++){
      idx[gg] = ixp[gg].second;
    }
    eles.resize(data_cell);
    for(int gg=0; gg<data_cell; gg++){
      eles[gg] = ixp[gg].first;
    }
    for(int hh=0; hh<nk; hh++){
      for(j=0; j<data_cell; j++){
        if(eles[j]==eles[hh]){
          kmat[i][idx[j]] = smat[i][idx[j]];
          kmat[idx[j]][i] = smat[idx[j]][i];
        }
      }
    }
  }

  zero_ct = 0;
  for(int gg=0; gg<data_cell; gg++){
    zero_ct = zero_ct + std::count(kmat[gg].begin(),kmat[gg].end(),0);
  }
  NumericMatrix matcoo((data_cell*data_cell-zero_ct), 3);
  count = 0;
  for(i=0; i<data_cell; i++) {
    for(j=0; j<data_cell; j++){
      if(kmat[i][j] != 0){
        matcoo(count,0) = i+1;
        matcoo(count,1) = j+1;
        matcoo(count,2) = kmat[i][j];
        count = count + 1;
      }
    }
  }
  return matcoo;
}




