#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
struct Match {
  int id; double distance, ri;
};
static bool better(const Match& a, const Match& b, double query_ri) {
  if(a.distance != b.distance) return a.distance < b.distance;
  double da=std::isfinite(a.ri)&&std::isfinite(query_ri)?std::abs(a.ri-query_ri):INFINITY;
  double db=std::isfinite(b.ri)&&std::isfinite(query_ri)?std::abs(b.ri-query_ri):INFINITY;
  return da!=db? da<db : a.id<b.id;
}
static double entropy_cross(double a, double b) {
  return a>0 && b>0 ? a*std::log1p(b/a)+b*std::log1p(a/b) : 0.;
}
// All R objects and work buffers are prepared before entering OpenMP workers.
// [[Rcpp::export]]
DataFrame library_search_cpp(NumericMatrix query, S4 library, NumericVector mass,
    List indices, NumericVector query_ri, double tolerance, int top,
    IntegerVector methods, int threads) {
  IntegerVector ptr=library.slot("p"), row=library.slot("i");
  NumericVector values=library.slot("x");
  const int d=query.nrow(), nq=query.ncol(), nl=ptr.size()-1, nm=methods.size();
  const double *raw=REAL(query), *v=REAL(values), *qr=REAL(query_ri);
  const int *ip=INTEGER(ptr), *ir=INTEGER(row), *which=INTEGER(methods);
  bool use[5]={false,false,false,false,false}; for(int z:methods) use[z-1]=true;
  std::vector<std::vector<double>> ri(nl);
  for(int j=0;j<nl;++j) {NumericVector r=indices[j];ri[j]=as<std::vector<double>>(r);}
  std::vector<double> n(nq,0), m(nl,0), qnorm(nq,0), lnorm(nl,0),
    qwnorm(nq,0), lwnorm(nl,0), qp((size_t)d*nq), qw(qp.size()), qs(qp.size()),
    lp(values.size()), lw(values.size()), ls(values.size()), weight(d);
  double largest=*std::max_element(mass.begin(),mass.end());
  for(int k=0;k<d;++k) weight[k]=std::pow(mass[k]/largest,3.);
  for(int i=0;i<nq;++i) {
    for(int k=0;k<d;++k) n[i]+=raw[k+(size_t)d*i];
    if(n[i]<=0) continue;
    for(int k=0;k<d;++k) {
      size_t at=k+(size_t)d*i;double z=raw[at]/n[i];
      qp[at]=z;qs[at]=std::sqrt(z);qw[at]=std::pow(z,.6)*weight[k];
      qnorm[i]+=z*z;qwnorm[i]+=qw[at]*qw[at];
    }
  }
  for(int j=0;j<nl;++j) {
    for(int a=ip[j];a<ip[j+1];++a) m[j]+=v[a];
    if(m[j]<=0) continue;
    for(int a=ip[j];a<ip[j+1];++a) {
      double z=v[a]/m[j];lp[a]=z;ls[a]=std::sqrt(z);lw[a]=std::pow(z,.6)*weight[ir[a]];
      lnorm[j]+=z*z;lwnorm[j]+=lw[a]*lw[a];
    }
  }
  std::vector<std::vector<Match>> output((size_t)nq*nm);
#ifdef _OPENMP
  threads=std::max(1,std::min(threads,omp_get_max_threads()));
#pragma omp parallel for num_threads(threads) schedule(dynamic)
#endif
  for(int i=0;i<nq;++i) {
    if(n[i]<=0) continue;
    for(int j=0;j<nl;++j) {
      if(m[j]<=0) continue;
      double near=NA_REAL, delta=INFINITY;
      for(double r:ri[j]) if(std::isfinite(r)&&std::isfinite(qr[i])&&std::abs(r-qr[i])<delta) {
        delta=std::abs(r-qr[i]);near=r;
      }
      if(std::isfinite(tolerance)&&delta>tolerance) continue;
      double dot=0,wdot=0,bc=0,cross=0,count_cross=0;
      for(int a=ip[j];a<ip[j+1];++a) {
        size_t at=ir[a]+(size_t)d*i;
        if(qp[at]<=0) continue;
        if(use[1]) dot+=qp[at]*lp[a]; if(use[0]) wdot+=qw[at]*lw[a]; if(use[2]) bc+=qs[at]*ls[a];
        if(use[3]) cross+=entropy_cross(qp[at],lp[a]); if(use[4]) count_cross+=entropy_cross(raw[at],v[a]);
      }
      double distance[5]={
        std::max(0.,1-wdot*wdot/(qwnorm[i]*lwnorm[j])),
        std::max(0.,1-dot*dot/(qnorm[i]*lnorm[j])),
        std::sqrt(std::max(0.,1-std::min(1.,bc))),
        std::sqrt(std::max(0.,1-cross/(2*std::log(2.)))),
        2*std::max(0.,entropy_cross(n[i],m[j])-count_cross)};
      for(int z=0;z<nm;++z) {
        auto &out=output[i+(size_t)nq*z];Match hit={j+1,distance[which[z]-1],near};
        if((int)out.size()==top&&!better(hit,out.back(),qr[i])) continue;
        auto pos=std::lower_bound(out.begin(),out.end(),hit,
          [i,qr](const Match& a,const Match& b){return better(a,b,qr[i]);});
        out.insert(pos,hit);if((int)out.size()>top) out.pop_back();
      }
    }
  }
  int length=0;for(auto &a:output) length+=a.size();
  IntegerVector qi(length),id(length),method(length);NumericVector dist(length),rval(length),total(length);
  int at=0;
  for(int z=0;z<nm;++z) for(int i=0;i<nq;++i) for(auto &a:output[i+(size_t)nq*z]) {
    qi[at]=i+1;id[at]=a.id;method[at]=which[z];dist[at]=a.distance;rval[at]=a.ri;total[at]=m[a.id-1];++at;
  }
  return DataFrame::create(_["query"]=qi,_["library_id"]=id,_["method"]=method,
    _["distance"]=dist,_["library_RI"]=rval,_["library_total"]=total);
}
