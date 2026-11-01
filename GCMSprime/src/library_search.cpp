#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
struct Match { int id; double distance, ri, ordering; };
static bool better(const Match& a, const Match& b, double query_ri) {
  if(a.ordering!=b.ordering) return a.ordering<b.ordering;
  double da=std::isfinite(a.ri)&&std::isfinite(query_ri)?std::abs(a.ri-query_ri):INFINITY;
  double db=std::isfinite(b.ri)&&std::isfinite(query_ri)?std::abs(b.ri-query_ri):INFINITY;
  return da!=db?da<db:a.id<b.id;
}
// Twice the scalar Poisson Bregman gap, stable near equality. Combining
// nonnegative gaps avoids cancellation amplified by large query count totals.
static double poisson_gap(double y,double mu) {
  if(y==mu) return 0.;
  if(y<=0) return 2*mu;
  double z=(y-mu)/mu;
  if(std::abs(z)<1e-4)
    return mu*z*z*(1-z/3+z*z/6-z*z*z/10+z*z*z*z/15);
  return 2*std::max(0.,y*std::log(y/mu)-(y-mu));
}
static int search_threads(int threads) {
#ifdef _OPENMP
  return std::max(1,std::min(threads,omp_get_max_threads()));
#else
  return 1;
#endif
}
// All R objects and work buffers are prepared before entering OpenMP workers.
// [[Rcpp::export]]
DataFrame library_search_cpp(NumericMatrix query, S4 library, NumericVector mass,
    List indices, NumericVector query_ri, double tolerance, int top,
    IntegerVector methods, int threads, NumericVector ion_weights,
    LogicalVector count_known, int reference_mode) {
  IntegerVector ptr=library.slot("p"),row=library.slot("i");NumericVector values=library.slot("x");
  const int d=query.nrow(),nq=query.ncol(),nl=ptr.size()-1,nm=methods.size();
  const double *raw=REAL(query),*v=REAL(values),*qr=REAL(query_ri),*w=REAL(ion_weights);
  const int *ip=INTEGER(ptr),*ir=INTEGER(row),*which=INTEGER(methods),*known=LOGICAL(count_known);
  bool use[7]={false,false,false,false,false,false,false};for(int z:methods) use[z-1]=true;
  bool need_cos=use[1]||use[5],need_stein=use[0]||use[6];
  double wmax=*std::max_element(ion_weights.begin(),ion_weights.end());
  std::vector<std::vector<double>> ri(nl);
  for(int j=0;j<nl;++j) {NumericVector r=indices[j];ri[j]=as<std::vector<double>>(r);}
  std::vector<double> n(nq),m(nl),qnorm(nq),lnorm(nl),qwnorm(nq),lwnorm(nl),
    qa(nq),la(nl),qp((size_t)d*nq),qw(qp.size()),qs(qp.size()),
    lp(values.size()),lw(values.size()),ls(values.size()),mass_weight(d);
  double largest=*std::max_element(mass.begin(),mass.end());
  for(int k=0;k<d;++k) mass_weight[k]=std::pow(mass[k]/largest,3.);
  for(int i=0;i<nq;++i) {
    for(int k=0;k<d;++k) n[i]+=raw[k+(size_t)d*i];
    if(n[i]<=0) continue;
    for(int k=0;k<d;++k) {
      size_t at=k+(size_t)d*i;double z=raw[at]/n[i];qp[at]=z;qa[i]+=w[k]*z;
      if(use[2]) qs[at]=std::sqrt(z);
      if(need_cos) qnorm[i]+=w[k]*z*z;
      if(need_stein) {qw[at]=std::pow(z,.6)*mass_weight[k];qwnorm[i]+=w[k]*qw[at]*qw[at];}
    }
  }
  for(int j=0;j<nl;++j) {
    for(int a=ip[j];a<ip[j+1];++a) m[j]+=v[a];
    if(m[j]<=0) continue;
    for(int a=ip[j];a<ip[j+1];++a) {
      int k=ir[a];double z=v[a]/m[j];lp[a]=z;la[j]+=w[k]*z;
      if(use[2]) ls[a]=std::sqrt(z);
      if(need_cos) lnorm[j]+=w[k]*z*z;
      if(need_stein) {lw[a]=std::pow(z,.6)*mass_weight[k];lwnorm[j]+=w[k]*lw[a]*lw[a];}
    }
  }
  std::vector<std::vector<Match>> output((size_t)nq*nm);threads=search_threads(threads);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic)
#endif
  for(int i=0;i<nq;++i) {
    if(n[i]<=0||qa[i]<=0) continue;
    for(int j=0;j<nl;++j) {
      if(m[j]<=0||la[j]<=0) continue;
      double near=NA_REAL,delta=INFINITY;
      for(double r:ri[j]) if(std::isfinite(r)&&std::isfinite(qr[i])&&std::abs(r-qr[i])<delta) {
        delta=std::abs(r-qr[i]);near=r;
      }
      if(std::isfinite(tolerance)&&delta>tolerance) continue;
      bool balanced=reference_mode==1||(reference_mode==0&&!known[j]);
      double dot=0,wdot=0,hgap=0,jgap=0,dgap=0,shared_q=0,shared_l=0;
      bool need_gap=use[2]||use[3]||use[4];
      double alpha=n[i]>=m[j]?1/(1+m[j]/n[i]):(n[i]/m[j])/(1+n[i]/m[j]);
      double beta=m[j]>=n[i]?1/(1+n[i]/m[j]):(m[j]/n[i])/(1+m[j]/n[i]);
      for(int a=ip[j];a<ip[j+1];++a) {
        int k=ir[a];size_t at=k+(size_t)d*i;
        if(qp[at]<=0||w[k]<=0) continue;
        if(need_cos) dot+=w[k]*qp[at]*lp[a];
        if(need_stein) wdot+=w[k]*qw[at]*lw[a];
        if(need_gap) {shared_q+=w[k]*qp[at];shared_l+=w[k]*lp[a];}
        if(use[2]) {double diff=qs[at]-ls[a];hgap+=.5*w[k]*diff*diff;}
        if(use[3]||(use[4]&&balanced)) {
          double r=.5*qp[at]+.5*lp[a];
          jgap+=.25*w[k]*(poisson_gap(qp[at],r)+poisson_gap(lp[a],r));
        }
        if(use[4]&&!balanced) {
          double r=alpha*qp[at]+beta*lp[a];
          dgap+=w[k]*(n[i]*poisson_gap(qp[at],r)+m[j]*poisson_gap(lp[a],r));
        }
      }
      double cosine=qnorm[i]>0&&lnorm[j]>0?std::min(1.,std::max(0.,(dot/std::sqrt(qnorm[i]))/std::sqrt(lnorm[j]))):NA_REAL;
      double stein=qwnorm[i]>0&&lwnorm[j]>0?std::min(1.,std::max(0.,(wdot/std::sqrt(qwnorm[i]))/std::sqrt(lwnorm[j]))):NA_REAL;
      double unmatched_q=std::max(0.,qa[i]-shared_q),unmatched_l=std::max(0.,la[j]-shared_l);
      double js=jgap+.5*std::log(2.)*(unmatched_q+unmatched_l);
      double distance[7]={std::isfinite(stein)?std::max(0.,1-stein*stein):INFINITY,
        std::isfinite(cosine)?std::max(0.,1-cosine*cosine):INFINITY,
        std::sqrt((hgap+.5*(unmatched_q+unmatched_l))/wmax),
        std::sqrt(js/(wmax*std::log(2.))),
        balanced?4*n[i]*js:dgap+2*(n[i]*unmatched_q*std::log1p(m[j]/n[i])+m[j]*unmatched_l*std::log1p(n[i]/m[j])),0,0};
      for(int z=0;z<nm;++z) {
        int method=which[z]-1;
        // Select angular hits with the cosine key to preserve ordering even at
        // machine precision; calculate acos only for retained candidates.
        double key=method==5?distance[1]:method==6?distance[0]:distance[method];
        if(!std::isfinite(key)) continue;
        auto &out=output[i+(size_t)nq*z];Match hit={j+1,distance[method],near,key};
        if((int)out.size()==top&&!better(hit,out.back(),qr[i])) continue;
        if(method==5) hit.distance=std::acos(cosine);
        if(method==6) hit.distance=std::acos(stein);
        auto pos=std::lower_bound(out.begin(),out.end(),hit,
          [i,qr](const Match& a,const Match& b){return better(a,b,qr[i]);});
        out.insert(pos,hit);if((int)out.size()>top) out.pop_back();
      }
    }
  }
  int length=0;for(auto &a:output) length+=a.size();
  IntegerVector qi(length),id(length),method(length);NumericVector dist(length),rval(length),total(length),effective(length),ordering(length);
  int at=0;
  for(int z=0;z<nm;++z) for(int i=0;i<nq;++i) for(auto &a:output[i+(size_t)nq*z]) {
    int j=a.id-1;bool balanced=reference_mode==1||(reference_mode==0&&!known[j]);
    qi[at]=i+1;id[at]=a.id;method[at]=which[z];dist[at]=a.distance;rval[at]=a.ri;
    total[at]=m[j];effective[at]=balanced?n[i]:m[j];ordering[at]=a.ordering;++at;
  }
  return DataFrame::create(_["query"]=qi,_["library_id"]=id,_["method"]=method,
    _["distance"]=dist,_["library_RI"]=rval,_["library_total"]=total,
    _["reference_total"]=effective,_["ordering"]=ordering);
}

// [[Rcpp::export]]
List ion_frequency_cpp(S4 library,double relative,int threads) {
  IntegerVector ptr=library.slot("p"),row=library.slot("i"),dim=library.slot("Dim");
  NumericVector val=library.slot("x");int d=dim[0],nl=dim[1];threads=search_threads(threads);
  const int *p=INTEGER(ptr),*r=INTEGER(row);const double *x=REAL(val);
  std::vector<int> freq((size_t)d*threads),used(threads);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for(int j=0;j<nl;++j) {
    int worker=0;
#ifdef _OPENMP
    worker=omp_get_thread_num();
#endif
    double peak=0;for(int a=p[j];a<p[j+1];++a) peak=std::max(peak,x[a]);
    if(peak<=0) continue;
    ++used[worker];for(int a=p[j];a<p[j+1];++a)
      if(x[a]>0&&x[a]>=relative*peak) ++freq[r[a]+(size_t)d*worker];
  }
  NumericVector count(d);int n=0;
  for(int worker=0;worker<threads;++worker) {
    n+=used[worker];for(int k=0;k<d;++k) count[k]+=freq[k+(size_t)d*worker];
  }
  return List::create(_["frequency"]=count,_["spectra"]=n);
}
