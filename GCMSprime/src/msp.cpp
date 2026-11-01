#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cstdlib>
struct MSPReader {
  std::ifstream in; std::string pending; size_t line=0; bool held=false;
  MSPReader(const std::string& path):in(path.c_str(),std::ios::binary) {
    if(!in) Rcpp::stop("Cannot open MSP file.");
  }
};
inline std::string trim(std::string s) {
  auto a=s.find_first_not_of(" \t\r\n");if(a==std::string::npos) return "";
  return s.substr(a,s.find_last_not_of(" \t\r\n")-a+1);
}
inline std::string key(std::string s) {
  std::string out;for(unsigned char c:s) if(std::isalnum(c)) out+=std::toupper(c);return out;
}
// [[Rcpp::export]]
SEXP msp_open_cpp(std::string path) { return Rcpp::XPtr<MSPReader>(new MSPReader(path),true); }
// [[Rcpp::export]]
Rcpp::List msp_read_cpp(SEXP handle,int count=1000,bool strict=true) {
  Rcpp::XPtr<MSPReader> reader(handle);Rcpp::List out(count);int filled=0;
  for(int record=0;record<count;++record) {
    std::vector<std::string> tags,values;std::vector<double> mass,intensity;
    std::string s;int declared=-1;size_t first=0;
    while(true) {
      if(reader->held) {s=reader->pending;reader->held=false;}
      else {if(!std::getline(reader->in,s)) break;++reader->line;}
      if(reader->line==1 && s.size()>=3 && (unsigned char)s[0]==239 && (unsigned char)s[1]==187 && (unsigned char)s[2]==191) s=s.substr(3);
      s=trim(s);if(s.empty()) {if(!tags.empty()&&(!mass.empty()||declared==0)) break;else continue;}
      if(!first) first=reader->line;
      size_t colon=s.find(':');
      if(colon!=std::string::npos && std::isalpha((unsigned char)s[0])) {
        std::string tag=trim(s.substr(0,colon)),value=trim(s.substr(colon+1)),k=key(tag);
        if(k=="NAME" && !tags.empty()) {reader->pending=s;reader->held=true;break;}
        tags.push_back(tag);values.push_back(value);
        if(k=="NUMPEAKS"||k=="NUMBEROFPEAKS") {try {declared=std::stoi(value);} catch(...) {Rcpp::stop("Invalid Num Peaks near line %d.",(int)reader->line);}}
      } else {
        const char* cursor=s.c_str();char* end;
        while(*cursor) {
          while(*cursor && (std::isspace((unsigned char)*cursor)||*cursor==';'||*cursor==',')) ++cursor;
          if(!*cursor) break;
          double mz=std::strtod(cursor,&end);if(end==cursor) break;cursor=end;
          while(*cursor&&(std::isspace((unsigned char)*cursor)||*cursor==':')) ++cursor;
          double value=std::strtod(cursor,&end);
          if(end==cursor||!R_FINITE(mz)||!R_FINITE(value)||mz<0||value<0) Rcpp::stop("Invalid peak pair near line %d.",(int)reader->line);
          mass.push_back(mz);intensity.push_back(value);cursor=end;
          while(*cursor&&std::isspace((unsigned char)*cursor)) ++cursor;
          if(*cursor=='"') {++cursor;while(*cursor&&*cursor!='"') ++cursor;if(*cursor) ++cursor;}
        }
      }
    }
    if(tags.empty()&&mass.empty()) break;
    if(strict && declared>=0 && declared!=(int)mass.size()) Rcpp::stop("Num Peaks mismatch for record beginning at line %d: declared %d, read %d.",(int)first,declared,(int)mass.size());
    if(mass.empty()&&declared!=0) Rcpp::stop("No peaks found for record beginning at line %d.",(int)first);
    out[filled++]=Rcpp::List::create(Rcpp::_["tags"]=tags,Rcpp::_["values"]=values,
      Rcpp::_["mass"]=mass,Rcpp::_["intensity"]=intensity,Rcpp::_["line"]=(double)first);
  }
  Rcpp::List answer(filled);for(int k=0;k<filled;++k) answer[k]=out[k];return answer;
}

// [[Rcpp::export]]
void msp_close_cpp(SEXP handle) {
  MSPReader* reader=(MSPReader*)R_ExternalPtrAddr(handle);
  if(reader) {delete reader;R_ClearExternalPtr(handle);}
}
