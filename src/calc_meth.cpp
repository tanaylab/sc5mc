// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include "CGDB.h"
using namespace Rcpp;


////////////////////////////////////////////////////////////////////////
RcppExport SEXP CGDB__new(SEXP db_dir_, SEXP CPG_NUM_){
    std::string db_dir = as<std::string>(db_dir_);
    int CPG_NUM = as<int>(CPG_NUM_);
    XPtr<CGDB> ptr(new CGDB(db_dir + "/data", CPG_NUM), true);
    return ptr;
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
DataFrame mean_meth(SEXP cgdb, const IntegerVector& idxs, const std::vector<std::string>& cells){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->mean_meth(idxs, cells));
}


////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
DataFrame bin_meth(SEXP cgdb, const IntegerVector& idxs, const IntegerVector& bins, const std::vector<std::string>& cells){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->bin_meth(idxs, bins, cells));    
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List extract_sc_data(SEXP cgdb, const IntegerVector& idxs, const std::vector<std::string>& cells){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->extract(idxs, cells));    
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
DataFrame extract_sc_data_sparse(SEXP cgdb, const IntegerVector& idxs, const std::vector<std::string>& cells){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->extract_sparse(idxs, cells));    
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List bin_meth_per_cell_cpp(SEXP cgdb, const IntegerVector& idxs, const IntegerVector& bins, const std::vector<std::string>& cells){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->bin_meth_per_cell_cpp(idxs, bins, cells));    
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
void freemem_cpp(SEXP cgdb){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->freemem());
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
std::vector<std::string> list_open_cells(SEXP cgdb){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->list_open_cells());   
}


