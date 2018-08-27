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
DataFrame extract_sc_data_sparse_all(SEXP cgdb, const std::vector<std::string>& cells){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->extract_sparse_all(cells));    
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List bin_meth_per_cell_cpp(SEXP cgdb, const IntegerVector& idxs, const IntegerVector& bins, const std::vector<std::string>& cells){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->bin_meth_per_cell_cpp(idxs, bins, cells));    
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
std::vector<int> count_pairs_cpp(SEXP cgdb, const IntegerVector& idxs, const std::string& cell1, const std::string& cell2){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    std::vector<int> res(4, 0);
    ptr->count_pairs(idxs, cell1, cell2, res);
    return(res);   
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericMatrix count_pairs_all_cpp(SEXP cgdb, const IntegerVector& idxs, const std::vector<std::string>& cells1, const std::vector<std::string>& cells2){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->count_pairs_all(idxs, cells1, cells2));   
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

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
unsigned get_cell_ncpgs(SEXP cgdb, const std::string& cell){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->get_cell_ncpgs(cell));      
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
std::string get_cell_filename(SEXP cgdb, const std::string& cell){
    Rcpp::XPtr<CGDB> ptr(cgdb);
    return(ptr->get_cell_filename(cell));      
}

