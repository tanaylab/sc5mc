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

// ////////////////////////////////////////////////////////////////////////
// // [[Rcpp::export]]
// DataFrame create_sparse_matrix(SEXP cgdb, const IntegerVector& idxs, const std::vector<std::string>& cells){
//     Rcpp::XPtr<CGDB> ptr(cgdb);
//     return(ptr->create_sparse_matrix(idxs, cells));
// }

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

// #include <mkl.h>
// // [[Rcpp::export]]
// void test(){
//     std::vector<float> values = {1,  -2,  -4,  -1,  5,   8,   4,   2,   -3,  6,   7,   4,   -5};
//     std::vector<int> rows =     {1,   2,   4,   1,  2,   5,   3,   4,    1,  3,   4,   3,    5};
//     std::vector<int> pointerB = {1,   4,   7,   9,   12};
//     std::vector<int> pointerE = {4,   7,   9,   12,  14};

//     sparse_matrix_t mat;

//     mkl_scsrmm (const char *transa , const MKL_INT *m , const MKL_INT *n , const MKL_INT *k , const float *alpha , const char *matdescra , const float *val , const MKL_INT *indx , const MKL_INT *pntrb , const MKL_INT *pntre , const float *b , const MKL_INT *ldb , const float *beta , float *c , const MKL_INT *ldc )
// }

    // mkl_sparse_s_create_csc (&mat, SPARSE_INDEX_BASE_ONE, 5, 5, &pointerE[0], &pointerB[0], &rows[0], &values[0]);
    // sparse_matrix_t res = NULL;
    // sparse_status_t status;
    // status = mkl_sparse_spmm (SPARSE_OPERATION_TRANSPOSE, mat, mat, &res);

    // std::cout << status << std::endl;

    // MKL_INT n_rows=0, n_cols=0;
    // sparse_index_base_t indexing;
    // MKL_INT *cols_start = NULL, *cols_end = NULL, *row_indx = NULL;
    // float  *values_C = NULL;
    // mkl_sparse_s_export_csc(res, &indexing, &n_rows, &n_cols, &cols_start, &cols_end, &row_indx, &values_C);

    // if (res == NULL){
    //     std::cout << "is null" << std::endl;
    // }
    // // sparse_status_t mkl_sparse_s_export_csc (const sparse_matrix_t source, sparse_index_base_t *indexing, MKL_INT *rows, MKL_INT *cols, MKL_INT **cols_start, MKL_INT **cols_end, MKL_INT **row_indx, float **values);
    // printf("\nrows = %i , cols = %i\n", n_rows, n_cols);


// sparse_status_t mkl_sparse_s_create_csc (sparse_matrix_t *A, sparse_index_base_t indexing, MKL_INT rows, MKL_INT cols, MKL_INT *cols_start, MKL_INT *cols_end, MKL_INT *row_indx, float *values)
// NumericMatrix test(IntegerVector i_v, IntegerVector j_v, std::vector<float> x_v){