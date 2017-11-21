// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// mean_meth
DataFrame mean_meth(SEXP cgdb, const IntegerVector& idxs, const std::vector<std::string>& cells);
RcppExport SEXP _sc5mc_mean_meth(SEXP cgdbSEXP, SEXP idxsSEXP, SEXP cellsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type cgdb(cgdbSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type cells(cellsSEXP);
    rcpp_result_gen = Rcpp::wrap(mean_meth(cgdb, idxs, cells));
    return rcpp_result_gen;
END_RCPP
}
// bin_meth
DataFrame bin_meth(SEXP cgdb, const IntegerVector& idxs, const IntegerVector& bins, const std::vector<std::string>& cells);
RcppExport SEXP _sc5mc_bin_meth(SEXP cgdbSEXP, SEXP idxsSEXP, SEXP binsSEXP, SEXP cellsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type cgdb(cgdbSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type bins(binsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type cells(cellsSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_meth(cgdb, idxs, bins, cells));
    return rcpp_result_gen;
END_RCPP
}
// extract_sc_data
List extract_sc_data(SEXP cgdb, const IntegerVector& idxs, const std::vector<std::string>& cells);
RcppExport SEXP _sc5mc_extract_sc_data(SEXP cgdbSEXP, SEXP idxsSEXP, SEXP cellsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type cgdb(cgdbSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type cells(cellsSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_sc_data(cgdb, idxs, cells));
    return rcpp_result_gen;
END_RCPP
}
// extract_sc_data_sparse
DataFrame extract_sc_data_sparse(SEXP cgdb, const IntegerVector& idxs, const std::vector<std::string>& cells);
RcppExport SEXP _sc5mc_extract_sc_data_sparse(SEXP cgdbSEXP, SEXP idxsSEXP, SEXP cellsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type cgdb(cgdbSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type cells(cellsSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_sc_data_sparse(cgdb, idxs, cells));
    return rcpp_result_gen;
END_RCPP
}
// bin_meth_per_cell_cpp
List bin_meth_per_cell_cpp(SEXP cgdb, const IntegerVector& idxs, const IntegerVector& bins, const std::vector<std::string>& cells);
RcppExport SEXP _sc5mc_bin_meth_per_cell_cpp(SEXP cgdbSEXP, SEXP idxsSEXP, SEXP binsSEXP, SEXP cellsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type cgdb(cgdbSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type bins(binsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type cells(cellsSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_meth_per_cell_cpp(cgdb, idxs, bins, cells));
    return rcpp_result_gen;
END_RCPP
}
// freemem_cpp
void freemem_cpp(SEXP cgdb);
RcppExport SEXP _sc5mc_freemem_cpp(SEXP cgdbSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type cgdb(cgdbSEXP);
    freemem_cpp(cgdb);
    return R_NilValue;
END_RCPP
}
// list_open_cells
std::vector<std::string> list_open_cells(SEXP cgdb);
RcppExport SEXP _sc5mc_list_open_cells(SEXP cgdbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type cgdb(cgdbSEXP);
    rcpp_result_gen = Rcpp::wrap(list_open_cells(cgdb));
    return rcpp_result_gen;
END_RCPP
}
// shuffle_mat_marginals
NumericMatrix shuffle_mat_marginals(const NumericMatrix& m_meth, const NumericMatrix& inds_mat, const int& n_shuff);
RcppExport SEXP _sc5mc_shuffle_mat_marginals(SEXP m_methSEXP, SEXP inds_matSEXP, SEXP n_shuffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type m_meth(m_methSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type inds_mat(inds_matSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_shuff(n_shuffSEXP);
    rcpp_result_gen = Rcpp::wrap(shuffle_mat_marginals(m_meth, inds_mat, n_shuff));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP CGDB__new(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_sc5mc_mean_meth", (DL_FUNC) &_sc5mc_mean_meth, 3},
    {"_sc5mc_bin_meth", (DL_FUNC) &_sc5mc_bin_meth, 4},
    {"_sc5mc_extract_sc_data", (DL_FUNC) &_sc5mc_extract_sc_data, 3},
    {"_sc5mc_extract_sc_data_sparse", (DL_FUNC) &_sc5mc_extract_sc_data_sparse, 3},
    {"_sc5mc_bin_meth_per_cell_cpp", (DL_FUNC) &_sc5mc_bin_meth_per_cell_cpp, 4},
    {"_sc5mc_freemem_cpp", (DL_FUNC) &_sc5mc_freemem_cpp, 1},
    {"_sc5mc_list_open_cells", (DL_FUNC) &_sc5mc_list_open_cells, 1},
    {"_sc5mc_shuffle_mat_marginals", (DL_FUNC) &_sc5mc_shuffle_mat_marginals, 3},
    {"CGDB__new",                     (DL_FUNC) &CGDB__new,                     2},
    {NULL, NULL, 0}
};

RcppExport void R_init_sc5mc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
