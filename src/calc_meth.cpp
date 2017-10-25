#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <Rcpp.h>
#include <mkl.h>
using namespace Rcpp;

unsigned get_file_length(std::string const& fname);
void* mmap_file(std::string const& fname, unsigned length);
void unmap_file(void* map, unsigned length);
unsigned get_cell_data(std::string const& cell, std::string const& db_dir, int*& cell_idx, float*& cell_met, float*& cell_cov);

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
DataFrame mean_meth(IntegerVector const& idxs, 
                    std::string const& db_dir, 
                    std::vector<std::string> const& cells, 
                    const int& CPG_NUM){

    std::vector<float> mask(CPG_NUM+1, 0);
    std::vector<float> ones(CPG_NUM+1, 1);    

    vsUnpackV(idxs.length(), &ones[0], &mask[0], idxs.begin());

    NumericVector meth(cells.size());
    NumericVector cov(cells.size());

    for (unsigned i=0; i<cells.size(); ++i) {
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;
        unsigned ncpgs = get_cell_data(cells[i], db_dir, cell_idx, cell_met, cell_cov);

        if ((cell_idx == NULL) || (cell_met == NULL) || (cell_cov == NULL) ) {
            meth[i] = NA_REAL;
            cov[i] = NA_REAL;
        }
        else {
            meth[i] = cblas_sdoti(ncpgs, cell_met, cell_idx, &mask[0]);            
            cov[i] = cblas_sdoti(ncpgs, cell_cov, cell_idx, &mask[0]);
        }

        unmap_file(cell_met, ncpgs*4);
        unmap_file(cell_cov, ncpgs*4);
        unmap_file(cell_idx, ncpgs*4);
    }    

    return DataFrame::create(_["cell"]=cells, _["cov"]=cov, _["meth"]=meth);
}

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
DataFrame bin_meth(IntegerVector const& idxs, 
                   IntegerVector const& bins, 
                   std::string const& db_dir, 
                   std::vector<std::string> const& cells, 
                   const int& CPG_NUM){

    std::vector<float> bins_full(CPG_NUM+1, 0);    
    
    vsUnpackV(idxs.length(), &as<std::vector<float> >(bins)[0], &bins_full[0], idxs.begin());
    
    unsigned int max_bin = *std::max_element(bins.begin(), bins.end());  
    NumericVector meth(max_bin + 1);
    NumericVector cov(max_bin + 1);

    for (unsigned i=0; i<cells.size(); ++i) {
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;
        unsigned ncpgs = get_cell_data(cells[i], db_dir, cell_idx, cell_met, cell_cov);
        
        if ((cell_idx != NULL) && (cell_met != NULL) && (cell_cov != NULL)) {

            // get bins for the covered CpGs        
            std::vector<float > cell_bins(ncpgs);
            
            vsPackV(ncpgs, &bins_full[0], cell_idx, &cell_bins[0]);
            
            float* cov_j = cell_cov;
            float* meth_j = cell_met;
            for (auto & j : cell_bins){
                cov[j]+= *(cov_j);
                meth[j]+= *(meth_j);
                ++cov_j;
                ++meth_j;
            }            
            
        } 
        
        unmap_file(cell_met, ncpgs*4);        
        unmap_file(cell_cov, ncpgs*4);        
        unmap_file(cell_idx, ncpgs*4);        
        
    }    
    
    // remove first element (sum of cpgs wihout a bin)
    cov.erase(0);
    meth.erase(0);
    
    std::vector<int> bin_id(max_bin);
    std::iota(bin_id.begin(), bin_id.end(), 1);

    return DataFrame::create(_["bin"] = bin_id, _["cov"] = cov, _["meth"] = meth);    
}


////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List bin_meth_per_cell_cpp(IntegerVector const& idxs, 
                   IntegerVector const& bins, 
                   std::string const& db_dir, 
                   std::vector<std::string> const& cells, 
                   const int& CPG_NUM){

    std::vector<float> bins_full(CPG_NUM+1, 0);    
    
    vsUnpackV(idxs.length(), &as<std::vector<float> >(bins)[0], &bins_full[0], idxs.begin());
    
    unsigned int max_bin = *std::max_element(bins.begin(), bins.end());  

    NumericMatrix meth(cells.size(), max_bin + 1);
    NumericMatrix cov(cells.size(), max_bin + 1);

    for (unsigned i=0; i<cells.size(); ++i) {
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;
        unsigned ncpgs = get_cell_data(cells[i], db_dir, cell_idx, cell_met, cell_cov);
   
        if ((cell_idx != NULL) && (cell_met != NULL) && (cell_cov != NULL)) {

            // get bins for the covered CpGs        
            std::vector<float > cell_bins(ncpgs);
            
            vsPackV(ncpgs, &bins_full[0], cell_idx, &cell_bins[0]);
            
            float* cov_j = cell_cov;
            float* meth_j = cell_met;
            for (auto & j : cell_bins){
                cov(i, j) += *(cov_j);
                meth(i, j) += *(meth_j);
                ++cov_j;
                ++meth_j;
            }            
            
        } else {
            Rcout << "asfasf" << std::endl;
        }
        
        unmap_file(cell_met, ncpgs*4);        
        unmap_file(cell_cov, ncpgs*4);        
        unmap_file(cell_idx, ncpgs*4);        
        
    }    
    
    std::vector<int> bin_id(max_bin);
    std::iota(bin_id.begin(), bin_id.end(), 1);

    List res;
    res["meth"] = meth;
    res["cov"] = cov;
    res["bin"] = bin_id;
    return res;     
}

////////////////////////////////////////////////////////////////////////
unsigned get_file_length(std::string const& fname){
    struct stat st;
    stat(fname.c_str(), &st);

    return st.st_size;
}


////////////////////////////////////////////////////////////////////////
void* mmap_file(std::string const& fname, unsigned length){
    int fd = open(fname.c_str(), O_RDONLY);
    if (fd < 0) {
        return NULL;
    }
    void* map = mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    return map;
}


////////////////////////////////////////////////////////////////////////
void unmap_file(void* map, unsigned length){
    if (map != NULL) {
        munmap(map, length);
        // map = NULL;
    }
}

////////////////////////////////////////////////////////////////////////
unsigned get_cell_data(std::string const& cell, std::string const& db_dir, int*& cell_idx, float*& cell_met, float*& cell_cov){
    unsigned point = cell.find(".");
    std::string fname = db_dir + "/" +  cell.substr(0, point) + "/" + cell.substr(point+1);

    std::string idx_fname = fname + ".idx.bin";
    std::string met_fname = fname + ".meth.bin";
    std::string cov_fname = fname + ".cov.bin";
    unsigned ncpgs = get_file_length(idx_fname) / 4;
        
    cell_idx = (int*)mmap_file(idx_fname, ncpgs*4);        
    cell_met = (float*)mmap_file(met_fname, ncpgs*4);        
    cell_cov = (float*)mmap_file(cov_fname, ncpgs*4);
    if (cell_idx == NULL || cell_met == NULL || cell_cov == NULL){
        Rcout << "yrs" << std::endl;
    }
    

    return(ncpgs);
}




// vsUnpackV( n, a, y, iy );
// n - Specifies the number of elements to be calculated.
// a (const float*) - Specifies the pointer to an array of size at least n that contains the input vector a.
// y (float*) - Specifies the pointer to an array that contains the output vector y.
// iy (const int*) - Specifies the pointer to an array of size at least n that contains the index vector for the elements of a.

// vsPackV( n, a, ia, y );
// n - Specifies the number of elements to be calculated.
// a (const float*) - Specifies pointer to an array that contains the input vector a. 
// ia (const int*) - Specifies the pointer to an array of size at least n that contains the index vector for the elements of a.
// y (float *) - Pointer to an array of size at least n that contains the output vector y.