#include "CGDB.h"
#include "ProgressReporter.h"
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////
void CGDB::add_cell_data(const std::string& cell){
    unsigned point = cell.find(".");
    std::string fname = m_db_dir + "/" +  cell.substr(0, point) + "/" + cell.substr(point+1);

    std::string idx_fname = fname + ".idx.bin";
    std::string met_fname = fname + ".meth.bin";
    std::string cov_fname = fname + ".cov.bin";
    unsigned ncpgs = get_file_length(idx_fname) / 4;

    m_cell_idx[cell] = (int*)mmap_file(idx_fname, ncpgs*4);        
    m_cell_met[cell] = (float*)mmap_file(met_fname, ncpgs*4);        
    m_cell_cov[cell] = (float*)mmap_file(cov_fname, ncpgs*4);
    m_ncpgs[cell] = ncpgs;
    // Rcout << "added "  << cell << std::endl;
}


////////////////////////////////////////////////////////////////////////
unsigned CGDB::get_cell_data(const std::string& cell, int*& cell_idx, float*& cell_met, float*& cell_cov){

    if (m_cell_idx.find(cell) == m_cell_idx.end()){        
        add_cell_data(cell);
    }

    cell_idx = m_cell_idx[cell];
    cell_met = m_cell_met[cell];
    cell_cov = m_cell_cov[cell];

    int ncpgs = m_ncpgs[cell];
    return(ncpgs);
}

////////////////////////////////////////////////////////////////////////
DataFrame CGDB::mean_meth(const IntegerVector& idxs, const std::vector<std::string>& cells){
	if (!valid_indexes(idxs)){
		stop("some indexes are out of scope");
	}

    std::vector<float> mask(m_CPG_NUM+1, 0);
    std::vector<float> ones(m_CPG_NUM+1, 1);    

    vsUnpackV(idxs.length(), &ones[0], &mask[0], idxs.begin());

    NumericVector meth(cells.size());
    NumericVector cov(cells.size());

    ProgressReporter progress;
    progress.init(cells.size(), 1);

    for (unsigned i=0; i<cells.size(); ++i) {
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;
        unsigned ncpgs = get_cell_data(cells[i], cell_idx, cell_met, cell_cov);

        

        if ((cell_idx == NULL) || (cell_met == NULL) || (cell_cov == NULL) ) {
            meth[i] = NA_REAL;
            cov[i] = NA_REAL;
        }
        else {
            meth[i] = cblas_sdoti(ncpgs, cell_met, cell_idx, &mask[0]);            
            cov[i] = cblas_sdoti(ncpgs, cell_cov, cell_idx, &mask[0]);
        }
        progress.report(1);
    }    
    progress.report_last();

    return DataFrame::create(_["cell"]=cells, _["cov"]=cov, _["meth"]=meth);
}

////////////////////////////////////////////////////////////////////////
DataFrame CGDB::bin_meth(const IntegerVector& idxs, const IntegerVector& bins, const std::vector<std::string>& cells){
	if (!valid_indexes(idxs)){
		stop("some indexes are out of scope");
	}

	if (idxs.length() != bins.length()){
		stop("indexes and bins have different lengths");
	}

    std::vector<float> bins_full(m_CPG_NUM+1, 0);    
    
    vsUnpackV(idxs.length(), &as<std::vector<float> >(bins)[0], &bins_full[0], idxs.begin());
    
    unsigned int max_bin = *std::max_element(bins.begin(), bins.end());  
    NumericVector meth(max_bin + 1);
    NumericVector cov(max_bin + 1);

    ProgressReporter progress;
    progress.init(cells.size(), 1);

    for (unsigned i=0; i<cells.size(); ++i) {
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;
        unsigned ncpgs = get_cell_data(cells[i], cell_idx, cell_met, cell_cov);
        
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
        progress.report(1);
    }    
    progress.report_last();
    
    // remove first element (sum of cpgs wihout a bin)
    cov.erase(0);
    meth.erase(0);
    
    std::vector<int> bin_id(max_bin);
    std::iota(bin_id.begin(), bin_id.end(), 1);

    return DataFrame::create(_["bin"] = bin_id, _["cov"] = cov, _["meth"] = meth);
}

////////////////////////////////////////////////////////////////////////
List CGDB::bin_meth_per_cell_cpp(const IntegerVector& idxs, const IntegerVector& bins, const std::vector<std::string>& cells){

	if (!valid_indexes(idxs)){
		stop("some indexes are out of scope");
	}

	if (idxs.length() != bins.length()){
		stop("indexes and bins have different lengths");
	}


    std::vector<float> bins_full(m_CPG_NUM+1, 0);    
    
    vsUnpackV(idxs.length(), &as<std::vector<float> >(bins)[0], &bins_full[0], idxs.begin());
    
    unsigned int max_bin = *std::max_element(bins.begin(), bins.end());  

    NumericMatrix meth(max_bin + 1, cells.size());
    NumericMatrix cov(max_bin + 1, cells.size());

    ProgressReporter progress;
    progress.init(cells.size(), 1);

    for (unsigned i=0; i<cells.size(); ++i) {
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;
        unsigned ncpgs = get_cell_data(cells[i], cell_idx, cell_met, cell_cov);
   
        if ((cell_idx != NULL) && (cell_met != NULL) && (cell_cov != NULL)) {

            // get bins for the covered CpGs        
            std::vector<float > cell_bins(ncpgs);
            
            vsPackV(ncpgs, &bins_full[0], cell_idx, &cell_bins[0]);
            
            float* cov_j = cell_cov;
            float* meth_j = cell_met;
            for (auto & j : cell_bins){
                cov(j, i) += *(cov_j);
                meth(j, i) += *(meth_j);
                ++cov_j;
                ++meth_j;
            }            
            
        }
        progress.report(1);         
    } 
    progress.report_last();   
    
    std::vector<int> bin_id(max_bin);
    std::iota(bin_id.begin(), bin_id.end(), 1);

    List res;
    res["meth"] = meth;
    res["cov"] = cov;
    res["bin"] = bin_id;
    return res;     
}

List CGDB::extract_sparse(const IntegerVector& idxs, const std::vector<std::string>& cells, const std::vector<std::string>& chrom, const std::vector<int>& start, const std::vector<int>& end){

	if (!valid_indexes(idxs)){
		stop("some indexes are out of scope");
	}

	List res;
	return res;
}

List CGDB::extract(const IntegerVector& idxs, const std::vector<std::string>& cells){
	if (!valid_indexes(idxs)){
		stop("some indexes are out of scope");
	}

	std::vector<float> all_cov(m_CPG_NUM+1, 0);	
	std::vector<float> all_meth(m_CPG_NUM+1, 0);
	std::vector<int> all_idxs = as<std::vector<int> >(idxs);	

	std::vector<std::vector<float> > cov(cells.size(), std::vector<float>(idxs.length(), 0));
    std::vector<std::vector<float> > meth(cells.size(), std::vector<float>(idxs.length(), 0));

    ProgressReporter progress;
    progress.init(cells.size(), 1);

    for (unsigned i=0; i < cells.size(); ++i) {
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;
        unsigned ncpgs = get_cell_data(cells[i], cell_idx, cell_met, cell_cov);
        if ((cell_idx != NULL) && (cell_met != NULL) && (cell_cov != NULL)) {
        	// scatter cell coverage to all_cov
        	vsUnpackV(ncpgs, cell_cov, &all_cov[0], cell_idx);

        	// gather cell coverage in idxs to result 2d array cov
        	vsPackV(idxs.length(), &all_cov[0], &all_idxs[0], &cov[i][0]);

        	// clean all_cov vector
        	cblas_sscal(all_cov.size(), 0, &all_cov[0], 1);
        	
        	// scatter cell methyaltion to all_cov
        	vsUnpackV(ncpgs, cell_met, &all_meth[0], cell_idx);

        	// gather cell methylation in idxs to result 2d array meth
        	vsPackV(idxs.length(), &all_meth[0], &all_idxs[0], &meth[i][0]);

        	// clean all_meth vector
        	cblas_sscal(all_meth.size(), 0, &all_meth[0], 1);        	
        }
        progress.report(1);
    }
    progress.report_last();

    Function cbind("cbind");
    Function do_call("do.call");

    NumericMatrix cov_mat = do_call(cbind, wrap(cov));
    NumericMatrix meth_mat = do_call(cbind, wrap(meth));

	List res;	
	res["cov"] = cov_mat;
	res["meth"] = meth_mat;
	return res;	
}

////////////////////////////////////////////////////////////////////////
unsigned get_file_length(const std::string& fname){
    struct stat st;
    stat(fname.c_str(), &st);

    return st.st_size;
}


////////////////////////////////////////////////////////////////////////
void* mmap_file(const std::string& fname, unsigned length){
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