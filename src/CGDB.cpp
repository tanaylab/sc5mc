#include "CGDB.h"
#include "ProgressReporter.h"
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////
void CGDB::add_cell_data(const std::string& cell){
    std::string fname = get_cell_filename(cell);
    // unsigned point = cell.find(".");
    // std::string fname = m_db_dir + "/" +  cell.substr(0, point) + "/" + cell.substr(point+1);

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
std::string CGDB::get_cell_filename(const std::string& cell){
    unsigned point = cell.find(".");
    std::string fname = m_db_dir + "/" +  cell.substr(0, point) + "/" + cell.substr(point+1);
    return(fname);
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
unsigned CGDB::get_cell_ncpgs(const std::string& cell){
    if (m_ncpgs.find(cell) != m_ncpgs.end()){        
        return(m_ncpgs[cell]);
    }

    unsigned point = cell.find(".");
    std::string fname = m_db_dir + "/" +  cell.substr(0, point) + "/" + cell.substr(point+1);

    std::string idx_fname = fname + ".idx.bin";    
    unsigned ncpgs = get_file_length(idx_fname) / 4;
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

    std::vector<double> bins_vec(bins.begin(), bins.end());
    std::vector<double> bins_full(m_CPG_NUM+1, 0);   
 
    vdUnpackV(idxs.length(), &bins_vec[0], &bins_full[0], idxs.begin());
 
    unsigned int max_bin = *std::max_element(bins.begin(), bins.end());  
    NumericVector meth(max_bin + 1);
    NumericVector cov(max_bin + 1);

    ProgressReporter progress;
    progress.init(cells.size(), 1);

    for (unsigned i=0; i < cells.size(); ++i) {
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;
        unsigned ncpgs = get_cell_data(cells[i], cell_idx, cell_met, cell_cov);
        
        if ((cell_idx != NULL) && (cell_met != NULL) && (cell_cov != NULL)) {

            // get bins for the covered CpGs        
            std::vector<double > cell_bins(ncpgs, 0);

            vdPackV(ncpgs, &bins_full[0], cell_idx, &cell_bins[0]);

            float* cov_j = cell_cov;
            float* meth_j = cell_met;

            for (auto & j : cell_bins){
                cov[(int)j]+= *(cov_j);
                meth[(int)j]+= *(meth_j);
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


    std::vector<double> bins_full(m_CPG_NUM+1, 0);    
    
    vdUnpackV(idxs.length(), &as<std::vector<double> >(bins)[0], &bins_full[0], idxs.begin());
    
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
            std::vector<double > cell_bins(ncpgs);
            
            vdPackV(ncpgs, &bins_full[0], cell_idx, &cell_bins[0]);
                        
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

////////////////////////////////////////////////////////////////////////
DataFrame CGDB::extract_sparse_all(const std::vector<std::string>& cells){
    unsigned tot_cpgs = 0;
    for (unsigned i = 0; i < cells.size(); ++i){
        tot_cpgs += get_cell_ncpgs(cells[i]);
    }

    std::vector<int> idxs_vec(tot_cpgs, 0);
    std::vector<float> cov_vec(tot_cpgs, 0);
    std::vector<float> meth_vec(tot_cpgs, 0);
    std::vector<std::string > cells_vec(tot_cpgs);

    auto idx_iter = idxs_vec.begin();
    auto cov_iter = cov_vec.begin();
    auto meth_iter = meth_vec.begin();
    auto cells_iter = cells_vec.begin();

    ProgressReporter progress;
    progress.init(cells.size(), 1);

    for (unsigned i=0; i < cells.size(); ++i) {        
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;        
        
        unsigned ncpgs = get_cell_data(cells[i], cell_idx, cell_met, cell_cov);

        if ((cell_idx != NULL) && (cell_met != NULL) && (cell_cov != NULL)) {

            std::copy(cell_idx, cell_idx + ncpgs, idx_iter);
            std::copy(cell_cov, cell_cov + ncpgs, cov_iter);
            std::copy(cell_met, cell_met + ncpgs, meth_iter);
            std::fill(cells_iter, cells_iter + ncpgs, cells[i]);
                        
            idx_iter += ncpgs;
            cov_iter += ncpgs;
            meth_iter += ncpgs;
            cells_iter += ncpgs;            
        }
        progress.report(1);
    }
    progress.report_last();   
    
    return DataFrame::create(_["cell_id"] = cells_vec, _["id"]=idxs_vec, _["cov"]=cov_vec, _["meth"]=meth_vec);        
}

////////////////////////////////////////////////////////////////////////
DataFrame CGDB::extract_sparse(const IntegerVector& idxs, const std::vector<std::string>& cells){

	if (!valid_indexes(idxs)){
		stop("some indexes are out of scope");
	}
    std::vector<float> all_cov(m_CPG_NUM+1, 0); 
    std::vector<float> all_meth(m_CPG_NUM+1, 0);
    std::vector<int> all_idxs = as<std::vector<int> >(idxs);    

    std::vector<std::vector<float> > cov_mat(cells.size(), std::vector<float>());
    std::vector<std::vector<float> > meth_mat(cells.size(), std::vector<float>());
    std::vector<std::vector<int> > idxs_mat(cells.size(), std::vector<int>());
    std::vector<std::vector<std::string > > cell_mat(cells.size(), std::vector<std::string>());

    ProgressReporter progress;
    progress.init(cells.size(), 1);

    std::vector<float> all_cell_cov(idxs.length(), 0);
    std::vector<float> all_cell_meth(idxs.length(), 0);    
    
    for (unsigned i=0; i < cells.size(); ++i) {        
        int* cell_idx = NULL;
        float* cell_met = NULL;
        float* cell_cov = NULL;        
        
        unsigned ncpgs = get_cell_data(cells[i], cell_idx, cell_met, cell_cov);
        if ((cell_idx != NULL) && (cell_met != NULL) && (cell_cov != NULL)) {
            
            // scatter cell coverage to all_cov
            vsUnpackV(ncpgs, cell_cov, &all_cov[0], cell_idx);        

            // gather cell coverage 
            vsPackV(idxs.length(), &all_cov[0], &all_idxs[0], &all_cell_cov[0]);
            
            // scatter cell methylation to all_meth
            vsUnpackV(ncpgs, cell_met, &all_meth[0], cell_idx);
            
            // gather cell methylation 
            vsPackV(idxs.length(), &all_meth[0], &all_idxs[0], &all_cell_meth[0]);
            
            for (unsigned j=0; j < all_cell_cov.size(); ++j){
                if (all_cell_cov[j] > 0){                
                    cov_mat[i].push_back(all_cell_cov[j]);
                    idxs_mat[i].push_back(idxs[j]);
                    meth_mat[i].push_back(all_cell_meth[j]);
                }
            }
            
            cell_mat[i].assign(cov_mat[i].size(), cells[i]);
            
            // // clean all_meth vector
            cblas_sscal(all_meth.size(), 0, &all_meth[0], 1);  
            cblas_sscal(all_cell_meth.size(), 0, &all_cell_meth[0], 1);    
            cblas_sscal(all_cov.size(), 0, &all_cov[0], 1);
            cblas_sscal(all_cell_cov.size(), 0, &all_cell_cov[0], 1);     
            
        }
        progress.report(1);
    }
    progress.report_last();
    
    Function c("c");
    Function do_call("do.call");
    
    NumericVector idxs_vector = do_call(c, wrap(idxs_mat));        
    NumericVector cov_vector = do_call(c, wrap(cov_mat));    
    NumericVector meth_vector = do_call(c, wrap(meth_mat));        
    StringVector cell_vector = do_call(c, wrap(cell_mat));    
    
    return DataFrame::create(_["cell_id"] = cell_vector, _["id"]=idxs_vector, _["cov"]=cov_vector, _["meth"]=meth_vector);
}

////////////////////////////////////////////////////////////////////////
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
        	
        	// scatter cell methylation to all_meth
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


NumericMatrix CGDB::count_pairs_all(const IntegerVector& idxs, const std::vector<std::string>& cells1, const std::vector<std::string>& cells2){

    std::vector<std::vector<int> > counts(cells1.size(), std::vector<int>(4, 0));

    std::vector<float> all_cov(m_CPG_NUM+1, 0); 
    std::vector<float> all_meth(m_CPG_NUM+1, 0);
    std::vector<int> all_idxs = as<std::vector<int> >(idxs);   

    std::vector<float> all_cell_cov(idxs.length(), 0);
    std::vector<float> all_cell_meth(idxs.length(), 0);  

    std::vector<float> cov1_vec;
    std::vector<float> meth1_vec;
    std::vector<int> idxs1_vec;  

    ProgressReporter progress;
    progress.init(cells1.size(), 1);
    std::string cell_id;
    for (unsigned i=0; i < cells1.size(); i++){
        if (i == 0 || cells1[i].compare(cell_id) != 0){
            cov1_vec.clear();
            meth1_vec.clear();
            idxs1_vec.clear();

            // get cell1 data
            int* cell1_idx = NULL;
            float* cell1_met = NULL;
            float* cell1_cov = NULL;
            unsigned ncpgs1 = get_cell_data(cells1[i], cell1_idx, cell1_met, cell1_cov);

            if ((cell1_idx != NULL) && (cell1_met != NULL) && (cell1_cov != NULL)) {

                cov1_vec.reserve(ncpgs1);
                meth1_vec.reserve(ncpgs1);
                idxs1_vec.reserve(ncpgs1);

                // scatter cell1 coverage to all_cov
                vsUnpackV(ncpgs1, cell1_cov, &all_cov[0], cell1_idx);        

                // gather cell1 coverage 
                vsPackV(idxs.length(), &all_cov[0], &all_idxs[0], &all_cell_cov[0]);
                
                // scatter cell1 methylation to all_meth
                vsUnpackV(ncpgs1, cell1_met, &all_meth[0], cell1_idx);
                
                // gather cell methylation 
                vsPackV(idxs.length(), &all_meth[0], &all_idxs[0], &all_cell_meth[0]);
                
                for (unsigned j=0; j < all_cell_cov.size(); ++j){
                    if (all_cell_cov[j] > 0){                                
                        idxs1_vec.push_back(idxs[j]);
                        meth1_vec.push_back(all_cell_meth[j]);
                    }
                }
            }
        }
          
        // clean all_meth vector
        cblas_sscal(all_meth.size(), 0, &all_meth[0], 1);  
        cblas_sscal(all_cell_meth.size(), 0, &all_cell_meth[0], 1);    
        cblas_sscal(all_cov.size(), 0, &all_cov[0], 1);
        cblas_sscal(all_cell_cov.size(), 0, &all_cell_cov[0], 1);    

        int* cell2_idx = NULL;
        float* cell2_met = NULL;
        float* cell2_cov = NULL;
        unsigned ncpgs2 = get_cell_data(cells2[i], cell2_idx, cell2_met, cell2_cov);  

        // scatter cell2 coverage to all_cov
        vsUnpackV(ncpgs2, cell2_cov, &all_cov[0], cell2_idx);      

        // scatter cell2 methylation to all_meth
        vsUnpackV(ncpgs2, cell2_met, &all_meth[0], cell2_idx);  

        std::vector<float> cov2_vec(meth1_vec.size(), 0);
        std::vector<float> meth2_vec(meth1_vec.size(), 0);

        // gather cell2 coverage according to cell1 indexes
        vsPackV(idxs1_vec.size(), &all_cov[0], &idxs1_vec[0], &cov2_vec[0]);

        // gather cell2 methylation according to cell1 indexes
        vsPackV(idxs1_vec.size(), &all_meth[0], &idxs1_vec[0], &meth2_vec[0]);

        // count pairs
        for (unsigned j=0; j < idxs1_vec.size(); ++j){
            if (cov2_vec[j] > 0){
                if (meth1_vec[j] == 0 && meth2_vec[j] == 0){
                    counts[i][0]++;
                } else if (meth1_vec[j] == 0 && meth2_vec[j] == 1){
                    counts[i][1]++;
                } else if (meth1_vec[j] == 1 && meth2_vec[j] == 0){
                    counts[i][2]++;
                } else if (meth1_vec[j] == 1 && meth2_vec[j] == 1){
                    counts[i][3]++;
                }
            }
        }
        progress.report(1);
    }      
    progress.report_last();

    Function rbind("rbind");
    Function do_call("do.call");
    NumericMatrix counts_mat = do_call(rbind, wrap(counts));
    return(counts_mat);
}

std::vector<int> CGDB::count_pairs(const IntegerVector& idxs, const std::string& cell1, const std::string& cell2, std::vector<int>& counts){

    // get cell1 data
    int* cell1_idx = NULL;
    float* cell1_met = NULL;
    float* cell1_cov = NULL;
    unsigned ncpgs1 = get_cell_data(cell1, cell1_idx, cell1_met, cell1_cov);

    // get cell2 data
    int* cell2_idx = NULL;
    float* cell2_met = NULL;
    float* cell2_cov = NULL;
    unsigned ncpgs2 = get_cell_data(cell2, cell2_idx, cell2_met, cell2_cov);

    std::vector<float> all_cov(m_CPG_NUM+1, 0); 
    std::vector<float> all_meth(m_CPG_NUM+1, 0);
    std::vector<int> all_idxs = as<std::vector<int> >(idxs);   

    std::vector<float> all_cell_cov(idxs.length(), 0);
    std::vector<float> all_cell_meth(idxs.length(), 0);    

    std::vector<float> cov1_vec;
    std::vector<float> meth1_vec;
    std::vector<int> idxs1_vec;

    if ((cell1_idx != NULL) && (cell1_met != NULL) && (cell1_cov != NULL) && (cell2_idx != NULL) && (cell2_met != NULL) && (cell2_cov != NULL)) {

        // scatter cell1 coverage to all_cov
        vsUnpackV(ncpgs1, cell1_cov, &all_cov[0], cell1_idx);        

        // gather cell1 coverage 
        vsPackV(idxs.length(), &all_cov[0], &all_idxs[0], &all_cell_cov[0]);
        
        // scatter cell1 methylation to all_meth
        vsUnpackV(ncpgs1, cell1_met, &all_meth[0], cell1_idx);
        
        // gather cell methylation 
        vsPackV(idxs.length(), &all_meth[0], &all_idxs[0], &all_cell_meth[0]);
        
        for (unsigned j=0; j < all_cell_cov.size(); ++j){
            if (all_cell_cov[j] > 0){                                
                idxs1_vec.push_back(idxs[j]);
                meth1_vec.push_back(all_cell_meth[j]);
            }
        }
        
        // clean all_meth vector
        cblas_sscal(all_meth.size(), 0, &all_meth[0], 1);  
        cblas_sscal(all_cell_meth.size(), 0, &all_cell_meth[0], 1);    
        cblas_sscal(all_cov.size(), 0, &all_cov[0], 1);
        cblas_sscal(all_cell_cov.size(), 0, &all_cell_cov[0], 1);  

        // scatter cell2 coverage to all_cov
        vsUnpackV(ncpgs2, cell2_cov, &all_cov[0], cell2_idx);      

        // scatter cell2 methylation to all_meth
        vsUnpackV(ncpgs2, cell2_met, &all_meth[0], cell2_idx);  

        std::vector<float> cov2_vec(meth1_vec.size(), 0);
        std::vector<float> meth2_vec(meth1_vec.size(), 0);

        // gather cell2 coverage according to cell1 indexes
        vsPackV(idxs1_vec.size(), &all_cov[0], &idxs1_vec[0], &cov2_vec[0]);

        // gather cell2 methylation according to cell1 indexes
        vsPackV(idxs1_vec.size(), &all_meth[0], &idxs1_vec[0], &meth2_vec[0]);

        // count pairs
        for (unsigned i=0; i < idxs1_vec.size(); ++i){
            if (cov2_vec[i] > 0){
                if (meth1_vec[i] == 0 && meth2_vec[i] == 0){
                    counts[0]++;
                } else if (meth1_vec[i] == 0 && meth2_vec[i] == 1){
                    counts[1]++;
                } else if (meth1_vec[i] == 1 && meth2_vec[i] == 0){
                    counts[2]++;
                } else if (meth1_vec[i] == 1 && meth2_vec[i] == 1){
                    counts[3]++;
                }
            }
        }      
    }

    return(counts);
}

////////////////////////////////////////////////////////////////////////
std::vector<std::string> CGDB::list_open_cells(){
	std::vector<std::string> cells;
	cells.reserve(m_cell_idx.size());	

	for(auto kv : m_cell_idx) {
		if (kv.second != NULL){
			cells.push_back(kv.first);    		
		}    	
	} 
	return cells;
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
        map = NULL;
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
