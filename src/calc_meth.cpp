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
void unmap_file(void*& map, unsigned length);

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
        std::string const& cell = cells[i];

        unsigned point = cell.find(".");
        std::string fname = db_dir + "/" +  cell.substr(0, point) + "/" + cell.substr(point+1);

        std::string idx_fname = fname + ".idx.bin";
        std::string met_fname = fname + ".meth.bin";
        std::string cov_fname = fname + ".cov.bin";
        unsigned ncpgs = get_file_length(idx_fname) / 4;

        int*   cell_idx = (int*)mmap_file(idx_fname, ncpgs*4);
        float* cell_met = (float*)mmap_file(met_fname, ncpgs*4);
        float* cell_cov = (float*)mmap_file(cov_fname, ncpgs*4);

        if ((cell_idx == NULL) || (cell_met == NULL)) {
            meth[i] = NA_REAL;
            cov[i] = NA_REAL;
        }
        else {
            meth[i] = cblas_sdoti(ncpgs, cell_met, cell_idx, &mask[0]);            
            cov[i] = cblas_sdoti(ncpgs, cell_cov, cell_idx, &mask[0]);
        }

        munmap(cell_met, ncpgs*4);
        munmap(cell_idx, ncpgs*4);
    }    

    return DataFrame::create(_["cell"]=cells, _["cov"]=cov, _["meth"]=meth);
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
void unmap_file(void*& map, unsigned length){
    if (map != NULL) {
        munmap(map, length);
        map = NULL;
    }
}
