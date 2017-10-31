#ifndef CGDB_H
#define CGDB_H

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

class CGDB {
    private: 
        const std::string m_db_dir;
        const int m_CPG_NUM;

        std::unordered_map<std::string, int*> m_cell_idx;        
        std::unordered_map<std::string, float*> m_cell_cov;
        std::unordered_map<std::string, float*> m_cell_met;
        std::unordered_map<std::string, int> m_ncpgs;

        void add_cell_data(std::string const& cell);

        unsigned get_cell_data(std::string const& cell, int*& cell_idx, float*& cell_met, float*& cell_cov);

        bool valid_indexes(const IntegerVector& idxs){
        	return (min(idxs) > 0 && max(idxs) <= m_CPG_NUM );
        }

    public:
        CGDB(const std::string& db_dir, const int& CPG_NUM): m_db_dir(db_dir), m_CPG_NUM(CPG_NUM){
        }

        void freemem(){
        	// unmap all files
            for (auto& n : m_cell_idx){                
                unmap_file(n.second, m_ncpgs[n.first]);                
            }
            m_cell_idx.clear();

            for (auto& n : m_cell_cov){                
                unmap_file(n.second, m_ncpgs[n.first]);                
            }
            m_cell_cov.clear();

            for (auto& n : m_cell_met){                
                unmap_file(n.second, m_ncpgs[n.first]);                
            }     
            m_cell_met.clear();
        }

        ~CGDB(){
       		freemem();
        }

        DataFrame mean_meth(const IntegerVector& idxs, const std::vector<std::string>& cells);
        DataFrame bin_meth(const IntegerVector& idxs, const IntegerVector& bins, const std::vector<std::string>& cells);        
        List bin_meth_per_cell_cpp(const IntegerVector& idxs, const IntegerVector& bins, const std::vector<std::string>& cells);
        List extract(const IntegerVector& idxs, const std::vector<std::string>& cells);
        List extract_sparse(const IntegerVector& idxs, const std::vector<std::string>& cells, const std::vector<std::string>& chrom, const std::vector<int>& start, const std::vector<int>& end);

        std::vector<std::string> list_open_cells();
};


#endif //CGDB_H
