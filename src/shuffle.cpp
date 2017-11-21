// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include <sstream>
#include <random>
#include "ProgressReporter.h"
using namespace Rcpp;

typedef std::pair<int, int> Point;

namespace std
{
    template<>
    struct hash<Point>
    {
        std::size_t operator () (const Point& p) const
        {
            return ((std::hash<int>()(p.first) >> 1 ) ^ (std::hash<int>()(p.second) >> 1));
        }
    };
}

// [[Rcpp::export]]
NumericMatrix shuffle_mat_marginals(const NumericMatrix& m_meth, const NumericMatrix& inds_mat, const int& n_shuff = 1000){	
	
	std::unordered_map< Point, int > m;  
	for (int i = 0; i < inds_mat.nrow(); ++i){
		int x = inds_mat(i, 0) - 1;
		int y = inds_mat(i, 1) - 1;

		Point p(x, y);
	
		m[p] = m_meth(x, y); 		
	}

	NumericVector r = {0, 0};
	int x1, x2, y1, y2;	
	ProgressReporter progress;
    progress.init(n_shuff, 1);
	for (int i = 0; i < n_shuff; i++){
		r = round(runif(2) * (m.size() - 1), 0);

		auto random_it1 = std::next(std::begin(m), r[0]);
		x1 = random_it1->first.first;
		y2 = random_it1->first.second;

		auto random_it2 = std::next(std::begin(m), r[1]);
		x2 = random_it2->first.first;
		y1 = random_it2->first.second;

		Point x1y1(x1, y1);
		Point x1y2(x1, y2);
		Point x2y2(x2, y2);
		Point x2y1(x2, y1);

		m[x1y1]++;
		m[x1y2]--;
		
		m[x2y2]++;
		m[x2y1]--;

		if (m.at(x1y2) == 0){
			m.erase(x1y2);
		}

		if (m.at(x2y1) == 0){
			m.erase(x2y1);
		}
		progress.report(1);
	}
	progress.report_last();
	
	NumericMatrix res(m_meth.nrow(), m_meth.ncol());		
	for (auto& i : m){	
		int x = i.first.first;
		int y = i.first.second;
		res(x, y) = i.second;		
	}
	
	return(res);
	
}

