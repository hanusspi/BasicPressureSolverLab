#pragma once
#include <functional>
#include <iostream>
#include <algorithm>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch.h>

namespace pressureSolver
{
	// DECLARATIONS ================================================================================
	class NeighborhoodSearch
	{
		static constexpr bool USE_CNS = true;

	public:
		/* Fields */
		CompactNSearch::NeighborhoodSearch cns = CompactNSearch::NeighborhoodSearch(0.1);
		double timing_cns = 0.0;
		double timing_tns = 0.0;
		int n_sets = 0;
		
		/* Methods */
		NeighborhoodSearch() = default;
		~NeighborhoodSearch() = default;

		void set_search_radius(const double r);
		int add_point_set(const double* data, const int n);
		int add_static_point_set(const double* data, const int n);
		void set_active_search(const int i, const int j);
		void run();
		void resize_point_set(const int set_i, const double* data, const int n);

		template<typename FUNC>
		void for_each_neighbor(const int set_i, const int set_j, const int i, FUNC f);
	};


	// DEFINITIONS ================================================================================
	inline void NeighborhoodSearch::set_search_radius(const double r)
	{
		if constexpr (USE_CNS) {
			this->cns = CompactNSearch::NeighborhoodSearch(r);
		}
	}
	inline int NeighborhoodSearch::add_point_set(const double* data, const int n)
	{
		this->n_sets++;
		int set_id = -1;
		if constexpr (USE_CNS) {
			set_id = (int)this->cns.add_point_set(data, n, true, false, false);
		}
		return set_id;
	}
	inline int NeighborhoodSearch::add_static_point_set(const double* data, const int n)
	{
		this->n_sets++;
		int set_id = -1;
		if constexpr (USE_CNS) {
			set_id = (int)this->cns.add_point_set(data, n, false);
		}
		return set_id;
	}
	inline void NeighborhoodSearch::set_active_search(const int i, const int j)
	{
		if constexpr (USE_CNS) {
			this->cns.set_active((unsigned int)i, (unsigned int)j, true);
		}
	}
	inline void NeighborhoodSearch::run()
	{
		if constexpr (USE_CNS) {
			const double t0 = omp_get_wtime();
			this->cns.find_neighbors();
			const double t1 = omp_get_wtime();
			this->timing_cns += t1 - t0;
		}
	}
	inline void NeighborhoodSearch::resize_point_set(const int set_i, const double* data, const int n)
	{
		if constexpr (USE_CNS) {
			this->cns.resize_point_set(set_i, data, n);
		}
	}
	template<typename FUNC>
	inline void NeighborhoodSearch::for_each_neighbor(const int set_i, const int set_j, const int i, FUNC f)
	{
		if (set_j >= this->n_sets) {
			return;
		}

		if constexpr (USE_CNS) {
			const auto& neighbors = this->cns.point_set((unsigned int)set_i);
			const int n_neighbors = (int)neighbors.n_neighbors((unsigned int)set_j, i);
			for (int loc_j = 0; loc_j < n_neighbors; loc_j++) {
				const int j = neighbors.neighbor((unsigned int)set_j, i, loc_j);
				f(j);
			}
		}
	}
}