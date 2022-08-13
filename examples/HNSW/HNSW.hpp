#ifndef _HNSW_HPP
#define _HNSW_HPP

#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <random>
#include <memory>
#include <atomic>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <set>
#include <iterator>
#include <type_traits>
#include <limits>
#include <thread>
// #include "parallelize.h"
#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/delayed_sequence.h>
#include <parlay/random.h>

#include <cpam/cpam.h>
#include <pam/pam.h>
#include "../graphs/aspen/aspen.h"
#define DEBUG_OUTPUT 0
#if DEBUG_OUTPUT
#define debug_output(...) fprintf(stderr, __VA_ARGS__)
#else
#define debug_output(...) do{[](...){}(__VA_ARGS__);}while(0)
#endif // DEBUG_OUTPUT

namespace ANN{

enum class type_metric{
	L2, ANGULAR, DOT
};

template<typename U, template<typename> class Allocator=std::allocator>
class HNSW
{
	using T = typename U::type_point;
	struct empty_weight{};
	using Graph = aspen::symmetric_graph<empty_weight>;
	using edge_tree = typename Graph::edge_tree;
	using edge_node = typename Graph::edge_node;
	using vertex_tree = typename Graph::vertex_tree;
	using vertex_node = typename Graph::vertex_node;
public:
	/*
		Construct from the vectors [begin, end).
		std::iterator_trait<Iter>::value_type ought to be convertible to T
		dim: 				vector dimension
		m_l: 				control the # of levels (larger m_l leads to more layer)
		m: 					max degree
		ef_construction:	beam size during the construction
		alpha:				parameter of the heuristic (similar to the one in vamana)
		batch_base: 		growth rate of the batch size (discarded because of two passes)
	*/
	template<typename Iter>
	HNSW(Iter begin, Iter end, uint32_t dim, float m_l=1, uint32_t m=100, uint32_t ef_construction=50, float alpha=5, float batch_base=2, bool do_fixing=false);

	// save the current model to a file
	void save(const std::string &filename_model) const;
public:
	typedef uint32_t node_id;

	
	struct node_orignal{
		// uint32_t id;
		uint32_t level;
		T data;
		std::vector<node_orignal*> *neighbors;
	};
	

	struct dist{
		float d;
		node_id u;
	};

	struct dist_ex : dist
	{
		uint32_t depth;
	};

	struct nearest{
		constexpr bool operator()(const dist &lhs, const dist &rhs) const{
			return lhs.d>rhs.d;
		}
	};

	struct farthest{
		constexpr bool operator()(const dist &lhs, const dist &rhs) const{
			return lhs.d<rhs.d;
		}
	};

	struct cmp_id{
		constexpr bool operator()(const dist &lhs, const dist &rhs) const{
			return lhs.u<rhs.u;
		}
	};

	std::vector<node_id> entrance; // To init
	// auto m, max_m0, m_L; // To init
	uint32_t dim;
	float m_l;
	uint32_t m;
	// uint32_t level_max = 30; // To init
	uint32_t ef_construction;
	float alpha;
	uint32_t n;
	// Allocator<node> allocator;
	// std::vector<node*> node_pool;
	std::vector<const T*> pointset;

	std::vector<Graph> G;
	// TODO: use uint8_t to save memory but
	// need to consider the compatibility of the saved models
	std::vector<uint32_t> level; 

	mutable std::atomic<size_t> total_visited = 0;
	mutable std::atomic<size_t> total_eval = 0;
	mutable std::atomic<size_t> total_size_C = 0;

	std::vector<node_id> neighbourhood(const node_id u, uint32_t level) const
	{
		std::vector<node_id> res;
		auto vtx = G[level].get_vertex(u);
		auto f = [&](node_id u, node_id v, empty_weight wgh){
			res.push_back(v);
			return true;
		};
		vtx.out_neighbors().foreach_cond(f);
		return res;
	}

	template<typename Iter>
	void insert(Iter begin, Iter end, bool from_blank);

	// To optimize
	auto select_neighbors_heuristic(const T *coord_target, 
		/*const std::priority_queue<dist,std::vector<dist>,farthest> &C*/
		const parlay::sequence<dist> &C, uint32_t M,
		uint32_t level, bool extendCandidate, bool keepPrunedConnections)
	{
		parlay::sequence<dist> W_d;
		std::set<dist,cmp_id> W_tmp;
		for(auto &e : C) // TODO: add const?
		{
			W_tmp.insert(e);
			if(extendCandidate)
			{
				for(node_id e_adj : neighbourhood(e.u,level))
				{
					// if(e_adj==nullptr) continue;
					if(W_tmp.find(dist{0,e_adj})==W_tmp.end())
						W_tmp.insert(dist{U::distance(pointset[e_adj],coord_target,dim),e_adj});
				}
			}
		}

		parlay::sequence<dist> W(W_tmp.begin(), W_tmp.end());
		std::sort(W.begin(), W.end(), farthest());
		W_tmp.clear();

		std::vector<node_id> R;
		for(const auto &e : W)
		{
			if(R.size()>=M) break;
			const auto d_q = e.d;

			bool is_good = true;
			for(const auto &r : R)
			{
				const auto d_r = U::distance(pointset[e.u], pointset[r], dim);
				//if(d_r*(level+1)<d_q*alpha*(entrance->level+1))
				if(d_r<d_q*alpha)
				{
					is_good = false;
					break;
				}
			}

			if(is_good)
				R.push_back(e.u);
			else
				W_d.push_back(e);
		}

		// elements in `W_d` are in order
		auto it = W_d.begin();
		if(keepPrunedConnections)
		{
			while(it!=W_d.end() && R.size()<M)
				R.push_back((it++)->u);
		}
		return R;
	}

	auto select_neighbors(const T *coord_target, 
		/*const std::priority_queue<dist,std::vector<dist>,farthest> &C,*/
		const parlay::sequence<dist> &C, uint32_t M,
		uint32_t level, bool extendCandidate=false, bool keepPrunedConnections=false)
	{
		/*
		(void)level, (void)extendCandidate, (void)keepPrunedConnections;
		return select_neighbors_simple(u,C,M);
		*/
		return select_neighbors_heuristic(coord_target, C, M, level, extendCandidate, keepPrunedConnections);
	}

	uint32_t get_level_random()
	{
		// static thread_local int32_t anchor;
		// uint32_t esp;
		// asm volatile("movl %0, %%esp":"=a"(esp));
		static thread_local std::hash<std::thread::id> h;
		static thread_local std::mt19937 gen{h(std::this_thread::get_id())};
		static thread_local std::uniform_real_distribution<> dis(std::numeric_limits<float>::min(), 1.0);
		const uint32_t res = uint32_t(-log(dis(gen))*m_l);
		return res;
	}

	auto search_layer(const node_id u, const std::vector<node_id> &eps, uint32_t ef, uint32_t l_c) const; // To static
	auto get_threshold_m(uint32_t level){
		return level==0? m*2: m;
	}
};

template<typename U, template<typename> class Allocator>
template<typename Iter>
HNSW<U,Allocator>::HNSW(Iter begin, Iter end, uint32_t dim_, float m_l_, uint32_t m_, uint32_t ef_construction_, float alpha_, float batch_base, bool do_fixing)
	: dim(dim_), m_l(m_l_), m(m_), ef_construction(ef_construction_), alpha(alpha_), n(std::distance(begin,end))
{
	static_assert(std::is_same_v<typename std::iterator_traits<Iter>::value_type, T>);
	static_assert(std::is_base_of_v<
		std::random_access_iterator_tag, typename std::iterator_traits<Iter>::iterator_category>);

	if(n==0) return;

	auto perm = parlay::random_permutation<uint32_t>(n, 1206);
	/*
	auto rand_seq = parlay::delayed_seq<const T&>(n, [&](uint32_t i){
		auto &t = *(begin+perm[i]);
		return (t);
	});
	*/
	auto rand_seq = parlay::make_slice(begin,end);

	level.resize(n); // TODO: decide the size

	const auto level_ep = get_level_random();
	// node *entrance_init = allocator.allocate(1);
	// new(entrance_init) node{level_ep, *rand_seq.begin(), new std::vector<node*>[level_ep+1]/*anything else*/};
	// node_pool.push_back(entrance_init);
	const auto id_ep = U::get_id(*rand_seq.begin());
	pointset.resize(id_ep+1);
	pointset[id_ep] = &*begin;

	level[id_ep] = level_ep;
	G.resize(level_ep+1);
	for(uint32_t l=0; l<=level_ep; ++l)
		G[l].insert_vertex_inplace(id_ep,nullptr);
	entrance.push_back(id_ep);

	uint32_t begin_batch=0, end_batch=1;
	const uint32_t batch_size_min = 1;
	const uint32_t batch_size_max = 20000;
	float progress = 0.0;
	while(end_batch<n)
	{
		begin_batch = end_batch;
		end_batch = std::min({
			n,
			(uint32_t)std::ceil(begin_batch*batch_base)+batch_size_min, 
			begin_batch+batch_size_max
		});

		insert(rand_seq.begin()+begin_batch, rand_seq.begin()+end_batch, true);
		// insert(rand_seq.begin()+begin_batch, rand_seq.begin()+end_batch, false);

		if(end_batch>n*(progress+0.1))
		{
			progress = float(end_batch)/n;
			fprintf(stderr, "Done: %.2f\n", progress);
			fprintf(stderr, "# visited: %lu\n", total_visited.load());
			fprintf(stderr, "# eval: %lu\n", total_eval.load());
			fprintf(stderr, "size of C: %lu\n", total_size_C.load());
		}
	}

	fprintf(stderr, "# visited: %lu\n", total_visited.load());
	fprintf(stderr, "# eval: %lu\n", total_eval.load());
	fprintf(stderr, "size of C: %lu\n", total_size_C.load());
	// if(do_fixing) symmetrize_edge();

	#if 0
		for(const auto *pu : node_pool)
		{
			fprintf(stderr, "[%u] (%.2f,%.2f)\n", U::get_id(pu->data), pu->data[0], pu->data[1]);
			for(int32_t l=pu->level; l>=0; --l)
			{
				fprintf(stderr, "\tlv. %d:", l);
				for(const auto *k : pu->neighbors[l])
					fprintf(stderr, " %u", U::get_id(k->data));
				fputs("\n", stderr);
			}
		}
	#endif
}

template<typename U, template<typename> class Allocator>
template<typename Iter>
void HNSW<U,Allocator>::insert(Iter begin, Iter end, bool from_blank)
{
	const auto level_ep = level[entrance[0]];
	const auto size_batch = std::distance(begin,end);
	auto node_new = std::make_unique<node_id[]>(size_batch+1);
	auto vertex_conn_forward = std::make_unique<std::tuple<node_id,edge_node*>[]>(size_batch);
	auto eps = std::make_unique<std::vector<node_id>[]>(size_batch);
	//const float factor_m = from_blank? 0.5: 1;
	const auto factor_m = 1;

	debug_output("Insert %lu elements; from blank? [%c]\n", size_batch, "NY"[from_blank]);

	// auto *pool = allocator.allocate(size_batch);
	// first, query the nearest point as the starting point for each node to insert
	parlay::parallel_for(0, size_batch, [&](uint32_t i){
		const node_id u = U::get_id(*(begin+i));
		node_new[i] = u;

		if(from_blank)
			level[u] = get_level_random();
	});

	// and add new nodes to the pool
	if(from_blank)
	{
		const node_id max_id = parlay::reduce(
			parlay::make_slice(node_new.get(),node_new.get()+size_batch),
			parlay::maxm<node_id>{}
		);
		if(max_id+1>pointset.size())
			pointset.resize(max_id+1);
		parlay::parallel_for(0, size_batch, [&](size_t i){
			// TODO: mind the dangling pointer in the real use case
			pointset[node_new[i]] = &*(begin+i);
		});
	}

	std::sort(node_new.get(), node_new.get()+size_batch, [&](node_id u, node_id v){
		return level[u] > level[v];
	}); // TODO: decide if using parallel sort
	auto index = parlay::delayed_seq<size_t>(size_batch, [&](size_t i){
		return i;
	});
	auto is_distinct = parlay::delayed_seq<size_t>(size_batch, [&](size_t i){
		return level[node_new[i]]!=level[node_new[i+1]];
	});
	node_new[size_batch] = node_new[size_batch-1];
	auto pos_split = parlay::pack(index, is_distinct);
	pos_split.push_back(size_batch-1);

	// insert new points into G
		/*
		
		*/

	// debug_output("Nodes are settled\n");

	parlay::parallel_for(0, size_batch, [&](size_t i){
		const node_id u = node_new[i];
		const auto level_u = level[u];
		auto &eps_u = eps[i]; 
		// eps_u.push_back(entrance);
		eps_u = entrance;
		for(uint32_t l=level_ep; l>level_u; --l)
		{
			const auto res = search_layer(u, eps_u, 1, l); // TODO: optimize
			eps_u.clear();
			eps_u.push_back(res[0].u);
		}
	});

	debug_output("Finish searching entrances\n");

	// then we process them layer by layer (from high to low)
	const uint32_t level_max = level[node_new[0]];
	if(level_max>level_ep)
		G.resize(level_max+1);
	parlay::parallel_for(0, size_batch, [&](size_t i){
		vertex_conn_forward[i] = {node_new[i], nullptr};
	});

	size_t j = 0;
	for(uint32_t l=level_max; l>level_ep; --l)
	{
		if(j<pos_split.size() && level[node_new[pos_split[j]]]==l)
			j++;
		assert(j>0);
		// inserted points in the range of [0, end_high) have their levels >= `l`
		const size_t end_high = pos_split[j-1]+1;
		G[l].insert_vertices_batch(end_high, vertex_conn_forward.get());
	}

	for(int32_t l=std::min(level_ep,level_max); l>=0; --l) // TODO: fix the type
	{
		debug_output("Finding neighbors on lev. %d\n", l);
		parlay::sequence<parlay::sequence<std::pair<node_id,node_id>>> edge_add(size_batch);

		if(j<pos_split.size() && level[node_new[pos_split[j]]]==l)
			j++;
		assert(j>0);
		const size_t end_high = pos_split[j-1]+1;

		parlay::parallel_for(0, end_high, [&](uint32_t i){
			const node_id u = node_new[i];

			auto &eps_u = eps[i];
			auto res = search_layer(u, eps_u, ef_construction, l);
			auto nbh_u = select_neighbors(pointset[u], res, get_threshold_m(l)*factor_m, l);
			// move the content from `nbh_u` to `u.neighbors[l]`
			auto &edge_u = edge_add[i];
			edge_u.clear();
			edge_u.reserve(nbh_u.size());

			parlay::sequence<std::tuple<node_id,empty_weight>> nbh_new;
			nbh_new.reserve(nbh_u.size());

			for(auto v : nbh_u)
			{
				edge_u.emplace_back(v, u);
				nbh_new.emplace_back(v, empty_weight{}); // TODO: if we can use `delay_seq` instead here
			}
			// insert new neighbors to `u`
			// nbh_new[i] = std::move(nbh_u);
			edge_tree edge_insert(nbh_new.begin(), nbh_new.end());
			vertex_conn_forward[i] = {u, edge_insert.root};
			edge_insert.root = nullptr;

			// prepare the entrance points for the next layer
			eps_u.clear();
			eps_u.reserve(res.size());
			for(const auto e : res)
				eps_u.push_back(e.u);
		});

		debug_output("Adding forward edges\n");
		G[l].insert_vertices_batch(end_high, vertex_conn_forward.get());

		debug_output("Adding reverse edges\n");
		// now we add edges in the other direction
		auto edge_add_flatten = parlay::flatten(edge_add);
		auto edge_add_grouped = parlay::group_by_key(edge_add_flatten);
		auto vertex_conn_backward = std::make_unique<std::tuple<node_id,edge_node*>[]>(edge_add_grouped.size());

		parlay::parallel_for(0, edge_add_grouped.size(), [&](size_t j){
			node_id v = edge_add_grouped[j].first;
			auto nbh_v = neighbourhood(v,l);
			auto &nbh_v_add = edge_add_grouped[j].second;

			for(auto it=nbh_v_add.begin(); it!=nbh_v_add.end();)
			{
				bool is_extant = *it==v||std::find_if(nbh_v.begin(), nbh_v.end(), [&](const node_id u_extant){
					return *it==u_extant;
				})!=nbh_v.end();
				it = is_extant? nbh_v_add.erase(it): std::next(it);
			}

			const uint32_t size_nbh_total = nbh_v.size()+nbh_v_add.size();

			const auto m_s = get_threshold_m(l)*factor_m;
			if(size_nbh_total>m_s) // TODO: optimize this part
			{
				auto dist_nbh = std::make_unique<dist[]>(size_nbh_total);
				for(size_t k=0; k<nbh_v.size(); ++k)
					dist_nbh[k] = dist{U::distance(pointset[nbh_v[k]],pointset[v],dim), nbh_v[k]};
				for(size_t k=0; k<nbh_v_add.size(); ++k)
					dist_nbh[k+nbh_v.size()] = dist{U::distance(pointset[nbh_v_add[k]],pointset[v],dim), nbh_v_add[k]};

				std::sort(dist_nbh.get(), dist_nbh.get()+size_nbh_total, farthest());

				nbh_v.resize(m_s);
				for(size_t k=0; k<m_s; ++k)
					nbh_v[k] = dist_nbh[k].u;
			}
			else nbh_v.insert(nbh_v.end(),nbh_v_add.begin(), nbh_v_add.end());

			parlay::sequence<std::tuple<node_id,empty_weight>> nbh_new;
			nbh_new.reserve(nbh_v.size());
			for(auto u : nbh_v)
				nbh_new.emplace_back(u, empty_weight{}); // TODO: if we can use `delay_seq` instead here

			// insert new neighbors to `v`
			edge_tree edge_insert(nbh_new.begin(), nbh_new.end());
			vertex_conn_backward[j] = {v, edge_insert.root};
			edge_insert.root = nullptr;
		});
		G[l].insert_vertices_batch(
			edge_add_grouped.size(),
			vertex_conn_backward.get()
		);
	} // for l

	debug_output("Updating entrance\n");
	// finally, update the entrances

	if(level_max>level_ep)
		entrance.clear();

	if(level_max>=level_ep)
		entrance.insert(entrance.end(), node_new.get(), node_new.get()+size_batch);
}

template<typename U, template<typename> class Allocator>
auto HNSW<U,Allocator>::search_layer(const node_id u, const std::vector<node_id> &eps, uint32_t ef, uint32_t l_c) const
{
	// std::vector<bool> visited(n);
	// TODO: Try hash to an array
	// TODO: monitor the size of `visited`
	std::set<uint32_t> visited;
	// std::priority_queue<dist,std::vector<dist>,nearest> C;
	// std::priority_queue<dist,std::vector<dist>,farthest> W;
	parlay::sequence<dist> /*C, W, */W_;
	std::set<dist,farthest> C, C_acc;
	uint32_t cnt_used = 0;

	for(node_id ep : eps)
	{
		// visited[U::get_id(ep->data)] = true;
		visited.insert(ep);
		const auto d = U::distance(pointset[u],pointset[ep],dim);
		C.insert({d,ep});
		C_acc.insert({d,ep});
		// C.push_back({d,ep,1});
		// W.push_back({d,ep,1});
	}
	// std::make_heap(C.begin(), C.end(), nearest());
	// std::make_heap(W.begin(), W.end(), farthest());

	while(C.size()>0)
	{
		// const auto &f = *(W[0].u);
		// if(U::distance(c.data,u.data,dim)>U::distance(f.data,u.data,dim))
		// if(C[0].d>W[0].d) break;
		if(C_acc.size()==cnt_used) break;
		// if(l_c==0) total_eval++;
		/*
		const auto &c = *(C[0].u);
		*/
		auto it = C.begin();
		const auto c = it->u;
		// W_.push_back(C[0]);
		W_.push_back(*it);
		// std::pop_heap(C.begin(), C.end(), nearest());
		// C.pop_back();
		C.erase(it);
		cnt_used++;
		for(node_id v: neighbourhood(c, l_c))
		{
			// if(visited[U::get_id(pv->data)]) continue;
			// visited[U::get_id(pv->data)] = true;
			if(!visited.insert(v).second) continue;
			// const auto &f = *(W[0].u);
			// if(W.size()<ef||U::distance(pv->data,u.data,dim)<U::distance(f.data,u.data,dim))
			const auto d = U::distance(pointset[u],pointset[v],dim);
			// if(W.size()<ef||d<W[0].d)
			// if(C.size()<ef||d<C.rend()->d)
			{
				// C.push_back({d,pv,dc+1});
				// std::push_heap(C.begin(), C.end(), nearest());
				/*
				W.push_back({d,pv,dc+1});
				std::push_heap(W.begin(), W.end(), farthest());
				if(W.size()>ef)
				{
					std::pop_heap(W.begin(), W.end(), farthest());
					W.pop_back();
				}
				*/
				if(C.size()<ef || d<C.rbegin()->d)
				{
				C.insert({d,v});
				if(C.size()>ef)
				{
					// std::pop_heap(C.begin(), C.end(), nearest());
					// C.pop_back();
					C.erase(std::prev(C.end()));
				}
				}
				if(C_acc.size()<ef || d<C_acc.rbegin()->d)
				{
				C_acc.insert({d,v});
				if(C_acc.size()>ef)
				{
					auto it = std::prev(C_acc.end());
					if(std::find_if(W_.begin(), W_.end(), [&](const dist &a){
						return a.u==it->u;
					})!=W_.end())
						cnt_used--;
					C_acc.erase(it);
				}
				}
			}
		}
	}
	if(l_c==0)
	{
		// total_visited += visited.size();
		// total_size_C += C.size();
	}
	/*
	std::sort(W.begin(), W.end(), farthest());
	if(W.size()>ef) W.resize(ef);
	*/
	return W_;
}

template<typename U, template<typename> class Allocator>
void HNSW<U,Allocator>::save(const std::string &filename_model) const
{
	std::ofstream model(filename_model, std::ios::binary);
	if(!model.is_open())
		throw std::runtime_error("Failed to create the model");

	const auto write = [&](const auto &data, auto ...args){
		auto write_impl = [&](auto &f, const auto &data, auto ...args){
			using T = std::remove_reference_t<decltype(data)>;
			if constexpr(std::is_pointer_v<std::decay_t<T>>)
			{
				auto write_array = [&](const auto &data, size_t size, auto ...args){
					for(size_t i=0; i<size; ++i)
						f(f, data[i], args...);
				};
				// use the array extent as the size
				if constexpr(sizeof...(args)==0 && std::is_array_v<T>)
				{
					write_array(data, std::extent_v<T>);
				}
				else
				{
					static_assert(sizeof...(args), "size was not provided");
					write_array(data, args...);
				}
			}
			else
			{
				static_assert(std::is_standard_layout_v<T>);
				model.write((const char*)&data, sizeof(data));
			}
		};
		write_impl(write_impl, data, args...);
	};
	// write header (version number, type info, etc)
	write("HNSW", 4);
	write(uint32_t(1));
	write(typeid(U).hash_code()^sizeof(U));
	// write(sizeof(node));
	write(sizeof(node_orignal));
	// write parameter configuration
	write(dim);
	write(m_l);
	write(m);
	write(ef_construction);
	write(alpha);
	write(n);
	// write indices
	for(const T *p : pointset)
	{
		const node_id u = U::get_id(*p);
		write(level[u]);
		write(u);
	}
	for(const T *p : pointset)
	{
		const node_id u = U::get_id(*p);
		for(uint32_t l=0; l<=level[u]; ++l)
		{
			const auto nbh_u = neighbourhood(u,l);
			write(nbh_u.size());
			for(const node_id v : nbh_u)
				write(v);
		}
	}
	// write entrances
	write(entrance.size());
	for(const node_id u : entrance)
		write(u);
}

} // namespace HNSW

#endif // _HNSW_HPP

