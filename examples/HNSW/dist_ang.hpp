#ifndef __DIST_ANG_HPP__
#define __DIST_ANG_HPP__

#include "type_point.hpp"

class descr_fvec
{
public:
	typedef fvec type_point;
	static float distance(const type_point &u, const type_point &v, uint32_t dim)
	{
		const auto &uc=u.coord, &vc=v.coord;
		float dot=0, nu=0, nv=0;
		for(uint32_t i=0; i<dim; ++i)
		{
			nu += float(uc[i])*uc[i];
			nv += float(vc[i])*vc[i];
			dot += float(uc[i])*vc[i];
		}
		return 1-dot/(sqrt(nu)*sqrt(nv));
	//	return acos(dot/(sqrt(nu)*sqrt(nv)));
	}

	static auto get_id(const type_point &u)
	{
		return u.id;
	}
};

#endif // _DIST_ANG_HPP_
