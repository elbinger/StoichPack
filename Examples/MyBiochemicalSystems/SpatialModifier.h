/* File:    SpatialModifier.h
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Create simple initial solutions that are not constant in space.
 */

#ifndef __H_SPATIAL_MODIFER__
#define __H_SPATIAL_MODIFER__

#include "StoichPackDefines.h"
#include <vector>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>

class SpatialModifier{
public:
	std::string name;
	StoichPack::sp_scalar value;
	std::vector<bool> mod;
	SpatialModifier(const std::string& n, StoichPack::sp_scalar v, const std::vector<bool>& m)
			     : name(n), value(v), mod(m) { assert(mod.size()==3 && value>=0); }
	bool operator==(const std::string& n) const { return name==n; }
	StoichPack::sp_scalar GetValue(const std::vector<StoichPack::sp_scalar>& X) const {
		StoichPack::sp_scalar result = 1;
		assert(X.size()<=3);
		for(size_t i=0;i<X.size();++i){
			if(mod[i]) result*=cos(M_PI*(1-X[i]));
			else result*=cos(M_PI*X[i]);
		}
		return 0.5*value*(1.+result);
	}
};

#endif
