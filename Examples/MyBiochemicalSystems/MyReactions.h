#ifndef __H_MYREACTIONS__
#define __H_MYREACTIONS__

#include "StoichPackIKineticReaction.h"

using namespace StoichPack;
using namespace std;

/*Michaelis-Menten reaction:
	educt --> product (irreversible)
 rate law description:
	r(educt,product) = vMax*educt/(Km+educt)*/

class MichaelisMenten : public IKineticReaction {
private:
	MichaelisMenten();
	const sp_scalar vMax, Km;
public:
	MichaelisMenten(string educt, string product, sp_scalar vmax, sp_scalar km);

	sp_scalar RateImpl(const vector<sp_scalar>& c) const;
	vector<sp_scalar> DiffRateImpl(const vector<sp_scalar>& c) const;
	vector<bool> Dependencies() const;

	string Name() const;
};

class BiBiPingPong : public IKineticReaction {
private:
	BiBiPingPong();
	const sp_scalar v1, v2, coeffA, coeffB, coeffP, coeffQ, coeffAB, coeffAP, coeffBQ, coeffPQ;
	vector<sp_scalar> diff_denominator(const vector<sp_scalar>& c) const;
	
public:
	BiBiPingPong(string A, string B, string P, string Q,
	             sp_scalar vfmax, sp_scalar vbmax, sp_scalar keq,
	             sp_scalar kma, sp_scalar kmb, sp_scalar kmp, sp_scalar kmq,
	             sp_scalar kia, sp_scalar kiq);

	sp_scalar RateImpl(const vector<sp_scalar>& c) const;
	vector<sp_scalar> DiffRateImpl(const vector<sp_scalar>& c) const;

	vector<bool> Dependencies() const;

	string Name() const;
};

class MassActionLaw : public IKineticReaction {
private:
	MassActionLaw();
	sp_scalar vf, vb;
public:
	MassActionLaw(const vector<string>& names, const vector<sp_scalar>& coefficients, sp_scalar ratef, sp_scalar rateb);
	sp_scalar RateImpl(const vector<sp_scalar>& c) const ;
	vector<sp_scalar> DiffRateImpl(const vector<sp_scalar>& c) const ;
	vector<bool> Dependencies() const;
	string Name() const;
};

#endif
