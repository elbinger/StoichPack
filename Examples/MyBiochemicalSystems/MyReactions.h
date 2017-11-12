/* File:    MyReactions.h
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Define some basic reactions
 */

#ifndef __H_MYREACTIONS__
#define __H_MYREACTIONS__

#include "StoichPackIKineticReaction.h"

using namespace StoichPack;
using namespace std;

/*Michaelis-Menten reaction:
	educt --> product (irreversible)
 rate law description:
	r(educt,product) = vMax*educt/(Km+educt)*/

class MichaelisMenten : public IKineticReaction { //IKineticReaction is the default interface for reactions in StoichPack
private:
	MichaelisMenten();
	const sp_scalar vMax, Km;
public:
	MichaelisMenten(string educt, string product, sp_scalar vmax, sp_scalar km);

	//functions inherited from interface
	sp_scalar RateImpl(const vector<sp_scalar>& c) const;
	vector<sp_scalar> DiffRateImpl(const vector<sp_scalar>& c) const;
	vector<bool> Dependencies() const;

	string Name() const;
};

/*Bi-Bi-Ping-Pong reaction:
	A + B <--> P + Q (reversible)
 rate law description:
	r(A,B,P,Q) = nom(A,B,P,Q)/denom(A,B,P,Q), where
	nom(A,B,P,Q) = vfmax*vbmax*(A*B-P*Q/keq) and
	denom(A,B,P,Q) = vbmax*kmb*A+vbmax*kma*B+vfmax*kmq*P/keq,vfmax*kmp*Q/keq
			+vbmax*A*B+vfmax*kmq*A*P/(keq*kia)+vfmax*kma*B*Q/kiq+vfmax*P*Q/keq
	all constants >= 0*/

class BiBiPingPong : public IKineticReaction {
private:
	BiBiPingPong();
	const sp_scalar v1, v2, coeffA, coeffB, coeffP, coeffQ, coeffAB, coeffAP, coeffBQ, coeffPQ;
	//v1= vfmax*vbmax, v2 = v1/keq, coeffA = vmax*kmb, coeffB = vbmax*kma, coeffP = vfmax*kmq/keq, coeffQ = vfmax*kmp/keq
	//coeffAB = vbmax, coeffAP = vfmax*kmq/(keq*kia), coeffBQ = vfmax*kma/kiq, coeffPQ = vfmax/keq
	// --> nom(A,B,P,Q) = v1*A*B-v2*P*Q
	// denom(A,B,P,Q) = coeffA*A+coeffB*B+coeffP*P+coeffQ*Q+coeffAB*A*B+coeffAP*A*P+coeffBQ*B*Q+coeffPQ*P*Q

	//claculate the derivative of denom
	vector<sp_scalar> diff_denominator(const vector<sp_scalar>& c) const;
	
public:
	BiBiPingPong(string A, string B, string P, string Q,
	             sp_scalar vfmax, sp_scalar vbmax, sp_scalar keq,
	             sp_scalar kma, sp_scalar kmb, sp_scalar kmp, sp_scalar kmq,
	             sp_scalar kia, sp_scalar kiq);

	//functions inherited from interface
	sp_scalar RateImpl(const vector<sp_scalar>& c) const;
	vector<sp_scalar> DiffRateImpl(const vector<sp_scalar>& c) const;

	vector<bool> Dependencies() const;

	string Name() const;
};

/*Reaction according to mass action law:
	a_1*X_1+a_2*X_2+ ... +a_n*X_n <--> b_1Y_1+b2Y_2+ ... + b_m*Y_m (reversible)
 rate law description:
	r(X,Y) = ratef*X_1^a1*X_2^a_2*...*X_n^a_n - rateb*Y_1^b_1*Y_2^b_2*...*Y_m^b_m
	nom(A,B,P,Q) = vfmax*vbmax*(A*B-P*Q/keq) and
	denom(A,B,P,Q) = vbmax*kmb*A+vbmax*kma*B+vfmax*kmq*P/keq,vfmax*kmp*Q/keq
			+vbmax*A*B+vfmax*kmq*A*P/(keq*kia)+vfmax*kma*B*Q/kiq+vfmax*P*Q/keq
	all constants >= 0*/

class MassActionLaw : public IKineticReaction {
private:
	MassActionLaw();
	sp_scalar vf, vb;
public:
	//Remark: the a_k are stored as negabtive values
	MassActionLaw(const vector<string>& names, const vector<sp_scalar>& coefficients, sp_scalar ratef, sp_scalar rateb);
	//Functions inherited from interface
	sp_scalar RateImpl(const vector<sp_scalar>& c) const ;
	vector<sp_scalar> DiffRateImpl(const vector<sp_scalar>& c) const ;
	vector<bool> Dependencies() const;
	string Name() const;
};

#endif
