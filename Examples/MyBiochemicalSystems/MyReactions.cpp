#include "MyReactions.h"

/* File:    MyReactions.cpp
 * Author:  Tobias Elbinger (elbinger@math.fau.de)
 * Purpose: Define some basic reactions, see MyReactions.h for more information about the reactions
 */


MichaelisMenten::MichaelisMenten(string educt, string product, sp_scalar vmax, sp_scalar km) : vMax(vmax), Km(km){
	IKineticReaction::AddSpecies(educt,-1);	//educts: always negative stoichiometric coefficient
	IKineticReaction::AddSpecies(product,1); //products: always positive stoichiometric coefficients
}

sp_scalar MichaelisMenten::RateImpl(const vector<sp_scalar>& c) const {
	return vMax*c[0]/(Km+c[0]);
}

vector<sp_scalar> MichaelisMenten::DiffRateImpl(const vector<sp_scalar>& c) const {
	return {vMax*Km/((Km+c[0])*(Km+c[0])),0.};
}

vector<bool> MichaelisMenten::Dependencies() const{
	// the rate function depends on the educt but not on the product
	return { true , false };
}

string MichaelisMenten::Name() const{
	return "Michaelis-Menten (irrev.)";
}

vector<sp_scalar> BiBiPingPong::diff_denominator(const vector<sp_scalar>& c) const{
	vector<sp_scalar> ret(4);
	ret[0]=coeffA+coeffAB*c[1]+coeffAP*c[2];
	ret[1]=coeffB+coeffAB*c[0]+coeffBQ*c[3];
	ret[2]=coeffP+coeffAP*c[0]+coeffPQ*c[3];
	ret[3]=coeffQ+coeffBQ*c[1]+coeffPQ*c[2];
	return ret;
}

BiBiPingPong::BiBiPingPong(string A, string B, string P, string Q,
             sp_scalar vfmax, sp_scalar vbmax, sp_scalar keq,
             sp_scalar kma, sp_scalar kmb, sp_scalar kmp, sp_scalar kmq,
             sp_scalar kia, sp_scalar kiq)
             	: v1(vfmax*vbmax), v2(vfmax*vbmax/keq), coeffA(vbmax*kmb), coeffB(vbmax*kma), coeffP(vfmax*kmq/keq), coeffQ(vfmax*kmp/keq),
	          coeffAB(vbmax), coeffAP(vfmax*kmq/(keq*kia)), coeffBQ(vfmax*kma/kiq), coeffPQ(vfmax/keq) {
	assert(keq!=0 && kia!=0 && kiq!=0);
	IKineticReaction::AddSpecies({A, B, P, Q},{-1, -1, 1, 1});
}
sp_scalar BiBiPingPong::RateImpl(const vector<sp_scalar>& c) const {
	const sp_scalar nom = v1*c[0]*c[1] - v2*c[2]*c[3];
	if(nom==0) return 0;
	const sp_scalar denom = coeffA*c[0]+coeffB*c[1]+coeffP*c[2]*coeffQ*c[3]+
	                        coeffAB*c[0]*c[1]+coeffAP*c[0]*c[2]+coeffBQ*c[1]*c[3]+coeffPQ*c[2]*c[3];
	assert(denom!=0);
	return nom/denom;
}

vector<sp_scalar> BiBiPingPong::DiffRateImpl(const vector<sp_scalar>& c) const {
	const sp_scalar nom = v1*c[0]*c[1] - v2*c[2]*c[3];
	const sp_scalar denom = coeffA*c[0]+coeffB*c[1]+coeffP*c[2]*coeffQ*c[3]+
	                        coeffAB*c[0]*c[1]+coeffAP*c[0]*c[2]+coeffBQ*c[1]*c[3]+coeffPQ*c[2]*c[3];

	if(denom==0){
		assert(c[0]==0 && c[1]==0 && c[2]==0 && c[3]==0);
		return vector<sp_scalar>(4,0);
	}

	vector<sp_scalar> ret(4);
	const vector<sp_scalar> diffdenom = diff_denominator(c);
	ret[0]=(denom * v1 *c[1] - nom*diffdenom[0])/(denom*denom);
	ret[1]=(denom * v1 *c[0] - nom*diffdenom[1])/(denom*denom);
	ret[2]=(denom * -v2 *c[3] - nom*diffdenom[2])/(denom*denom);
	ret[3]=(denom * -v2* c[2] - nom*diffdenom[3])/(denom*denom);

	return ret;
}

vector<bool> BiBiPingPong::Dependencies() const{
	return vector<bool>(4,true);
}

string BiBiPingPong::Name() const{
	return "Bi-Bi-Ping-Pong";
}

MassActionLaw::MassActionLaw(const vector<string>& names, const vector<sp_scalar>& coefficients, sp_scalar ratef, sp_scalar rateb)
                             : vf(ratef), vb(rateb) {
	if(vf*vb>0) throw StoichPackException("reaction rates have equal sign");

	vector<sp_scalar>::const_iterator it = find(coefficients.begin(),coefficients.end(),0);
	if(it!=coefficients.end()) throw StoichPackException("Cannot handle species with coefficient 0");

	class predicate{
	public:
		bool operator()(sp_scalar x) { return x>0; }
	};
	predicate p;
	const size_t npos=count_if(coefficients.begin(),coefficients.end(),p);
	if(npos==0 || npos==coefficients.size()) throw StoichPackException("Both sides must have species!");

	IKineticReaction::AddSpecies(names,coefficients);
}

sp_scalar MassActionLaw::RateImpl(const vector<sp_scalar>& c) const {
	sp_scalar forward=0, backward=0;
	const vector<sp_scalar>& coeff = this->Coefficients();
	const size_t s=c.size();
	for(size_t i=0;i<s;++i){
		if(coeff[i]>0) forward*=pow(c[i],coeff[i]);
		else backward*=pow(c[i],abs(coeff[i]));
	}
	return vf*forward-vb*backward;
}
vector<sp_scalar> MassActionLaw::DiffRateImpl(const vector<sp_scalar>& c) const {
	vector<sp_scalar> result(c.size());
	const vector<sp_scalar>& coeff = this->Coefficients();
	const size_t s=c.size();
	for(size_t i=0;i<s;++i) {
		if(coeff[i]>0) result[i]=vf*coeff[i]*pow(c[i],coeff[i]-1);
		else result[i]=vb*coeff[i]*pow(c[i],abs(coeff[i])-1);
	}
	return result;
}

vector<bool> MassActionLaw::Dependencies() const {
	vector<bool> result;
	for(auto c : this->Coefficients()){
		if(c>0 && vf==0) result.push_back(false);
		else if(c<0 && vb==0) result.push_back(false);
		else result.push_back(true);
	}
	return result;
}
string MassActionLaw::Name() const {
	if(vf==0 || vb==0) return "Mass Action Law (irrev.)";
	else return "Mass Action Law (rev.)";
}

