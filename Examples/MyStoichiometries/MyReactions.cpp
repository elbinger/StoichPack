#include "MyReactions.h"

MichaelisMenten::MichaelisMenten(string educt, string product, sp_scalar vmax, sp_scalar km) : vMax(vmax), Km(km){
	IKineticReaction::AddSpecies(educt,-1);
	IKineticReaction::AddSpecies(product,1);
}

sp_scalar MichaelisMenten::RateImpl(const vector<sp_scalar>& c) const {
	return vMax*c[0]/(Km+c[0]);
}

std::vector<sp_scalar> MichaelisMenten::DiffRateImpl(const vector<sp_scalar>& c) const {
	vector<sp_scalar> ret(2,0);
	ret[0]=vMax*Km/((Km+c[0])*(Km+c[0]));
	return ret;
}

std::vector<sp_scalar> MichaelisMenten::Dependencies() const{
	vector<sp_scalar> ret(2,0);
	ret[0]=1;
	return ret;
}

std::string MichaelisMenten::Name() const{
	return string("Michaelis-Menten kinetic: ")+this->Participants()[0].Name()+string(" --> ")+this->Participants()[1].Name();
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
	IKineticReaction::AddSpecies(A,-1);
	IKineticReaction::AddSpecies(B,-1);
	IKineticReaction::AddSpecies(P,-1);
	IKineticReaction::AddSpecies(Q,-1);
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

vector<sp_scalar> BiBiPingPong::Dependencies() const{
	return vector<sp_scalar>(4,1);
}

std::string BiBiPingPong::Name() const{
	return string("BiBiPingPong kinetic: ")+this->Participants()[0].Name()+string(" + ")+this->Participants()[1].Name()+string(" <--> ")
	              +this->Participants()[2].Name()+string(" + ")+this->Participants()[3].Name();
}

