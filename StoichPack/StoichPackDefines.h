#ifndef __H_STOICHPACK_DEFINES__
#define __H_STOICHPACK_DEFINES__

#include <exception>
#include <string>

#define USE1
//#define USE2

namespace StoichPack{

 //type for scalars, change if you want to. default = double.
 typedef double sp_scalar;

 class StoichPackException : public std::exception {
	private:
		std::string description;
		StoichPackException(); // FORBID
	public:
		// create StoichPackException: desc: description of exception
		StoichPackException(const std::string& desc) : description("StoichPackException: "+desc) {}

		const char* what() const throw() { return description.c_str(); } // inherited
 };

}

#endif
