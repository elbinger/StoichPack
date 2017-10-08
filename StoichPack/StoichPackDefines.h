/* File: StoichPackDefines.h
 * Purpose: Provide basic definitions for StoichPack.
 * Author: Tobias Elbinger <elbinger@math.fau.de> */

#ifndef __H_STOICHPACK_DEFINES__
#define __H_STOICHPACK_DEFINES__

#include <exception>
#include <string>

namespace StoichPack{

 //type for scalars, change if you want to. default = double.
 typedef double sp_scalar;

 // exception class
 class StoichPackException : public std::exception {
	private:
		std::string description; // description of the error
		StoichPackException(); // forbid empty construction
	public:
		// create StoichPackException: desc: description of exception
		StoichPackException(const std::string& desc) : description("StoichPackException: "+desc) {}

		const char* what() const throw() { return description.c_str(); } // inherited

		virtual ~StoichPackException() throw() {}
 };

}

#endif
