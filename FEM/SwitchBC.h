/*
 * Switch.h
 *
 *  Created on: Oct 21, 2015
 *      Author: waltherm
 */

#ifndef FEM_SWITCHBC_H_
#define FEM_SWITCHBC_H_


#include "FEMEnums.h"

struct SwitchBC {
	double switchValue;
	double switchOnValue;
	double switchOffValue;
	FiniteElement::ProcessType switchProcessType;
	FiniteElement::PrimaryVariable switchPrimVar;


	SwitchBC () :
		switchValue(0.0),
		switchOnValue(1.0),
		switchOffValue(0.0),
		switchProcessType(FiniteElement::INVALID_PROCESS),
		switchPrimVar(FiniteElement::INVALID_PV)
	{}

};



#endif /* FEM_SWITCHBC_H_ */
