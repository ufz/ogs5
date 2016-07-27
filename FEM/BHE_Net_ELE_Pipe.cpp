#include "BHE_Net_ELE_Pipe.h"

using namespace BHE;

BHE_Net_ELE_Pipe::BHE_Net_ELE_Pipe(std::string & name, BHE_NET_ELE::type type )
    : BHE_Net_ELE_Abstract(name, type, 1, 1)
{
	this->set_penalty_factor(1.0e6);
}

double BHE_Net_ELE_Pipe::get_RHS_value()
{
    return 0; 
}