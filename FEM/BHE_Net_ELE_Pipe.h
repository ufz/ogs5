/**
* \file BHE_Net_ELE_Pipe.h
* 2015/08/28 HS inital implementation
* borehole heat exchanger network element pipeline class
*
*
*/

#ifndef BHE_NET_ELE_PIPE_H
#define BHE_NET_ELE_PIPE_H

#include "BHE_Net_ELE_Abstract.h"

namespace BHE  // namespace of borehole heat exchanger
{
    class BHE_Net_ELE_Pipe : public BHE_Net_ELE_Abstract {

    public:
        /**
          * constructor
          */
        BHE_Net_ELE_Pipe(std::string & name, BHE_NET_ELE::type type = BHE_NET_ELE::BHE_NET_PIPE);

        double get_RHS_value();

		double set_BC(double T_in, double current_time)
		{
			return 0;
		}

		double get_flowrate()
		{
			return 0;
		}

    };

}

#endif