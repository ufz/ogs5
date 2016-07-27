/**
* \file BHE_Net_ELE_Pipe_Inner_1U.h
* 2016/04/14 HS inital implementation
* borehole heat exchanger network element pipeline class
* particularly designed for 1U type of BHE
* used to connect the bottom of BHE, linking the inlet and 
* outlet pipelines
*/

#ifndef BHE_NET_ELE_PIPE_INNER_1U_H
#define BHE_NET_ELE_PIPE_INNER_1U_H

#include "BHE_Net_ELE_Abstract.h"
#include "BHEAbstract.h"
#include "BHE_Net_ELE_Pipe.h"

namespace BHE  // namespace of borehole heat exchanger
{
    class BHE_Net_ELE_Pipe_Inner_1U : public BHE_Net_ELE_Pipe {

    public:
        /**
        * constructor
        */
        BHE_Net_ELE_Pipe_Inner_1U(std::string & name, BHE::BHEAbstract * m_BHE);

    protected:
        /**
          * obtain the global index at the pipeline inlet
          */
        std::size_t get_global_idx_in(); 

        /**
          * obtain the global index at the pipeline outlet
          */
        std::size_t get_global_idx_out();

    private:
        /**
          * the global index at the pipeline inlet
          */
        std::size_t _global_idx_in; 

        /**
          * the global index at the pipeline outlet
          */
        std::size_t _global_idx_out;

        /**
          * the BHE which this pipeline is applied on
          */
        const BHE::BHEAbstract * _m_BHE; 
    };

}

#endif