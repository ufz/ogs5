/**
* \file BHE_Net.h
* 2015/08/28 HS inital implementation
* borehole heat exchanger network class
*
*
*/

#ifndef BHE_NET_H
#define BHE_NET_H

#include <map>
#include "BHE_Net_ELE_Abstract.h"
#include "BHE_Net_ELE_Pipe.h"
#include "BHE_Net_ELE_Pipe_Inner_1U.h"
#include "BHE_Net_ELE_HeatPump.h"

namespace BHE  // namespace of borehole heat exchanger
{
    typedef std::map<std::string, BHE_Net_ELE_Abstract*> bhe_map;

    class BHE_Net {
    public:
        /**
          * constructor
          */
        BHE_Net(); 

        void add_bhe_net_elem(BHE_Net_ELE_Abstract* element);

        void add_bhe_net_pipe(BHE_Net_ELE_Pipe* pipe,
                              std::string & from,
                              int from_ele_which_port,
                              std::string & to,
                              int to_ele_which_port);

        /**
          * get the number of unknowns
          */
        int get_n_unknowns(); 

        /**
          * get the number of elements in the network
          */
        int get_n_elems(); 

        /**
          * set the global and local indices for all elements in the network
          */
        void set_network_elem_idx(long n_nodes, long n_dofs_BHE);

        bhe_map get_network()
        {
            return _bhe_net;
        }

        long get_global_start_idx()
        {
            return _global_start_idx; 
        }

    private:

        void count_n_unknowns();

        /**
          * a map including all bhes, distributors, and pipelines
          */
        bhe_map _bhe_net;

        /**
          * number of unknown temperatures in the network
          */
        int n_unknowns; 

        /**
          * the starting index in the global linear equation system
          */
        long _global_start_idx; 
    };
}

#endif