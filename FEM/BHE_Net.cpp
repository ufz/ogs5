#include <iostream>
#include "BHE_Net.h"

using namespace BHE; 

BHE_Net::BHE_Net()
{
    n_unknowns = 0; 
}

void BHE_Net::add_bhe_net_elem(BHE_Net_ELE_Abstract* element)
{
    bhe_map::iterator itr = _bhe_net.begin();
    std::string name = element->get_ele_name();
    itr = _bhe_net.find(name);
    if (itr == _bhe_net.end())
    {
        _bhe_net[name] = element; 
        return;
    }
    std::cout << "BHE net element already exists!\n";

    _global_start_idx = 0;
}

void BHE_Net::add_bhe_net_pipe(BHE_Net_ELE_Pipe* pipe,
                               std::string & from,
                               int from_ele_which_port,
                               std::string & to,
                               int to_ele_which_port)
{
    bhe_map::iterator itr      = _bhe_net.begin();
    bhe_map::iterator itr_from = _bhe_net.begin();
    bhe_map::iterator itr_to   = _bhe_net.begin();
    std::string name = pipe->get_ele_name();
    itr = _bhe_net.find(name);
    itr_from = _bhe_net.find(from); 
    itr_to   = _bhe_net.find(to);

    if (itr == _bhe_net.end())
    {
        // not having it in the map yet
        _bhe_net[name] = pipe;

        // check inlet link
        if (itr_from == _bhe_net.end())
        {
            // not existing
            std::cout << "BHE net pipeline inlet link does not exist!\n";
            exit(1);
        }
        else
        {
            pipe->add_inlet_connet(itr_from->second);
            pipe->add_inlet_connet_port(from_ele_which_port);
            itr_from->second->add_outlet_connet(pipe);
        }

        // check outlet link
        if (itr_to == _bhe_net.end())
        {
            // not existing
            std::cout << "BHE net pipeline outlet link does not exist!\n";
            exit(1);
        }
        else
        {
            pipe->add_outlet_connet(itr_to->second);
            pipe->add_outlet_connet_port(to_ele_which_port);
            itr_to->second->add_inlet_connet(pipe);
        }

        return;
    }
    else
    {
        std::cout << "BHE net pipeline already exists!\n";
        exit(1);
    }
}

void BHE_Net::count_n_unknowns()
{
    n_unknowns = 0;
    // loop over all elements in the map 
    typedef bhe_map::iterator it_type;
    for (it_type iterator = _bhe_net.begin(); iterator != _bhe_net.end(); iterator++) {
        // not counting the BHE, not counting the pipe
        if (iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_BOREHOLE       || 
            iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE           || 
            iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE_INNER_1U  ||
            iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE_INNER_2U  || 
            iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE_INNER_CXC ||
            iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE_INNER_CXA )
            continue; 
        else
        {
            n_unknowns += iterator->second->get_n_T_in();
            n_unknowns += iterator->second->get_n_T_out();
        }
    }
}

int BHE_Net::get_n_unknowns()
{
    // first count
    count_n_unknowns();

    // then return
    return n_unknowns; 
}

int BHE_Net::get_n_elems()
{
    // return the number of elements in the network
    return _bhe_net.size(); 
}

void BHE_Net::set_network_elem_idx(long n_nodes, long n_dofs_BHE)
{
    int i;
    long idx = n_nodes + n_dofs_BHE;
    int local_idx = 0; 

    _global_start_idx = idx;

    // loop over all elements in the map 
    typedef bhe_map::iterator it_type;
    for (it_type iterator = _bhe_net.begin(); iterator != _bhe_net.end(); iterator++) {
        // BHE
        if ( iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_BOREHOLE )
        {
            // do nothing on the global side, 
            // since it has already been assigned with global index before

            // onlye the local index needs to be added
            for (i = 0; i < iterator->second->get_n_T_in(); i++)
            {
                iterator->second->set_T_in_local_index(local_idx, i);
                local_idx++;
            }
            for (i = 0; i < iterator->second->get_n_T_out(); i++)
            {
                iterator->second->set_T_out_local_index(local_idx, i);
                local_idx++;
            }

        }
        // now DISTRIBUTOR or heat pump
        else if (iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_DISTRIBUTOR || iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_HEATPUMP)
        {
			long pipe_global_index;
			int pipe_local_index;
			int connected_port;

			// assgin index to T_in
            for (i = 0; i < iterator->second->get_n_T_in(); i++)
            {
                iterator->second->set_T_in_local_index(local_idx, i);
                local_idx++; 
                iterator->second->set_T_in_global_index(idx, i);
                idx++;
            }

            // assgin index to T_out
            for (i = 0; i < iterator->second->get_n_T_out(); i++)
            {
                iterator->second->set_T_out_local_index(local_idx, i);
                local_idx++;
                iterator->second->set_T_out_global_index(idx, i);
                idx++;
            }
        }
        else
        {
            continue; 
        }
    }

    // second loop, only deal with the pipelines
    for (it_type iterator = _bhe_net.begin(); iterator != _bhe_net.end(); iterator++) 
    {
        if ( iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE )
        {
            long pipe_global_index;
            int pipe_local_index;
            int connected_port;

            // the pipeline T_in T_out index is obtained from BHE, heat pump and distributors

            // assgin index to T_in
            connected_port = iterator->second->get_inlet_connet_port();
            pipe_global_index = iterator->second->get_inlet_connect()->get_T_out_global_index(connected_port);
            pipe_local_index = iterator->second->get_inlet_connect()->get_T_out_local_index(connected_port);
            iterator->second->set_T_in_global_index(pipe_global_index);
            iterator->second->set_T_in_local_index(pipe_local_index);

            // assgin index to T_out
            connected_port = iterator->second->get_outlet_connet_port();
            pipe_global_index = iterator->second->get_outlet_connect()->get_T_in_global_index(connected_port);
            pipe_local_index = iterator->second->get_outlet_connect()->get_T_in_local_index(connected_port);
            iterator->second->set_T_out_global_index(pipe_global_index);
            iterator->second->set_T_out_local_index(pipe_local_index);
        }
        else if ( iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE_INNER_1U )
        {
            // do not need the local index anymore

            long pipe_global_index;
            int connected_port;

            // the pipeline T_in T_out index is obtained from BHE, heat pump and distributors

            // assgin index to T_in
            connected_port = iterator->second->get_inlet_connet_port();
            // notice the inlet of this pipe is from the bottom of the BHE
            pipe_global_index = iterator->second->get_inlet_connect()->get_T_in_bottom_global_index(connected_port);
            iterator->second->set_T_in_global_index(pipe_global_index);

            // assgin index to T_out
            connected_port = iterator->second->get_outlet_connet_port();
            // notice the outlet of this pipe is from the bottom of the BHE
            pipe_global_index = iterator->second->get_outlet_connect()->get_T_out_bottom_global_index(connected_port);
            iterator->second->set_T_out_global_index(pipe_global_index);
            
        }
        else if (iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE_INNER_2U)
        {
            // TODO
        }
        else if (iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE_INNER_CXC)
        {
            // TODO
        }
        else if (iterator->second->get_net_ele_type() == BHE_NET_ELE::BHE_NET_PIPE_INNER_CXA)
        {
            // TODO
        }
        else
        {
            continue;
        }
    }

    
}