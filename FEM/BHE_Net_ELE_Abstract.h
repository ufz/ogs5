/**
* \file BHE_Net_ELE_Abstract.h
* 2015/08/27 HS inital implementation
* borehole heat exchanger network element abstract class
*
*
*/

#ifndef BHE_NET_ELE_ABSTRACT_H
#define BHE_NET_ELE_ABSTRACT_H

#include "Eigen/Eigen"
#include <vector>

using namespace std;

namespace BHE  // namespace of borehole heat exchanger
{
    /**
    * list the types of BHE net element
    */
    namespace BHE_NET_ELE{
      enum type {
        BHE_NET_PIPE,          // pipeline, used to connect different BHEs
        BHE_NET_PIPE_INNER_1U, // pipeline at the bottom of a 1U BHE
        BHE_NET_PIPE_INNER_2U, // pipeline at the bottom of a 2U BHE
        BHE_NET_PIPE_INNER_CXC,// pipeline at the bottom of a CXC BHE
        BHE_NET_PIPE_INNER_CXA,// pipeline at the bottom of a CXA BHE
        BHE_NET_DISTRIBUTOR,   // distributor
        BHE_NET_HEATPUMP,      // heat pump
        BHE_NET_BOREHOLE	   // borehole
      };
    }

    class BHE_Net_ELE_Abstract {
    
    public:
        /**
          * constructor
          */
      BHE_Net_ELE_Abstract(std::string name, BHE_NET_ELE::type type, int n_inlet = 1, int n_outlet = 1)
            : N_IN(n_inlet), N_OUT(n_outlet), _name(name), _ele_type(type)
        {
            int i; 
            // initialize T_in
            T_in = new double[N_IN];
            global_idx_T_in = new long[N_IN];
            global_idx_T_in_bottom = new long[N_IN];
            local_idx_T_in = new int[N_IN];
            for (i = 0; i < N_IN; i++)
            {
                T_in[i] = 0.0;
                global_idx_T_in[i] = -1; // uninitialized index value
                global_idx_T_in_bottom[i] = -1; // uninitialized index value
                local_idx_T_in[i] = -1;
            }

            // initialize T_out
            T_out = new double[N_OUT];
            global_idx_T_out = new long[N_OUT];
            global_idx_T_out_bottom = new long[N_OUT];
            local_idx_T_out = new int[N_OUT];
            for (i = 0; i < N_OUT; i++)
            {
                T_out[i] = 0.0;
                global_idx_T_out[i] = -1; // uninitialized index value
                global_idx_T_out_bottom[i] = -1; // uninitialized index value
                local_idx_T_out[i] = -1;
            }

            // initialize penalty factor
            _penalty_factor = 0.0;

        }

        /**
          * destructor
          */
        ~BHE_Net_ELE_Abstract(){
            delete [] T_in; 
            delete [] T_out; 
        }

        /**
          * return the net element type
          */
        BHE_NET_ELE::type get_net_ele_type()
        {
            return _ele_type; 
        }

        std::string get_ele_name()
        {
            return _name; 
        }

        /**
          * get inlet temperature
          */
        double get_T_in(int idx = 0) {
            return T_in[idx];
        }

        /**
          * get the global index of T_in
          */
        long get_T_in_global_index(int idx_T_in = 0) {
            return global_idx_T_in[idx_T_in];
        }

        /**
          * set the global index of T_in
          */
        void set_T_in_global_index(long new_idx, int idx_T_in = 0) {
            global_idx_T_in[idx_T_in] = new_idx;
        }

        /**
        * set the local index of T_in
        */
        void set_T_in_local_index(long new_idx, int idx_T_in = 0) {
            local_idx_T_in[idx_T_in] = new_idx;
        }

        /**
        * get the local index of T_out
        */
        long get_T_in_local_index(int idx_T_in = 0) {
            return local_idx_T_in[idx_T_in];
        }

        /**
        * set the global index of T_in at the bottom of BHE
        */
        void set_T_in_bottom_global_index(long new_idx, int idx_T_in = 0) {
            global_idx_T_in_bottom[idx_T_in] = new_idx;
        }

        /**
        * get the global index of T_in at the bottom of BHE
        */
        long get_T_in_bottom_global_index(int idx_T_in = 0) {
            return global_idx_T_in_bottom[idx_T_in];
        }

        /**
          * get outlet temperature
          */
        double get_T_out(int idx = 0) {
            return T_out[idx];
        }

        /**
          * set the global index of T_out
          */
        void set_T_out_global_index(long new_idx, int idx_T_out = 0) {
            global_idx_T_out[idx_T_out] = new_idx;
        }

        /**
          * get the global index of T_out
          */
        long get_T_out_global_index(int idx_T_out = 0) {
            return global_idx_T_out[idx_T_out];
        }

        /**
        * set the global index of T_out at the bottom of BHE
        */
        void set_T_out_bottom_global_index(long new_idx, int idx_T_out = 0) {
            global_idx_T_out_bottom[idx_T_out] = new_idx;
        }

        /**
        * get the global index of T_out at the bottom of BHE
        */
        long get_T_out_bottom_global_index(int idx_T_out = 0) {
            return global_idx_T_out_bottom[idx_T_out];
        }

        /**
        * set the local index of T_out
        */
        void set_T_out_local_index(long new_idx, int idx_T_out = 0) {
            local_idx_T_out[idx_T_out] = new_idx;
        }

        /**
        * get the local index of T_out
        */
        int get_T_out_local_index(int idx_T_out = 0) {
            return local_idx_T_out[idx_T_out];
        }

        /**
          * set inlet temperature
          */
        void set_T_in(double val, int idx = 0){
            T_in[idx] = val;
        }

        /**
          * set outlet temperature
          */
        void set_T_out(double val, int idx = 0){
            T_out[idx] = val; 
        }

        /**
          * return the number of inlet temperatures
          */
        int get_n_T_in()
        {
            return N_IN; 
        }

        /**
          * return the number of outlet temperatures
          */
        int get_n_T_out()
        {
            return N_OUT;
        }

        void add_inlet_connet(BHE_Net_ELE_Abstract* connect)
        {
            _vec_ele_inlet.push_back(connect);
        }

        void add_inlet_connet_port(int port)
        {
            _vec_ele_inlet_port.push_back(port);
        }

        int get_inlet_connet_port(int idx = 0)
        {
            return _vec_ele_inlet_port[idx];
        }

        void add_outlet_connet(BHE_Net_ELE_Abstract* connect)
        {
            _vec_ele_outlet.push_back(connect);
        }

        void add_outlet_connet_port(int port)
        {
            _vec_ele_outlet_port.push_back(port);
        }

        int get_outlet_connet_port(int idx=0)
        {
            return _vec_ele_outlet_port[idx];
        }

        BHE_Net_ELE_Abstract* get_inlet_connect(int idx = 0)
        {
            return _vec_ele_inlet[idx];
        }

        BHE_Net_ELE_Abstract* get_outlet_connect(int idx = 0)
        {
            return _vec_ele_outlet[idx];
        }



        double get_inlet_ratio(int idx = 0)
        {
            return _vec_inlet_ratio(idx);
        }

        double get_outlet_ratio(int idx = 0)
        {
            return _vec_outlet_ratio(idx);
        }

        void set_penalty_factor(double val)
        {
            _penalty_factor = val; 
        }

        double get_penalty_factor()
        {
            return _penalty_factor;
        }

		// this one is needed for the setting the heat pump B.C.
		virtual double set_BC(double T_in, double current_time) = 0;

		virtual double get_flowrate() = 0;

        /**
          * return the RHS value, needs to be implemented.
          */
        // virtual double get_RHS_value() = 0;

    protected:
        /**
        * how the inlet flow rate is determined.
        */
        Eigen::VectorXd _vec_inlet_ratio;

        /**
        * how the outlet flow rate is determined.
        */
        Eigen::VectorXd _vec_outlet_ratio;

    private:

        /**
          * array of inlet temperature 
          */
        double * T_in; 

        /**
          * array of outlet temperature
          */
        double * T_out; 

        /**
          * array of inlet temperature global index
          */
        long * global_idx_T_in;

        /**
          * array of outlet temperature global index
          */
        long * global_idx_T_out;

        /**
        * array of inlet temperature global index at the bottom of the BHE
        */
        long * global_idx_T_in_bottom;

        /**
        * array of outlet temperature global index at the bottom of the BHE
        */
        long * global_idx_T_out_bottom;

        /**
          * array of penalty factor
          */
        double _penalty_factor; 

        /**
          * array of inlet temperature local index
          */
        int * local_idx_T_in;

        /**
          * array of outlet temperature local index
          */
        int * local_idx_T_out;

        const int N_IN; 

        const int N_OUT;

        const BHE_NET_ELE::type _ele_type;

        std::string _name; 

        std::vector<BHE_Net_ELE_Abstract*> _vec_ele_inlet; 

        std::vector<BHE_Net_ELE_Abstract*> _vec_ele_outlet;

        std::vector<int> _vec_ele_inlet_port; 

        std::vector<int> _vec_ele_outlet_port; 
       
    };





}

#endif
