#include "BHE_Net_ELE_Distributor.h"

using namespace BHE;

BHE_Net_ELE_Distributor::BHE_Net_ELE_Distributor(std::string & name, Eigen::VectorXd & vec_Inlet_Ratio, Eigen::VectorXd & vec_Outlet_Ratio)
: BHE_Net_ELE_Abstract(name, BHE_NET_ELE::BHE_NET_DISTRIBUTOR, vec_Inlet_Ratio.size(), vec_Outlet_Ratio.size())
{
    _vec_inlet_ratio = vec_Inlet_Ratio;
    _vec_outlet_ratio = vec_Outlet_Ratio;

}

double BHE_Net_ELE_Distributor::get_RHS_value()
{
    return 0; 
}