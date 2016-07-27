#include "BHE_Net_ELE_Pipe_Inner_1U.h"

using namespace BHE;

BHE_Net_ELE_Pipe_Inner_1U::BHE_Net_ELE_Pipe_Inner_1U(std::string & name, BHE::BHEAbstract * m_BHE)
    : BHE_Net_ELE_Pipe(name, BHE_NET_ELE::BHE_NET_PIPE_INNER_1U),
    _m_BHE(m_BHE)
{

    // configure the penalty factor
    this->set_penalty_factor(1.0e6);

}

std::size_t BHE::BHE_Net_ELE_Pipe_Inner_1U::get_global_idx_in()
{
    return _global_idx_in;
}

std::size_t BHE::BHE_Net_ELE_Pipe_Inner_1U::get_global_idx_out()
{
    return _global_idx_out;
}

