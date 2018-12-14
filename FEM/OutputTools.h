/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OUTPUTTOOLS_H_
#define OUTPUTTOOLS_H_

#include <iostream>

#include "rf_mmp_new.h"
#include "fem_ele_std.h"

using namespace std;

struct ELEMENT_MMP_VALUES
{
    static double getValue(CMediumProperties* mmp, int mmp_id, long i_e,
                           double* gp, double theta)
    {
        double mat_value = .0;
        switch (mmp_id)
        {
            case 0:
                mat_value = mmp->Porosity(i_e, theta);
                break;
            case 1:
                mat_value = mmp->PermeabilityTensor(i_e)[0];
                break;
            case 2:
                mat_value = mmp->StorageFunction(i_e, gp, theta);
                break;
            default:
                cout
                    << "ELEMENT_MMP_VALUES::getValue(): no MMP values specified"
                    << endl;
                break;
        }
        return mat_value;
    }

    static int getMMPIndex(const std::string& mmp_name)
    {
        int mmp_id = -1;
        if (mmp_name.compare("POROSITY") == 0)
        {
            mmp_id = 0;
        }
        else if (mmp_name.compare("PERMEABILITY") == 0)
        {
            mmp_id = 1;
        }
        else if (mmp_name.compare("STORAGE") == 0)
        {
            mmp_id = 2;
        }
        else
        {
            cout << "ELEMENT_MMP_VALUES::getMMPIndex(): no valid MMP values "
                    "specified. "
                 << mmp_name << endl;
        }
        return mmp_id;
    }
};

struct ELEMENT_MFP_VALUES
{
    static double getValue(CFluidProperties* mfp, int mfp_id)
    {
        double mat_value = .0;
        switch (mfp_id)
        {
            case 0:
                mat_value = mfp->Density();
                break;
            case 1:
                mat_value = mfp->Viscosity();
                break;
            default:
                cout << "ELEMENT_MFP_VALUES: no MFP values specified" << endl;
                break;
        }
        return mat_value;
    }

    static int getMFPIndex(const std::string& mfp_name)
    {
        int mfp_id = -1;
        if (mfp_name.compare("DENSITY") == 0)
        {
            mfp_id = 0;
        }
        else if (mfp_name.compare("VISCOSITY") == 0)
        {
            mfp_id = 1;
        }
        else
        {
            cout << "ELEMENT_MFP_VALUES: no valid MFP values specified. "
                 << mfp_name << endl;
        }
        return mfp_id;
    }
};

inline double getElementMMP(int mmp_id, MeshLib::CElem* ele, CRFProcess* m_pcs)
{
    double gp[3] = {.0, .0, .0};
    double theta = 1.0;
    ele->SetOrder(false);
    CFiniteElementStd* fem = m_pcs->GetAssember();
    fem->ConfigElement(ele, false);
    fem->Config();
    fem->getShapeFunctionCentroid();
    CMediumProperties* mmp = mmp_vector[ele->GetPatchIndex()];
    double val =
        ELEMENT_MMP_VALUES::getValue(mmp, mmp_id, ele->GetIndex(), gp, theta);
    return val;
}

inline double getNodeMMP(int mmp_id, MeshLib::CFEMesh* m_msh,
                         MeshLib::CNode* node, CRFProcess* m_pcs)
{
    const std::vector<size_t>& connected_ele_ids =
        node->getConnectedElementIDs();
    double ele_avg = .0;
    for (long i_e = 0; i_e < (long)connected_ele_ids.size(); i_e++)
    {
        MeshLib::CElem* ele = m_msh->ele_vector[connected_ele_ids[i_e]];
        ele_avg += getElementMMP(mmp_id, ele, m_pcs);
    }
    ele_avg /= connected_ele_ids.size();
    return ele_avg;
}

inline double getNodeElementValue(int ele_value_id, MeshLib::CFEMesh* m_msh,
                                  MeshLib::CNode* node, CRFProcess* m_pcs)
{
    const std::vector<size_t>& connected_ele_ids =
        node->getConnectedElementIDs();
    double ele_avg = .0;
    for (long i_e = 0; i_e < (long)connected_ele_ids.size(); i_e++)
    {
        MeshLib::CElem* ele = m_msh->ele_vector[connected_ele_ids[i_e]];
        ele_avg += m_pcs->GetElementValue(ele->GetIndex(), ele_value_id);
    }
    ele_avg /= connected_ele_ids.size();
    return ele_avg;
}

#endif  // OUTPUTTOOLS_H_
