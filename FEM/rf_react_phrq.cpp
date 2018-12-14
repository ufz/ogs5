/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
rf_react_phrq.cpp phrq
Reaction package to go with PHREEQC
Dedong Li
Kiel, 01/2011
*/

//#include "stdafx.h" /* MFC */
#include <signal.h>
#include "display.h"
#include "makros.h"
#include "rf_pcs.h"
//#include "nodes.h"
#include "rf_pcs.h"
#include "mathlib.h"
#include "rfmat_cp.h"
#include "rf_react_cap.h"
#include "rf_react_phrq.h"
#include "rf_react.h"
#include "rf_react_int.h"
#include "rf_ic_new.h"
#include "stdio.h"
//#include "geo_strings.h"
#include "rf_pcs.h"  //OK_MOD"
//#include "files.h"
#include "files0.h"
//#include "nodes.h"
//#include "elements.h"           /* für ElGetElementGroupNumber */
#include <iostream>
#include <iomanip>
#include <fstream>

#include <sstream>
#include <string>

#include <vector>
#include "rf_tim_new.h"
#include "rf_mmp_new.h"
#include "rf_kinreact.h"
// Elem object
#include "fem_ele_std.h"
#ifdef OGS_FEM_CAP   // CAP_REACT
#include "cacint.h"  //DL
#endif
#include "gs_project.h"  // CB

using namespace std;

#ifdef OGS_FEM_CAP  // CAP_REACT

vector<vector<std::string> > PHREEQC_TEMPLATE;
vector<int> PHREEQC_replace_position;
vector<REACT_PRQ*> REACT_PRQ_vec;

/**************************************************************************/
/* Constructor */
REACT_PRQ::REACT_PRQ(void)
{
    flag_prq = false;                 /* DL 28,10,08*/
    check_no_reaction_nodes = false;  // ToDo
    nodenumber = 0;
}
/* Destructor */
REACT_PRQ::~REACT_PRQ(void) {}
/**************************************************************************/

// split string line to pieces, and store in a vector
vector<std::string> REACT_PRQ::string2vector(std::string line)
{
    stringstream in;
    std::string sp;
    vector<std::string> pies;
    pies.clear();
    in.str(line);
    while (1)
    {
        if (in.eof())
            break;
        sp = "";
        in >> sp;
        if (sp != "")
            pies.push_back(sp);
    }
    return pies;
}

void REACT_PRQ::CreateREACT(void)
{
    int i, vector_size;
    CRFProcess* m_pcs = NULL;

    vector_size = (int)pcs_vector.size();
    for (i = 0; i < vector_size; i++)
    {
        m_pcs = pcs_vector[i];
        // if(m_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0){
        if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
        {
            nodenumber = (long)m_pcs->m_msh->GetNodesNumber(false);
            // elenumber = (long) m_pcs->m_msh->ele_vector.size();
        }
    }
    rateflag = (int*)Malloc(sizeof(int) * nodenumber);
}

vector<vector<std::string> > REACT_PRQ::file2vec(std::string phrq_file)
{
    int i;
    ifstream in;
    std::string instr;
    vector<std::string> pies, piesr;
    vector<vector<std::string> > res;
    in.open(phrq_file.c_str(), ios::in);

    res.clear();
    if (in.good())
    {
        while (1)
        {
            if (in.eof())
                break;
            getline(in, instr);
            pies = REACT_PRQ::string2vector(instr);
            piesr.clear();
            if ((int)pies.size() > 0 && pies[0].substr(0, 1) != "#")
            {
                for (i = 0; i < (int)pies.size(); i++)
                {
                    if (pies[i].substr(0, 1) == "#")
                        break;
                    else
                        piesr.push_back(pies[i]);
                }
                res.push_back(piesr);
            }
        }
    }
    return res;
}

int REACT_PRQ::isKey(std::string pies0)
{
    int i, res = -1;
    const std::string PHREEQC_KEYWORD[8] = {
        "SOLUTION",        "EQUILIBRIUM_PHASES", "KINETICS", "PRINT",
        "SELECTED_OUTPUT", "USER_PUNCH",         "END",      "RATES"};
    for (i = 0; i < 8; i++)
    {
        if (pies0 == PHREEQC_KEYWORD[i])
            res = i;
    }
    return res;
}

void REACT_PRQ::SetInterface(void)
{
    int i, j, ids, ips, no_pcs;  // number of process and transport species
    vector<int> idx, ipx, ipcs, id;
    CRFProcess* m_pcs = NULL;
    no_pcs = (int)pcs_vector.size();

    this->kin_no_steps = 1;
    this->pcs_name.clear();
    for (i = 0; i < no_pcs; i++)
    {
        m_pcs = pcs_vector[i];
        pcs_name.push_back(m_pcs->pcs_primary_function_name[0]);
    }

    id_key.clear();
    for (i = 0; i < (int)PHREEQC_TEMPLATE.size(); i++)
    {
        id_key.push_back(isKey(PHREEQC_TEMPLATE[i][0]));
    }

    ids = -1;
    idx.clear();
    for (i = 0; i < (int)id_key.size(); i++)
    {
        if (id_key[i] != -1)
            ids = id_key[i];
        idx.push_back(ids);
    }

    ipcs.clear();
    for (i = 0; i < (int)idx.size(); i++)
    {
        ips = -1;
        for (j = 0; j < no_pcs; j++)
        {
            if (pcs_name[j] == PHREEQC_TEMPLATE[i][0])
                ips = j;
        }
        ipcs.push_back(ips);
    }

    ips = -1;
    for (i = 0; i < (int)idx.size(); i++)
    {
        if (idx[i] == 2)
        {
            if (ipcs[i] != -1)
                ips = ipcs[i];
            ipcs[i] = ips;
        }
    }

    ipx.clear();
    for (i = 0; i < (int)idx.size(); i++)
    {
        ips = -1;
        if (ipcs[i] != -1)
        {
            if (idx[i] == 0)
                ips = 1;
            else if (idx[i] == 1)
                ips = 2;
            else if (idx[i] == 2 && PHREEQC_TEMPLATE[i][0] == "-m")
                ips = 1;
            else if (idx[i] == 2 && PHREEQC_TEMPLATE[i][0] == "-steps")
            {
                ips = 1;
                ipcs[i] = -2;  // mark for time step
                this->kin_no_steps = (int)atoi(PHREEQC_TEMPLATE[i][3].c_str());
            }
        }
        ipx.push_back(ips);
        // cout << " " << id_key[i] << " " << idx[i] << " " << ipcs[i] << " " <<
        // ipx[i] << "\n";
        //  idx-id_KEY, ipcs process, ipx position
    }

    this->phrq_id.clear();
    for (i = 0; i < (int)idx.size(); i++)
    {
        if (ipx[i] != -1)
            this->phrq_id.push_back(2);
        else if (idx[i] == 0 || idx[i] == 1 || idx[i] == 2 || idx[i] == 6)
            this->phrq_id.push_back(1);
        else
            this->phrq_id.push_back(0);
    }
    this->phrq_id_pos = ipx;
    this->phrq_id_pcs = ipcs;
    this->idx_key = idx;

    cout << " id_key idx_key id id_pos id_pcs "
         << "\n";
    for (i = 0; i < (int)phrq_id.size(); i++)
    {
        cout << right << setw(7) << id_key[i] << setw(8) << idx_key[i]
             << setw(3) << phrq_id[i] << setw(7) << phrq_id_pos[i] << setw(7)
             << phrq_id_pcs[i] << "\n";
    }
}

void REACT_PRQ::vec2file(int f)
{
    int i, j, ii, n;
    double pcs_node_value;
    ofstream out;
    CRFProcess* m_pcs = NULL;
    vector<int> phrq_id, phrq_id_pcs, phrq_id_pos;

    double unitfactor_l = 1, unitfactor_s = 1;
    REACTINT* m_rei = NULL;
    if (REACTINT_vec.size() > 0)
        m_rei = REACTINT_vec[0];

    double dt;
    CTimeDiscretization* m_tim = NULL;
    if (time_vector.size() > 0)
    {
        m_tim = time_vector[0];
        dt = m_tim->this_stepsize;
        // cout << " dt 0 " << dt << "\n";
    }

    out.clear();
    out.open(in_file.c_str(), ios::out);
    for (ii = 0; ii < this->nodenumber; ii++)
    {  // ii==0 as boundary point without reaction calc

        // CB 19.1.2011
        // calculate unit conversion factors for phreeqc molarity-->molality
        if (m_rei)
        {
            if (m_rei->unitconversion)
            {
                m_rei->CalcUnitConversionFactors(ii, &unitfactor_l,
                                                 &unitfactor_s, true);
                // unitfactor_l =  MOLH2OPERKG / m_rei->water_conc[ii];
                // unitfactor_s = (1 - m_rei->node_porosity[ii]) * MOLH2OPERKG /
                // (m_rei->water_conc[ii] * m_rei->node_porosity[ii] *
                // m_rei->GetWaterSaturation(ii)); if(unitfactor_s ==0)
                // unitfactor_s = (1 - m_rei->node_porosity[ii]) * MOLH2OPERKG /
                // (m_rei->water_conc[ii] * m_rei->node_porosity[ii] * 1);
                // cout << ii << " unitfactor_l " << unitfactor_l  << "
                // unitfactor_s " << unitfactor_s << "\n";
            }
        }

        for (i = 0; i < (int)this->phrq_id.size(); i++)
        {
            if (this->id_key[i] == -1)
                if (this->phrq_id[i] != 0 || ii == 0)
                    out << "  ";

            n = (int)PHREEQC_TEMPLATE[i].size();
            if (this->phrq_id[i] == 0 && ii == 0)
            {
                for (j = 0; j < n; j++)
                    out << PHREEQC_TEMPLATE[i][j] << " ";
                out << "\n";
            }

            if (this->phrq_id[i] == 1)
            {
                for (j = 0; j < n; j++)
                    out << PHREEQC_TEMPLATE[i][j] << " ";
                out << "\n";
            }

            if (this->phrq_id[i] == 2)
            {
                if (idx_key[i] != 1 || f != -1)
                {
                    for (j = 0; j < n; j++)
                    {
                        if (this->phrq_id_pos[i] == j)
                        {
                            if (this->phrq_id_pcs[i] >= 0)
                            {
                                m_pcs = pcs_vector[this->phrq_id_pcs[i]];
                                pcs_node_value = m_pcs->GetNodeValue(ii, 1);
                                if (ii == 1)
                                {  // DL to test
                                    cout << this->phrq_id_pcs[i] << " value "
                                         << pcs_node_value << "\n";
                                }
                            }
                            else if (this->phrq_id_pcs[i] ==
                                     -2)  // mark for current time step
                                pcs_node_value = dt;

                            // idx =
                            // pcs_vector[pqc_process[j]]->GetProcessComponentNumber();
                            ////mi,w = Ci,w * n *55.5 / CH2O
                            ////mi,s = Ci,w * (1-n) *55.5 / CH2O
                            // if(cp_vec[idx]->transport_phase==0) // liquid
                            // phase dval /= unitfactor_l; else
                            // if(cp_vec[idx]->transport_phase==1) // solid
                            // phase dval /= unitfactor_s;

                            if (this->phrq_id_pcs[i] >= 0)
                            {
                                if (idx_key[i] == 0)
                                {
                                    pcs_node_value *= unitfactor_l;
                                    if (pcs_node_value < 0)
                                        pcs_node_value = 0;  // 1.0e-16;

                                    //----for PHREEQC precision----
                                    if (pcs_node_value < 1.0e-9)
                                        pcs_node_value = 1.0e-9;
                                    //-----------------------------
                                }
                                else if (idx_key[i] == 1 || idx_key[i] == 2)
                                {
                                    pcs_node_value *= unitfactor_s;
                                    if (pcs_node_value < 0)  // 1.0e-16)
                                        pcs_node_value = 0;
                                }
                            }
                            out << pcs_node_value << " ";
                        }
                        else
                            out << PHREEQC_TEMPLATE[i][j] << " ";
                    }
                }
                out << "\n";
            }
        }
        out << "\n";
    }
    out.close();
}

int REACT_PRQ::Call_Phreeqc(void)
{
    std::string m_phrq;
    m_phrq = exe_file + " " + in_file + " " + out_file + " " + database_file;
    if (!system(m_phrq.c_str()))
    {
        DisplayMsgLn("Phreeqc runs succesfully! ");
        return 1;
    }
    else
    {
        DisplayMsgLn("Warnung: Phreeqc doesn't run properly!!! ");
        exit(0);
    }
}

void REACT_PRQ::file2pcs(int f)
{
    int i, ii, j, idx, iout_length = 0, iout_step;
    double pcs_node_value;
    ifstream in, in0;
    std::string instr;
    vector<int> pcs_id;
    vector<std::string> pies;
    CRFProcess* m_pcs = NULL;

    double unitfactor_l = 1, unitfactor_s = 1;
    REACTINT* m_rei = NULL;
    if (REACTINT_vec.size() > 0)
        m_rei = REACTINT_vec[0];

    in0.open(value_file.c_str(), ios::in);
    if (in0.good())
        while (1)
            if (in0.eof())
                break;
            else
            {
                getline(in0, instr);
                iout_length++;
            }
    else
        cout << " Warning ! Could not open value file "
             << "\n";
    in0.close();
    iout_step = (iout_length - 2) / nodenumber - 1;

    in.open(value_file.c_str(), ios::in);
    if (in.good())
    {
        pcs_id.clear();
        getline(in, instr);
        pies = string2vector(instr);
        for (i = 0; i < (int)pies.size(); i++)
        {
            pcs_id.push_back(-1);
            for (j = 0; j < (int)pcs_name.size(); j++)
                if (pies[i] == pcs_name[j])
                    pcs_id[i] = j;
        }

        for (ii = 0; ii < this->nodenumber; ii++)
        {
            // calculate unit conversion factors for phreeqc molarity-->molality
            if (m_rei)
            {
                if (m_rei->unitconversion)
                {
                    m_rei->CalcUnitConversionFactors(ii, &unitfactor_l,
                                                     &unitfactor_s, true);
                    // unitfactor_l =  MOLH2OPERKG / m_rei->water_conc[ii];
                    // unitfactor_s = (1 - m_rei->node_porosity[ii]) *
                    // MOLH2OPERKG / (m_rei->water_conc[ii] *
                    // m_rei->node_porosity[ii] *
                    // m_rei->GetWaterSaturation(ii)); if(unitfactor_s ==0)
                    // unitfactor_s = (1 - m_rei->node_porosity[ii]) *
                    // MOLH2OPERKG / (m_rei->water_conc[ii] *
                    // m_rei->node_porosity[ii] * 1);
                }
            }

            for (i = 0; i < iout_step; /*i<this->kin_no_steps;*/ i++)
                getline(in, instr);
            getline(in, instr);
            pies = string2vector(instr);
            for (i = 0; i < (int)pies.size(); i++)
            {
                if (pcs_id[i] != -1)
                {
                    m_pcs = pcs_vector[pcs_id[i]];
                    pcs_node_value = atof(pies[i].c_str());
                    // if(pcs_node_value< 1.0e-12) pcs_node_value=0;

                    idx = pcs_vector[pcs_id[i]]->GetProcessComponentNumber();
                    // mi,w = Ci,w * n *55.5 / CH2O
                    // mi,s = Ci,w * (1-n) *55.5 / CH2O
                    if (cp_vec[idx]->transport_phase == 0)  // liquid phase
                        pcs_node_value /= unitfactor_l;
                    else if (cp_vec[idx]->transport_phase == 1)  // solid phase
                        pcs_node_value /= unitfactor_s;

                    if (cp_vec[idx]->transport_phase == 0)  // liquid phase
                        m_pcs->SetNodeValue(ii, 1, pcs_node_value);
                    else if (cp_vec[idx]->transport_phase == 1)
                    {                // solid phase
                        if (f == 1)  // for full system
                            m_pcs->SetNodeValue(ii, 1, pcs_node_value);
                        else if (f == 0 && m_rei->icSolidUpdate)
                            m_pcs->SetNodeValue(ii, 1, pcs_node_value);
                    }
                }
            }
        }
    }
    else
        cout << " Warning ! Could not open value file "
             << "\n";
    in.close();
}

void REACT_PRQ::ExecuteReactionsPHRQ_new(int f)
{
    if (f == 0)
    {
        this->CreateREACT();
        PHREEQC_TEMPLATE = this->file2vec(phrq_file);
        this->SetInterface();
        this->vec2file(f);
        this->Call_Phreeqc();
        // cin.get();
        this->file2pcs(f);
    }
    if (f == 1 || f == -1)
    {
        this->vec2file(f);
        this->Call_Phreeqc();
        // cin.get();
        this->file2pcs(f);
    }
}

/**************************************************************************/
/*  Global function  */
/**************************************************************************/
bool REACT_PRQ_Read(std::string file_base_name)
{
    char line[MAX_ZEILE];
    std::string file_name_prq, line_string;
    ios::pos_type position;

    REACT_PRQ* rc_prq = new REACT_PRQ();
    // look if file is there
    file_name_prq = file_base_name + REACTION_EXTENSION_PHRQ;
    ifstream prq_file(file_name_prq.data(), ios::in);
    if (!prq_file.good())
    {
        delete rc_prq;
        rc_prq = NULL;
    }
    else
    {  // file is there - use NEW PHREEQC
        if (REACT_vec.capacity() > 0)
        {  // Test, if OLD PHREEQC is used also
            cout << "\n"
                 << " Warning!  PHREEQC is actived. NEW PHREEQC will NOT be "
                    "used ! "
                 << "\n";
            delete rc_prq;
            rc_prq = NULL;
        }
        else
        {
            rc_prq->flag_prq = true;
            // Read input file *.prq
            prq_file.clear();
            prq_file.seekg(0, ios::beg);
            cout << "PRQ_Read"
                 << "\n";
            prq_file.getline(line, MAX_ZEILE);  // first line
            prq_file.getline(line, MAX_ZEILE);  // second line ToDo
            line_string = line;
            if (line_string.find("#STOP") != string::npos)
                return true;

            // Call the object read function
            position = rc_prq->Read(&prq_file);
            // store data in global vector
            REACT_PRQ_vec.push_back(rc_prq);
            // close file
            prq_file.close();
        }
    }
    return true;
}

ios::pos_type REACT_PRQ::Read(ifstream* prq_file)
{
    bool new_keyword = false;
    std::string hash("#");
    std::string line_string;
    std::stringstream in;
    ios::pos_type position;

    while (!new_keyword)
    {
        line_string = GetLineFromFile1(prq_file);
        // if(line_string.size() < 1) break;
        if (line_string.find(hash) != string::npos)
        {
            new_keyword = true;
            break;
        }
        /* read keywords */
        //....................................................................
        if (line_string.find("$PHRQ_FILE") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(prq_file));
            in >> phrq_file;  // sub_line
            in.clear();
        }
        if (line_string.find("$EXE_FILE") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(prq_file));
            in >> exe_file;  // sub_line
            in.clear();
        }
        if (line_string.find("$DATABASE_FILE") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(prq_file));
            in >> database_file;  // sub_line
            in.clear();
        }
        if (line_string.find("$IN_FILE") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(prq_file));
            in >> in_file;  // sub_line
            in.clear();
        }
        if (line_string.find("$OUT_FILE") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(prq_file));
            in >> out_file;  // sub_line
            in.clear();
        }
        if (line_string.find("$VALUE_FILE") != string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(prq_file));
            in >> value_file;  // sub_line
            in.clear();
        }
    }  // end while
    return position;
}

#endif
