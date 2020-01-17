/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>
// using namespace std;

#ifdef OGS_FEM_CAP  // CAP_REACT

class REACT_PRQ
{
private:
public:
    REACT_PRQ(void);   // Initialize
    ~REACT_PRQ(void);  // Clean up

    // Data

    std::string phrq_file, exe_file, database_file, in_file, out_file,
        value_file;
    int* rateflag;   /* flag used for determining if reaction are calculated */
    long nodenumber; /* number of nodes, on which reactions are calculated */
    bool flag_prq;   /* flag if *.prq file exists   DL 28,10,08*/
    bool check_no_reaction_nodes; /* flag if CheckNoReactionNodes has been
                                     performed */
    int kin_no_steps;

    std::vector<int> id_key, idx_key, phrq_id, phrq_id_pcs, phrq_id_pos;

    std::vector<std::string> pcs_name;

    // Member functions

    std::ios::pos_type Read(std::ifstream*);

    std::vector<std::vector<std::string> > file2vec(std::string);

    static std::vector<std::string> string2vector(
        std::string
            line);  // split string line to pieces, and store in a vector
    int isKey(std::string);  // return the no. in Key words list

    void CreateREACT(void);
    void SetInterface(void);

    void vec2file(int f);
    int Call_Phreeqc(void);
    void file2pcs(int f);

    void ExecuteReactionsPHRQ_new(int f);
};

extern std::vector<std::vector<std::string> > PHREEQC_TEMPLATE;
extern std::vector<REACT_PRQ*> REACT_PRQ_vec;
// extern bool REACT_PRQ_Read2(std::string);
extern bool REACT_PRQ_Read(std::string);

#endif
