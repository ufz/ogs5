/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object: Volume
   Task:
   Programing:
   07/2003 OK Implementation
   09/2005 CC GeoLib2
**************************************************************************/

#include <stdlib.h>
// GEOLib
#include "files0.h"
#include "geo_lib.h"
#include "geo_pnt.h"
#include "geo_sfc.h"
#include "geo_vol.h"

using namespace std;  // 11.08.2011. WW
//------------------------------------------------------------------------
std::vector<CGLVolume*> volume_vector;  // CC
//////////////////////////////////////////////////////////////////////////
// Construction
CGLVolume::CGLVolume(void)
{
    name = "VOLUME";
    selected = false;
    mat_group = 0;
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
**************************************************************************/
CGLVolume::~CGLVolume(void)
{
    surface_vector.clear();
    // surface_name_list.clear();
}
//////////////////////////////////////////////////////////////////////////
// I/O
/**************************************************************************
   GEOLib-Method:
   Task: Volumes read function
   Programing:
   03/2004 OK Implementation
**************************************************************************/
void GEOReadVolumes(std::string file_name_path_base)
{
    CGLVolume* m_volume = NULL;
    char line[MAX_ZEILEN];
    std::string sub_line;
    std::string line_string;
    std::string gli_file_name;
    std::ios::pos_type position;
    //========================================================================
    // File handling
    gli_file_name = file_name_path_base + ".gli";
    std::ifstream gli_file(gli_file_name.data(), std::ios::in);
    if (!gli_file.good())
        return;
    gli_file.seekg(0L, std::ios::beg);  // rewind?
    //========================================================================
    // Keyword loop
    while (!gli_file.eof())
    {
        gli_file.getline(line, MAX_ZEILEN);
        line_string = line;
        if (line_string.find("#STOP") != string::npos)  // 11.08.2011. WW
            break;

        //----------------------------------------------------------------------
        if (line_string.find("#VOLUME") != std::string::npos)  // keyword found
        {
            m_volume = new CGLVolume();
            position = m_volume->Read(&gli_file);
            volume_vector.push_back(m_volume);
            gli_file.seekg(position, std::ios::beg);
        }  // keyword found
    }      // eof
    gli_file.close();
}

/**************************************************************************
   GeoLib-Method: Read
   Task: Read volume data from file
   Programing:
   03/2004 OK Implementation
   04/2005 OK/PCH Reading surfaces
   04/2005 PCH Modified to read according to the change made in Write function
   09/2005 OK $LAYER
   09/2005 OK MAT group name
   ToDo: CC streaming
**************************************************************************/
std::ios::pos_type CGLVolume::Read(std::ifstream* gli_file)
{
    char line[MAX_ZEILEN];
    std::string sub_line;
    std::string line_string;
    std::string delimiter(",");
    bool new_keyword = false;
    std::string hash("#");
    std::ios::pos_type position;
    std::string sub_string;
    Surface* m_surface = NULL;
    std::stringstream in;
    std::string sfc_name;
    bool ok_true = true;
    //========================================================================
    // Schleife ueber alle Phasen bzw. Komponenten
    while (!new_keyword)
    {
        position = gli_file->tellg();
        gli_file->getline(line, MAX_ZEILEN);
        line_string = line;
        if (line_string.find(hash) != std::string::npos)
        {
            new_keyword = true;
            break;
        }
        //....................................................................
        if (line_string.find("$NAME") != std::string::npos)  // subkeyword found
        {
            gli_file->getline(line, MAX_ZEILEN);
            line_string = line;
            remove_white_space(&line_string);
            name = line_string.substr(0);
        }  // subkeyword found
        //....................................................................
        if (line_string.find("$TYPE") != std::string::npos)  // subkeyword found
        {
            gli_file->getline(line, MAX_ZEILEN);
            line_string = line;
            remove_white_space(&line_string);
            type_name = line_string;
            type = strtol(line_string.data(), NULL, 0);
        }  // subkeyword found
        //....................................................................
        if (line_string.find("$SURFACES") !=
            std::string::npos)  // subkeyword found
            while (ok_true)
            {
                line_string = GetLineFromFile1(gli_file);
                in.str(line_string);
                if (line_string.find("$") != std::string::npos)
                {
                    in.clear();
                    break;
                }
                if (line_string.find("#") != std::string::npos)
                {
                    in.clear();
                    return position;
                }
                in >> sfc_name;
                m_surface = GEOGetSFCByName(sfc_name);  // CC
                if (m_surface)
                    surface_vector.push_back(m_surface);
                in.clear();
            }
        // subkeyword found
        //....................................................................
        if (line_string.find("$MAT_GROUP") !=
            std::string::npos)  // subkeyword found
        {
            gli_file->getline(line, MAX_ZEILEN);
            line_string = line;
            // OK mat_group = strtol(line_string.data(),NULL,0);
            remove_white_space(&line_string);
            mat_group_name = line_string;
        }  // subkeyword found
        //....................................................................
        if (line_string.find("$LAYER") !=
            std::string::npos)  // subkeyword found
        {
            gli_file->getline(line, MAX_ZEILEN);
            line_string = line;
            layer = strtol(line_string.data(), NULL, 0);
            if (layer > 0)
                type_name = "LAYER";
        }  // subkeyword found
    }
    return position;
}

/**************************************************************************
   GeoLib-Method:
   Task: Write volume data
   Programing:
   01/2005 OK Implementation
   04/2005 PCH Modified to include surface_list in .gli files
   09/2005 OK $LAYER
   09/2005 OK MAT group name
   11/05 CC Write function
**************************************************************************/
void CGLVolume::Write(std::string path_name)
{
    const char* char_string;
    Surface* m_sfc = NULL;
    //-----------------------------------------------------------------------
    // File handling
    FILE* gli_file = NULL;
    std::string gli_file_name;
    gli_file_name = path_name + ".gli";
    const char* gli_file_name_char = 0;
    gli_file_name_char = gli_file_name.data();
    gli_file = fopen(gli_file_name_char, "a");
    if (gli_file == NULL)
        return;
    else
    {
        fprintf(gli_file, "%s\n", "#VOLUME");
        //-----------------------------------------------------------------
        fprintf(gli_file, " %s\n", "$NAME");
        char_string = name.data();
        fprintf(gli_file, "  %s\n", char_string);

        // Instead of the lines above, the following lines are added by PCH
        int surface_list_size = (int)surface_vector.size();  // CC
        if (surface_list_size > 0)
        {
            fprintf(gli_file, " %s\n", "$SURFACES");
            std::vector<Surface*>::iterator p_sfc =
                surface_vector.begin();            // CC
            while (p_sfc != surface_vector.end())  // CC
            {
                m_sfc = *p_sfc;
                fprintf(gli_file, "  %s\n", m_sfc->name.c_str());
                ++p_sfc;
            }
        }
    }  // dat_in
    fclose(gli_file);
    //-----------------------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////
// Access
/**************************************************************************
   GeoLib-Method: GEOGetVolume
   Task: select volume instance from list by name
   Programing:
   07/2003 OK Implementation
   09/2005 CC m_volume
**************************************************************************/
CGLVolume* GEOGetVOL(std::string name)
{
    CGLVolume* m_volume = NULL;
    std::vector<CGLVolume*>::iterator p = volume_vector.begin();  // CC
    while (p != volume_vector.end())
    {
        m_volume = *p;
        if (m_volume->name == name)
            return m_volume;
        ++p;
    }
    return NULL;
}

//////////////////////////////////////////////////////////////////////////
// GEO-Methods
/**************************************************************************
   GeoLib-Method: PointInVolume
   Task:
   Programing:
   03/2004 OK Implementation
   08/2005 CC
**************************************************************************/
bool CGLVolume::PointInVolume(CGLPoint* m_point, int dim_type)
{
    std::vector<Surface*>::iterator p_sfc;
    std::string surface_name;
    Surface* m_surface = NULL;
    bool ok = false;
    long i;
    CGLPoint m_element_point;
    std::vector<CTriangle*> prism_triangles;
    long surface_TIN_length;
    double xp[6], yp[6], zp[6];
    CTriangle* m_triangle = NULL;
    bool found = false;

    //=======================================================================
    switch (dim_type)
    {
        //---------------------------------------------------------------------
        case 0:  // 2.5z-D (vertical extension)
            //...................................................................
            // check x-y
            p_sfc = surface_vector.begin();  // CC
            while (p_sfc != surface_vector.end())
            {
                m_surface = *p_sfc;
                surface_TIN_length = (long)m_surface->TIN->Triangles.size();
                for (i = 0; i < surface_TIN_length; i++)
                {
                    m_triangle = m_surface->TIN->Triangles[i];
                    xp[0] = m_triangle->x[0];
                    yp[0] = m_triangle->y[0];
                    zp[0] = m_triangle->z[0];
                    xp[1] = m_triangle->x[1];
                    yp[1] = m_triangle->y[1];
                    zp[1] = m_triangle->z[1];
                    xp[2] = m_triangle->x[2];
                    yp[2] = m_triangle->y[2];
                    zp[2] = m_triangle->z[2];
                    if (m_point->IsInTriangleXYProjection(xp, yp))
                    {
                        found = true;
                        prism_triangles.push_back(m_triangle);
                        continue;
                    }
                }
                ++p_sfc;
            }
            if (!found)
                return false;
            // should be 2 triangles
            int triangles_vector_length = (int)prism_triangles.size();
            for (i = 0; i < triangles_vector_length; i++)
            {
                xp[0 + 3 * i] = prism_triangles[i]->x[0];
                yp[0 + 3 * i] = prism_triangles[i]->y[0];
                zp[0 + 3 * i] = prism_triangles[i]->z[0];
                xp[1 + 3 * i] = prism_triangles[i]->x[1];
                yp[1 + 3 * i] = prism_triangles[i]->y[1];
                zp[1 + 3 * i] = prism_triangles[i]->z[1];
                xp[2 + 3 * i] = prism_triangles[i]->x[2];
                yp[2 + 3 * i] = prism_triangles[i]->y[2];
                zp[2 + 3 * i] = prism_triangles[i]->z[2];
            }
            ok = m_point->IsInsidePrism(xp, yp, zp);
            break;
    }
    // Memory
    prism_triangles.clear();

    return ok;
}

//////////////////////////////////////////////////////////////////////////
// Old stuff
/**************************************************************************
   GeoLib-Method: GEOInsertVolume2List
   Task: construct volume instance
   Programing:
   07/2003 OK Implementation
**************************************************************************/
long CGLVolume::GEOInsertVolume2List(CGLVolume* m_volume)
{
    volume_vector.push_back(m_volume);
    return (long)volume_vector.size();
}
void CGLVolume::AddSurface(Surface* P_Surface)
{
    surface_vector.push_back(P_Surface);  // CC
}
/**************************************************************************
   GeoLib-Method: GEOGetVolumes
   Task:
   Programing:
   07/2003 OK Implementation
   09/2005 CC
**************************************************************************/
std::vector<CGLVolume*> GEOGetVolumes(void)
{
    return volume_vector;
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   07/2003 OK Implementation
   09/2005 CC
**************************************************************************/
CGLVolume* GetVolume(long nsel)
{
    CGLVolume* m_volume = NULL;
    std::vector<CGLVolume*>::iterator p = volume_vector.begin();  // CC
    long iter = 0;
    while (p != volume_vector.end())
    {
        m_volume = *p;

        ++p;

        if (iter == nsel)
        {
            return m_volume;
            break;
        }
        ++iter;
    }
    return NULL;
}
/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   08/2005 CC
**************************************************************************/
void GEOWriteVolumes2TecplotV1(std::string file_base_name)
{
    (void)file_base_name;
    long i;
    CGLVolume* m_volume = NULL;
    std::vector<CGLVolume*>::iterator p_vol = volume_vector.begin();  // CC
    std::string name;
    Surface* m_surface = NULL;
    std::vector<Surface*>::iterator p_vol_sfc;
    long no_triangles;
    long no_nodes;
    std::string delimiter(", ");

    while (p_vol != volume_vector.end())
    {
        m_volume = *p_vol;
        //--------------------------------------------------------------------
        // file handling
        std::string vol_file_name =
            "VOL_" + m_volume->name + TEC_FILE_EXTENSIONS;
        std::fstream vol_file(vol_file_name.data(),
                              std::ios::trunc | std::ios::out);
        vol_file.setf(std::ios::scientific, std::ios::floatfield);
        vol_file.precision(12);
        if (!vol_file.good())
            return;
        vol_file.seekg(0L, std::ios::beg);
        //--------------------------------------------------------------------
        vol_file << "VARIABLES = X,Y,Z"
                 << "\n";
        //--------------------------------------------------------------------
        p_vol_sfc = m_volume->surface_vector.begin();        // CC
        while (p_vol_sfc != m_volume->surface_vector.end())  // CC
        {
            m_surface = *p_vol_sfc;
            //..................................................................
            if (m_surface->TIN && (m_surface->TIN->Triangles.size() > 0))
            {
                no_triangles = (long)m_surface->TIN->Triangles.size();
                no_nodes = 3 * no_triangles;
                vol_file << "ZONE T = " << m_surface->TIN->name << delimiter
                         << "N = " << no_nodes << delimiter
                         << "E = " << no_triangles << delimiter << "F = FEPOINT"
                         << delimiter << "ET = TRIANGLE"
                         << "\n";
                for (i = 0; i < no_triangles; i++)
                {
                    vol_file << m_surface->TIN->Triangles[i]->x[0] << " "
                             << m_surface->TIN->Triangles[i]->y[0] << " "
                             << m_surface->TIN->Triangles[i]->z[0] << "\n";
                    vol_file << m_surface->TIN->Triangles[i]->x[1] << " "
                             << m_surface->TIN->Triangles[i]->y[1] << " "
                             << m_surface->TIN->Triangles[i]->z[1] << "\n";
                    vol_file << m_surface->TIN->Triangles[i]->x[2] << " "
                             << m_surface->TIN->Triangles[i]->y[2] << " "
                             << m_surface->TIN->Triangles[i]->z[2] << "\n";
                }
                for (i = 0; i < no_triangles; i++)
                    vol_file << 3 * i + 1 << " " << 3 * i + 2 << " "
                             << 3 * i + 3 << "\n";
            }
            ++p_vol_sfc;
        }
        //--------------------------------------------------------------------
        // bounding surface
        // 1 - write header
        p_vol_sfc = m_volume->surface_vector.begin();  // CC
        while (p_vol_sfc != m_volume->surface_vector.end())
        {
            m_surface = *p_vol_sfc;
            //..................................................................
            no_nodes = (long)m_surface->polygon_point_vector.size();
            long no_elements = (long)m_surface->polygon_point_vector.size() - 1;
            vol_file << "ZONE T = " << m_surface->TIN->name << delimiter
                     << "N = " << 2 * no_nodes << delimiter
                     << "E = " << no_elements << delimiter << "F = FEPOINT"
                     << delimiter << "ET = QUADRILATERAL"
                     << "\n";
            //++p_vol_sfc;
            break;
        }
        //--------------------------------------------------------------------
        // 2 - write nodes
        p_vol_sfc = m_volume->surface_vector.begin();
        while (p_vol_sfc != m_volume->surface_vector.end())
        {
            m_surface = *p_vol_sfc;
            //..................................................................
            no_nodes = (long)m_surface->polygon_point_vector.size();
            for (i = 0; i < no_nodes; i++)
                vol_file << m_surface->polygon_point_vector[i]->x << " "
                         << m_surface->polygon_point_vector[i]->y << " "
                         << m_surface->polygon_point_vector[i]->z << "\n";
            //++p_vol_sfc;
            for (i = 0; i < no_nodes; i++)
                vol_file << m_surface->polygon_point_vector[i]->x << " "
                         << m_surface->polygon_point_vector[i]->y << " "
                         << m_surface->polygon_point_vector[i]->z - 1000.
                         << "\n";
            break;
        }
        //--------------------------------------------------------------------
        // 3 - write elements
        p_vol_sfc = m_volume->surface_vector.begin();
        while (p_vol_sfc != m_volume->surface_vector.end())
        {
            m_surface = *p_vol_sfc;
            //..................................................................
            no_nodes = (long)m_surface->polygon_point_vector.size();
            long no_elements = (long)m_surface->polygon_point_vector.size() - 1;
            for (i = 0; i < no_elements; i++)
                vol_file << i + 1 << " " << i + 2 << " " << i + no_nodes + 2
                         << " " << i + no_nodes + 1 << "\n";
            //++ps;
            break;
        }
        //--------------------------------------------------------------------
        ++p_vol;
    }
}

/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   03/2006 CC
**************************************************************************/
void GEOWriteVolumes(std::string path_name)
{
    CGLVolume* m_vol = NULL;
    for (int i = 0; i < (int)volume_vector.size(); i++)
    {
        m_vol = volume_vector[i];
        m_vol->Write(path_name);
    }
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation, list destructor
   09/2005 CC clear
   11/2005 TK "delete" loops
**************************************************************************/
void GEORemoveAllVolumes()
{
    int i = 0;
    for (i = 0; i < (int)volume_vector.size(); i++)
    {
        delete volume_vector[i];
        volume_vector[i] = NULL;
    }
    volume_vector.clear();
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   10/2005 OK Implementation
**************************************************************************/
void GEORemoveVOL(CGLVolume* m_vol)
{
    CGLVolume* m_vol_this = NULL;
    for (int i = 0; i < (int)volume_vector.size(); i++)
    {
        m_vol_this = volume_vector[i];
        if (m_vol_this->name.compare(m_vol->name) == 0)
        {
            delete m_vol_this;
            volume_vector.erase(volume_vector.begin() + i);
            return;
        }
    }
}
