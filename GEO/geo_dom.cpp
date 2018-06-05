/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object: Domain
   Task:
   Programing:
   09/2003 OK Implementation
   09/2005 CC GeoLib2
**************************************************************************/
#include <stdio.h>
#include <string.h>
// C++ STL
// GeoLib
#include "files0.h"
#include "geo_dom.h"
// Vector
std::vector<CGLDomain*> domain_vector; // CC

/**************************************************************************
   GeoLib-Method: CGLDomain::Insert
   Task: insert to list
   Programing:
   09/2003 OK Implementation
**************************************************************************/
long CGLDomain::Insert(CGLDomain* m_domain)
{
	domain_vector.push_back(m_domain);
	return (long)domain_vector.size();
}

/**************************************************************************
   GeoLib-Method: CGLDomain::GetList
   Task:
   Programing:
   09/2003 OK Implementation
**************************************************************************/
std::vector<CGLDomain*> CGLDomain::GetVector(void)
{
	return domain_vector;
}
/*----------------------------------------------------------------------*/
// constructor
CGLDomain::CGLDomain(void)
{
	name = "DOMAIN";
}
// deconstructor
CGLDomain::~CGLDomain(void)
{
}
/**************************************************************************
   GeoLib-Method: GEOReadVolume
   Task: Read volume data from file
   Programing:
   07/2003 OK Implementation
**************************************************************************/
int CGLDomain::Read(char* data, FILE* f)
{
	int pos = 0, pos_s = 0;
	int p = 0;
	char* sub;
	int begin;
	int ok = 1;
	int p_sub = 0;
	char name[80];
	double ddummy;

	LineFeed(f);
	FilePrintString(f, "; ------------------------------------------");
	LineFeed(f);
	FilePrintString(f, "; GeoLib - Domain");
	LineFeed(f);

	//---------------------------------------------------------------------
	// Loop over all volumes
	while (StrTestHash(&data[p += pos], &pos))
	{
		CGLDomain* m_domain = NULL;
		m_domain = new CGLDomain;
		/* Write keyword */
		LineFeed(f);
		FilePrintString(f, "#DOMAIN");
		LineFeed(f);
		//-------------------------------------------------------------------
		// Check sub keywords
		sub = new char[(int)strlen(data) + 2];
		while (StrReadSubKeyword(sub, data, p += pos, &begin, &p))
		{
			ok = StrReadStr(name, sub, f, /*TFString,*/ &p_sub) && ok;
			pos = 0;
			//-----------------------------------------------------------------
			if (!strcmp(name, "$NAME"))
			{
				ok = (StrReadStr(name, &sub[p_sub], f, /*TFString,*/ &pos) && ok);
				LineFeed(f);
				m_domain->name = name;
			}
			//-----------------------------------------------------------------
			if (!strcmp(name, "$COORDINATES"))
			{
				pos_s = 0;
				ok = (StrReadDouble(&ddummy, &sub[p_sub += pos_s], f, &pos_s) && ok);
				m_domain->x_min = ddummy;
				ok = (StrReadDouble(&ddummy, &sub[p_sub += pos_s], f, &pos_s) && ok);
				m_domain->x_max = ddummy;
				ok = (StrReadDouble(&ddummy, &sub[p_sub += pos_s], f, &pos_s) && ok);
				m_domain->y_min = ddummy;
				ok = (StrReadDouble(&ddummy, &sub[p_sub += pos_s], f, &pos_s) && ok);
				m_domain->y_max = ddummy;
				LineFeed(f);
			}
			pos = 0;
		} // sub-keyword
		delete (sub);
		// insert into list
		Insert(m_domain);
	}

	return 1;
}

int GEOReadDomain(char* data, int found, FILE* f)
{
	CGLDomain* m_domain = NULL;
	m_domain->Read(data, f);
	found = found;
	// delete(m_domain);
	return 1;
}
/**************************************************************************
   GeoLib-Method: GEOGetVolume
   Task: select volume instance from list by name
   Programing:
   07/2003 OK Implementation
**************************************************************************/
CGLDomain* CGLDomain::Get(std::string name)
{
	CGLDomain* m_domain;
	std::vector<CGLDomain*>::iterator p = domain_vector.begin(); // CC
	while (p != domain_vector.end())
	{
		m_domain = *p;
		if (m_domain->name == name)
			return m_domain;
		++p;
	}
	return NULL;
}
