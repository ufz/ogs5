/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   ROCKFLOW - Modul: timer.c

   Aufgabe:
   Funktionen zur Laufzeitermittlung im Testbetrieb.
   Es werden beliebig viele Zeitspeicher bereitgestellt.

   Programmaenderungen:
   07/1994     MSR        Erste Version
   6/1997      C.Thorenz  Komplett neue zweite Version
   1/1999      C.Thorenz  Dritte Version: CPU-Zeit auf POSIX-Rechner
   09/1999     AH         Funktionen: TGetTime und TGetTicksPerSecond global.
   11/1999     C.Thorenz  Beliebige Anzahl Zeitspeicher
**************************************************************************/
#include "makros.h"
#include "timer.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
//#include<windows.h>
using namespace std;

/* Auf POSIX-Rechern ist exaktere Zeitmessung vorhanden */
#ifdef _POSIX_SOURCE
#ifndef WIN32
#include <sys/times.h>
#include <unistd.h>
#endif // WIN32
#include <time.h>
#endif

#include "makros.h"
#include "memory.h"

/* Zeitspeicher */
static int max_zeitspeicher = -1;
static long* zeit = NULL;
static int* running = NULL;
vector<CClockTime*> ClockTimeVec;

/*************************************************************************
   ROCKFLOW - Funktion: TInitTimer

   Aufgabe:
   Setzt Zeitspeicher auf 0 aber startet den Timer noch nicht. Die
   Funktion muss vor TGetTimer aufgerufen werden ! Der Timer0 ist
   bei ROCKFLOW fuer die Gesamtlaufzeit reserviert.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int speicher: Nummer (0..9) des Zeitspeichers

   Ergebnis:
   - void -

   Programmaenderungen:
   07/1994     MSR        Erste Version
   6/1997      C.Thorenz  Komplett neue zweite Version
   1/1999      C.Thorenz  Dritte Version: CPU-Zeit auf POSIX-Rechner
 **************************************************************************/
void TInitTimer(int speicher)
{
	/* Ggf. Liste der Speicher erweitern */
	if (speicher > max_zeitspeicher)
	{
		zeit = (long*)Realloc(zeit, sizeof(long) * (speicher + 1));
		running = (int*)Realloc(running, sizeof(int) * (speicher + 1));
		max_zeitspeicher = speicher;
	}

	/* Der Timer wird "genullt", aber nicht gestartet. */
	zeit[speicher] = 0;
	running[speicher] = 0;
}

/*************************************************************************
   ROCKFLOW - Funktion: TStartTimer

   Aufgabe:
   Setzt Timer auf 0 und startet den Timer neu. Der Timer0 ist bei
   ROCKFLOW fuer die Gesamtlaufzeit reserviert.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int speicher: Nummer (0..x) des Zeitspeichers

   Ergebnis:
   - void -

   Programmaenderungen:
   07/1994     MSR        Erste Version
   6/1997      C.Thorenz  Komplett neue zweite Version
   1/1999      C.Thorenz  Dritte Version: CPU-Zeit auf POSIX-Rechner
 **************************************************************************/
void TStartTimer(int speicher)
{
	/* Ggf. Liste der Speicher erweitern */
	if (speicher > max_zeitspeicher)
	{
		zeit = (long*)Realloc(zeit, sizeof(long) * (speicher + 1));
		running = (int*)Realloc(running, sizeof(int) * (speicher + 1));
		max_zeitspeicher = speicher;
	}

	/* Der Timer wird "genullt" indem die aktuelle Zeit im Speicher abgelegt
	   wird. */
	zeit[speicher] = TGetTime();
	running[speicher] = 1;
}

/**************************************************************************
   ROCKFLOW - Funktion: TGetTimerDouble

   Aufgabe:
   Liefert Zeitdifferenz seit dem letzten Aufruf von TStartTimer
   als Double mit der Aufloesung, die zur Verfuegung steht.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int speicher       : Nummer (0..9) des Zeitspeichers
   R double time_gone_by: Laufzeit

   Ergebnis:
   (CPU-)Zeitdifferenz seit dem letzten Aufruf von TStartTimer

   Programmaenderungen:
   07/1994     MSR        Erste Version
   6/1997      C.Thorenz  Komplett neue zweite Version
   1/1999      C.Thorenz  Dritte Version: CPU-Zeit auf POSIX-Rechner
**************************************************************************/
double TGetTimerDouble(int speicher)
{
	double time_gone_by;

	if (!running[speicher])
		/* Der Timer war angehalten */
		time_gone_by = (double)zeit[speicher] / (double)TGetTicksPerSecond();
	else
		/* Der Timer lief */
		time_gone_by = (double)(TGetTime() - zeit[speicher]) / (double)TGetTicksPerSecond();

	return time_gone_by;
}

/**************************************************************************
   ROCKFLOW - Funktion: TGetTimer

   Aufgabe:
   Liefert Zeitdifferenz seit dem letzten Aufruf von TStartTimer
   als long mit Sekunden-Aufloesung. Wenn eine hoehere Aufloesung
   zur Verfuegung steht, wird diese intern genutzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int speicher: Nummer (0..9) des Zeitspeichers
   R long time_gone_by: Laufzeit

   Ergebnis:
   (CPU-)Zeitdifferenz seit dem letzten Aufruf von TStartTimer

   Programmaenderungen:
   07/1994     MSR        Erste Version
   6/1997      C.Thorenz  Komplett neue zweite Version
   1/1999      C.Thorenz  Dritte Version: CPU-Zeit auf POSIX-Rechner
   08/2007 OK Test
**************************************************************************/
long TGetTimer(int speicher)
{
	if (!running) // OK
		return -1;
	long time_gone_by;
	if (!running[speicher])
		/* Der Timer war angehalten */
		time_gone_by = (long)(zeit[speicher] / TGetTicksPerSecond());
	else
		/* Der Timer lief */
		time_gone_by = (long)((TGetTime() - zeit[speicher]) / TGetTicksPerSecond());
	return time_gone_by;
}

/**************************************************************************
   ROCKFLOW - Funktion: TStopTimer

   Aufgabe:
   Haelt den zug. Timer an.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int speicher: Nummer (0..9) des Zeitspeichers

   Ergebnis:

   Programmaenderungen:
   6/1997   C.Thorenz        Erste Version
**************************************************************************/
void TStopTimer(int speicher)
{
	/* Im Speicher wird die bisher verstrichene Zeit abgelegt.
	   Wird nur bei laufendem Timer ausgefuehrt. */

	if (running[speicher])
	{
		zeit[speicher] = TGetTime() - zeit[speicher];
		running[speicher] = 0;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: TRestartTimer

   Aufgabe:
   Laesst den zugehoerigen Timer wieder weiterlaufen.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int speicher: Nummer (0..9) des Zeitspeichers

   Ergebnis:

   Programmaenderungen:
   6/1997   C.Thorenz        Erste Version

**************************************************************************/
void TRestartTimer(int speicher)
{
	/* Im Speicher liegt die bisher vom Timer gezaehlte Zeit. Mit der
	   aktuellen Zeit wird die Startzeit errechnet.
	   Wird nur bei angehaltenem Timer ausgefuehrt. */

	if (!running[speicher])

	{
		zeit[speicher] = TGetTime() - zeit[speicher];
		running[speicher] = 1;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: TGetTime

   Aufgabe:
   Probiert eine "Zeit" zu erhalten. In ANSI-C ist das die
   absolute Uhrzeit in Sekunden, bei POSIX-Systemen (UNIX) die
   echte Prozessorzeit in Ticks.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R long zeit: Zeit

   Ergebnis:

   Programmaenderungen:
   1/1999   C.Thorenz        Erste Version

**************************************************************************/
long TGetTime(void)
{
	long runtime;

#ifdef _POSIX_SOURCE
	runtime = clock();
#else
	runtime = (long)time(NULL);
#endif

	return runtime;
}

/**************************************************************************
   ROCKFLOW - Funktion: TGetTicksPerSecond

   Aufgabe:
   Liefert die Aufloesung der mit TGetTime erhaltenen
   "ticks" in 1/Sekunde.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   R long ticks: Ticks pro Sekunde

   Ergebnis:

   Programmaenderungen:
   1/1999   C.Thorenz        Erste Version

**************************************************************************/
long TGetTicksPerSecond(void)
{
	long TicksPerSecond;

	TicksPerSecond = 1;

#ifdef _POSIX_SOURCE
	TicksPerSecond = CLOCKS_PER_SEC;
#endif

	return TicksPerSecond;
}

/*************************************************************************
   ROCKFLOW - Funktion: TDestroyTimers

   Aufgabe:
   Zerstoert alle Zeitspeicher

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)

   Ergebnis:
   - void -

   Programmaenderungen:
   10/1999      C.Thorenz  Erste Version
 **************************************************************************/
void TDestroyTimers(void)
{
	/* Speicherfreigaben */
	zeit = (long*)Free(zeit);
	running = (int*)Free(running);
}

/*************************************************************************
   ROCKFLOW - Funktion: ctime_

   Aufgabe:
   Interface zu Timing des AMG-Loesers
   R double* time: CPU-Zeit

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)

   Ergebnis:
   - void -

   Programmaenderungen:
   10/2001      C.Thorenz  Erste Version
 **************************************************************************/
void ctime_(float* time)
{
	*time = (float)TGetTime() / (float)TGetTicksPerSecond();
}

// New SB time

CClockTime::CClockTime(void)
{
	delta_clocktime = 0.0;
	time_flow.clear();
	time_transport.clear();
	time_kinreact.clear();
	time_equireact.clear();
	time_reactdeact.clear();
	start = 0;
	end = 0;
	time_total_flow = 0.0;
	time_total_transport = 0.0;
	time_total_kinreact = 0.0;
	time_total_equireact = 0.0;
	time_total_reactdeact = 0.0;
	time1 = 0;
	time2 = 0;
	difftime = 0;
}

CClockTime::~CClockTime(void)
{
	delta_clocktime = 0.0;
	time_flow.clear();
	time_transport.clear();
	time_kinreact.clear();
	time_equireact.clear();
	time_reactdeact.clear();
}

void CClockTime::StartTime(void)
{
	start = clock();
	// WW	time1=GetTickCount();
}

void CClockTime::StopTime(const std::string& name)
{
	char name1;
	name1 = name[0];

	end = clock();
	this->delta_clocktime = (double)(end - start) / CLOCKS_PER_SEC;

	// WW time2=GetTickCount();
	difftime = (time2 - time1) / 1000.0;
	// cout << " ClockTime: " << delta_clocktime << ", TickTime: " << difftime << "\n";

	switch (name1)
	{
		default:
			break;
		case ('F'):
			time_flow.push_back(delta_clocktime);
			time_total_flow += delta_clocktime;
			break;
		case ('T'):
			this->time_transport.push_back(delta_clocktime);
			time_total_transport += delta_clocktime;
			break;
		case ('K'):
			this->time_kinreact.push_back(delta_clocktime);
			time_total_kinreact += delta_clocktime;
			break;
		case ('E'):
			this->time_equireact.push_back(delta_clocktime);
			time_total_equireact += delta_clocktime;
			break;
		case ('R'):
			this->time_reactdeact.push_back(delta_clocktime);
			time_total_reactdeact += delta_clocktime;
			break;
	}
}

void CClockTime::PrintTimes(void)
{
#ifndef OGS_OUTPUT_TIME_ALLPCS
	if (time_kinreact.size() == 0 || time_equireact.size() == 0 || time_reactdeact.size() == 0) // WW
		return;
#endif

	int i, length;
	double tot = 0., help = 0.0, tot_zeitschritt = 0.;
	string outname = FileName + "_ClockTimes.txt";

	cout.precision(2);
	tot = time_total_flow + time_total_transport + time_total_kinreact + time_total_equireact + time_total_reactdeact;
	cout << "ClockTimes: "
	     << "\n"
	     << "Unit   Flow:  Transport:  KinReactions:  EquiReactions:  ReactDeact:  total: "
	     << "\n";
	cout << "[sec] " << setw(6) << time_total_flow << "  " << setw(10) << time_total_transport << "  " << setw(13)
	     << time_total_kinreact << "  " << setw(14) << time_total_equireact << "  " << setw(14) << time_total_reactdeact
	     << "  " << setw(6) << tot << "\n";
	cout << "[%]   " << setw(6) << time_total_flow / tot * 100 << "  " << setw(10) << time_total_transport / tot * 100
	     << "  " << setw(13) << time_total_kinreact / tot * 100 << "  " << setw(14) << time_total_equireact / tot * 100
	     << "  " << setw(14) << time_total_reactdeact / tot * 100 << "  " << setw(6) << tot / tot * 100 << "\n";

	length = (int)this->time_flow.size();
	if ((int)this->time_transport.size() > length)
		length = (int)this->time_transport.size();
	if ((int)this->time_kinreact.size() > length)
		length = (int)this->time_kinreact.size();
	if ((int)this->time_equireact.size() > length)
		length = (int)this->time_equireact.size();
	if ((int)this->time_reactdeact.size() > length)
		length = (int)this->time_reactdeact.size();

	ofstream out_file(outname.data(), ios::out);
	out_file.precision(6);

	out_file << "Step   Flow   Transport  KinReactions  EquiReactions ReactDeact Total"
	         << "\n";
	for (i = 0; i < length; i++)
	{
		out_file << i << "  ";
		// flow
		if ((int)time_flow.size() > i)
			help = time_flow[i];
		else
			help = 0.0;
		out_file << help << "  ";
		tot_zeitschritt = help;
		// transport
		if ((int)time_transport.size() > i)
			help = time_transport[i];
		else
			help = 0.0;
		out_file << help << "  ";
		tot_zeitschritt += help;
		// kinetic reactions
		if ((int)time_kinreact.size() > i)
			help = time_kinreact[i];
		else
			help = 0.0;
		out_file << help << "  ";
		tot_zeitschritt += help;
		// equilibrium reactions
		if ((int)time_equireact.size() > i)
			help = time_equireact[i];
		else
			help = 0.0;
		out_file << help << "  ";

		// reaction deactivation
		if ((int)time_reactdeact.size() > i)
			help = time_reactdeact[i];
		else
			help = 0.0;
		out_file << help << "  ";

		tot_zeitschritt += help;

		out_file << tot_zeitschritt << "\n";
	}
	out_file << "\n";
	out_file << "TotalFlow   TotalTransport  TotalKinReactions  TotalEquiReactions  TotalReactDeact  Total"
	         << "\n";
	out_file << time_total_flow << "  " << time_total_transport << "  " << time_total_kinreact << "  "
	         << time_total_equireact << "  " << time_total_reactdeact << "  " << tot << " "
	         << "\n";
	out_file << time_total_flow / tot * 100 << "  " << time_total_transport / tot * 100 << "  "
	         << time_total_kinreact / tot * 100 << "  " << time_total_equireact / tot * 100 << "  "
	         << time_total_reactdeact / tot * 100 << "  " << tot / tot * 100 << "\n";
	out_file.close();
}

void CreateClockTime(void)
{
	CClockTime* m_ct = NULL;
	m_ct = new CClockTime();
	m_ct->delta_clocktime = 0.0;
	ClockTimeVec.push_back(m_ct);
}

void DestroyClockTime(void)
{
	for (size_t i = 0; i < ClockTimeVec.size(); i++)
		delete ClockTimeVec[i];
	ClockTimeVec.clear();
}
