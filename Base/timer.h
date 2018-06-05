/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   ROCKFLOW - Modul: timer.h

   Aufgabe:
   Funktionen zur Laufzeitermittlung im Testbetrieb.
**************************************************************************/

#ifndef timer_INC

/* Schutz gegen mehrfaches Einfuegen */
#define timer_INC

/* Andere oeffentlich benutzte Module */
#include <time.h>
#include <vector>
//#include <windows.h>

/* Deklarationen */

/* Setzt Timer auf 0, aber startet noch nicht */
extern void TInitTimer(int speicher);

/* Setzt Timer auf 0 und startet ihn */
extern void TStartTimer(int speicher);

/* Liefert Laufzeit des gew. Timers seit TStartTimer als Long */
extern long TGetTimer(int speicher);

/* Liefert Laufzeit des gew. Timers seit TStartTimer als Double */
extern double TGetTimerDouble(int speicher);

/*  Haelt Timer an */
extern void TStopTimer(int speicher);

/*  Laesst Timer weiterlaufen */
extern void TRestartTimer(int speicher);

/* Liefert Laufzeit/CPU-Zeit in "Ticks" */
extern long TGetTime(void);

/* Liefert Laenge eines "Ticks" */
extern long TGetTicksPerSecond(void);

/* Zerstoert alle Zeitspeicher */
extern void TDestroyTimers(void);

/* Weitere externe Objekte */

/* Interface fuer AMG-Loeser */
extern void ctime_(float* time);

// SB_time

class CClockTime
{
private:
public:
	CClockTime(void);
	~CClockTime(void);
	std::vector<double> time_flow;
	std::vector<double> time_transport;
	std::vector<double> time_kinreact;
	std::vector<double> time_equireact;
	std::vector<double> time_reactdeact;
	double time_total_flow;
	double time_total_transport;
	double time_total_kinreact;
	double time_total_equireact;
	double time_total_reactdeact;
	double delta_clocktime;
	clock_t start;
	clock_t end;

	void StopTime(const std::string& name = "");
	void StartTime(void);
	void PrintTimes(void);

	long time1;
	long time2;
	double difftime;
};
extern void CreateClockTime(void);
extern void DestroyClockTime(void);
extern std::vector<CClockTime*> ClockTimeVec;
#endif
