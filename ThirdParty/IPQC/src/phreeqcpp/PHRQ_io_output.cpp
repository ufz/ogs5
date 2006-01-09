#include <assert.h>
#include "Phreeqc.h"
#include "phqalloc.h"

/* ---------------------------------------------------------------------- */
int Phreeqc::
warning_msg(const char *err_str)
/* ---------------------------------------------------------------------- */
{
	if (state == TRANSPORT && transport_warnings == FALSE)
		return (OK);
	if (state == ADVECTION && advection_warnings == FALSE)
		return (OK);
	count_warnings++;
	if (pr.warnings >= 0)
	{
		if (count_warnings > pr.warnings)
			return (OK);
	}
	if (phrq_io)
	{
		if (status_on)
		{
			phrq_io->screen_msg("\n");
		}
		std::ostringstream msg;
		msg << "WARNING: " << err_str;
		phrq_io->warning_msg(msg.str().c_str());
		status_on = false;
	}
	
	return OK;
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
echo_msg(const char *str)
/* ---------------------------------------------------------------------- */
{
	if (pr.echo_input == TRUE)
	{
		if (phrq_io) phrq_io->echo_msg(str);
	}
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
set_forward_output_to_log(int value)
/* ---------------------------------------------------------------------- */
{
	forward_output_to_log = value;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_forward_output_to_log(void)
/* ---------------------------------------------------------------------- */
{
	return forward_output_to_log;
}

void Phreeqc::
fpunchf_heading(const char *name)
{
	if (pr.punch == TRUE && current_selected_output != NULL)
	{
		punch_msg(name);
	}
}
void Phreeqc::
fpunchf(const char *name, const char *format, double d)
{
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, d);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}
void Phreeqc::
fpunchf(const char *name, const char *format, char * s)
{
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, s);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}
void Phreeqc::
fpunchf(const char *name, const char *format, int d)
{
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, d);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}

void Phreeqc::
fpunchf_user(int user_index, const char *format, double d)
{
	const char *name;
	
	if (current_user_punch == NULL)
		return;
	// check headings
	//if (user_index < user_punch_count_headings)
	int user_punch_count_headings = (int) current_user_punch->Get_headings().size();
	if (user_index < user_punch_count_headings)
	{
		//name = user_punch_headings[user_index];
		name = current_user_punch->Get_headings()[user_index].c_str();
	}
	else
	{
		if (fpunchf_user_s_warning == 0)
		{
			error_string = sformatf(
					"USER_PUNCH: Headings count does not match number of calls to PUNCH.\n");
			warning_msg(error_string);
			fpunchf_user_s_warning = 1;
		}
		sprintf(fpunchf_user_buffer, "no_heading_%d",
				(user_index - user_punch_count_headings) + 1);
		name = fpunchf_user_buffer;
	}
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, (double) d);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}

void Phreeqc::
fpunchf_user(int user_index, const char *format, char * d)
{
	const char *name;
	
	if (current_user_punch == NULL)
		return;
	int user_punch_count_headings = (int) current_user_punch->Get_headings().size();
	// check headings
	if (user_index < user_punch_count_headings)
	{
		//name = user_punch_headings[user_index];
		name = current_user_punch->Get_headings()[user_index].c_str();
	}
	else
	{
		if (fpunchf_user_s_warning == 0)
		{
			error_string = sformatf(
					"USER_PUNCH: Headings count does not match number of calls to PUNCH.\n");
			warning_msg(error_string);
			fpunchf_user_s_warning = 1;
		}
		sprintf(fpunchf_user_buffer, "no_heading_%d",
				(user_index - user_punch_count_headings) + 1);
		name = fpunchf_user_buffer;
	}
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, d);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}

int Phreeqc::
fpunchf_end_row(const char *format)
{
	if (phrq_io) 
	{
		phrq_io->fpunchf_end_row(format);
	}
	return OK;
}
#ifdef ERROR_OSTREAM
/* ---------------------------------------------------------------------- */
int Phreeqc::
process_file_names(int argc, char *argv[], std::istream **db_cookie,
				   std::istream **input_cookie, int log)
/* ---------------------------------------------------------------------- */
{
	int l;
	char token[2 * MAX_LENGTH], default_name[2 * MAX_LENGTH];
	char query[2 * MAX_LENGTH];
	char in_file[2 * MAX_LENGTH], out_file[2 * MAX_LENGTH], db_file[2 * MAX_LENGTH];
	char *env_ptr;
	char *ptr;
/*
 *   Prepare error handling
 */
	try {
		if (phrq_io == NULL) 
		{
			std::cerr << "No PHRQ_io output handler defined in process_file_names" << "\n";
		}
/*
 *   Prep for get_line
 */
		max_line = MAX_LINE;
		space((void **) ((void *) &line), INIT, &max_line, sizeof(char));
		space((void **) ((void *) &line_save), INIT, &max_line, sizeof(char));
/*
 *   Open error ostream
 */
		if (argc > 4)
		{
			if (!phrq_io->error_open(argv[4]))
			{
				error_string = sformatf( "Error opening file, %s.", argv[4]);
				warning_msg(error_string);
			}
		}
		else
		{
			phrq_io->error_open(NULL);
		}
/*
 *   Open user-input file
 */
		strcpy(query, "Name of input file?");
		std::ifstream * local_input_stream = NULL;
		if (argc <= 1)
		{
			default_name[0] = '\0';
			local_input_stream = open_input_stream(query, default_name, std::ios_base::in, false);
		}
		else
		{
			strcpy(default_name, argv[1]);
			local_input_stream = open_input_stream(query, default_name, std::ios_base::in, true);
		}
		screen_msg(sformatf("Input file: %s\n\n", default_name));
		strcpy(in_file, default_name);
/*
 *   Open file for output
 */
		strcpy(query, "Name of output file?");
		ptr = default_name;
		copy_token(token, &ptr, &l);
		strcpy(token, default_name);
		strcat(token, ".out");
		std::ofstream * local_output_stream;
		if (argc <= 1)
		{
			local_output_stream = open_output_stream(query, token, std::ios_base::out, false);
		}
		else if (argc == 2)
		{
			local_output_stream = open_output_stream(query, token, std::ios_base::out, true);
		}
		else if (argc >= 3)
		{
			strcpy(token, argv[2]);
			local_output_stream = open_output_stream(query, token, std::ios_base::out, true);
		}
		screen_msg(sformatf("Output file: %s\n\n", token));
		strcpy(out_file, token);
		phrq_io->Set_output_ostream(local_output_stream);
/*
 *   Open log file
 */
		if (log == TRUE)
		{
			if (!phrq_io->log_open("phreeqc.log"))
			{
				error_msg("Cannot open log file, phreeqc.log.", STOP);
			}
		}
/*
 *  Read input file for DATABASE keyword
 */
		if (local_input_stream->is_open())
		{
			phrq_io->push_istream(local_input_stream);
			if (get_line() == KEYWORD)
			{
				ptr = line;
				copy_token(token, &ptr, &l);
				if (strcmp_nocase(token, "database") == 0)
				{
#ifdef PHREEQ98
					user_database = string_duplicate(prefix_database_dir(ptr));
#else
					user_database = string_duplicate(ptr);
#endif
					if (string_trim(user_database) == EMPTY)
					{
						warning_msg("DATABASE file name is missing; default database will be used.");
						user_database = (char *) free_check_null(user_database);
					}
				}
			}
			phrq_io->pop_istream();
		}
		else
		{
			delete local_input_stream;
			error_string = sformatf( "Error opening file, %s.", in_file);
			error_msg(error_string, STOP);
		}
		
/*
 *   Open data base
 */
		strcpy(query, "Name of database file?");
		env_ptr = getenv("PHREEQC_DATABASE");
		if (user_database != NULL)
		{
			strcpy(token, user_database);
		}
		else if (env_ptr != NULL)
		{
			strcpy(token, env_ptr);
		}
		else
		{
			strcpy(token, default_data_base);
		}

		std::ifstream * local_database_file;
		if (argc <= 1)
		{
			local_database_file = open_input_stream(query, token, std::ios_base::in, false);
		}
		else if (argc < 4)
		{
			local_database_file = open_input_stream(query, token, std::ios_base::in, true);
		}
		else if (argc >= 4)
		{
			if (user_database == NULL)
			{
				strcpy(token, argv[3]);
			}
			else
			{
#ifndef PHREEQCI_GUI
				warning_msg	("Database file from DATABASE keyword is used; command line argument ignored.");
#endif
			}
			local_database_file = open_input_stream(query, token, std::ios_base::in, true);
		}
		local_database_file->close();
		delete local_database_file;
		screen_msg(sformatf("Database file: %s\n\n", token));
		strcpy(db_file, token);
#ifdef NPP
		output_msg(sformatf("Using PHREEQC: version 3.0.5, compiled on May 29, 2013\n"));
#endif
		output_msg(sformatf("   Input file: %s\n", in_file));
		output_msg(sformatf("  Output file: %s\n", out_file));
		output_msg(sformatf("Database file: %s\n\n", token));
#ifdef NPP
		output_flush();
#endif
		/*
		*   local cleanup
		*/
		user_database = (char *) free_check_null(user_database);
		line = (char *) free_check_null(line);
		line_save = (char *) free_check_null(line_save);

		*db_cookie = new std::ifstream(db_file, std::ios_base::in);
		*input_cookie = new std::ifstream(in_file, std::ios_base::in);
	}
	catch (PhreeqcStop e)
	{
		return get_input_errors();
	}
	return 0;
}
#else
/* ---------------------------------------------------------------------- */
int Phreeqc::
process_file_names(int argc, char *argv[], std::istream **db_cookie,
				   std::istream **input_cookie, int log)
/* ---------------------------------------------------------------------- */
{
	int l;
	char token[2 * MAX_LENGTH], default_name[2 * MAX_LENGTH];
	char query[2 * MAX_LENGTH];
	char in_file[2 * MAX_LENGTH], out_file[2 * MAX_LENGTH], db_file[2 * MAX_LENGTH];
	char *env_ptr;
	char *ptr;
/*
 *   Prepare error handling
 */
	try {
		if (phrq_io == NULL) 
		{
			std::cerr << "No PHRQ_io output handler defined in process_file_names" << "\n";
		}
/*
 *   Prep for get_line
 */
		max_line = MAX_LINE;
		space((void **) ((void *) &line), INIT, &max_line, sizeof(char));
		space((void **) ((void *) &line_save), INIT, &max_line, sizeof(char));
/*
 *   Open error ostream
 */
		if (argc > 4)
		{
			if (!phrq_io->error_open(argv[4]))
			{
				error_string = sformatf( "Error opening file, %s.", argv[4]);
				warning_msg(error_string);
			}
		}
		else
		{
			phrq_io->error_open(NULL);
		}
/*
 *   Open user-input file
 */
		strcpy(query, "Name of input file?");
		std::ifstream * local_input_stream = NULL;
		if (argc <= 1)
		{
			default_name[0] = '\0';
			local_input_stream = open_input_stream(query, default_name, std::ios_base::in, false);
		}
		else
		{
			strcpy(default_name, argv[1]);
			local_input_stream = open_input_stream(query, default_name, std::ios_base::in, true);
		}
		screen_msg(sformatf("Input file: %s\n\n", default_name));
		strcpy(in_file, default_name);
/*
 *   Open file for output
 */
		strcpy(query, "Name of output file?");
		ptr = default_name;
		copy_token(token, &ptr, &l);
		strcat(token, ".out");
		std::ofstream * local_output_stream;
		if (argc <= 1)
		{
			local_output_stream = open_output_stream(query, token, std::ios_base::out, false);
		}
		else if (argc == 2)
		{
			local_output_stream = open_output_stream(query, token, std::ios_base::out, true);
		}
		else if (argc >= 3)
		{
			strcpy(token, argv[2]);
			local_output_stream = open_output_stream(query, token, std::ios_base::out, true);
		}
		screen_msg(sformatf("Output file: %s\n\n", token));
		strcpy(out_file, token);
		phrq_io->Set_output_ostream(local_output_stream);
/*
 *   Open log file
 */
		if (log == TRUE)
		{
			if (!phrq_io->log_open("phreeqc.log"))
			{
				error_msg("Cannot open log file, phreeqc.log.", STOP);
			}
		}
/*
 *  Read input file for DATABASE keyword
 */
		if (local_input_stream->is_open())
		{
			phrq_io->push_istream(local_input_stream);
			if (get_line() == KEYWORD)
			{
				ptr = line;
				copy_token(token, &ptr, &l);
				if (strcmp_nocase(token, "database") == 0)
				{
#ifdef PHREEQ98
					user_database = string_duplicate(prefix_database_dir(ptr));
#else
					user_database = string_duplicate(ptr);
#endif
					if (string_trim(user_database) == EMPTY)
					{
						warning_msg("DATABASE file name is missing; default database will be used.");
						user_database = (char *) free_check_null(user_database);
					}
				}
			}
			phrq_io->pop_istream();
		}
		else
		{
			delete local_input_stream;
			error_string = sformatf( "Error opening file, %s.", in_file);
			error_msg(error_string, STOP);
		}
		
/*
 *   Open data base
 */
		strcpy(query, "Name of database file?");
		env_ptr = getenv("PHREEQC_DATABASE");
		if (user_database != NULL)
		{
			strcpy(token, user_database);
		}
		else if (env_ptr != NULL)
		{
			strcpy(token, env_ptr);
		}
		else
		{
			strcpy(token, default_data_base);
		}

		std::ifstream * local_database_file;
		if (argc <= 1)
		{
			local_database_file = open_input_stream(query, token, std::ios_base::in, false);
		}
		else if (argc < 4)
		{
			local_database_file = open_input_stream(query, token, std::ios_base::in, true);
		}
		else if (argc >= 4)
		{
			if (user_database == NULL)
			{
				strcpy(token, argv[3]);
			}
			else
			{
#ifndef PHREEQCI_GUI
				warning_msg	("Database file from DATABASE keyword is used; command line argument ignored.");
#endif
			}
			local_database_file = open_input_stream(query, token, std::ios_base::in, true);
		}
		local_database_file->close();
		delete local_database_file;
		screen_msg(sformatf("Database file: %s\n\n", token));
		strcpy(db_file, token);

		output_msg(sformatf("   Input file: %s\n", in_file));
		output_msg(sformatf("  Output file: %s\n", out_file));
		output_msg(sformatf("Database file: %s\n\n", token));
		/*
		*   local cleanup
		*/
		user_database = (char *) free_check_null(user_database);
		line = (char *) free_check_null(line);
		line_save = (char *) free_check_null(line_save);

		*db_cookie = new std::ifstream(db_file, std::ios_base::in);
		*input_cookie = new std::ifstream(in_file, std::ios_base::in);
	}
	catch (PhreeqcStop e)
	{
		return get_input_errors();
	}
	return 0;
}
#endif
/* ---------------------------------------------------------------------- */
std::ifstream * Phreeqc::
open_input_stream(char *query, char *default_name, std::ios_base::openmode mode, bool batch)
/* ---------------------------------------------------------------------- */
{
	char name[MAX_LENGTH];
	std::ifstream *new_stream;
	int l;
#ifdef ERROR_OSTREAM
	std::ostream * error_ostream_save = phrq_io->Get_error_ostream();
#else
	FILE * error_file_save = phrq_io->Get_error_file();
#endif

	for (;;)
	{
/*
 *   Get file name
 */
		strcpy(name, default_name);
		if (!batch )
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			screen_msg(sformatf("%s\n", query));
			if (default_name[0] != '\0')
			{
				screen_msg(sformatf("Default: %s\n", default_name));
			}
			fgets(name, MAX_LENGTH, stdin);
			l = (int) strlen(name);
			name[l - 1] = '\0';
			if (name[0] == '\0')
			{
				strcpy(name, default_name);
			}
		}
/*
 *   Open existing file to read
 */
		new_stream = new std::ifstream(name, mode);
		if (new_stream == NULL || !new_stream->is_open())
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			error_string = sformatf( "\nERROR: Cannot open file, %s.\n", name);
			screen_msg(error_string);
#ifdef NPP
			error_msg(sformatf( "\nERROR: Cannot open file, %s.\n       Please check, and give the correct, full path + name.\n", name), STOP);
			break;
#endif
			error_flush();
			batch = FALSE;
			continue;		
		}
		break;
	}
	strncpy(default_name, name, MAX_LENGTH);
	if (!batch )
	{
		//phrq_io->Set_error_ostream(error_file_save);
#ifdef ERROR_OSTREAM
		phrq_io->Set_error_ostream(error_ostream_save);
#else
		phrq_io->Set_error_file(error_file_save);
#endif
	}
	return (new_stream);
}
/* ---------------------------------------------------------------------- */
std::ofstream * Phreeqc::
open_output_stream(char *query, char *default_name, std::ios_base::openmode mode, bool batch)
/* ---------------------------------------------------------------------- */
{
	char name[MAX_LENGTH];
	std::ofstream *new_stream;
	int l;
#ifdef ERROR_OSTREAM
	std::ostream * error_ostream_save = phrq_io->Get_error_ostream();
#else
	FILE * error_file_save = phrq_io->Get_error_file();
#endif
	
	for (;;)
	{
/*
 *   Get file name
 */
		strcpy(name, default_name);
		if (!batch )
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			
			screen_msg(sformatf("%s\n", query));
			if (default_name[0] != '\0')
			{
				screen_msg(sformatf("Default: %s\n", default_name));
			}
			fgets(name, MAX_LENGTH, stdin);
			l = (int) strlen(name);
			name[l - 1] = '\0';
			if (name[0] == '\0')
			{
				strcpy(name, default_name);
			}
		}
/*
 *   Open existing file to read
 */
		new_stream = new std::ofstream(name, mode);
		if (new_stream == NULL || !new_stream->is_open())
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			error_string = sformatf( "\nERROR: Cannot open file, %s.\n", name);
			screen_msg(error_string);
			error_flush();
			batch = FALSE;
			continue;		
		}
		break;
	}
	strncpy(default_name, name, MAX_LENGTH);
	if (!batch )
	{
#ifdef ERROR_OSTREAM
		phrq_io->Set_error_ostream(error_ostream_save);
#else
		phrq_io->Set_error_file(error_file_save);
#endif
	}
	return (new_stream);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
std::ofstream * Phreeqc::
open_output_file(char *query, char *default_name, std::ios_base::openmode mode, bool batch)
/* ---------------------------------------------------------------------- */
{
	char name[MAX_LENGTH];
	std::ofstream *new_stream;
	int l;
#ifdef ERROR_OSTREAM
		std::ostream * error_ostream_save = phrq_io->Get_error_ostream();
#else
		FILE * error_file_save = phrq_io->Get_error_file();
#endif
	

	for (;;)
	{
/*
 *   Get file name
 */
		strcpy(name, default_name);
		if (!batch )
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			screen_msg(sformatf("%s\n", query));
			if (default_name[0] != '\0')
			{
				screen_msg(sformatf("Default: %s\n", default_name));
			}
			fgets(name, MAX_LENGTH, stdin);
			l = (int) strlen(name);
			name[l - 1] = '\0';
			if (name[0] == '\0')
			{
				strcpy(name, default_name);
			}
		}
/*
 *   Open existing file to read
 */
		new_stream = new std::ofstream(name, mode);
		if (new_stream == NULL || !new_stream->is_open())
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			
			error_string = sformatf( "\nERROR: Cannot open file, %s.\n", name);
			screen_msg(error_string);
			error_flush();
			batch = FALSE;
			continue;		
		}
		break;
	}
	strncpy(default_name, name, MAX_LENGTH);
	if (!batch )
	{
#ifdef ERROR_OSTREAM
		phrq_io->Set_error_ostream(error_ostream_save);
#else
		phrq_io->Set_error_file(error_file_save);
#endif
	}
	return (new_stream);
}
#endif
/* ---------------------------------------------------------------------- */
void Phreeqc::
screen_msg(const char *err_str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) phrq_io->screen_msg(err_str);
}
// ---------------------------------------------------------------------- */
// dump file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
dump_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->dump_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_close();
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_msg(str);
}
// ---------------------------------------------------------------------- */
// error file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
error_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->error_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
error_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->error_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
error_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->error_close();
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
error_msg(const char *err_str, bool stop)
/* ---------------------------------------------------------------------- */
{
	if (get_input_errors() <= 0)
		input_error = 1;
	if (phrq_io)
	{
		std::ostringstream msg;
		msg << "ERROR: " << err_str << "\n";

		phrq_io->output_msg(msg.str().c_str());
		phrq_io->log_msg(msg.str().c_str());

		if (status_on)
		{
			phrq_io->screen_msg("\n");
		}
		status_on = false;
		phrq_io->error_msg(msg.str().c_str(), stop);
	}

	if (stop)
	{
		throw PhreeqcStop();
	}
}
// ---------------------------------------------------------------------- */
// log file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
log_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->log_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
log_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
log_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_close();
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
log_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_msg(str);
}
// ---------------------------------------------------------------------- */
// output_temp file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
output_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) 
		return this->phrq_io->output_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->output_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->output_close();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
	{
		if (get_forward_output_to_log())
		{
			phrq_io->log_msg(str);
		}
		else
		{
			phrq_io->output_msg(str);
		}
	}
}
// ---------------------------------------------------------------------- */
// punch file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
punch_open(const char *file_name, int n_user)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->punch_open(file_name, std::ios_base::out, n_user);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_close();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_msg(str);
}
