#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

#include "cir_parser.h"
#include "../hashtable/hashtable.h"
#include "../spicy.h"
#include "../lists/lists.h"


// converts a string to uppercase string
void strtoupper(char *str) {
	int i;

	for (i=0; i < (int)strlen(str); i++) {
		str[i] = toupper(str[i]);
	}
}


// returns 1 if character c exists in components string
byte component_type_is_valid(char c) {
	if (strchr("VIRCLDMQ", toupper(c)) == NULL)
		return 0;
	else
		return 1;
}


// Splits plot commands with multiple nodes into multiple commands
// for example .PLOT V(1) V(2) will be split into .PLOT V(1) and .PLOT V(2)
void parse_plot(char *command){					/* .PLOT V(1) V(2) -> .PLOT V(1) / .PLOT V(2)  */

	const char delim[5] = " \r\t\n";
	char *token = NULL;
	char *plot_command = NULL;
	unsigned int size_cmd;

	token = strtok(command, delim);       // First token is the command .PLOT/.PRINT
	plot_command = strdup(token);
	size_cmd = strlen(token);

	token = strtok(NULL,delim);           // Parse each node for the current command
	while ( token!= NULL){

		plot_command = (char *)realloc(plot_command,( size_cmd + strlen(token) + 3)*sizeof(char));

		sprintf(&plot_command[size_cmd]," %s%c",token,'\0');
		add_command_to_list(plot_command);

		token = strtok(NULL,delim);
	}

	free(plot_command);
	plot_command = NULL;
}


// parses the command and stores it into the command list
// Unknown commands are bypassed. Prints a Warning message in that case
void parse_command(char *command) {
	char *str_ptr;

	strtoupper(command);

	// each one of those commands are supposed to have arguments
	if ( (strncmp(command, ".DC ", 4) == 0) || \
		 (strncmp(command, ".OPTIONS ", 9) == 0) || \
		 (strcmp(command, ".OPTIONS\n") == 0) || \
		 (strcmp(command, ".OPTIONS\t") == 0) || \
		 (strcmp(command, ".OPTIONS") == 0) || \
		 (strncmp(command, ".PRINT ", 7) == 0) || \
		 (strncmp(command, ".PLOT ", 6) == 0)  || \
		 (strncmp(command, ".TRAN ", 6) == 0) || \
		 (strncmp(command, ".DC ", 4) == 0) ) {

		// erase possible \n at the end of the command (probably not necessary)
		str_ptr = strchr(command, '\n');
		if (str_ptr != NULL)
			*str_ptr = '\0';

		// adds the command to list
		// Additional checks for the commands will be performed during command execution
		if((strncmp(command, ".PRINT ", 7) == 0) || \
		 (strncmp(command, ".PLOT ", 6) == 0)){
		 	parse_plot(command);

		} else {
			add_command_to_list(command);
		}
	}
	else {
		printf(YEL "Warning" NRM ": Unknown command "
			  "or command with no arguments. Ignoring %s\n", command);
	}

}


// return 1 if str can be converted to a double, otherwise 0
byte parse_double(double *d, char *str) {
	char *end_ptr = str;

	errno = 0;	// errno.h's errno

	*d = strtod(str, &end_ptr);

	// overflow or underflow
	if (errno != 0)
		return 0;

	// no convertion applied
	if (str == end_ptr)
		return 0;

	// trailing characters
	if (*end_ptr != 0)
		return 0;

	return 1;
}



// prints parsed information
void print_parsed_information(char type, char *name, char *node1_name, \
							  char *node2_name, char *node3_name, char *node4_name, \
							  double val, byte has_G2, char *model_name, double l, double w) {

	printf("---\n" BLU "*PARSED DATA*\n" NRM\
		   "<type> \t\t%c\n<name>  \t%s\n<node1_name> \t%s\n<node2_name> \t%s\n"
		   "<node3_name> \t%s\n<node4_name> \t%s\n<val> \t\t%lf\n<has_G2> \t%u\n"
		   "<model_name> \t%s\n<l_val> \t%lf\n<w_val> \t%lf\n" \
		   BLU "*END PARSED DATA*" NRM "\n---\n", \
								type, (name==NULL)?"(null)":name, \
								(node1_name==NULL)?"(null)":node1_name, \
								(node2_name==NULL)?"(null)":node2_name, \
								(node3_name==NULL)?"(null)":node3_name, \
								(node4_name==NULL)?"(null)":node4_name, val, has_G2, \
								(model_name==NULL)?"(null)":model_name, l, w);
}



// returns 1 when line has only whitespace characters, otherwise 0
int whitespaces_only(char *line) {
	unsigned int i;

	// stops when something other than whitespace is found
	for(i=0; (i < (unsigned int)strlen(line)) && (isspace(line[i]) != 0); i++);

	if (i != (unsigned int)strlen(line))
		return 0;	// line does not contain only whitespaces
	else
		return 1;	// line contains only whitespaces
}



unsigned long get_components_num(char *filename) {
	FILE *fp = NULL;
	char *line = NULL;
	size_t len = 0;
	unsigned long components_num = 0;

	fp = fopen(filename, "r");
	if (fp == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}

	while(getline(&line, &len, fp) != -1) {
		if (component_type_is_valid(line[0]))
				components_num++;
	}

	free(line);
	fclose(fp);
	return components_num;
}



// parse the Circuit file
void parse_cir(char *filename) {

	FILE *fp = NULL;

	byte has_G2 = 0;			// "G2" Optional Variable of I, R, C
	char type = 'X';			// Component Type: V, I, R, C, L, D, M, Q (unknown = X)
	char *name = NULL;			// Component Name: String right after type
	char *node1_name = NULL;	// for node types: + or D (MOS) or C (BJT)
	char *node2_name = NULL;	// for node types: - or G (MOS) or B (BJT)
	char *node3_name = NULL;	// for node types:		S (MOS) or E (BJT)
	char *node4_name = NULL;	// for node type :		B (MOS)

	char *model_name = NULL;	// model name of D, M, Q
	double val = -1;			// holds optional area value for D and Q (-1 if not existent)
	double l = -1;				// holds L value of MOS transistor
	double w = -1;				// holds W value of MOS transistor


	// hashtable entry field variables
	unsigned long id = 0;
	element_h *node1 = NULL;
	element_h *node2 = NULL;
	element_h *node3 = NULL;
	element_h *node4 = NULL;

	// used for parsing
	char *line = NULL;		// the line parsed
	char *line_pos;			// used to point at a line's byte
	size_t len = 0;			// the length of the line parsed


	// variables used with strtok
	const char delim[5] = " \r\t\n";
	char *token = NULL;
	unsigned int tok_count = 0;
	unsigned int min_tok_count = 0;
	byte rest_line_commented = 0;

	// offset to pass possible non ascci invisible characters from line
	unsigned int line_offset = 0;

	// variables used for transient and AC analysis (we have THAT many for optical reasons)
	int tr_type;
	int ac_usability;
	void *tran_spec_data;
	void *ac_spec;
	double v1, v2, v3, v4, v5, v6, v7;
	double mag, phase;
	double *times, *values;
	unsigned int total_tuples;

	int i;

	fp = fopen(filename, "r");
	if (fp == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}


	// read the file line-by-line
	while((getline(&line, &len, fp)) != -1) {

		// bypass possible non ascii characters (for some reason they exist in a benchmark)
		line_offset = 0;
		while ((line[line_offset] < 0) || (line[line_offset] > 127)) { line_offset++; }


		// bypass lines containing only whitespaces
		if (whitespaces_only(&line[line_offset]))
			continue;

		// bypass comment lines
		if ((line[line_offset] == '%') || (line[line_offset] == '*'))
			continue;


		/*#ifdef DEBUG*/
		printf(GRN "<line read> %s" NRM, &line[line_offset]);
		/*#endif*/



		/* ************** *
		 * PARSE THE LINE *
		 * ************** */

		/* re-initialise variables and parse line */
		has_G2 = 0;
		type = 'X';
		name = NULL;
		node1_name = NULL;
		node2_name = NULL;
		node3_name = NULL;
		node4_name = NULL;
		model_name = NULL;
		// Transient
		tr_type = TR_TYPE_NONE;
		tran_spec_data = NULL;
		ac_spec = NULL;
		v1 = 0; v2 = 0; v3 = 0;
		v4 = 0; v5 = 0; v6 = 0; v7 = 0;
		times = NULL; values = NULL;
		total_tuples = 0;
		// AC
		mag = 0; phase = 0;
		ac_usability = 0;

		// -1 if not present or part of the component parsed
		val = -1;
		l = -1;
		w = -1;

		rest_line_commented = 0;
		min_tok_count = 0;
		tok_count = 0;




		// check if the first line character is '.'
		if (line[line_offset] == '.') {
			parse_command(&line[line_offset]);
			continue; // nothing more for this line
		}


		// not a command. parse the component information

		token = strtok(&line[line_offset], delim);

		while (token != NULL) {
			/*printf("token = %s\n", token);*/
			// a comment might start with an '*' or a '%'
			line_pos = strpbrk(token, "%*");
			if (line_pos != NULL) {

				// the line comment starts at the beggining of the token
				if (*line_pos == token[0])
					break;

				// comment starts somewhere in the token
				*line_pos = '\0';
				rest_line_commented = 1;
			}


			/* ************************ *
			 * TOKEN PROCESSING SECTION *
			 * ************************ */

			tok_count++;
			/*printf("Token #%u: %s\n", tok_count, token);*/

			if (tok_count == 1) {

				type = token[0];
				printf("type = %c\n", type);
				if (!component_type_is_valid(type)) {
					printf("Syntax Error: Invalid component type (%s)\n", token);
					exit(EXIT_FAILURE);
				}
				name = strdup(&token[1]);
				if (name == NULL) {
					printf("Error. Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}

				if (strchr("VIRCLD", toupper(type)) != NULL)
					min_tok_count = 4;
				else if (toupper(type) == 'Q')
					min_tok_count = 5;
				else // if (toupper(type) == 'M')
					min_tok_count = 8;

			}
			else if (tok_count == 2) {

				node1_name = strdup(token);
				if (node1_name == NULL) {
					printf("Error. Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}
			}
			else if (tok_count == 3) {

				node2_name = strdup(token);
				if (node2_name == NULL) {
					printf("Error. Memory allocation problems. Exitin..\n");
					exit(EXIT_FAILURE);
				}
			}
			else if (tok_count == 4) {

				if (toupper(type) == 'D') {
					model_name = strdup(token);
					if (model_name == NULL) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}
				}
				else if ((toupper(type) == 'M') || (toupper(type) == 'Q')) {
					node3_name = strdup(token);
					if (node3_name == NULL) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}
				}
				else {	// V, I, R, C, L
					if (parse_double(&val, token) == 0) {
						printf("Syntax error. Value (%s) cannot be converted to double\n", token);
						exit(EXIT_FAILURE);
					}
				}
			}
			else if (tok_count == 5) {

				if (((toupper(type) == 'I') || (toupper(type) == 'R') || (toupper(type) == 'C')) &&
				   ((strcmp(token, "G2") == 0) || (strcmp(token, "g2") == 0))) {
					has_G2 = 1;
				}
				else if ((toupper(type) == 'D')) {
					if (parse_double(&val, token) == 0) {
						printf("Syntax error. Value (%s) cannot be converted to double\n", token);
						exit(EXIT_FAILURE);
					}
				}
				else if ((toupper(type) == 'M')) {
					node4_name = strdup(token);
					if (node4_name == NULL) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}
				}
				else if ((toupper(type) == 'Q')) {
					model_name = strdup(token);
					if (model_name == NULL) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}
				}
				else if ((toupper(type) == 'I') || (toupper(type) == 'V')) { // transient analysis

				    /* ********************************** *
					 * START OF TRANSIENT/AC SPEC PARSING *
					 * ********************************** */

					if ((strstr("EXP SIN PULSE PWL", token) != NULL) \
					 || (strstr("exp sin pulse pwl", token) != NULL)) {
						printf("Transient Spec Found\n");

						// The second letter of the words EXP, SIN, PULSE and PWL is unique
						// Take advantage of this for some easy and lazy checks
						if (toupper(token[1] == 'W')) { // PWL
							printf("PWL function\n");
							tr_type = TR_TYPE_PWL;

							// iterate through tuples
							total_tuples = 0;
							while (1) {

								// read time
								token = strtok(NULL, delim);

								// no more tuples
								if (token == NULL)
									break;

								// temporarily store time into v1
								if (parse_double(&v1, &token[1]) == 0) {
									printf("Error. Value (%s) cannot be converted to double\n", token);
									free(times);
									free(values);
									exit(EXIT_FAILURE);
								}

								// each time should be bigger than the previous one
								if (total_tuples > 0) {
									if (v1 <= times[total_tuples-1]) {
										printf("Error. Tuple time shouldn't be smaller "
												"than previous time\n");
										free(times);
										free(values);
										exit(EXIT_FAILURE);
									}
								}

								// read value
								token = strtok(NULL, delim);
								if (token == NULL) {
									printf("Error. Incomplete PWL tuple.\n");
									free(times);
									free(values);
									exit(EXIT_FAILURE);
								}


								// check for ')'
								if (strchr(token, ')') == NULL) {
									printf("Error. Parentheses must close directly after value (%s)\n",\
									token);
									free(times);
									free(values);
									exit(EXIT_FAILURE);
								}

								token[strlen(token)-1] = '\0';  // get rid of ')'
								if (parse_double(&v2, token) == 0) {
									printf("Error. Value (%s) cannot be converted to double\n", token);
									free(times);
									free(values);
									exit(EXIT_FAILURE);
								}
								token[strlen(token)] = ')';     // recover previously deleted ')'


								// store the (time value) tuple
								total_tuples++;
								times  = (double *) realloc(times, total_tuples*sizeof(double));
								values = (double *) realloc(values, total_tuples*sizeof(double));
								if ((times == NULL) || (values == NULL)) {
									printf("Error. Memory allocation problems. Exiting..\n");
									exit(EXIT_FAILURE);
								}

								times[total_tuples-1] = v1;
								values[total_tuples-1] = v2;

							}

							printf("tuples (%d): (times, values) = ", total_tuples);
							for (i = 0; i < total_tuples; i++) {
								printf("(%lf, %lf) ", times[i], values[i]);
							}
							printf("\n");

						}
						else if (strchr("XIU", toupper(token[1])) != NULL) { // eXp sIn pUlse

							if (toupper(token[1]) == 'X')
								tr_type = TR_TYPE_EXP;
							else if (toupper(token[1]) == 'I')
								tr_type = TR_TYPE_SIN;
							else
								tr_type = TR_TYPE_PULSE;

							// parse i1
							token = strtok(NULL, delim);
							if (parse_double(&v1, &token[1]) == 0) {
								printf("Syntax error. Value (%s) cannot be converted to double\n", token);
								exit(EXIT_FAILURE);
							}

							// parse i2 (for exp and pulse) or ia (for sin)
							token = strtok(NULL, delim);
							if (parse_double(&v2, token) == 0) {
								printf("Syntax error. Value (%s) cannot be converted to double\n", token);
								exit(EXIT_FAILURE);
							}

							// parse td1 (exp) or fr (sin) or td (pulse)
							token = strtok(NULL, delim);
							if (parse_double(&v3, token) == 0) {
								printf("Syntax error. Value (%s) cannot be converted to double\n", token);
								exit(EXIT_FAILURE);
							}

							// parse tc1 (exp) or td (sin) or tr (pulse)
							token = strtok(NULL, delim);
							if (parse_double(&v4, token) == 0) {
								printf("Syntax error. Value (%s) cannot be converted to double\n", token);
								exit(EXIT_FAILURE);
							}

							// parse td2 (exp) or df (sin) or tf (pulse)
							token = strtok(NULL, delim);
							if (parse_double(&v5, token) == 0) {
								printf("Syntax error. Value (%s) cannot be converted to double\n", token);
								exit(EXIT_FAILURE);
							}

							// parse the last arguments
							if (tr_type == TR_TYPE_EXP) {  // exp function
								printf("EXP function\n");

								// parse tc2 (this is the 6th and last argument of exp function)
								token = strtok(NULL, delim);

								// check for ')'
								if (strchr(token, ')') == NULL) {
									printf("Error. Parenthesis must close directly after value (%s)\n",\
									token);
									exit(EXIT_FAILURE);
								}

								token[strlen(token)-1] = '\0'; // get rid of ')'
								if (parse_double(&v6, token) == 0) {
									printf("Error. Value (%s) cannot be converted to double\n", token);
									exit(EXIT_FAILURE);
								}
								token[strlen(token)] = ')';    // recover previously deleted ')'


								// check if there is another argument (syntax error)
								token = strtok(NULL, delim);
								if (token != NULL) {
									printf("Error. Too many fields in transient function\n");
									exit(EXIT_FAILURE);
								}
							}
							else if (tr_type == TR_TYPE_SIN) {
								printf("SIN function\n");

								// parse tc2 (this is the 6th and last argument of sin function)
								token = strtok(NULL, delim);

								// check for ')'
								if (strchr(token, ')') == NULL) {
									printf("Error. Parenthesis must close directly after value (%s)\n",\
									token);
									exit(EXIT_FAILURE);
								}

								token[strlen(token)-1] = '\0'; // get rid of ')'
								if (parse_double(&v6, token) == 0) {
									printf("Error. Value (%s) cannot be converted to double\n", token);
									exit(EXIT_FAILURE);
								}
								token[strlen(token)] = ')';    // recover previously deleted ')'


								// check if there is another argument (syntax error)
								token = strtok(NULL, delim);
								if (token != NULL) {
									printf("Error. Too many fields in transient function\n");
									exit(EXIT_FAILURE);
								}

							}
							else  { // if (toupper(token[1]) == 'U')
								printf("PULSE function\n");

								// parse pw (this is the 6th argument of pulse function)
								token = strtok(NULL, delim);
								if (parse_double(&v6, token) == 0) {
									printf("Error. Value (%s) cannot be converted to double\n", token);
									exit(EXIT_FAILURE);
								}

								// parse per (this is the 7th and last argument of pulse function)
								token = strtok(NULL, delim);

								// check for ')'
								if (strchr(token, ')') == NULL) {
									printf("Error. Parenthesis must close directly after value (%s)\n",\
									token);
									exit(EXIT_FAILURE);
								}

								token[strlen(token)-1] = '\0';	// get rid of ')'
								if (parse_double(&v7, token) == 0) {
									printf("Error. Value (%s) cannot be converted to double\n", token);
									exit(EXIT_FAILURE);
								}
								token[strlen(token)] = ')';		// recover previously deleted ')'


								// check if there is another argument (syntax error)
								token = strtok(NULL, delim);
								if (token != NULL) {
									printf("Error. Too many fields in transient function\n");
									exit(EXIT_FAILURE);
								}

								// check if the period of the pulse if not big enough
								// td + per must be bigger than td + tr + tf + pw
								if ((v3 + v7) < (v3 + v4 + v5 + v6)) {
									printf(RED "Error: " NRM "Pulse period not big enough\n");
									exit(EXIT_FAILURE);
								}
							}


							// print the parsed arguments
							printf("arg1=%lf, arg2=%lf, arg3=%lf, arg4=%lf, "
								   "arg5=%lf, arg6=%lf, arg7=%lf\n", \
									v1, v2, v3, v4, v5, v6, v7);

						}
						else {
							printf("Error. Unknown transient function (%s)\n", token);
							exit(EXIT_FAILURE);
						}

						break;
					}
					else if ((strncmp("AC", token, 2) == 0) || (strncmp("ac", token, 2) == 0)) {
						printf("AC Spec Found\n");

						// parse mag
						token = strtok(NULL, delim);
						if (parse_double(&mag, token) == 0) {
							printf("Syntax error. Value (%s) cannot be converted to double\n", token);
							exit(EXIT_FAILURE);
						}

						// parse phase
						token = strtok(NULL, delim);
						if (parse_double(&phase, token) == 0) {
							printf("Syntax error. Value (%s) cannot be converted to double\n", token);
							exit(EXIT_FAILURE);
						}

						// check if there is another argument (syntax error)
						token = strtok(NULL, delim);
						if (token != NULL) {
							printf("Error. Too many fields in AC analysys\n");
							exit(EXIT_FAILURE);
						}

						ac_usability = 1;
					}
					else {
						printf("syntax error. Type '%c' has unknown fifth field (%s)\n", type, token);
						exit(EXIT_FAILURE);
					}

					/* ******************************** *
					 * END OF TRANSIENT/AC SPEC PARSING *
					 * ******************************** */

				}
				else {
					printf("Syntax error. Type '%c' has unknown fifth field (%s)\n", type, token);
					exit(EXIT_FAILURE);
				}

			}
			else if (tok_count == 6) {

				if (toupper(type) == 'M') {
					model_name = strdup(token);
					if (model_name == NULL) {
						printf("Error. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}
				}
				else if (toupper(type) == 'Q') {
					if (parse_double(&val, token) == 0) {
						printf("Syntax error. Value (%s) cannot be converted to double\n", token);
						exit(EXIT_FAILURE);
					}
				}
				else if (toupper(type) == 'I') {  // transient analysis (I had G2 field)

					/* ******************************* *
					 * START OF TRANSIENT SPEC PARSING *
					 * ******************************* */

					// is G2 was found at 'I' component than transient spec might start
					// at the sixth field/token instead of the fifth
					if (has_G2 == 1) {
						if ((strstr("EXP PULSE SIN PWL", token) != NULL) \
						 || (strstr("exp pulse sin pwl", token) != NULL)) {
							printf("Transient Spec Found\n");

							// The second letter of the words EXP, SIN, PULSE and PWL is unique
							// Take advantage of this for some easy and lazy checks
							if (toupper(token[1] == 'W')) { // PWL
								printf("PWL function\n");
								tr_type = TR_TYPE_PWL;

								// iterate through tuples
								total_tuples = 0;
								while (1) {

									// read time
									token = strtok(NULL, delim);

									// no more tuples
									if (token == NULL)
										break;

									// temporarily store time into v1
									if (parse_double(&v1, &token[1]) == 0) {
										printf("Error. Value (%s) cannot be converted to double\n", token);
										free(times);
										free(values);
										exit(EXIT_FAILURE);
									}

									// each time should be bigger than the previous one
									if (total_tuples > 0) {
										if (v1 <= times[total_tuples-1]) {
											printf("Error. Tuple time shouldn't be smaller "
													"than previous time\n");
											free(times);
											free(values);
											exit(EXIT_FAILURE);
										}
									}

									// read value
									token = strtok(NULL, delim);
									if (token == NULL) {
										printf("Error. Incomplete PWL tuple.\n");
										free(times);
										free(values);
										exit(EXIT_FAILURE);
									}


									// check for ')'
									if (strchr(token, ')') == NULL) {
										printf("Error. Parentheses must close directly after value (%s)\n",\
										token);
										free(times);
										free(values);
										exit(EXIT_FAILURE);
									}

									token[strlen(token)-1] = '\0';  // get rid of ')'
									if (parse_double(&v2, token) == 0) {
										printf("Error. Value (%s) cannot be converted to double\n", token);
										free(times);
										free(values);
										exit(EXIT_FAILURE);
									}
									token[strlen(token)] = ')';     // recover previously deleted ')'


									// store the (time value) tuple
									total_tuples++;
									times  = (double *) realloc(times, total_tuples*sizeof(double));
									values = (double *) realloc(values, total_tuples*sizeof(double));
									if ((times == NULL) || (values == NULL)) {
										printf("Error. Memory allocation problems. Exiting..\n");
										exit(EXIT_FAILURE);
									}

									times[total_tuples-1] = v1;
									values[total_tuples-1] = v2;

								}


								if (total_tuples == 0) {
									printf(RED "Error: " NRM "No tuples given in PWL\n");
									exit(EXIT_FAILURE);
								}


								printf("tuples (%d): (times, values) = ", total_tuples);
								for (i = 0; i < total_tuples; i++) {
									printf("(%lf, %lf) ", times[i], values[i]);
								}
								printf("\n");

							}
							else if (strchr("XIU", toupper(token[1])) != NULL) { // eXp sIn pUlse

								if (toupper(token[1]) == 'X')
									tr_type = TR_TYPE_EXP;
								else if (toupper(token[1]) == 'I')
									tr_type = TR_TYPE_SIN;
								else
									tr_type = TR_TYPE_PULSE;

								// parse i1
								token = strtok(NULL, delim);
								if (parse_double(&v1, &token[1]) == 0) {
									printf("Syntax error. Value (%s) cannot be converted to double\n", token);
									exit(EXIT_FAILURE);
								}

								// parse i2 (for exp and pulse) or ia (for sin)
								token = strtok(NULL, delim);
								if (parse_double(&v2, token) == 0) {
									printf("Syntax error. Value (%s) cannot be converted to double\n", token);
									exit(EXIT_FAILURE);
								}

								// parse td1 (exp) or fr (sin) or td (pulse)
								token = strtok(NULL, delim);
								if (parse_double(&v3, token) == 0) {
									printf("Syntax error. Value (%s) cannot be converted to double\n", token);
									exit(EXIT_FAILURE);
								}

								// parse tc1 (exp) or td (sin) or tr (pulse)
								token = strtok(NULL, delim);
								if (parse_double(&v4, token) == 0) {
									printf("Syntax error. Value (%s) cannot be converted to double\n", token);
									exit(EXIT_FAILURE);
								}

								// parse td2 (exp) or df (sin) or tf (pulse)
								token = strtok(NULL, delim);
								if (parse_double(&v5, token) == 0) {
									printf("Syntax error. Value (%s) cannot be converted to double\n", token);
									exit(EXIT_FAILURE);
								}

								// parse the last arguments
								if (tr_type == TR_TYPE_EXP) {  // exp function
									printf("EXP function\n");

									// parse tc2 (this is the 6th and last argument of exp function)
									token = strtok(NULL, delim);

									// check for ')'
									if (strchr(token, ')') == NULL) {
										printf("Error. Parenthesis must close directly after value (%s)\n",\
										token);
										exit(EXIT_FAILURE);
									}

									token[strlen(token)-1] = '\0'; // get rid of ')'
									if (parse_double(&v6, token) == 0) {
										printf("Error. Value (%s) cannot be converted to double\n", token);
										exit(EXIT_FAILURE);
									}
									token[strlen(token)] = ')';    // recover previously deleted ')'


									// check if there is another argument (syntax error)
									token = strtok(NULL, delim);
									if (token != NULL) {
										printf("Error. Too many fields in transient function\n");
										exit(EXIT_FAILURE);
									}
								}
								else if (tr_type == TR_TYPE_SIN) {
									printf("SIN function\n");

									// parse tc2 (this is the 6th and last argument of sin function)
									token = strtok(NULL, delim);

									// check for ')'
									if (strchr(token, ')') == NULL) {
										printf("Error. Parenthesis must close directly after value (%s)\n",\
										token);
										exit(EXIT_FAILURE);
									}

									token[strlen(token)-1] = '\0'; // get rid of ')'
									if (parse_double(&v6, token) == 0) {
										printf("Error. Value (%s) cannot be converted to double\n", token);
										exit(EXIT_FAILURE);
									}
									token[strlen(token)] = ')';    // recover previously deleted ')'


									// check if there is another argument (syntax error)
									token = strtok(NULL, delim);
									if (token != NULL) {
										printf("Error. Too many fields in transient function\n");
										exit(EXIT_FAILURE);
									}

								}
								else  { // if (toupper(token[1]) == 'U')
									printf("PULSE function\n");

									// parse pw (this is the 6th argument of pulse function)
									token = strtok(NULL, delim);
									if (parse_double(&v6, token) == 0) {
										printf("Error. Value (%s) cannot be converted to double\n", token);
										exit(EXIT_FAILURE);
									}

									// parse per (this is the 7th and last argument of pulse function)
									token = strtok(NULL, delim);

									// check for ')'
									if (strchr(token, ')') == NULL) {
										printf("Error. Parenthesis must close directly after value (%s)\n",\
										token);
										exit(EXIT_FAILURE);
									}

									token[strlen(token)-1] = '\0';	// get rid of ')'
									if (parse_double(&v7, token) == 0) {
										printf("Error. Value (%s) cannot be converted to double\n", token);
										exit(EXIT_FAILURE);
									}
									token[strlen(token)] = ')';		// recover previously deleted ')'


									// check if there is another argument (syntax error)
									token = strtok(NULL, delim);
									if (token != NULL) {
										printf("Error. Too many fields in transient function\n");
										exit(EXIT_FAILURE);
									}
								}


								// print the parsed arguments
								printf("arg1=%lf, arg2=%lf, arg3=%lf, arg4=%lf, "
									   "arg5=%lf, arg6=%lf, arg7=%lf\n", \
										v1, v2, v3, v4, v5, v6, v7);

							}
							else {
								printf("Error. Unknown transient function (%s)\n", token);
								exit(EXIT_FAILURE);
							}

							break;
						}
						else {
							printf("syntax error. Type '%c' has unknown sixth field (%s)\n", type, token);
							exit(EXIT_FAILURE);
						}
					}
					else {
						printf("Syntax error. Type %c has unknown sixth field (%s)\n", type, token);
						exit(EXIT_FAILURE);
					}

					/* ***************************** *
					 * END OF TRANSIENT SPEC PARSING *
					 * ***************************** */

				}
				else {
					printf("Syntax error. Type '%c' cannot have a 6th field\n", type);
					exit(EXIT_FAILURE);
				}

			}
			else if (tok_count == 7) {
				if (toupper(type) == 'M') {
					if ((toupper(token[0]) != 'L') && (token[1] != '=')) {
						printf("Syntax error: Invalid field value (%s)\n", token);
						exit(EXIT_FAILURE);
					}
					else {
						if (parse_double(&l, &token[2]) == 0) {
							printf("Syntax error. Value (%s) cannot be converted to double\n", token);
							exit(EXIT_FAILURE);
						}
					}
				}
				else {
					printf("Syntax error. Type '%c' doesn't have a 7th field\n", type);
					exit(EXIT_FAILURE);
				}
			}
			else if (tok_count == 8) {
				if (toupper(type) == 'M') {
					if ((toupper(token[0]) != 'W') && (token[1] != '=')) {
						printf("Syntax error: Invalid field value (%s)\n", token);
						exit(EXIT_FAILURE);
					}
					else {
						if (parse_double(&w, &token[2]) == 0) {
							printf("Syntax error. Value (%s) cannot be converted to double\n", token);
							exit(EXIT_FAILURE);
						}
					}
				}
				else {
					printf("Syntax error. Type '%c' doesn't have a 8th field\n", type);
					exit(EXIT_FAILURE);
				}
			}


			// token contains start of comment line. continue parsing next line
			if (rest_line_commented)
				break;


			token = strtok(NULL, delim);
		}

		// check if all the necessary fields were parsed
		if (tok_count < min_tok_count) {
			printf("Syntax error. Missing field\n");
			exit(EXIT_FAILURE);
		}

		print_parsed_information(type, name, node1_name, node2_name, node3_name, \
								 node4_name, val, has_G2, model_name, l, w);


		// Search for node and if it doesn't exist add it to the hash table and
		// give assign: node->id = ++id;

		// V<name> <node1_name> <node2_name> <val> [transient_spec]
		// I<name> <node1_name> <node2_name> <val> [G2] [transient_spec]
		// R<name> <node1_name> <node2_name> <val> [G2]
		// C<name> <node1_name> <node2_name> <val> [G2]
		// L<name> <node1_name> <node2_name> <val>
		// D<name> <node1_name> <node2_name> <model_name> [<val>]
		// M<name> <node1_name> <node2_name> <node3_name> <node4_name> <model_name> L=<val> W=<val>
		// Q<name> <node1_name> <node2_name> <node3_name> <model_name> [<val>]



		// Create structures for AC analysys info
		if (ac_usability) {
			ac_spec = (void *) malloc(sizeof(ACSpecT));
			if (ac_spec == NULL) {
				printf("Error Memory allocation problems. Exiting..\n");
				exit(EXIT_FAILURE);
			}

			((ACInfoT *)ac_spec)->mag = mag;
			((ACInfoT *)ac_spec)->phase = phase;
		}
		else {
			ac_spec = NULL;
		}

		// Create structures for transient analysys info
		// Note: update the free functions of the lists and delete the frees below (done)
		/*free(times);*/
		/*free(values);*/
		switch (tr_type) {
			case TR_TYPE_NONE:
				tran_spec_data = NULL;
				break;
			case TR_TYPE_EXP:
				tran_spec_data = (void *) malloc(sizeof(ExpInfoT));
				if (tran_spec_data == NULL) {
					printf("Error Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}

				((ExpInfoT *)tran_spec_data)->i1 = v1;
				((ExpInfoT *)tran_spec_data)->i2 = v2;
				((ExpInfoT *)tran_spec_data)->td1 = v3;
				((ExpInfoT *)tran_spec_data)->tc1 = v4;
				((ExpInfoT *)tran_spec_data)->td2 = v5;
				((ExpInfoT *)tran_spec_data)->tc2 = v6;

				break;
			case TR_TYPE_SIN:
				tran_spec_data = (void *) malloc(sizeof(SinInfoT));
				if (tran_spec_data == NULL) {
					printf("Error Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}

				((SinInfoT *)tran_spec_data)->i1 = v1;
				((SinInfoT *)tran_spec_data)->ia = v2;
				((SinInfoT *)tran_spec_data)->fr = v3;
				((SinInfoT *)tran_spec_data)->td = v4;
				((SinInfoT *)tran_spec_data)->df = v5;
				((SinInfoT *)tran_spec_data)->ph = v6;

				break;
			case TR_TYPE_PULSE:
				tran_spec_data = (void *) malloc(sizeof(PulseInfoT));
				if (tran_spec_data == NULL) {
					printf("Error Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}

				((PulseInfoT *)tran_spec_data)->i1 = v1;
				((PulseInfoT *)tran_spec_data)->i2 = v2;
				((PulseInfoT *)tran_spec_data)->td = v3;
				((PulseInfoT *)tran_spec_data)->tr = v4;
				((PulseInfoT *)tran_spec_data)->tf = v5;
				((PulseInfoT *)tran_spec_data)->pw = v6;
				((PulseInfoT *)tran_spec_data)->per = v7;

				break;
			case TR_TYPE_PWL:
				tran_spec_data = (void *) malloc(sizeof(PwlInfoT));
				if (tran_spec_data == NULL) {
					printf("Error Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}

				((PwlInfoT *)tran_spec_data)->times = times;
				((PwlInfoT *)tran_spec_data)->values = values;
				((PwlInfoT *)tran_spec_data)->total_tuples = total_tuples;

				break;
			default:
				printf("Error. Unknown transient spec function. Exiting..\n");
				exit(EXIT_FAILURE);
		}



		switch (toupper(type)) {
			case 'I':
			case 'R':
			case 'C':
				node1 = ht_get(node1_name);
				if (node1 == NULL) {
					id++;
					node1 = ht_put(node1_name, id);
				}


				node2 = ht_get(node2_name);
				if (node2 == NULL) {
					id++;
					node2 = ht_put(node2_name, id);
				}

				// add component to node component list (hashtable field)
				// update components list (lists)
				if (has_G2 == 0){
					if (insert_element(&team1_list, (toupper(type)=='I'?S_I:(toupper(type)=='R'?R:C)), \
										name, node1, node2, val, tr_type, tran_spec_data, ac_usability, ac_spec) == -1) {
						printf("insert_element. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}
				} else {
					if (insert_element(&team2_list, (toupper(type)=='I'?S_I:(toupper(type)=='R'?R:C)), \
										name, node1, node2, val, tr_type, tran_spec_data, ac_usability, ac_spec) == -1) {
						printf("insert_element. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}
				}

				break;

			case 'V':
			case 'L':

				node1 = ht_get(node1_name);
				if (node1 == NULL) {
					id++;
					node1 = ht_put(node1_name, id);
				}


				node2 = ht_get(node2_name);
				if (node2 == NULL) {
					id++;
					node2 = ht_put(node2_name, id);
				}

				// add component to node component list (hashtable field)
				// update components list (lists)
				if (insert_element(&team2_list, (toupper(type)=='V'?V:L), name, node1, node2, val, \
							       tr_type, tran_spec_data, ac_usability, ac_spec) == -1) {
					printf("insert_element. Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}

				break;
			case 'D':

				node1 = ht_get(node1_name);
				if (node1 == NULL) {
					id++;
					node1 = ht_put(node1_name, id);
				}


				node2 = ht_get(node2_name);
				if (node2 == NULL) {
					id++;
					node2 = ht_put(node2_name, id);
				}

				if (insert_diode(name, node1, node2, model_name) == -1){
					printf("insert_element. Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}

				break;
			case 'M':

				node1 = ht_get(node1_name);
				if (node1 == NULL) {
					id++;
					node1 = ht_put(node1_name, id);
				}


				node2 = ht_get(node2_name);
				if (node2 == NULL) {
					id++;
					node2 = ht_put(node2_name, id);
				}


				node3 = ht_get(node3_name);
				if (node3 == NULL) {
					id++;
					node3 = ht_put(node3_name, id);
				}


				node4 = ht_get(node4_name);
				if (node4 == NULL) {
					id++;
					node4 = ht_put(node4_name, id);
				}

				if (insert_mos(name, node1, node2, node3, node4, l, w, model_name) == -1){
					printf("insert_element. Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}


				break;
			case 'Q':

				node1 = ht_get(node1_name);
				if (node1 == NULL) {
					id++;
					node1 = ht_put(node1_name, id);
				}


				node2 = ht_get(node2_name);
				if (node2 == NULL) {
					id++;
					node2 = ht_put(node2_name, id);
				}


				node3 = ht_get(node3_name);
				if (node3 == NULL) {
					id++;
					node3 = ht_put(node3_name, id);
				}

				if (insert_bjt(name, node1, node2, node3, model_name) == -1){
					printf("insert_element. Memory allocation problems. Exiting..\n");
					exit(EXIT_FAILURE);
				}

				break;
			default :
				printf("Unknown type (%c)\n", toupper(type));
				exit(EXIT_FAILURE);
		}


		// free varialbes as space is being allocated iteratively for each line
		free(name);
		free(node1_name);
		free(node2_name);
		free(node3_name);
		free(node4_name);
		free(model_name);

		printf("\n");
		/*printf("parsed %s\n", token);*/

	}

	fclose(fp);
	free(line);
}



