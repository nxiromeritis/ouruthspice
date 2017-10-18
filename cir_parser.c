#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

#include "cir_parser.h"
#include "hashtable.h"
#include "spicy.h"
#include "lists.h"



// returns 1 if character c exists in components string
byte component_type_is_valid(char c) {
	if (strchr("VIRCLDMQ", toupper(c)) == NULL)
		return 0;
	else
		return 1;
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



	fp = fopen(filename, "r");
	if (fp == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}


	// read the file line-by-line
	while((getline(&line, &len, fp)) != -1) {

		// bypass comments lines or lines containing only whitespaces
		if (whitespaces_only(line))
			continue;

		if ((line[0] == '%') || (line[0] == '*'))
			continue;


		/*#ifdef DEBUG*/
		printf(GRN "<line read> %s" NRM, line);
		/*#endif*/


		/* ************** *
		 * PARSE THE LINE *
		 * ************** */

		/* reinitialise variables and parse line */
		has_G2 = 0;
		type = 'X';
		name = NULL;
		node1_name = NULL;
		node2_name = NULL;
		node3_name = NULL;
		node4_name = NULL;
		model_name = NULL;

		// -1 if not present or part of the component parsed
		val = -1;
		l = -1;
		w = -1;

		rest_line_commented = 0;
		min_tok_count = 0;
		tok_count = 0;


		// bypass possible non ascii characters (for some reason they exist in a benchmark)
		line_offset = 0;
		while ((line[line_offset] < 0) || (line[line_offset] > 127)) { line_offset++; }


		// check if the first line character is '.'
		// this indicates a possible spice command which will be ignored
		if (line[line_offset] == '.') {
			continue;
		}


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

				if ((toupper(type) == 'I') || (toupper(type) == 'R') || (toupper(type) == 'C')) {
					if ((strcmp(token, "G2") == 0) || (strcmp(token, "g2") == 0)) {
						has_G2 = 1;
					}
					else {
						printf("Syntax error. Unknown field (%s)\n", token);
						exit(EXIT_FAILURE);
					}
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
				else {
					printf("Syntax error. Type '%c' doesn't have a fifth field\n", type);
					exit(EXIT_FAILURE);
				}

			}
			else if (tok_count == 6) {

				if (toupper(type) == 'M') {
					model_name = strdup(token);
					if (model_name == NULL) {
						printf("Error. Memory allocation problem.s Exiting..\n");
						exit(EXIT_FAILURE);
					}
				}
				else if (toupper(type) == 'Q') {
					if (parse_double(&val, token) == 0) {
						printf("Syntax error. Value (%s) cannot be converted to double\n", token);
						exit(EXIT_FAILURE);
					}
				}
				else {
					printf("Syntax error. Type '%c' doesn't have a 6th field\n", type);
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

		// V<name> <node1_name> <node2_name> <val>
		// I<name> <node1_name> <node2_name> <val> [G2]
		// R<name> <node1_name> <node2_name> <val> [G2]
		// C<name> <node1_name> <node2_name> <val> [G2]
		// L<name> <node1_name> <node2_name> <val>
		// D<name> <node1_name> <node2_name> <model_name> [<val>]
		// M<name> <node1_name> <node2_name> <node3_name> <node4_name> <model_name> L=<val> W=<val>
		// Q<name> <node1_name> <node2_name> <node3_name> <model_name> [<val>]

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
					if (insert_element(&team1_list, (toupper(type)=='I'?I:(toupper(type)=='R'?R:C)), \
										name, node1, node2, val) == -1) {
						printf("insert_element. Memory allocation problems. Exiting..\n");
						exit(EXIT_FAILURE);
					}
				}else{
					if (insert_element(&team2_list, (toupper(type)=='I'?I:(toupper(type)=='R'?R:C)), \
										name, node1, node2, val) == -1) {
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
				if (insert_element(&team2_list, (toupper(type)=='V'?V:L), name, node1, node2, val) == -1) {
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
				printf("Unknown type\n");
				exit(EXIT_FAILURE);
		}

		// TODO These are temporary. Actual pointers will be passed to
		// data structure insert functions and the free will be performed by the
		// another data structure function
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



