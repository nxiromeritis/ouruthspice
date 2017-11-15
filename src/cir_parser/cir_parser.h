#ifndef _CIR_PARSER_H_
#define _CIR_PARSER_H_

extern unsigned long get_components_num(char *filename);
extern void parse_cir(char *filename);
extern void parse_command(char *command);
extern unsigned char parse_double(double *d, char *str);
extern void strtoupper(char *str);
#endif
