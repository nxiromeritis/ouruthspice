#ifndef _SPICY_H_
#define _SPICY_H_

#define BUF_MAX 128



#ifdef COLORS_ON
#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"
#define YEL  "\x1B[33m"
#define BLU  "\x1B[34m"
#define MAG  "\x1B[35m"
#define CYN  "\x1B[36m"
#define WHT  "\x1B[37m"
#else
#define NRM  ""
#define RED  ""
#define GRN  ""
#define YEL  ""
#define BLU  ""
#define MAG  ""
#define CYN  ""
#define WHT  ""
#endif


typedef unsigned char byte;
#endif
