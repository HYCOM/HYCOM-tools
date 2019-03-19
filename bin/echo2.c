
/* echo2 - writes characters to standard error
           Option:    -n   do not append newline
*/
#include <stdio.h>

main(argc,argv)
	short argc;
	char *argv[];
{
	short nl = 1;			/* set newline flag on */

	if (strcmp(argv[1],"-n") == 0) {
		nl = 0;			/* request no newline */
		argc--;			/* decrease arg count */
		argv++;			/* move to next arg */
	}

	while (--argc) {		/* loop through the args */
		fputs(*++argv,stderr);	/* put the word out to std err */
		if (argc > 1)		/* put a space between words */
		     putc(' ',stderr);	/* except for the last word */
	}
	if (nl)				/* append newline? */
	     putc('\n',stderr);
}
