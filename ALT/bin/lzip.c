
/* lzip - Usage:  lzip file
              Run length encodes 1.e10 sequences in file to produce
              a shorter file.lz.  Reverse with lunzip file.lz.

   This version only works if float is IEEE 32-bit.
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

main(argc,argv)
	short argc;
	char *argv[];
{
	FILE *infile, *outfile;
	char filelz[256];
	struct stat  status;
	unsigned long i, iout, j, n, ns, nitems;
	long k;
        float r4[16384];
        const float spval  =  1.e10;
        const float spvalm = -1.e10;

	if (argc != 2) {
		fprintf(stderr, "Usage: lzip file\n");
		exit(1);
	}

	if ( (infile = fopen(argv[1], "r")) == (FILE *) NULL ) {
		fprintf(stderr, "Can't read %s\n", argv[1]);
		exit(2);
	}
        if ( stat(argv[1],&status) != -1 ) {
                if ( status.st_size == 0 ) {
                        fprintf(stderr, "input file %s is zero length", argv[1]);
                        exit(5);
                }
        }

	strcpy(filelz,argv[1]);
	strcat(filelz,".lz");
	if ( (outfile = fopen(filelz, "w")) == (FILE *) NULL ) {
		fprintf(stderr, "Can't write %s\n", filelz);
		exit(6);
	}

	iout = 0;
	for ( i = 1; i < status.st_size; i += 65536 ) {
		nitems = status.st_size+1 - i;
		nitems = (nitems > 65536) ? 16384 : nitems/4;
		n = fread(r4,sizeof(float),nitems,infile);
		k  = -1;
		ns =  0;
		for ( j=0; j < n; j++ ) {
			if (r4[j] == spval) {
				ns = ns + 1;
			} else {
				if (ns == 0) {
					/* no spval */
				} else if (ns == 1) {
					k++;
					r4[k] = spval;
					ns = 0;
				} else if (ns == 2) {
					k++;
					r4[k] = spval;
					k++;
					r4[k] = spval;
					ns = 0;
				} else if (ns > 2) {
					/* align RLE on 8-byte boundary */
					if ((iout+k+1)%2 == 0) {
						k++;
						r4[k] = spvalm;
						k++;
						r4[k] = ns;
						ns = 0;
					} else {
						k++;
						r4[k] = spval;
						k++;
						r4[k] = spvalm;
						k++;
						r4[k] = ns-1;
						ns = 0;
					}
				}
				k++;
				r4[k] = r4[j];
			}
		}
		if (ns == 1) {
			k++;
			r4[k] = spval;
		} else if (ns == 2) {
			k++;
			r4[k] = spval;
			k++;
			r4[k] = spval;
		} else if (ns > 2) {
			if ((iout+k+1)%2 == 0) {
				k++;
				r4[k] = spvalm;
				k++;
				r4[k] = ns;
			} else {
				k++;
				r4[k] = spval;
				k++;
				r4[k] = spvalm;
				k++;
				r4[k] = ns-1;
			}
		}
		nitems = fwrite(r4,sizeof(float),k+1,outfile);
		iout += nitems;
	}
}
