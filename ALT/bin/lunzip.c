
/* lunzip - Usage:  lunzip file.lz
              Input is file.lz with 1.e10 sequences run length encoded by lzip.
              Output is the original file.

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
	char filelz[256], fileorig[256];
	struct stat  status;
	unsigned long i, j, l, n, ns, nitems;
	long k;
        float r4[16384], s4[16384];
        const float spval  =  1.e10;
        const float spvalm = -1.e10;

	if (argc != 2) {
		fprintf(stderr, "Usage: lunzip file.lz\n");
		exit(1);
	}

	if ( (infile = fopen(argv[1], "r")) == (FILE *) NULL ) {
		fprintf(stderr, "Can't read %s\n", argv[1]);
		exit(2);
	}
        if ( stat(argv[1],&status) != -1 ) {
                if ( status.st_size == 0 ) {
                        fprintf(stderr, "input file %s is zero length\n", argv[1]);
                        exit(3);
                }
        }

	strcpy(filelz,argv[1]);
	n = strlen(filelz);
	if (filelz[n-3] != '.' || filelz[n-2] != 'l' || filelz[n-1] != 'z') {
		fprintf(stderr, "input file %s is not .lz\n", argv[1]);
		exit(4);
	}
	strncpy(fileorig,filelz,n-3);
	if ( (outfile = fopen(fileorig, "w")) == (FILE *) NULL ) {
		fprintf(stderr, "Can't write %s\n", fileorig);
		exit(6);
	}

	for ( i = 1; i < status.st_size; i += 65536 ) {
		nitems = status.st_size+1 - i;
		nitems = (nitems > 65536) ? 16384 : nitems/4;
		n = fread(r4,sizeof(float),nitems,infile);
		k = -1;
		j =  0;
		while ( j < n ) {
			if (r4[j] == spvalm) {
				ns = r4[j+1];
				for ( l = 0; l < ns; l++) {
					k++;
					s4[k] = spval;
					if (k == 16384 - 1) {
						nitems = fwrite(s4,sizeof(float),16384,outfile);
						k = -1;
					}
				}
				j = j + 2;
			} else {
				k++;
				s4[k] = r4[j];
				if (k == 16384 - 1) {
					nitems = fwrite(s4,sizeof(float),16384,outfile);
					k = -1;
				}
				j++;
			}
		}
		if (k != -1) {
			nitems = fwrite(s4,sizeof(float),k+1,outfile);
		}
	}
}
