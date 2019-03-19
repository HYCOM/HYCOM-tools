 
/* 
   Convert files consisting of 32-bit words between big and little endian.
 
   Usage:  endian <input-file> <output-file>
 
   Will work on f77 unformatted files, providing the record length tags
   are 32-bits.
 
   Norman Francis Beekwilder,  June 1997.
   BUFSIZ logic added by Alan J. Wallcraft.
*/
 
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<errno.h>

int main(int argc, char *argv[])
{
  char buf[BUFSIZ];
  char tmp;
  int j, n;
  FILE *in;
  FILE *out;

  if(argc != 3) {
    fprintf(stderr, "usage: endian <input-file> <output-file>\n");
    return 1;
  }

  if((in = fopen(argv[1], "rb")) == (FILE *) NULL) {
    fprintf(stderr, "endian unable to open input file for reading\n");
    return 2;
  }

  if((out = fopen(argv[2], "wb")) == (FILE *) NULL) {
    fprintf(stderr, "endian unable to open output file for writing\n");
    return 3;
  }

  for (;;) {
    n = fread(buf, 1, BUFSIZ, in);
    if (n < 4) {
      if(feof(in) != 0)
      {
        break;
      }
      fprintf(stderr, "error reading input file %d\n", errno);
      fclose(in);
      fclose(out);
      return 4;
    }
    for ( j = 0; j < n; j += 4 ) {
      tmp      = buf[j+0];
      buf[j+0] = buf[j+3];
      buf[j+3] = tmp;
      tmp      = buf[j+1];
      buf[j+1] = buf[j+2];
      buf[j+2] = tmp;
    }
    if ( fwrite(buf, 1, n, out) != n ) {
      fprintf(stderr, "error writing output file %d\n", errno);
      fclose(in);
      fclose(out);
      return 6;
    }
  }
  fclose(in);
  fclose(out);
  return 0;
}
