/* SCWRLBINARY.C 

   by Roland L. Dunbrack, Jr.
   Fox Chase Cancer Center
   7701 Burholme Avenue
   Philadelphia PA 19111
   
   Copyright 1998 Roland Dunbrack.
   Last modified: March 17, 1998.
   Email: RL_Dunbrack@fccc.edu

   Makes a binary library for use with SCWRL from the 
   full, sorted backbone-dependent rotamer library.
   The library is available via anonymous ftp at

        fccc.edu

   file 

        dunbrack/pub/compressed_files/bbdep98.Feb.sortlib.gz.

   The library will be updated periodically, so the date within
   the name of the file will change.

   To make the binary file, first compile this program in this
   directory with the command "make".

   Gunzip the library with

        gunzip bbdep98.Feb.sortlib.gz

   Then run the program as follows:

        scwrlbinary < bbdep98.Feb.sortlib 

   You should get two files as output:

        scwrlbin34.lib
	README.TO.CHECK.LIBRARY

   The second file should have the first library entries for phi=-60, psi=-60
   for all amino acids.  This is just to check the binary form
   of the library.  A comparison version is presented in the directory,
   COMPARE.TO.README.TO.CHECK.LIBRARY.  They should be the same,
   or similar (if the library has been updated).

   The binary library, scwrlbin34.lib, is now ready for use.  Make
   sure its location is correctly determined in scwrl.h before compiling
   scwrl.

   */
   
#include   <stdio.h>
#include   <stdlib.h>
#include   <ctype.h>
#include   <math.h>
#include   <string.h>
#include   <time.h>

typedef struct {
  char res[4];
  int phi;
  int psi;
  int N;
  int r1;
  int r2;
  int r3;
  int r4;
  float p;
  float chi1;
  float chi2;
  float chi3;
  float chi4;
} 
librecord;

int main(int argc, char *argv[]) {
  
  int option;
  FILE *depfile;
  FILE *outfile;
  FILE *datafile;
  char res[4];
  int phi;
  int psi;
  long int iphi,ipsi;
  int N;
  int r1;
  int r2;
  int r3;
  int r4;
  float p;
  float chi1;
  float chi2;
  float chi3;
  float chi4;
  librecord librec;
  librecord record;
  long int ramasize;
  long int argseek;
  long int argsize=81L;
  long int asnseek;
  long int asnsize=9L;
  long int aspseek;
  long int aspsize=9L;
  long int cysseek;
  long int cyssize=3L;
  long int glnseek;
  long int glnsize=27L;
  long int gluseek;
  long int glusize=27L;
  long int hisseek;
  long int hissize=6L;
  long int ileseek;
  long int ilesize=9L;
  long int leuseek;
  long int leusize=9L;
  long int lysseek;
  long int lyssize=81L;
  long int metseek;
  long int metsize=27L;
  long int pheseek;
  long int phesize=6L;
  long int proseek;
  long int prosize=2L;
  long int serseek;
  long int sersize=3L;
  long int thrseek;
  long int thrsize=3L;
  long int trpseek;
  long int trpsize=9L;
  long int tyrseek;
  long int tyrsize=6L;
  long int valseek;
  long int valsize=3L;
  
  outfile=fopen("README.TO.CHECK.LIBRARY","w");
  depfile= fopen("scwrlbin34.lib","w+b");
  
  while (scanf("%s %d %d %d %d %d %d %d %f %f %f %f %f\n",
		 res,&phi,&psi,&N,&r1,&r2,&r3,&r4,&p,&chi1,&chi2,&chi3,&chi4)!=EOF) {
    
    strcpy(librec.res,res);
    librec.phi=phi;
    librec.psi=psi;
    librec.N=N;
    librec.r1=r1;
    librec.r2=r2;
    librec.r3=r3;
    librec.r4=r4;
    librec.p=p;
    librec.chi1=chi1;
    librec.chi2=chi2;
    librec.chi3=chi3;
    librec.chi4=chi4;
    
    fwrite(&librec, sizeof(librecord), 1, depfile);
  }
  
  rewind(depfile);

  ramasize=1369L;
  argseek=0L;
  asnseek=argseek+argsize*ramasize*sizeof(librecord);
  aspseek=asnseek+asnsize*ramasize*sizeof(librecord);
  cysseek=aspseek+aspsize*ramasize*sizeof(librecord);
  glnseek=cysseek+cyssize*ramasize*sizeof(librecord);
  gluseek=glnseek+glnsize*ramasize*sizeof(librecord);
  hisseek=gluseek+glusize*ramasize*sizeof(librecord);
  ileseek=hisseek+hissize*ramasize*sizeof(librecord);
  leuseek=ileseek+ilesize*ramasize*sizeof(librecord);
  lysseek=leuseek+leusize*ramasize*sizeof(librecord);
  metseek=lysseek+lyssize*ramasize*sizeof(librecord);
  pheseek=metseek+metsize*ramasize*sizeof(librecord);
  proseek=pheseek+phesize*ramasize*sizeof(librecord);
  serseek=proseek+prosize*ramasize*sizeof(librecord);
  thrseek=serseek+sersize*ramasize*sizeof(librecord);
  trpseek=thrseek+thrsize*ramasize*sizeof(librecord);
  tyrseek=trpseek+trpsize*ramasize*sizeof(librecord);
  valseek=tyrseek+tyrsize*ramasize*sizeof(librecord);

  iphi=12L*sizeof(librec);
  ipsi=12L*sizeof(librec);

  fseek(depfile, argseek + 37L*argsize*iphi +argsize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, asnseek + 37L*asnsize*iphi +asnsize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, aspseek + 37L*aspsize*iphi +aspsize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, cysseek + 37L*cyssize*iphi +cyssize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, glnseek + 37L*glnsize*iphi +glnsize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, gluseek + 37L*glusize*iphi +glusize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, hisseek + 37L*hissize*iphi +hissize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, ileseek + 37L*ilesize*iphi +ilesize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, leuseek + 37L*leusize*iphi +leusize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, lysseek + 37L*lyssize*iphi +lyssize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, metseek + 37L*metsize*iphi +metsize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, pheseek + 37L*phesize*iphi +phesize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, proseek + 37L*prosize*iphi +prosize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, serseek + 37L*sersize*iphi +sersize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, thrseek + 37L*thrsize*iphi +thrsize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, trpseek + 37L*trpsize*iphi +trpsize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, tyrseek + 37L*tyrsize*iphi +tyrsize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);
  
  fseek(depfile, valseek + 37L*valsize*iphi +valsize*ipsi, SEEK_SET);
  fread(&record, sizeof(librec),1,depfile);
  fprintf(outfile,"%3s %5d%5d%5d%5d%5d%5d%3d%10.6f%7.1f%7.1f%7.1f%7.1f\n",
	  record.res, 
	  record.phi, 
	  record.psi, 
	  record.N, 
	  record.r1, 
	  record.r2, 
	  record.r3, 
	  record.r4, 
	  record.p, 
	  record.chi1, 
	  record.chi2, 
	  record.chi3, 
	  record.chi4);

  exit(0);
}
