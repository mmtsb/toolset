#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

int first=1;
int connection=-1;

double rgval=0.0;
double rgforce=0.0;
double rhoval=0.0;
double rhoforce=0.0;

void parseString(char *buf, double *target, double *force) {
  char *cptr=buf;
  *target=-1;
  *force=-1;
  for (;*cptr!=0 && *cptr!=','; cptr++);
    if (!strncmp(buf,"force=",6)) 
      *force=atof(buf+6);
    if (!strncmp(buf,"target=",7))
      *target=atof(buf+7);
    if (!strncmp(cptr+1,"force=",6)) 
      *force=atof(cptr+7);
    if (!strncmp(cptr+1,"target=",7))
      *target=atof(cptr+8);
}

void copyFortranString(char *dest, char *src, int len) {
  int i;
  for (i=0; i<len && src[i]!=0 && src[i]!=' '; i++) 
    dest[i]=src[i];
  dest[i]=0;
}

int getSocket(char *sname, char *sport) {
  struct hostent *servhost;
  struct sockaddr_in server;
  struct in_addr in;

  int sockfd;
  int ret;

  servhost=gethostbyname(sname);
  if (servhost==0) {
    fprintf(stderr,"cannot resolve host name\n");
    exit(1);
  }

  sockfd=socket(AF_INET,SOCK_STREAM,0);
  if (sockfd<0) {
    fprintf(stderr,"cannot create socket\n");
    exit(1);
  }
  
  memset((void *)&server,0,sizeof(server));
  server.sin_family=AF_INET;

#ifdef CRAY
  memcpy(&server.sin_addr.s_da,servhost->h_addr,SIZEOF_in_addr);
#else
  memcpy(&server.sin_addr.s_addr,servhost->h_addr,sizeof(in.s_addr));
#endif
  server.sin_port=htons(atoi(sport));

  ret=connect(sockfd,(struct sockaddr *)&server,sizeof(server));

  if (ret<0) {
    fprintf(stderr,"cannot connect to server %s:%s\n",sname,sport);
    exit(1);
  }

  return sockfd;
}

int sockWriteN(int sock, char *buf, int len) {
  int lefttowrite=len;
  int nwrite;
  
  do {
    nwrite=write(sock,buf,lefttowrite);
    if (nwrite<=0) return 0;
    lefttowrite-=nwrite;
    buf+=nwrite;
  } while (lefttowrite>0);

  return 1;
}

int sockWrite(int sock, char *buf) {
  return sockWriteN(sock,buf,strlen(buf));
}

int sockRead(int sock, char *buf) {
  int nread;
  int totread=0;

  do {
    nread=read(sock,buf,4095-totread);
    if (nread==0) return 0;
    totread+=nread;
    buf+=nread;
  } while (*(buf-1)!='\n');

  *buf=0;
  return nread;
}

void connectToServer(char *svname, char *svport, char *svid, char *jid) {
  char cmd[512];
  char answer[1024];

  connection=getSocket(svname,svport);

  sprintf(cmd,"CONNECT %s %s\n",svid,jid);
  
  sockWrite(connection,cmd);
  sockRead(connection,answer);
  
  if (!strstr(answer,"connected")) {
    fprintf(stderr,"cannot establish connection to %s at port %s\n%s\n",
	    svname,svport,answer);
    exit(1);
  }
}

int checkFile(char *datadir, char *fname, char *xname) {
  sprintf(xname,"%s",fname);
  if (!access(xname,R_OK)) 
    return 1;
  sprintf(xname,"%s/%s",datadir,fname);
  if (!access(xname,R_OK))
    return 1;
  return 0;
}

void sendFile(char *datadir, char *fname, char *buffer) {
  char line[512];
  char xname[512];
  FILE *inptr;
  int len;
  int nread;
  int totread;

  if (!checkFile(datadir,fname,xname)) 
    return;

  /*  fprintf(stderr,"%s file %s ok\n",datadir,xname); */

  sprintf(line,"FILE %s/%s\n",datadir,fname);
  sockWrite(connection,line);
  sockRead(connection,line);

  /*  fprintf(stderr,"%s response: %s\n",datadir,line); */
  if (strstr(line,"accepted")) {
    inptr=fopen(xname,"r");
    fseek(inptr,0,SEEK_END);
    len=ftell(inptr);
    fseek(inptr,0,SEEK_SET);
    /*    fprintf(stderr,"%s length: %d\n",datadir,len); */
    sprintf(line,"%8d",len);
    sockWrite(connection,line);
    
    totread=0;
    do {
      nread=fread(buffer,1,50000,inptr);
      totread+=nread;
      sockWriteN(connection,buffer,nread);
    } while(totread<len);
    fclose(inptr);

    /*    fprintf(stderr,"%s finished\n",datadir);  */
  }
}


#if defined CRAY
void GETBIAS(double *rg, double *krg, double *rho, double *krho) {
#else
#if defined IBM
void getbias(double *rg, double *krg, double *rho, double *krho) {
#else
void getbias_(double *rg, double *krg, double *rho, double *krho) {
#endif
#endif
  *rg=rgval;
  *krg=rgforce;
  *rho=rhoval;
  *krho=rhoforce;
}

#if defined CRAY
double NEWTEMP(char *svname, char *svport, char *svid, char *jid, char *dir,
	       double *energy, double *rg, double *rho, int *sendfiles,
	       int lensvname, int lensvport, int lensvid, int lenjid, int lendir) {
#else
#if defined IBM
double newtemp(char *svname, char *svport, char *svid, char *jid, char *dir,
	       double *energy, double *rg, double *rho, int *sendfiles,
	       int lensvname, int lensvport, int lensvid, int lenjid, int lendir) {
#else
double newtemp_(char *svname, char *svport, char *svid, char *jid, char *dir,
		double *energy, double *rg, double *rho, int *sendfiles,
		int lensvname, int lensvport, int lensvid, int lenjid, int lendir) {

#endif
#endif
  char servername[256];
  char serverport[256];
  char serverid[256];
  char jobid[256];
  char datadir[256];

  int len;

  char line[4096];
  char *buffer=(char *)malloc(50000);
  FILE *inptr;

  char mode[20];
  int run;
  int trun;
  double tempval;

  char buf1[100];
  char buf2[100];

  int nst;

  *buf1=0;
  *buf2=0;

  copyFortranString(servername,svname,lensvname);
  copyFortranString(serverport,svport,lensvport);
  copyFortranString(serverid,svid,lensvid);
  copyFortranString(jobid,jid,lenjid);
 
  if (connection<-1)
    return -99999.0;

  if (connection<0) 
    connectToServer(servername,serverport,serverid,jobid);

  if (first) {
    sprintf(line,"DATA ener N/A val N/A:N/A\n");
    first=0;
  } else {
    if (*sendfiles==1) {
    /*    fprintf(stderr,"pclient: sending files\n"); */
      copyFortranString(datadir,dir,lendir);
    /*    fprintf(stderr,"pclient: datadir: %s\n",datadir); */
      sendFile(datadir,"monsster.final.chain",buffer);
      sendFile(datadir,"restart",buffer);
      sendFile(datadir,"traj.x",buffer);
    }

    if (*rg>0 && *rho>0) {
      sprintf(line,"DATA ener %lf val %1.5lf:%1.5lf\n",*energy,*rg,*rho);
    } else if (*rho<-999) {
      sprintf(line,"DATA ener %lf val %1.5lf\n",*energy,*rg);
    } else {
      sprintf(line,"DATA ener %lf\n",*energy);
    }
  }

  sockWrite(connection,line);  

  sockRead(connection,line);
  
  nst=sscanf(line,"%s%d%d%lf%s%s",mode,&run,&trun,&tempval,buf1,buf2);

  rgval=-1.0;
  rgforce=-1.0;
  rhoval=-1.0;
  rhoforce=-1.0;

  if (nst>4) {
    parseString(buf1,&rgval,&rgforce);
    if (nst>5) {
      parseString(buf2,&rhoval,&rhoforce);
    }
  }

  free(buffer);

  /*  fprintf(stderr,"pclient: nst: %d, >%s< >%s< %lf/%lf %lf/%lf %lf\n",
      nst,buf1,buf2,rgval,rgforce,rhoval,rhoforce,tempval); */

  if (!strncmp(mode,"DONE",4)) {
    close(connection);
    connection=-99;
    return -999999.0;
  } else 
    return tempval;
}

 
#ifdef CRAY
double NEXTCONF(char *svname, char *svport, char *svid, char *jid, char *dir,
		double *energy, double *ecent, double *etemp, int *sendfiles,
		int lensvname, int lensvport, int lensvid, int lenjid, int lendir) {
#else
#if defined IBM
double nextconf(char *svname, char *svport, char *svid, char *jid, char *dir,
		double *energy, double *ecent, double *etemp, int *sendfiles,
		int lensvname, int lensvport, int lensvid, int lenjid, int lendir) {
#else
double nextconf_(char *svname, char *svport, char *svid, char *jid, char *dir,
		double *energy, double *ecent, double *etemp, int *sendfiles,
		int lensvname, int lensvport, int lensvid, int lenjid, int lendir) {
#endif
#endif
  char servername[256];
  char serverport[256];
  char serverid[256];
  char jobid[256];
  char datadir[256];

  int len;

  char line[4096];
  char *buffer=(char *)malloc(100000);
  FILE *inptr;

  int mode;

  copyFortranString(servername,svname,lensvname);
  copyFortranString(serverport,svport,lensvport);
  copyFortranString(serverid,svid,lensvid);
  copyFortranString(jobid,jid,lenjid);
 
  if (connection<-1)
    return -99999.0;

  if (connection<0) 
    connectToServer(servername,serverport,serverid,jobid);

  if (*sendfiles==1) {
    copyFortranString(datadir,dir,lendir);
    sprintf(line,"FILE %s/monsster.final.chain\n",datadir);
    sockWrite(connection,line);
    sockRead(connection,line);

    if (strstr(line,"accepted")) {
      inptr=fopen("monsster.final.chain","r");
      len=fread(buffer,1,200000,inptr);
      fclose(inptr);
      buffer[len]=0;
      sprintf(line,"%8d",len);
      sockWrite(connection,line);
      sockWrite(connection,buffer);
    }
    fclose(inptr);
  }

  sprintf(line,"DATA ener %lf ecent %lf etemp %lf\n",*energy,*ecent,*etemp);
  sockWrite(connection,line);  

  sockRead(connection,line);

  sscanf(line,"%d",&mode);

  free(buffer);

  if (mode<0) {
    close(connection);
    connection=-99;
    return -1.0;
  } else 
    return 1.0;
}
 
 
