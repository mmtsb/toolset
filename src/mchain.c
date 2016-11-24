// http://mmtsb.scripps.edu/doc/mchain.html
// 2000, Michael Feig, Brooks group, TSRI
//
// *** C++ code ***
//
// requires:
//   vector.h
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "field.h"
#include "vector.h"
#include "pdb.h"

void usage() {
  fprintf(stderr,"usage:    mchain [options] { -p pdb | -n residues }\n");
  fprintf(stderr,"options:  [-o ox oy oz] [-g gridsize] [-center]\n");
  fprintf(stderr,"          [-l min:max[=min:max ...]]\n");
  fprintf(stderr,"          [-seed value]\n");
  exit(1);
}

const int MAXRESIDUES=2000;

class IVector {
public:
  int x,y,z;
};

Vector offset;
double resolution=1.45;
int gridsize;
int center;

int outofbounds=0;

class UnSuccessful {};

class Stuck {
 private:
  int num;

 public:
  Stuck(int n) : num(n) {}
  int number() { return num; }
};

void readPDB(char *pdbName, Vector *sc, int *valid) {
  FILE *fptr;

  fptr= (!strcmp(pdbName,"-")) ? stdin : fopen(pdbName,"r");
  
  if (fptr==0) {
    fprintf(stderr,"Cannot open PDB file %s\n",pdbName);
    exit(1);
  }

  int lastnum=-1;
  
  Vector sum(0.0,0.0,0.0);
  int natoms=0;

  PDBEntry pdbEntry;

  while (!feof(fptr)) {
    if (pdbEntry.read(fptr)>0) {
      if (pdbEntry.residueNumber()!=lastnum) {
	if (lastnum>=0) {
	  if (valid[lastnum]==0 && lastnum<MAXRESIDUES) {
	    valid[lastnum]=1;
	    sc[lastnum]=sum/((double)natoms);
	    sc[lastnum]/=resolution;
	  }
	  sum=Vector(0.0,0.0,0.0);
	  natoms=0;
	}
	lastnum=pdbEntry.residueNumber();
      }

      if (strcmp(pdbEntry.atomName(),"O") && 
	  strcmp(pdbEntry.atomName(),"N") && 
	  strcmp(pdbEntry.atomName(),"C") && 
	  pdbEntry.type()!='H') {
	sum+=pdbEntry.coordinates();
	natoms++;
      }
    }
  }

  fclose(fptr);

  if (valid[lastnum]==0 && lastnum<MAXRESIDUES) {
    valid[lastnum]=1;
    sc[lastnum]=sum/(double)natoms;
    sc[lastnum]/=resolution;
  }

  if (center) {
    Vector cofm(0.0,0.0,0.0);
    int ncofm=0;
    int i;

    for (i=0; i<MAXRESIDUES; i++) {
      if (valid[i]) {
	cofm+=sc[i];
	ncofm++;
      }
    }
    cofm/=ncofm;

    //    fprintf(stderr,"ncofm=%d, cofm: %lf %lf %lf\n",ncofm,cofm.x(),cofm.y(),cofm.z());
    for (i=0; i<MAXRESIDUES; i++) {
      if (valid[i]) {
	sc[i]-=cofm;
      }
    }
  }
}

void genint(IVector *icofm, Vector *scofm, int i) {
  icofm[i].x=int(scofm[i].x()+offset.x()+0.5000001);
  icofm[i].y=int(scofm[i].y()+offset.y()+0.5000001);
  icofm[i].z=int(scofm[i].z()+offset.z()+0.5000001);
}

int colinear(IVector* icofm, int i, int j, int k) {
  register int ix=icofm[i].x-icofm[j].x;
  register int iy=icofm[i].y-icofm[j].y;
  register int iz=icofm[i].z-icofm[j].z;

  register int jx=icofm[j].x-icofm[k].x;
  register int jy=icofm[j].y-icofm[k].y;
  register int jz=icofm[j].z-icofm[k].z;

  register int kx=iy*jz-iz*jy;
  register int ky=jx*iz-ix*jz;
  register int kz=ix*jy-iy*jx;

  return (kx*kx+ky*ky+kz*kz)==0;
}

inline int idiff(IVector *icofm, int i, int j) {
  register int ix=icofm[i].x-icofm[j].x;
  register int iy=icofm[i].y-icofm[j].y;
  register int iz=icofm[i].z-icofm[j].z;
  
  return ix*ix+iy*iy+iz*iz;
}

int gridQ(int* grid, int x, int y, int z) {
  if (x<10 || y<10 || z<10 || x>=gridsize-10 || y>=gridsize-10 || z>=gridsize-10 ) {
    outofbounds=1;
    return 0;
  } else
    return !grid[x*gridsize*gridsize+y*gridsize+z];
}

int checkCList(int inx, IVector *cl, int *ic, int nc, int x, int y, int z) {
  for (int i=0; i<nc; i++) {
    if (inx < ic[i]) {
      int dx=cl[i].x-x;
      int dy=cl[i].y-y;
      int dz=cl[i].z-z;
      if ((dx*dx+dy*dy+dz*dz)<9) 
	return 0;
    }
  }
  return 1;
}

void setGrid(int *grid, int x, int y, int z) {

  for (int ix=-2; ix<=2; ix++) {
    for (int iy=-2; iy<=2; iy++) {
      for (int iz=-2; iz<=2; iz++) {
	int sq=ix*ix+iy*iy+iz*iz;
	if (sq<9) {
	  int inx=(x+ix)*gridsize*gridsize+(y+iy)*gridsize+(z+iz);
	  grid[inx]=1;
	}
      }
    }
  }
}

void setGrid(int *grid, int x, int y, int z, int *slist, int& nslist) {
  for (int ix=-2; ix<=2; ix++) {
    for (int iy=-2; iy<=2; iy++) {
      for (int iz=-2; iz<=2; iz++) {
	int sq=ix*ix+iy*iy+iz*iz;
	if (sq<9) {
	  int inx=(x+ix)*gridsize*gridsize+(y+iy)*gridsize+(z+iz);
	  slist[nslist++]=inx;
	  slist[nslist++]=grid[inx];
	  grid[inx]=1;
	}
      }
    }
  }
}

void resetGrid(int *grid, int *slist, int start, int end) {
  //  fprintf(stderr,"resetGrid: %d %d\n",start,end);
  for (int i=start; i<end; i++) {
    int inx=slist[i++];
    grid[inx]=slist[i];
  }
}

class Location {
 public:
  int x,y,z;
  double distance;
};

void showGrid(int *grid, int x, int y, int z) {
  fprintf(stderr,"grid %d / %d / %d\n",x,y,z);
  for (int ix=-2; ix<=2; ix++) {
    for (int iy=-2; iy<=2; iy++) {
      for (int iz=-2; iz<=2; iz++) {
        fprintf(stderr,"%d ",grid[(x+ix)*gridsize*gridsize+(y+iy)*gridsize+(z+iz)]);
      }
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
  }
}

  
int add(int *have, int i, Vector *scofm, IVector *icofm, int *grid, 
	IVector *clist, int *iclist, int nclist,
	int *slist, int& nslist, double lim, int jstart=0) {

  int j,k;

  Location loc[25];
  int nloc=0;
  int iloc[25];

  genint(icofm,scofm,i);

  int x=icofm[i].x;
  int y=icofm[i].y;
  int z=icofm[i].z;
  for (int ix=-2; ix<=2; ix++) {
    int tx=x+ix;
    for (int iy=-2; iy<=2; iy++) {
      int ty=y+iy;
      for (int iz=-2; iz<=2; iz++) {
	int tz=z+iz;
	if (gridQ(grid,tx,ty,tz) && (nclist==0 || checkCList(i,clist,iclist,nclist,tx,ty,tz))) {
	  double dx=scofm[i].x()-(double)(tx-offset.x());
	  double dy=scofm[i].y()-(double)(ty-offset.y());
	  double dz=scofm[i].z()-(double)(tz-offset.z());
  
	  double d=dx*dx+dy*dy+dz*dz;

	  if (d<lim && nloc<25) {
	    loc[nloc].x=tx;
	    loc[nloc].y=ty;
	    loc[nloc].z=tz;

	    loc[nloc].distance=d;
	    nloc++;
	  }
	}
      }
    }
  }

  //    fprintf (stderr, " %2d ",nloc);

  for (j=0; j<nloc; j++) 
    iloc[j]=j;
      
  for (j=0; j<nloc-1; j++)
    for (k=j+1; k<nloc; k++)
      if (loc[iloc[j]].distance>loc[iloc[k]].distance) {
	int t=iloc[j];
	iloc[j]=iloc[k];
	iloc[k]=t;
      }

  // fprintf(stderr,"  i=%d, lim=%lf, jstart=%d, found:\n",i,lim,jstart);

  //  for (j=0; j<nloc; j++) {
  //    fprintf(stderr,"    %d %d %d %lf\n",loc[iloc[j]].x,loc[iloc[j]].y,loc[iloc[j]].z,
  //  	    loc[iloc[j]].distance);
  //  }

  for (j=jstart; j<nloc; j++) {
    icofm[i].x=loc[iloc[j]].x;
    icofm[i].y=loc[iloc[j]].y;
    icofm[i].z=loc[iloc[j]].z;

    if (i>0 && have[i-1] && idiff(icofm,i-1,i)>30) ;
    else if (i>1 && have[i-2] && idiff(icofm,i-2,i)>68) ;
    else if (i>1 && have[i-1] && have[i-2] && colinear(icofm,i-2,i-1,i)) ;
    else {
      //            fprintf(stderr,"  %d: set %d, set grid around %d / %d / %d, Q:%d\n",
      //      	      i,j,icofm[i].x,icofm[i].y,icofm[i].z,
      //      	      gridQ(grid,icofm[i].x,icofm[i].y,icofm[i].z));

      //            fprintf(stderr," s %d ",j);
      
      setGrid(grid,icofm[i].x,icofm[i].y,icofm[i].z,slist,nslist);
      
      break;
    }
  }

  if (j>=nloc) {
    icofm[i].x=x;
    icofm[i].y=y;
    icofm[i].z=z;
    return -1;
  }

  return j;
}

void validStart(IVector& ip2, IVector& ip1, IVector& iv, int *grid, 
		int tx, int ty, int tz) {

  int nvalid=0;
  int validx[1000];
  int validy[1000];
  int validz[1000];

  int jx=ip1.x-ip2.x;
  int jy=ip1.y-ip2.y;
  int jz=ip1.z-ip2.z;

  for (int m=1; m<=4; m++) {
    for (int ix=-m; ix<=m; ix++) 
      for (int iy=-m; iy<=m; iy++) 
	for (int iz=-m; iz<=m; iz++) {
	  int vx=iv.x+ix;
	  int vy=iv.y+iy;
	  int vz=iv.z+iz;

	  if (gridQ(grid,vx,vy,vz)) {
	    int dx=vx-ip1.x;
	    int dy=vy-ip1.y;
	    int dz=vz-ip1.z;

	    int id=dx*dx+dy*dy+dz*dz;
	    if (id<=30) {
	      int pdx=vx-ip2.x;
	      int pdy=vy-ip2.y;
	      int pdz=vz-ip2.z;
	    
	      int pd=pdx*pdx+pdy*pdy+pdz*pdz;
	      if (pd<=68) {
		int kx=dy*jz-dz*jy;
		int ky=jx*dz-dx*jz;
		int kz=dx*jy-dy*jx;
		if ((kx*kx+ky*ky+kz*kz)!=0) {
                  if (nvalid<1000) {
		    validx[nvalid]=iv.x+ix;
		    validy[nvalid]=iv.y+iy;
		    validz[nvalid]=iv.z+iz;
		    ++nvalid;
		  }
		}
	      }
	    }
	  }
	}
  }

  if (nvalid>0) {
    double r=(drand48()*(double)nvalid);
    int ir=(int)r;
    iv.x=validx[ir];
    iv.y=validy[ir];
    iv.z=validz[ir];
    return;
  } else {
    fprintf(stderr,"cannot find valid starting point\n");
    exit(1);
  }
}


int findClosest(IVector& v, int x, int y, int z, int *grid, int maxval) {
  for (int ff=0; ff<=maxval; ff++) 
    for (int ix=-ff; ix<=ff; ix++) 
      for (int iy=-ff; iy<=ff; iy++) 
	for (int iz=-ff; iz<=ff; iz++) {
	  if (gridQ(grid,x+ix,y+iy,z+iz)) {
	    v.x=x+ix;
	    v.y=y+iy;
	    v.z=z+iz;
	    return 1;
	  }
	}
  return 0;
}

int randomwalk(IVector *icofm, int *grid, int first, int last,
	       IVector *walk, int maxwalk) {

  int nwalk=1;

  walk[0].x=icofm[first-1].x;
  walk[0].y=icofm[first-1].y;
  walk[0].z=icofm[first-1].z;

  int tx=icofm[last+1].x;
  int ty=icofm[last+1].y;
  int tz=icofm[last+1].z;

  validStart(icofm[first-2],icofm[first-1],walk[0],grid,tx,ty,tz);
  
  double dist;
  do {
    int x=walk[nwalk-1].x;
    int y=walk[nwalk-1].y;
    int z=walk[nwalk-1].z;
    
    Vector d=Vector((double)tx,(double)ty,(double)tz)-
      Vector((double)x,(double)y,(double)z);
    
    dist=d.norm();
    d/=d.norm();
    d*=2.0;

    if (dist>4) {
      int nx,ny,nz;
      int ecnt=0;
      do {
	double rx=drand48()*3.0-1.0+d.x();
	double ry=drand48()*3.0-1.0+d.y();
	double rz=drand48()*3.0-1.0+d.z();
	
	nx=x+(int)rx;
	ny=y+(int)ry;
	nz=z+(int)rz;
        
        ecnt++;
      } while (!findClosest(walk[nwalk],nx,ny,nz,grid,5) && ecnt<10);

      if (ecnt>=10) {
	//        fprintf(stderr,"randomwalk stuck\n");
	throw UnSuccessful();
      }

      if (++nwalk>=maxwalk) {
	//        fprintf(stderr,"reached maxwalk\n");
        throw UnSuccessful();
      }
    }
  } while(dist>4);
  
  return nwalk;
}

void initialdist(IVector *icofm, IVector *walk, int nwalk,
		 int first, int last) {
  icofm[first].x=walk[0].x;
  icofm[first].y=walk[0].y;
  icofm[first].z=walk[0].z;

  icofm[last].x=walk[nwalk-1].x;
  icofm[last].y=walk[nwalk-1].y;
  icofm[last].z=walk[nwalk-1].z;

  double del=(double)(nwalk-2)/(double)(last-first-1);
  double v=-0.0001;
  for (int i=1; i<=last-first-1; i++) {
    v+=del;
    icofm[first+i].x=walk[(int)v].x;
    icofm[first+i].y=walk[(int)v].y;
    icofm[first+i].z=walk[(int)v].z;
  }    
}

void move(IVector *icofm, int i, int *grid, int sx, int sy, int sz) {
  Vector sdi=Vector((double)icofm[i].x,(double)icofm[i].y,
		   (double)icofm[i].z)-
             Vector((double)sx,(double)sy,(double)sz);
  double dist=sdi.norm();
  if (dist>0.0001) {
    sdi/=dist;
    sdi*=1.5;
  }

  int nx,ny,nz;
  int cnt=0;
  do {
    double rx=drand48()*2.0-1.0+sdi.x();
    double ry=drand48()*2.0-1.0+sdi.y();
    double rz=drand48()*2.0-1.0+sdi.z();
    nx=icofm[i].x+((rx<-0.5)?-1:((rx>0.5)?1:0));
    ny=icofm[i].y+((ry<-0.5)?-1:((ry>0.5)?1:0));
    nz=icofm[i].z+((rz<-0.5)?-1:((rz>0.5)?1:0));
  } while(!gridQ(grid,nx,ny,nz) && ++cnt<1000);

  if (gridQ(grid,nx,ny,nz)) {
    icofm[i].x=nx;
    icofm[i].y=ny;
    icofm[i].z=nz;
  }
}


void stretch(IVector *icofm, int i, int j, int *grid, int sx, int sy, int sz) {
  Vector d=Vector((double)icofm[i].x,(double)icofm[i].y,
		  (double)icofm[i].z)-
           Vector((double)icofm[j].x,(double)icofm[j].y,
		  (double)icofm[j].z);

  double dist=d.norm();
  if (dist>0.0001) {
    d/=dist;
  }

  Vector sdi=Vector((double)icofm[i].x,(double)icofm[i].y,
		   (double)icofm[i].z)-
             Vector((double)sx,(double)sy,(double)sz);
  dist=sdi.norm();
  if (dist>0.0001) {
    sdi/=dist;
    sdi*=1.8;
  }

  Vector sdj=Vector((double)icofm[j].x,(double)icofm[j].y,
		    (double)icofm[j].z)-
             Vector((double)sx,(double)sy,(double)sz);
  dist=sdj.norm();
  if (dist>0.0001) {
    sdj/=dist;
    sdj*=1.8;
  }

  int nx,ny,nz;
  int cnt=0;
  do {
    double rx=drand48()*2.0+d.x()-1.0+sdi.x();
    double ry=drand48()*2.0+d.y()-1.0+sdi.y();
    double rz=drand48()*2.0+d.z()-1.0+sdi.z();
    nx=icofm[i].x+((rx<-0.5)?-1:((rx>0.5)?1:0));
    ny=icofm[i].y+((ry<-0.5)?-1:((ry>0.5)?1:0));
    nz=icofm[i].z+((rz<-0.5)?-1:((rz>0.5)?1:0));
  } while(!gridQ(grid,nx,ny,nz) && ++cnt<1000);

  if (gridQ(grid,nx,ny,nz)) {
    icofm[i].x=nx;
    icofm[i].y=ny;
    icofm[i].z=nz;
  }

  cnt=0;
  do {
    double rx=drand48()*2.0+d.x()-1.0+sdj.x();
    double ry=drand48()*2.0+d.y()-1.0+sdj.y();
    double rz=drand48()*2.0+d.z()-1.0+sdj.z();
    nx=icofm[i].x-((rx<-0.5)?-1:((rx>0.5)?1:0));
    ny=icofm[i].y-((ry<-0.5)?-1:((ry>0.5)?1:0));
    nz=icofm[i].z-((rz<-0.5)?-1:((rz>0.5)?1:0));
  } while(!gridQ(grid,nx,ny,nz) && ++cnt<1000);

  if (gridQ(grid,nx,ny,nz)) {
    icofm[i].x=nx;
    icofm[i].y=ny;
    icofm[i].z=nz;
  }
}

void contract(IVector *icofm, int i, int j, int *grid, int fixed) {
  Vector d=Vector((double)icofm[i].x,(double)icofm[i].y,
		  (double)icofm[i].z)-
           Vector((double)icofm[j].x,(double)icofm[j].y,
		  (double)icofm[j].z);
  double dist=d.norm();
  if (dist>0.0001) {
    d/=d.norm();
  }

  int cnt=0;

  int nx,ny,nz;
  do {
    double rx=drand48()*2.0+d.x()-1.0;
    double ry=drand48()*2.0+d.y()-1.0;
    double rz=drand48()*2.0+d.z()-1.0;
    nx=icofm[i].x-((rx<-0.5)?-1:((rx>0.5)?1:0));
    ny=icofm[i].y-((ry<-0.5)?-1:((ry>0.5)?1:0));
    nz=icofm[i].z-((rz<-0.5)?-1:((rz>0.5)?1:0));
  } while(!gridQ(grid,nx,ny,nz) && ++cnt<1000);

  if (gridQ(grid,nx,ny,nz)) {
    icofm[i].x=nx;
    icofm[i].y=ny;
    icofm[i].z=nz;
  }

  if (!fixed) {
    cnt=0;
    do {
      double rx=drand48()*2.0+d.x()-1.0;
      double ry=drand48()*2.0+d.y()-1.0;
      double rz=drand48()*2.0+d.z()-1.0;
      nx=icofm[j].x+((rx<-0.5)?-1:((rx>0.5)?1:0));
      ny=icofm[j].y+((ry<-0.5)?-1:((ry>0.5)?1:0));
      nz=icofm[j].z+((rz<-0.5)?-1:((rz>0.5)?1:0));
    } while(!gridQ(grid,nx,ny,nz) && ++cnt<1000);
   
    if (gridQ(grid,nx,ny,nz)) {
      icofm[j].x=nx;
      icofm[j].y=ny;
      icofm[j].z=nz;
    } 
  }
}

int adjust(IVector *icofm, int *grid, int first, int last, 
	   int sx, int sy, int sz, int *valid) {
  int i,j;
  int id;

  int cnt=0;
  int violate;
  do {
    violate=0;
    for (i=first; i<last; i++) {
      for (j=i+1; j<=last; j++) {
	id=idiff(icofm,i,j);
	if (id<9) {
	  stretch(icofm,i,j,grid,sx,sy,sz);
	  violate++;
	}
      }
      if ((id=idiff(icofm,i,i+1))>30) {
	contract(icofm,i,i+1,grid,0);
	violate++;
      }
      if (i>first && (id=idiff(icofm,i-1,i+1))>68) {
	contract(icofm,i-1,i+1,grid,0);
	violate++;
      }
      if (i>first && (colinear(icofm,i-1,i,i+1))) {
	move(icofm,i,grid,sx,sy,sz);
	violate++;
      }
    }

    if ((id=idiff(icofm,first,first-1))>30) {
      contract(icofm,first,first-1,grid,1);
      violate++;
    }

    if ((id=idiff(icofm,last,last+1))>30) {
      contract(icofm,last,last+1,grid,1);
      violate++;
    }

    if ((id=idiff(icofm,first,first-2))>68) {
      contract(icofm,first,first-2,grid,1);
      violate++;
    }

    if ((id=idiff(icofm,first+1,first-1))>68) {
      contract(icofm,first+1,first-1,grid,1);
      violate++;
    }
    
    if ((id=idiff(icofm,last,last+2))>68) {
      contract(icofm,last,last+2,grid,1);
      violate++;
    }

    if ((id=idiff(icofm,last-1,last+1))>68) {
      contract(icofm,last-1,last+1,grid,1);
      violate++;
    }

    if ((colinear(icofm,first-1,first,first+1))) {
      move(icofm,first,grid,sx,sy,sz);
      violate++;
    }

    if ((colinear(icofm,last-1,last,last+1))) {
      move(icofm,last,grid,sx,sy,sz);
      violate++;
    }

    if (valid[first-2] && (colinear(icofm,first-2,first-1,first))) {
      move(icofm,first,grid,sx,sy,sz);
      violate++;
    }

    if (valid[last+2] && (colinear(icofm,last,last+1,last+2))) {
      move(icofm,last,grid,sx,sy,sz);
      violate++;
    }
  } while (violate>0 && ++cnt<1000);

  return (violate==0);
}

void addbead(int inx, int *grid, IVector *icofm, int a1, int a2) {
  int dx=icofm[a1].x-icofm[a2].x;
  int dy=icofm[a1].y-icofm[a2].y;
  int dz=icofm[a1].z-icofm[a2].z;

  icofm[inx].x=icofm[a1].x+dx;
  icofm[inx].y=icofm[a1].y+dy;
  icofm[inx].z=icofm[a1].z+dz;
  
  double rx=drand48()*3.0;
  double ry=drand48()*3.0;
  double rz=drand48()*3.0;
  icofm[inx].x+=((rx<1.0)?-1:((rx<2.0)?0:1));
  icofm[inx].y+=((ry<1.0)?-1:((ry<2.0)?0:1));
  icofm[inx].z+=((rz<1.0)?-1:((rz<2.0)?0:1));

  int px=icofm[a2].x;
  int py=icofm[a2].y;
  int pz=icofm[a2].z;
  
  int nx=icofm[a1].x;
  int ny=icofm[a1].y;
  int nz=icofm[a1].z;

  int x=icofm[inx].x;
  int y=icofm[inx].y;
  int z=icofm[inx].z;

  int mincd=9999999;
  int savx,savy,savz;
  int found=0;

  int best=(drand48()<0.2)?1:0;

  for (int ff=0; ff<=5; ff++) 
    for (int ix=-ff; ix<=ff; ix++) 
      for (int iy=-ff; iy<=ff; iy++) 
	for (int iz=-ff; iz<=ff; iz++) {
          if (gridQ(grid,x+ix,y+iy,z+iz)) { 
	    int ndx=nx-(x+ix);
	    int ndy=ny-(y+iy);
	    int ndz=nz-(z+iz);
	    int nd=ndx*ndx+ndy*ndy+ndz*ndz;
	    
	    int pdx=px-(x+ix);
	    int pdy=py-(y+iy);
	    int pdz=pz-(z+iz);
	    int pd=pdx*pdx+pdy*pdy+pdz*pdz;

	    int kx=ndy*dz-ndz*dy;
	    int ky=ndz*dx-ndx*dz;
	    int kz=ndx*dy-ndy*dx;
	    int kd=kx*kx+ky*ky+kz*kz;

	    if (nd<=30 && pd<=68 && kd!=0) {
	      if (!best) {
		icofm[inx].x=x+ix;
		icofm[inx].y=y+iy;
		icofm[inx].z=z+iz;
		setGrid(grid,icofm[inx].x,icofm[inx].y,icofm[inx].z);      
		return;
	      } else {
		int cdx=(x+ix-gridsize/2);
		int cdy=(y+iy-gridsize/2);
		int cdz=(z+iz-gridsize/2);
		int cd=cdx*cdx+cdy*cdy+cdz*cdz;
		if (cd<mincd) {
		  savx=x+ix;
		  savy=y+iy;
		  savz=z+iz;
		  mincd=cd;
		  found=1;
		}
	      }
	    }
	  }
	}

  if (found) {
    icofm[inx].x=savx;
    icofm[inx].y=savy;
    icofm[inx].z=savz;
    setGrid(grid,icofm[inx].x,icofm[inx].y,icofm[inx].z);      
  } else 
    throw UnSuccessful();
}


void generateChain(Vector *scofm, IVector *icofm, int *grid, int *valid) {
  int *list=new int[MAXRESIDUES];
  int nlist=0;

  int *slist=new int[MAXRESIDUES*200];
  int nslist=0;

  int *islist=new int[MAXRESIDUES];

  double *limit=new double[MAXRESIDUES];
  int start[MAXRESIDUES];
  int i;

  int *cflag=new int[MAXRESIDUES];
  for (i=0; i<MAXRESIDUES; cflag[i++]=0);

  IVector *clist=new IVector[MAXRESIDUES];
  int nclist=0;
  int *iclist=new int[MAXRESIDUES];

  int *have=new int[MAXRESIDUES];

  for (i=0; i<MAXRESIDUES; i++) {
    have[i]=0;
    start[i]=0;
    limit[i]=0.75;
    if (valid[i]) {
      list[nlist++]=i;
    }
  }
  limit[0]=8.0;

  int redo=0;
  int il=0;
  while(il<nlist) {
    i=list[il];
    //    fprintf(stderr,"*** RES %d %d %d\n",i,start[il],nslist);

    if (redo)
      start[il]++;

    //        fprintf(stderr,"%3d RES %3d %4s %2d %6d %8.2lf %6d :: ",
    //    	    il,i,(redo)?"REDO":"",start[il],islist[il],limit[il],nslist);

    islist[il]=nslist;
    
    int gotinx=add(have,i,scofm,icofm,grid,clist,iclist,nclist,slist,nslist,limit[il],start[il]);
    if (gotinx<0) {
      if (il==0) {
	if (outofbounds) {
	  fprintf(stderr,"Sorry. Cannot generate chain\n");
	  fprintf(stderr,"Structure outside grid?\n");
	} else {
	  fprintf(stderr,"beginning of chain, should not have trouble here\n");
	}
	exit(1);
      }
      
      //            fprintf(stderr," :: restart\n");

      limit[il]+=0.1;
      //                  fprintf(stderr,"---- restart ----\n");

      if (limit[il]>8.0) {

	limit[il]=0.75;

	int imin=1;
	int imin3=-1;
	int imin4=-1;
	int imin5=-1;
	//	double dmin=99999.0;
	for (int is=0; is<il-1; is++) {
	  double dx=(double)icofm[list[is]].x-scofm[list[il]].x()-offset.x();
	  double dy=(double)icofm[list[is]].y-scofm[list[il]].y()-offset.y();
	  double dz=(double)icofm[list[is]].z-scofm[list[il]].z()-offset.z();
	  double dd=sqrt(dx*dx+dy*dy+dz*dz);
	  if (dd<3.0 && imin3<0) {
	    imin3=is;
	  }
	  if (dd<4.0 && imin4<0) {
	    imin4=is;
	  } 
	  if (dd<5.0 && imin5<0) {
	    imin5=is;
	  }
	}

	if (imin3>=0) 
	  imin=imin3;
	else if (imin4>=0) 
	  imin=imin4;
	else if (imin5>=0) 
	  imin=imin5;

	if (cflag[i]==0) {
	  iclist[nclist]=i;
	  genint(icofm,scofm,i);
	  clist[nclist].x=icofm[i].x;
	  clist[nclist].y=icofm[i].y;
	  clist[nclist].z=icofm[i].z;
	  nclist++;
	}

	cflag[i]++;

	//	fprintf(stderr,"set %d (%d) going back to %d\n",i,cflag[i],list[imin]);

	//	for (j=0; j<nclist; j++) {
	//	  fprintf(stderr,"set %d %d : %d %d %d\n",j,iclist[j],clist[j].x,clist[j].y,clist[j].z);
	//	}
	
	//	fprintf(stderr,"  %d %d %d\n",
	//		icofm[list[imin]].x,icofm[list[imin]].y,icofm[list[imin]].z);
	for (;il>imin; il--) {
	  start[il]=0;
	  //	  limit[il]=0.75;
	  resetGrid(grid,slist,islist[il-1],islist[il]);
	  nslist=islist[il-1];
	}
      } else {
	start[il]=0;
	resetGrid(grid,slist,islist[il-1],islist[il]);
	nslist=islist[il-1];
	il--;
      }
      redo=1;
    } else {
      start[il]=gotinx;
      //            fprintf(stderr," :: ok\n");
	    //      fprintf(stderr,"set %d  %d %d %d\n",
      //	      i,icofm[i].x,icofm[i].y,icofm[i].z);
      
      redo=0;
      il++;
    }
    
    have[i]=1;
  }

  delete slist;
  delete list;
  delete limit;
  delete have;
}

int finalcheck(IVector* icofm, int minres, int maxres) {
  int i,j,id;

  int violate=0;
  for (i=minres; i<maxres; i++) {
    for (j=i+1; j<=maxres; j++) {
      id=idiff(icofm,i,j);
      if (id<9) {
	fprintf(stderr,"%d <-> %d  too short (%d)\n",i,j,id);
	violate++;
      }
    }
    if ((id=idiff(icofm,i,i+1))>30) {
      fprintf(stderr,"%d <-> %d too large (%d)\n",i,i+1,id);
      violate++;
    }
    if (i>minres && (id=idiff(icofm,i-1,i+1))>68) {
      fprintf(stderr,"%d <-> %d too large (%d)\n",i-1,i+1,id);
      violate++;
    }
    if (i>minres && (colinear(icofm,i-1,i,i+1))) {
      fprintf(stderr,"%d <-> %d colinear \n",i-1,i+1);
      violate++;
    }
  }
  return violate;
}

void addrandomchain(int anchor, int tores, IVector* icofm, int *grid) {
  int dx=(tores<anchor)?-1:1;
  for (int i=anchor+dx; i!=tores+dx; i+=dx) 
    addbead(i,grid,icofm,i-dx,i-dx-dx);
}

int main(int argc, char **argv) {
  int i,j;
  char *initpdb=0;

  gridsize=100;
  offset=50.0;
  resolution=1.45;
  center=0;

  if (argc<2) 
    usage();

  int *valid=new int[MAXRESIDUES];
  for (i=0; i<MAXRESIDUES; valid[i++]=0);

  int *resmin=new int[MAXRESIDUES];
  int *resmax=new int[MAXRESIDUES];
  int nres=0;

  long seed=(long)time(0)+(long)getpid();

  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i],"-p")) {
      initpdb=argv[++i];
    } else if (!strcmp(argv[i],"-l")) {
      Field f(argv[++i],"=");
      for (j=0; j<f.number(); j++) {
	Field g(f[j],":");
	resmin[nres]=atoi(g[0]);
	resmax[nres++]=atoi(g[1]);
	if (resmax<resmin) {
	  fprintf(stderr,"invalid interval\n");
	  exit(1);
	}
      }
    } else if (!strcmp(argv[i],"-g")) {
      gridsize=atoi(argv[++i]);
    } else if (!strcmp(argv[i],"-n")) {
      resmin[0]=1;
      resmax[0]=atoi(argv[++i]);
      nres=1;
    } else if (!strcmp(argv[i],"-center")) {
      center=1;
    } else if (!strcmp(argv[i],"-o")) {
      offset.x()=atof(argv[++i]);
      offset.y()=atof(argv[++i]);
      offset.z()=atof(argv[++i]);
    } else if (!strcmp(argv[i],"-seed")) {
      seed=atol(argv[++i]);
    } else if (!strcmp(argv[i],"-help")) {
      usage();
    }
  }

  srand48(seed);

  int *grid;
  grid=new int[gridsize*gridsize*gridsize];
  for (i=0; i<gridsize*gridsize*gridsize; grid[i++]=0);

  int *tgrid;
  tgrid=new int[gridsize*gridsize*gridsize];
  for (i=0; i<gridsize*gridsize*gridsize; tgrid[i++]=0);

  IVector *icofm;

  int templ=0;

  if (initpdb!=0) {
    Vector *scofm=new Vector[MAXRESIDUES];
    readPDB(initpdb,scofm,valid);
    for (i=0; i<nres; i++) 
      for (j=resmin[i]; j<=resmax[i]; j++) 
	valid[j]=0;
    icofm=new IVector[MAXRESIDUES];
    generateChain(scofm,icofm,grid,valid);
    templ=1;
  } else 
    icofm=new IVector[MAXRESIDUES];
  
  int redocnt,redo;
  
  for (int ires=0; ires<nres; ires++) {
    for (i=0; i<gridsize*gridsize*gridsize; i++)
      tgrid[i]=grid[i];
  
    int sx=0;
    int sy=0;
    int sz=0;
    int ns=0;

    redocnt=0;
    do {
      redo=0;
      try {
	if (!templ) {
	  int midres=(resmax[ires]-resmin[ires]+1)/2;
	  icofm[midres].x=gridsize/2;
	  icofm[midres].y=gridsize/2;
	  icofm[midres].z=gridsize/2;
	  icofm[midres+1].x=gridsize/2+2;
	  icofm[midres+1].y=gridsize/2+2;
	  icofm[midres+1].z=gridsize/2+2;
	  setGrid(grid,icofm[midres].x,icofm[midres].y,icofm[midres].z);
	  setGrid(grid,icofm[midres+1].x,icofm[midres+1].y,icofm[midres+1].z);
	  addrandomchain(midres,resmin[ires],icofm,grid);
	  addrandomchain(midres+1,resmax[ires],icofm,grid);
	} else if (!valid[resmin[ires]-1]) {
	  if (!valid[resmax[ires]+1]) {
	    fprintf(stderr,"ERROR! Gapped fragments requested\n");
	    exit(1);
	  } 
	  addrandomchain(resmax[ires]+1,resmin[ires],icofm,grid);
	} else if (valid[resmax[ires]+1]==0) {
	  addrandomchain(resmin[ires]-1,resmax[ires],icofm,grid);
	} else {
	  sx=sy=sz=0;
	  ns=0;
	  for (i=0; i<MAXRESIDUES; i++) {
	    if (valid[i]) {
	      sx+=icofm[i].x;
	      sy+=icofm[i].y;
	      sz+=icofm[i].z;
	      ns++;
	    }
	  }
	  sx/=ns;
	  sy/=ns;
	  sz/=ns;

	  const int maxwalk=5000;
	  IVector *walk=new IVector[maxwalk];
	  int nwalk;
	  
	  int cnt=0;
	  do {
	    nwalk=randomwalk(icofm,grid,resmin[ires],resmax[ires],walk,maxwalk);
	    initialdist(icofm,walk,nwalk,resmin[ires],resmax[ires]);
	  } while (!adjust(icofm,grid,resmin[ires],resmax[ires],sx,sy,sz,valid) 
		   && ++cnt<20);
	  if (cnt>=20) {
	    throw UnSuccessful();
	  } else 
	    for (i=resmin[ires]; i<=resmax[ires]; i++)
	      setGrid(grid,icofm[i].x,icofm[i].y,icofm[i].z);
	  delete[] walk;
	}

	ns=0;
	sx=sy=sz=0;
	for (i=0; i<MAXRESIDUES; i++) {
	  if (valid[i] || (i>=resmin[ires] && i<=resmax[ires])) {
	    sx+=icofm[i].x;
	    sy+=icofm[i].y;
	    sz+=icofm[i].z;
	    ns++;
	  }
	}
	sx/=ns;
	sy/=ns;
	sz/=ns;
	sx-=gridsize/2;
	sy-=gridsize/2;
	sz-=gridsize/2;

        for (i=0; i<MAXRESIDUES; i++) {
	  if (valid[i] || (i>=resmin[ires] && i<=resmax[ires])) {
	    int tx=icofm[i].x-sx;
	    int ty=icofm[i].y-sy;
	    int tz=icofm[i].z-sz;
	    if (tx<5 || ty<5 || tz<5 ||
		tx>=gridsize-5 || ty>=gridsize-5 || tz>=gridsize-5) {
	      throw UnSuccessful();
	    }
	  }
	}
      } catch (UnSuccessful&) {
	for (i=0; i<gridsize*gridsize*gridsize; i++)
	  grid[i]=tgrid[i];
	if (++redocnt>100) {
	  printf("Sorry. Cannot generate chain. Dimensions too tight?\n");
	  exit(1);
	}
	redo=1;
      }
    } while (redo);
    for (i=resmin[ires]; i<=resmax[ires]; valid[i++]=1);
  }

  for (i=0; i<MAXRESIDUES && !valid[i]; i++);
  int minres=i;
  for (; i<MAXRESIDUES && valid[i]; i++);
  int maxres=i-1;
  for (; i<MAXRESIDUES; i++) 
    if (valid[i]) { 
      fprintf(stderr,"ERROR! Structure has gaps\n");
      exit(1);
    }

  srand48((long)0);
  redocnt=0;
  do {
    redo=0;
    try {
      addbead(minres-1,grid,icofm,minres,minres+1);
    } catch (UnSuccessful&) {
      redo=1;
      if (++redocnt>100) {
	fprintf(stderr,"cannot add first bead\n");
	exit(1);
      }
    }
  } while(redo);

  redo=redocnt=0;
  do {
    redo=0;
    try {
      addbead(maxres+1,grid,icofm,maxres,maxres-1);
    } catch (UnSuccessful&) {
      redo=1;
      if (++redocnt>100) {
	fprintf(stderr,"cannot add last bead\n");
	exit(1);
      }
    }
  } while(redo);
  

  if (finalcheck(icofm,minres-1,maxres+1)) {
    fprintf(stderr,"chain is invalid\n");
    exit(1);
  }

  fprintf(stdout,"%5d\n",maxres-minres+3);
  for (i=minres-1; i<=maxres+1; i++) {
    fprintf(stdout,"%5d%5d%5d\n",
	   icofm[i].x,icofm[i].y,icofm[i].z);
  }

  delete grid;
  delete tgrid;
  delete valid;
  delete resmin;
  delete resmax;

  return 0;
}

