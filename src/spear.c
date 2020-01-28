#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double erfcc(double x) {
  double t,z,ans;
  
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
	    t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
	    t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}

double gammln(double xx) {
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double betacf(double a, double b, double x) {
  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;

  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT) {
    fprintf(stderr,"a or b too big, or MAXIT too small in betacf\n");
    exit(1);
  }
  return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN

double betai(double a, double b, double x) {
  double bt;

  if (x < 0.0 || x > 1.0) {
    fprintf(stderr,"x: %lf invalid in betai\n",x);
    exit(1);
  }
  if (x == 0.0 || x == 1.0) bt=0.0;
  else
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}



void crank(int n, double *w, double& s) {
  unsigned long j=1,ji,jt;
  double t,rank;
  
  s=0.0;
  while (j < n) {
    if (w[j+1] != w[j]) {
      w[j]=j;
      ++j;
    } else {
      for (jt=j+1;jt<=n && w[jt]==w[j];jt++);
      rank=0.5*(j+jt-1);
      for (ji=j;ji<=(jt-1);ji++) w[ji]=rank;
      t=jt-j;
      s += t*t*t-t;
      j=jt;
    }
  }
  if (j == n) w[n]=n;
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort2(int n, double *arr, double *brr) {
  int i,ir=n,j,k,l=1;
  int *istack,jstack=0;
  double a,b,temp;

  istack=new int[NSTACK+1];
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=1;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (!jstack) {
	delete istack;
	return;
      }
      ir=istack[jstack];
      l=istack[jstack-1];
      jstack -= 2;
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
	SWAP(brr[k],brr[l+1])
	if (arr[l+1] > arr[ir]) {
	  SWAP(arr[l+1],arr[ir])
	    SWAP(brr[l+1],brr[ir])
	    }
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	  SWAP(brr[l],brr[ir])
	  }
      if (arr[l+1] > arr[l]) {
	SWAP(arr[l+1],arr[l])
	  SWAP(brr[l+1],brr[l])
	  }
      i=l+1;
      j=ir;
      a=arr[l];
      b=brr[l];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
	  SWAP(brr[i],brr[j])
	  }
      arr[l]=arr[j];
      arr[j]=a;
      brr[l]=brr[j];
      brr[j]=b;
      jstack += 2;
      if (jstack > NSTACK) {
	fprintf(stderr,"NSTACK too small in sort2\n");
	exit(1);
      }
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}
#undef M
#undef NSTACK
#undef SWAP

void spear(double *data1, double *data2, 
	   int n, double& d, double& zd, 
	   double& probd, double& rs, double& probrs) {
  int j;
  double vard,t,sg,sf,fac,en3n,en,df,aved;
  
  double *wksp1=new double[n+1];
  double *wksp2=new double[n+1];

  for (j=0;j<n;j++) {
    wksp1[j+1]=data1[j];
    wksp2[j+1]=data2[j];
  }

  sort2(n,wksp1,wksp2);
  crank(n,wksp1,sf);
  sort2(n,wksp2,wksp1);
  crank(n,wksp2,sg);

  d=0.0;
  for (j=1;j<=n;j++)
    d+=(wksp1[j]-wksp2[j])*(wksp1[j]-wksp2[j]);
  en=n;
  en3n=en*en*en-en;
  aved=en3n/6.0-(sf+sg)/12.0;
  fac=(1.0-sf/en3n)*(1.0-sg/en3n);
  vard=((en-1.0)*en*en*(en+1.0)*(en+1.0)/36.0)*fac;
  zd=(d-aved)/sqrt(vard);
  probd=erfcc(fabs(zd)/1.4142136);
  rs=(1.0-(6.0/en3n)*(d+(sf+sg)/12.0))/sqrt(fac);
  fac=(rs+1.0)*(1.0-(rs));
  if (fac > 0.0) {
    t=(rs)*sqrt((en-2.0)/fac);
    df=en-2.0;
    probrs=betai(0.5*df,0.5,df/(df+t*t));
  } else
    probrs=0.0;

  delete wksp1;
  delete wksp2;
}

#define MAXN 100000
int main(int argc, char **argv) {
  double *dx=new double[MAXN];
  double *dy=new double[MAXN];
  
  int n=0;
  while (!feof(stdin)) {
    int v=fscanf(stdin,"%lf%lf",&dx[n],&dy[n]);
    if (v>0) n++;
  }
  fclose(stdin);

  double d,zd,probd,rs,probrs;
  
  spear(dx,dy,n,d,zd,probd,rs,probrs);

  printf("%lf %lf %lf %lf %lf\n",rs,zd,probd,probrs,d);
}
