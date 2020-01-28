// Michael Feig, 1998

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>

// ### Vector ######

class Vector {
 private:
  double xval,yval,zval;

  double anint(double x) {
    if (x>=0) {
     return floor(x+0.5);
    } else {
     return -floor(0.5-x);
    }
  }

 public:
  Vector(double xi=0.0, double yi=0.0, double zi=0.0) : 
    xval(xi), yval(yi), zval(zi) {}

  Vector(const Vector& from) {
    xval=from.xval;
    yval=from.yval;
    zval=from.zval;
  }

  inline Vector& operator=(const Vector& from) {
    xval=from.xval;
    yval=from.yval;
    zval=from.zval;
    return(*this);
  }

  inline Vector& operator=(const double v) {
    xval=v;
    yval=v;
    zval=v;
    return (*this);
  }

  inline Vector& operator+=(const Vector& from) {
    xval+=from.xval;
    yval+=from.yval;
    zval+=from.zval;
    return(*this);
  }

  inline Vector& operator-=(const Vector& from) {
    xval-=from.xval;
    yval-=from.yval;
    zval-=from.zval;
    return(*this);
  }

  inline Vector& operator*=(const Vector& from) {
    xval*=from.xval;
    yval*=from.yval;
    zval*=from.zval;
    return(*this);
  }

  inline Vector& operator*=(const double r) {
    xval*=r;
    yval*=r;
    zval*=r;
    return(*this);
  }

  inline Vector& operator/=(const double r) {
    xval/=r;
    yval/=r;
    zval/=r;
    return(*this);
  }

  inline Vector& operator+=(const double r) {
    xval+=r;
    yval+=r;
    zval+=r;
    return(*this);
  }

  inline Vector& operator-=(const double r) {
    xval-=r;
    yval-=r;
    zval-=r;
    return(*this);
  }

  inline Vector operator-() const {
    return Vector(-xval,-yval,-zval);
  }

  inline Vector operator+(const Vector& v) const {
    return Vector(xval+v.xval,yval+v.yval,zval+v.zval);
  }

  inline Vector operator-(const Vector& v) const {
    return Vector(xval-v.xval,yval-v.yval,zval-v.zval);
  }

  inline double operator*(const Vector& v) const {
    return xval*v.xval+yval*v.yval+zval*v.zval;
  }

  friend inline Vector operator*(double r, const Vector& v) {
    return Vector(r*v.xval,r*v.yval,r*v.zval);
  }

  friend inline Vector operator*(const Vector& v, double r) {
    return Vector(r*v.xval,r*v.yval,r*v.zval);
  }

  friend inline Vector operator/(const Vector& v, double r) {
    return Vector(v.xval/r,v.yval/r,v.zval/r);
  }

  inline double& x() { return xval; }
  inline double& y() { return yval; }
  inline double& z() { return zval; }

  double operator[](int inx) const {
    return (inx==0?xval:inx==1?yval:zval);
  }

  inline double norm() {
    return sqrt(xval*xval+yval*yval+zval*zval);
  }

  inline Vector cross(const Vector& v) const {
    return Vector(yval*v.zval - zval*v.yval,
		  zval*v.xval - xval*v.zval,
		  xval*v.yval - yval*v.xval);
  }

  inline Vector rotateX(double phi) const {
    return Vector(xval,
		  yval*cos(phi) + zval*sin(phi),
		 -yval*sin(phi) + zval*cos(phi));
  }
  
  inline Vector rotateY(double phi) const {
      return Vector(xval*cos(phi) - zval*sin(phi),
		    yval,
		    xval*sin(phi) + zval*cos(phi));
  }

  inline Vector rotateZ(double phi) const {
    return Vector(xval*cos(phi) + yval*sin(phi),
		 -xval*sin(phi) + yval*cos(phi),
		  zval);  
  }

  inline Vector transformB(Vector ax, Vector ay, Vector az) {
    return Vector(xval*ax.x()+yval*ay.x()+zval*az.x(),
		  xval*ax.y()+yval*ay.y()+zval*az.y(),
		  xval*ax.z()+yval*ay.z()+zval*az.z());
  }

  inline Vector transformF(Vector ax, Vector ay, Vector az) {
    return Vector((*this)*ax,(*this)*ay,(*this)*az);
  }

  inline double dot(const Vector& v) {
     return (xval*v.xval + yval*v.yval + zval*v.zval);
  }

  inline void pbc(Vector box) {
     xval-=anint(xval/box.x())*box.x();
     yval-=anint(yval/box.y())*box.y();
     zval-=anint(zval/box.z())*box.z();
  }
};

#endif /* VECTOR_H */
