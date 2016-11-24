#ifndef FIELD_H
#define FIELD_H

class Field {
 private:
  static const int maxfields;
  
  char *tstr;
  char *fstr[300];
  
  int nfields;

  int first;

  int contains(char c, char *set);

 public:
  Field(char *str, char *sep=" \t");
  ~Field();

  void shift();

  int number() { return nfields-first; }
  
  char *operator[](int inx);
};

#endif /* FIELD_H */








