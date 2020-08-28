#ifndef __posizione__
#define __posizione__

class posizione {

private:
  double x ,y, z;

protected:

public:
  // constructors
  posizione();
  posizione( double , double , double );
  // destructor
  ~posizione();
  // methods
  void singolo_passo(Random&);
  double Getx(){return x;}
  double Gety(){return y;}
  double Getz(){return z;}
  double dist2(){return x*x + y*y + z*z;}
  void Setx(double a){x=a;}
  void Sety(double a){y=a;}
  void Setz(double a){z=a;}
  void singolo_passo_cont(Random&);
  void singolo_passo_cont2(Random&);

};

#endif 