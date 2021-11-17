#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

#include "physconst.h"
#include "random.h"

#include "vecmat.h"
#include "solver.h"
#include "leastsquare.h"

enum { NDATA = 28, NPARAM = 3 };

static double repdistdata[NDATA] = {
  0.1012,
  0.076600000000000001,
  0.1406,
  0.1268,
  0.152,
  0.17860000000000001,
  0.14500000000000002,
  0.14580000000000001,
  0.20900000000000002,
  0.2054,
  0.17519999999999999,
  0.17500000000000002,
  0.193,
  0.28360000000000002,
  0.23280000000000001,
  0.24680000000000002,
  0.20960000000000001,
  0.17780000000000001,
  0.2142,
  0.26240000000000002,
  0.2384,
  0.18860000000000002,
  0.1734,
  0.2772,
  0.4572,
  0.62080000000000002,
  0.70140000000000002,
  0.42800000000000005};

const double highestene = 0.020086056243873936, 
  lowestene = -0.02053671650329171;

// initial parameters for model
const double oldavrpre = 1.0, oldavrbeta = 3383.284140381932,
  oldavrq = 2.1;

const double avrpopwidth = ( - 1.0 + pow( 1.0 - oldavrbeta*(1.0-oldavrq)*
					( highestene - lowestene ), 
					(2.0-oldavrq)/(1.0-oldavrq)) )/ 
double(NDATA);


class TestData {
private:
  // True Parameters for Model
  //double avrpre, avrbeta, avrq;

  double data[NDATA];

public:
  TestData() { for (int i=0; i<NDATA; i++) data[i] = repdistdata[i]; }
  int dataNum() const { return NDATA; }
  double operator() (int i) const { return data[i]; }
  double sigma(int i) const { return 1.0; }

  void printParam() const { }
};

class TestModel {
private:
  enum { iavrpre = 0, iavrbeta = 1, iavrq = 2 };
  Vec<double> tmpparam, param;

  double model[NDATA], modelgrad[NDATA][NPARAM];

  bool dampingflag;

  bool areParamValid() const {
    if ( param(iavrpre) <= 0.0 ) return false;
    if ( param(iavrbeta) < 1.0/(KB*1000.0) ) return false;
    if ( param(iavrq) <= 2.0 || param(iavrq) >= 200.0 ) return false;
    return true;
  }

  void update() {
    double avrpre = param(iavrpre), avrbeta = param(iavrbeta),
      avrq = param(iavrq);

    for (int i=0; i<NDATA; i++) {
      double oldt0 = pow( 1.0 + avrpopwidth*i, (1.0-oldavrq)/(2.0-oldavrq) ),
	oldt1 = pow( 1.0 + avrpopwidth*(i+1), (1.0-oldavrq)/(2.0-oldavrq) ),

	ene0 = (1.0-oldt0)/(oldavrbeta*(1.0-oldavrq)),
	ene1 = (1.0-oldt1)/(oldavrbeta*(1.0-oldavrq)),
	
	t0 = 1.0 - (1.0-avrq)*avrbeta*ene0,
	t1 = 1.0 - (1.0-avrq)*avrbeta*ene1,
	deltapow = pow(t1,(2.0-avrq)/(1.0-avrq)) - 
	pow(t0,(2.0-avrq)/(1.0-avrq));

      model[i] = - avrpre / ( avrbeta*(2.0-avrq) ) * deltapow;

      modelgrad[i][iavrpre] = deltapow/(avrbeta*(-2.0 + avrq));
      modelgrad[i][iavrbeta] =  (avrpre*(- deltapow - 
       avrbeta*(-2.0 + avrq)*
        (ene0*pow(t0, 1.0/(1.0 - avrq)) - ene1*pow(t1, 1.0/(1.0 - avrq)))))/
	(avrbeta*avrbeta*(-2.0 + avrq));
      modelgrad[i][iavrq] = (avrpre*(-deltapow - 
       ((-2.0 + avrq)*pow(t0, 1.0/(1.0 - avrq))*
          (avrbeta*(2.0 - 3.0*avrq + avrq*avrq)*ene0 + 
            t0*log(t0)))/
        (-1.0 + avrq)*(-1.0 + avrq) + 
       ((-2.0 + avrq)*pow(t1, 1.0/(1.0 - avrq))*
          (avrbeta*(2.0 - 3.0*avrq + avrq*avrq)*ene1 + 
            t1*log(t1)))/
        (-1.0 + avrq)*(-1.0 + avrq)))/
	(avrbeta*(-2.0 + avrq)*(-2.0 + avrq));
    }
  }

public:
  // Just allocating memory, DON'T FORGET TO CALL init() 
  // before improving param.
  TestModel() : tmpparam(NPARAM), param(NPARAM), dampingflag(false) {
    param(iavrpre) = oldavrpre;
    param(iavrbeta) = 1.0/(KB*3.0); //oldavrbeta;
    param(iavrq) = 5.0; //oldavrq;
  }
  void init() { update(); }

  int dataNum() const { return NDATA; }
  int paramNum() const { return NPARAM; }

  double operator() (int i) const { return model[i]; }
  double grad(int idat, int iparam) const { return modelgrad[idat][iparam]; }
  bool paramAreDamped() const { return dampingflag; }

  void increaseParamBy(const Vec<double>& deltaparam) { 
    tmpparam = param;
    param += deltaparam; 

    if ( ! areParamValid() ) {
      dampingflag = true;
      param = tmpparam;
      double damping = 1.0, dampsum = 0.0;
      for (int i=0; i<70; i++) {
	damping *= 0.5;
	param += deltaparam * damping;
	if ( ! areParamValid() ) {
	  param = tmpparam;
	} else {
	  tmpparam = param;
	  dampsum += damping;
	}
      }
    } else {
      dampingflag = false;
    }
     
    update();
  }

  void printParam() const { 
    cout << "# Improved Param : " << param(iavrpre) << " " << 
      param(iavrbeta) << " " << param(iavrq) << endl;
  }
};


int main(int argc, char * argv[]) {
  TestData data;
  TestModel model;
  model.init();

  const double initiallambda = 1.0e-2;
  LevenbergMarquardtFitter fitter(data.dataNum(), model.paramNum(),
				  initiallambda);

  const double smallestfreq = 1.0e-60;
  fitter.update(data,model);
  fitter.improveParam(model,smallestfreq);
  double prechisqr = 0.0, chisqr =fitter.updateChiSqr(data,model); 
  cout << "# chisqr " << chisqr << endl;
  for (int iter=0; fabs(chisqr - prechisqr) > 1.0e-14 && iter<1000; iter++) {
    fitter.update(data,model);
    fitter.improveParam(model,smallestfreq);
    prechisqr = chisqr;
    chisqr = fitter.updateChiSqr(data,model);
    cout << "# chisqr " << chisqr << endl;
  }

  cout << endl;

  data.printParam();
  model.printParam();
  cout << endl;

  cout << "# coordinate   data      model" << endl;
  for (int i=0; i<data.dataNum(); i++) 
    cout << i << " " << data(i) << " " << model(i) << endl;

  return 0;
}


