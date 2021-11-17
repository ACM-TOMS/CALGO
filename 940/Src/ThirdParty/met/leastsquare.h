//-*- C++ -*-
//
//  Lesst Square Fitting Class Declaration
//
//     Copyright (C) 2000 Masakatsu Ito
//     Institute for Molecular Science
//     Myodaiji  Okazaki  444-0867 Japan
//     E-mail mito@ims.ac.jp

//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 

#ifndef LEASTSQUARE_H_
#define LEASTSQUARE_H_

#include "vecmat.h"
#include "solver.h"

class BaseLeastSquareFitter {
protected:
  int ndata, nparam;
  double chisqr, prechisqr;
  //double *data, *model,
  //  *modelparadata, **modelpara;
  
  Vec<double> mgrad, deltaparam;
  Mat<double> hess;
  MatDiagonalized<double> diaghess;

  double *mgradelemdata, **mgradelem;

  virtual void error(const char* msg="") const {
    cerr << "BaseLeastSquareFitter error : " << msg << endl;
    exit(1);
  }

  template <class Data, class Model>
  void updateGradHessian(const Data& data, const Model& model) {
    int i,j,k;
    for (j=0; j<nparam; j++) {
      mgrad(j) = 0.0;
      for (k=0; k<nparam; k++) hess(j,k) = 0.0;
    }

    for (i=0; i<ndata; i++) {
      for (j=0; j<nparam; j++) {
	mgradelem[i][j] = 
	  - ( model(i) - data(i) ) / data.sigma(i) * model.grad(i,j);
	mgrad(j) += mgradelem[i][j];
	for (k=0; k<nparam; k++) 
	  hess(j,k) += model.grad(i,j) * model.grad(i,k) / data.sigma(i);
      }
    }
  }

  explicit BaseLeastSquareFitter(int ndata_, int nparam_) :  
    ndata(ndata_), nparam(nparam_), chisqr(0.0), prechisqr(chisqr),
    mgrad(nparam), deltaparam(nparam), hess(nparam), diaghess(hess),
    mgradelemdata(new double[ndata*nparam]), mgradelem(new double*[ndata]) {
    for (int i=0; i<ndata; i++) mgradelem[i] = mgradelemdata + i*nparam;
  }

  virtual ~BaseLeastSquareFitter() {
    delete [] mgradelem;
    delete [] mgradelemdata;
  }

public:
  template <class Data, class Model>
  double updateChiSqr(const Data& data, const Model& model){
    prechisqr = chisqr;
    chisqr = 0.0;
    
    for (int i=0; i<ndata; i++) 
      chisqr += ( model(i) - data(i) ) * ( model(i) - data(i) ) / 
	data.sigma(i);

    return chisqr;
  }
};



// ref. "Scientific and Engineerging C++", Section 19.6
//      "Numerical Recipes", Section 15.5
// Marquardt technique starts with the empirical observation that
// solution by linear approximation of nonlinear equation can lead
// to too large parameter changes. Hence the hessian diagonal elements are
// augumented with a peralty for large corrections.
//
// "Numerical Recipes" suggests beggining with lambda = 0.001 and 
// changing by factors of 10.

class LevenbergMarquardtFitter : public BaseLeastSquareFitter {
private:
  double lambda;
  const double maxlambda, lambdafactor;

public:
  LevenbergMarquardtFitter(int ndata_, int nparam_, double initlambda = 0.0001,
			   double maxlambda_ = 1.0e+3,
			   double lambdafactor_ = 10.0) :
    BaseLeastSquareFitter(ndata_, nparam_), lambda(initlambda),
    maxlambda(maxlambda_), lambdafactor(lambdafactor_) {}

  void initLambda(double initlambda) { lambda = initlambda; }

  template <class Data, class Model>
  void update(const Data& data, const Model& model) {
    if ( chisqr > prechisqr || model.paramAreDamped() ) { 
      lambda *= lambdafactor; 
      if ( lambda > maxlambda ) 
	lambda /= lambdafactor;
    } else { 
      lambda /= lambdafactor; 
    }

    updateGradHessian(data, model);

    for (int j=0; j<nparam; j++) hess(j,j) *= 1.0 + lambda;

    diaghess.update();
  }

  // DON'T FORGET to call updateChiSqr(data,model) after using this.
  template <class Model>
  void improveParam(Model& model) {
    diaghess.solve(mgrad, deltaparam);
    model.increaseParamBy(deltaparam);
  }

  // DON'T FORGET to call updateChiSqr(data,model) after using this.
  template <class Model>
  void improveParam(Model& model, 
		    double smallest_eigenvalue) {
    diaghess.solve(mgrad, deltaparam, smallest_eigenvalue);
    model.increaseParamBy(deltaparam);
  }

  template <class Data, class Model>
  double fit(const Data& data, Model& model,
	     double initiallambda = 1000.0, double deltachisqrthres = 1.0e-14,
	     double smallestfreq = 1.0e-80) {
    initLambda(initiallambda);
    update(data, model);

    improveParam(model, smallestfreq);
    chisqr = updateChiSqr(data, model);
    //cout << "# chisqr " << chisqr << endl;

    for (int iter=0; 
	 ( fabs(chisqr-prechisqr) > deltachisqrthres || 
	   model.paramAreDamped() ) && iter<2000; 
	 iter++) {
      update(data, model);
      improveParam(model, smallestfreq);
      prechisqr = chisqr;
      chisqr = updateChiSqr(data, model);
      //cout << "# chisqr " << chisqr << endl;
    }
    return chisqr;
  }
};


class BaseLeastDiffFitter {
protected:
  int ndata, nparam;
  double chisqr, prechisqr;
  //double *data, *model,
  //  *modelparadata, **modelpara;
  
  Vec<double> mgrad, deltaparam;
  Mat<double> hess;
  MatDiagonalized<double> diaghess;

  virtual void error(const char* msg="") const {
    cerr << "BaseLeastDiffFitter error : " << msg << endl;
    exit(1);
  }

  template <class Data, class Model>
  void updateGradHessian(const Data& data, const Model& model) {
    int i,j,k;
    for (j=0; j<nparam; j++) {
      mgrad(j) = 0.0;
      for (k=0; k<nparam; k++) hess(j,k) = 0.0;
    }

    for (i=0; i<ndata-1; i++) {
      for (j=0; j<nparam; j++) {
	mgrad(j) -= ( model(i+1) - data(i+1) - model(i) + data(i) ) / 
	  data.sigma(i) * ( model.grad(i+1,j) - model.grad(i,j) );
	for (k=0; k<nparam; k++) 
	  hess(j,k) += ( model.grad(i+1,j) - model.grad(i,j) ) * 
	    ( model.grad(i+1,k) - model.grad(i,k) ) / data.sigma(i);
      }
    }
  }

  explicit BaseLeastDiffFitter(int ndata_, int nparam_) :  
    ndata(ndata_), nparam(nparam_), chisqr(0.0), prechisqr(chisqr),
    mgrad(nparam), deltaparam(nparam), hess(nparam), diaghess(hess) {}

public:
  template <class Data, class Model>
  double updateChiSqr(const Data& data, const Model& model){
    prechisqr = chisqr;
    chisqr = 0.0;
    
    for (int i=0; i<ndata; i++) {
      double df = model(i+1) - data(i+1) - model(i) + data(i);
      chisqr += df * df / data.sigma(i);
    }

    return chisqr;
  }
};



// ref. "Scientific and Engineerging C++", Section 19.6
//      "Numerical Recipes", Section 15.5
// Marquardt technique starts with the empirical observation that
// solution by linear approximation of nonlinear equation can lead
// to too large parameter changes. Hence the hessian diagonal elements are
// augumented with a peralty for large corrections.
//
// "Numerical Recipes" suggests beggining with lambda = 0.001 and 
// changing by factors of 10.

class LevenbergMarquardtDiffFitter : public BaseLeastDiffFitter {
private:
  double lambda;
  const double maxlambda, lambdafactor;

public:
  LevenbergMarquardtDiffFitter(int ndata_, int nparam_, 
			       double initlambda = 0.0001,
			       double maxlambda_ = 1.0e+3,
			       double lambdafactor_ = 10.0) :
    BaseLeastDiffFitter(ndata_, nparam_), lambda(initlambda),
    maxlambda(maxlambda_), lambdafactor(lambdafactor_) {}

  void initLambda(double initlambda) { lambda = initlambda; }

  template <class Data, class Model>
  void update(const Data& data, const Model& model) {
    if ( chisqr > prechisqr || model.paramAreDamped() ) { 
      lambda *= lambdafactor; 
      if ( lambda > maxlambda ) 
	lambda /= lambdafactor;
    } else { 
      lambda /= lambdafactor; 
    }

    updateGradHessian(data, model);

    for (int j=0; j<nparam; j++) hess(j,j) *= 1.0 + lambda;

    diaghess.update();
  }

  // DON'T FORGET to call updateChiSqr(data,model) after using this.
  template <class Model>
  void improveParam(Model& model) {
    diaghess.solve(mgrad, deltaparam);
    model.increaseParamBy(deltaparam);
  }

  // DON'T FORGET to call updateChiSqr(data,model) after using this.
  template <class Model>
  void improveParam(Model& model, 
		    double smallest_eigenvalue) {
    diaghess.solve(mgrad, deltaparam, smallest_eigenvalue);
    model.increaseParamBy(deltaparam);
  }
};

#endif // LEASTSQUARE_H_
