%
% This MATLAB routine tests the accuracy of the algorithms in 
%
% Algorithm AS 241: The Percentage Points of the Normal Distribution,
% Wichura, M.J., Applied Statistics, 37(8):477-484, 1988
%

function as241

close all
clear all

%
% test main range
%

x1 = linspace(0,0.5,1000);
x1 = x1(2:end-1);

y1 = [];

for p = x1
  y1 = [y1 norminv_as241(p)];
end

%
% test tail
%

x2 = 10.^linspace(-306,-1,1000);
x2 = x2(1:end-1);

y2 = [];

for p = x2
  y2 = [y2 norminv_as241(p)];
end

%
% plot relative errors
%

figure
subplot(2,1,1)
plot(x1,(y1-norminv(x1))./y1)
title('relative error in main interval')
xlabel('p')

subplot(2,1,2)
semilogx(x2,(y2-norminv(x2))./y2)
title('relative error in tail')
xlabel('p')

end

%
% version 1: max relative error of about 1e-7
%

function res = norminv_as241_1(p)

  q = p - 0.5;
  if (abs(q) <= 0.425)
    r = 0.180625 - q*q;

    num =         5.9109374720e+01;
    num = r*num + 1.5929113202e+02;
    num = r*num + 5.0434271938e+01;
    num = r*num + 3.3871327179e+00;

    den =         6.7187563600e+01;
    den = r*den + 7.8757757664e+01;
    den = r*den + 1.7895169469e+01;
    den = r*den + 1.0000000000e+00;

    res = q * num / den;
    return

  else

    if (q < 0.0)
      r = p;
    else
      r = 1.0 - p;
    end

    r = sqrt(-log(r));

    if (r <= 5.0)
      r = r - 1.6;

      num =         1.7023821103e-01;
      num = r*num + 1.3067284816e+00;
      num = r*num + 2.7568153900e+00;
      num = r*num + 1.4234372777e+00;

      den =         1.2021132975e-01;
      den = r*den + 7.3700164250e-01;
      den = r*den + 1.0000000000e+00;

      res = num / den;

    else
      r = r - 5.0;

      num =         1.7337203997e-02;
      num = r*num + 4.2868294337e-01;
      num = r*num + 3.0812263860e+00;
      num = r*num + 6.6579051150e+00;

      den =         1.2258202635e-02;
      den = r*den + 2.4197894225e-01;
      den = r*den + 1.0000000000e+00;

      res = num / den;
    end

    if (q < 0.0)
      res = - res;
    end

    return
  end
end

%
% version 2: max relative error of less than 1e-15
%

function res = norminv_as241(p)

  q = p - 0.5;
  if (abs(q) <= 0.425)
    r = 0.180625 - q*q;

    num =         2.5090809287301226727e+3;
    num = r*num + 3.3430575583588128105e+4;
    num = r*num + 6.7265770927008700853e+4;
    num = r*num + 4.5921953931549871457e+4;
    num = r*num + 1.3731693765509461125e+4;
    num = r*num + 1.9715909503065514427e+3;
    num = r*num + 1.3314166789178437745e+2;
    num = r*num + 3.3871328727963666080e0;

    den =         5.2264952788528545610e+3;
    den = r*den + 2.8729085735721942674e+4;
    den = r*den + 3.9307895800092710610e+4;
    den = r*den + 2.1213794301586595867e+4;
    den = r*den + 5.3941960214247511077e+3;
    den = r*den + 6.8718700749205790830e+2;
    den = r*den + 4.2313330701600911252e+1;
    den = r*den + 1.0000000000e+00;

    res = q * num / den;
    return

  else

    if (q < 0.0)
      r = p;
    else
      r = 1.0 - p;
    end

    r = sqrt(-log(r));

    if (r <= 5.0)
      r = r - 1.6;

      num =         7.74545014278341407640e-4;
      num = r*num + 2.27238449892691845833e-2;
      num = r*num + 2.41780725177450611770e-1;
      num = r*num + 1.27045825245236838258e0;
      num = r*num + 3.64784832476320460504e0;
      num = r*num + 5.76949722146069140550e0;
      num = r*num + 4.63033784615654529590e0;
      num = r*num + 1.42343711074968357734e0;

      den =         1.05075007164441684324e-9;
      den = r*den + 5.47593808499534494600e-4;
      den = r*den + 1.51986665636164571966e-2;
      den = r*den + 1.48103976427480074590e-1;
      den = r*den + 6.89767334985100004550e-1;
      den = r*den + 1.67638483018380384940e0;
      den = r*den + 2.05319162663775882187e0;
      den = r*den + 1.0000000000e+00;

      res = num / den;

    else
      r = r - 5.0;

      num =         2.01033439929228813265e-7;
      num = r*num + 2.71155556874348757815e-5;
      num = r*num + 1.24266094738807843860e-3;
      num = r*num + 2.65321895265761230930e-2;
      num = r*num + 2.96560571828504891230e-1;
      num = r*num + 1.78482653991729133580e0;
      num = r*num + 5.46378491116411436990e0;
      num = r*num + 6.65790464350110377720e0;

      den =         2.04426310338993978564e-15;
      den = r*den + 1.42151175831644588870e-7;
      den = r*den + 1.84631831751005468180e-5;
      den = r*den + 7.86869131145613259100e-4;
      den = r*den + 1.48753612908506148525e-2;
      den = r*den + 1.36929880922735805310e-1;
      den = r*den + 5.99832206555887937690e-1;
      den = r*den + 1.0000000000e+00;

      res = num / den;
    end

    if (q < 0.0)
      res = - res;
    end

    return
  end
end

