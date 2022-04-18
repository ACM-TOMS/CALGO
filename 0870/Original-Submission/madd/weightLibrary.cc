
using namespace std;

#include <math.h>

#define max(x, y) ((x) > (y)? (x): (y))


extern "C" double  weightFunction1(double x, double y)
{

  return 5;

}

extern "C" double  nohole1(double x, double y)
{
  if (y <= 2)
	return 0;

  if ((y + 6 * x - 20) < 0)
	return 0;

  if ((y - 6 * x + 16) < 0)
	return 0;

  return 5;
}


extern "C" double  nohole2(double x, double y)
{
  // the triangle area defined by (3,2), (2.5, 5), (3.5 ,5)
  if ((y + 6 * x - 20) >= 0 && (y - 6 * x + 16) >= 0)
	return 4;

  // the area up of (x, 1.5) right of (2.5, 1.5), (2,5)
  if (y >= 1.5 && (y  + 7 * x - 19) >= 0)
	return 2;

  return 0;
}

extern "C" double  nohole3(double x, double y)
{
  return (5-sqrt((x-3) * (x-3) + (y - 4) * (y - 4)));
}

extern "C" double  pipe1(double x, double y)
{
  return (53 - sqrt(x * x + y * y ));
}

extern "C" double  pipe2(double x, double y)
{
  double value = (53 - sqrt(x * x + y * y ));
  return value * value;
}


extern "C" double  keyweight(double x, double y)
{
  return (100 - y);
}


extern "C" double  pipearea1(double x, double y)
{
  double value = sqrt(x*x + y*y);
  value = value * value * value / 12000;
  return value * value * value * value;
}

extern "C" double  pipearea2(double x, double y)
{
  double value = sqrt(x*x + y*y);
  value = value * value * value / 20000;
  return value * value * value * value;
}


extern "C" double ff(double x, double y)
{

  return 35.29;

}
