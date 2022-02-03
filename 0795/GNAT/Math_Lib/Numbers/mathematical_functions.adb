with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics;                     use Ada.Numerics;

package body Mathematical_Functions is

  package Double_Elementary_Functions is new
                  Ada.Numerics.Generic_Elementary_Functions (double_float);

  function "**" ( x,y : double_float ) return double_float is
  begin
    return Double_Elementary_Functions."**"(x,y);
  end "**";

  function LOG10 ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.LOG(x,10.0);
  end LOG10;

  function SQRT ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.SQRT(x);
  end SQRT;

  function SIN ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.SIN(x);
  end SIN;

  function COS ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.COS(x);
  end COS;

  function TAN ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.TAN(x);
  end TAN;

  function ARCSIN ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.ARCSIN(x);
  end ARCSIN;

  function ARCCOS ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.ARCCOS(x);
  end ARCCOS;

  function ARCTAN ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.ARCTAN(x);
  end ARCTAN;

  function Radius ( x,y : double_float ) return double_float is
  begin
    return Sqrt(x**2+y**2);
  end Radius;

  function Angle ( x,y : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.ARCTAN(x,y);
  end Angle;

end Mathematical_Functions;
