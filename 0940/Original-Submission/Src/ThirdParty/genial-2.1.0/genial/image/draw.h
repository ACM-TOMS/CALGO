
//GENIAL - GENeric Image & Array Library
//Copyright (C) 2005  IENT - RWTH Aachen
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#ifndef DRAW_H
#define DRAW_H

#include "functional.h"
#include "image.h"

//namespace genial
//{


template<class Image>
class image_drawing
{
  public:
    typedef Image image_type;
    typedef typename image_type::value_type value_type;
    typedef typename image_type::int_type   int_type;
    typedef typename image_type::index_type index_type;
    
    typedef int_type argument_type;

  protected:
    mutable image_type &X;

  public:
    inline explicit image_drawing(image_type &x) : X(x) {}
    
    inline image_type       &image()       { return X; }
    inline const image_type &image() const { return X; }
    
    inline argument_type left  () const { return image().col_lower_bound(); }
    inline argument_type right () const { return image().col_upper_bound(); }
    inline argument_type top   () const { return image().row_lower_bound(); }
    inline argument_type bottom() const { return image().row_upper_bound(); }
};


template<class Image>
class value_drawing : public image_drawing<Image>
{
  public:
    typedef value_drawing self;
    typedef image_drawing<Image> base;
    
    typedef typename base::image_type    image_type;
    typedef typename base::value_type    value_type;
    typedef typename base::int_type      int_type;
    typedef typename base::index_type    index_type;
    typedef typename base::argument_type argument_type;

  protected:
    value_type val;

  public:
    inline explicit value_drawing(image_type &x, const value_type &v=value_type()) : base(x), val(v) {}
    
    using base::image;
    
    inline value_type value() const { return val; }
};

template<class Image>
class point_drawing : public value_drawing<Image>
{
  public:
    typedef point_drawing self;
    typedef value_drawing<Image> base;

    typedef typename base::image_type    image_type;
    typedef typename base::value_type    value_type;
    typedef typename base::int_type      int_type;
    typedef typename base::index_type    index_type;
    typedef typename base::argument_type argument_type;

  public:
    inline explicit point_drawing(image_type &x, const value_type &v=value_type()) : base(x,v) {}
    
    using base::X;
    using base::image;
    using base::value;

    inline void operator()(const index_type &p) const { if (image().withinbounds(p)) X[p]=value(); }
    inline void draw      (const index_type &p) const { (*this)(p); }
    
    inline void operator()(argument_type x, argument_type y) const { (*this)(index_type(y,x)); }
    inline void draw      (argument_type x, argument_type y) const { (*this)(x,y); }
    
    inline const self &line(int_type x1, int_type y1, int_type x2, int_type y2) const;
    inline const self &line(const index_type &p1, const index_type &p2) const { return line(p1.x,p1.y,p2.x,p2.y); }    
};

template<class Image>
const typename point_drawing<Image>::self &point_drawing<Image>::line(int_type x1, int_type y1, int_type x2, int_type y2) const
{
  typedef PROMOTE2(float,int_type) float_type;

  int n = __max(abs(y2-y1),abs(x2-x1));
  if (n==0) { draw(x1,y1); return *this; }

  float_type dx = (float_type(x2)-float_type(x1))/n;
  float_type dy = (float_type(y2)-float_type(y1))/n;
          
  for (int i=0; i<=n; ++i)
    draw(x1+round(i*dx),y1+round(i*dy));
    
  return *this;
}

static boolImage point3_brush(3,3,"0 1 0  1 1 1  0 1 0");
static boolImage point5_brush(5,5,"0 1 1 1 0  1 1 1 1 1  1 1 1 1 1  1 1 1 1 1  0 1 1 1 0");
static boolImage point7_brush(7,7,"0 0 1 1 1 0 0  0 1 1 1 1 1 0  1 1 1 1 1 1 1  1 1 1 1 1 1 1  1 1 1 1 1 1 1  0 1 1 1 1 1 0  0 0 1 1 1 0 0");

template<class Image,class Brush=boolImage >
class brush_drawing : public point_drawing<Image>
{
  public:
    typedef brush_drawing self;
    typedef point_drawing<Image> base;
    typedef Brush brush_type;

    typedef typename base::image_type    image_type;
    typedef typename base::value_type    value_type;
    typedef typename base::int_type      int_type;
    typedef typename base::index_type    index_type;
    typedef typename base::argument_type argument_type;
    
  protected:
    brush_type mask;
  
  public:
    inline explicit brush_drawing(image_type &x, const brush_type &brush, const value_type &v=value_type()) : base(x,v), mask(brush) {}
    
    using base::image;
    using base::value;

    inline void operator()(argument_type x, argument_type y) const
    {
      int m=mask.nrows(), n=mask.ncols();
      int m2=m/2, n2=n/2; 
      for (int i=0; i<m; ++i)
        for (int j=0; j<n; ++j)
          if (mask(i,j)) base::draw(x+i-m2,y+j-n2);
    }
    inline void draw      (argument_type x, argument_type y) const { (*this)(x,y); }
    
    inline void operator()(const index_type &p) const { (*this)(p.x,p.y); }
    inline void draw      (const index_type &p) const { (*this)(p); }    
};

template<class Draw>
struct drawing_traits
{
  typedef Draw drawing_type;
  typedef typename drawing_type::image_type image_type;
  typedef typename drawing_type::int_type   int_type;
  typedef typename drawing_type::index_type index_type;
  typedef typename drawing_type::argument_type argument_type;
};

template<class Draw>
class drawing_drawing : public drawing_traits<Draw>
{
  public:
    typedef drawing_drawing self;
    typedef drawing_traits<Draw> base;

    typedef typename base::drawing_type  drawing_type;
    typedef typename base::image_type    image_type;
    typedef typename base::int_type      int_type;
    typedef typename base::index_type    index_type;
    typedef typename base::argument_type argument_type;    

  protected:
    drawing_type dr;
    
  public:
    inline drawing_drawing(const self &x) : dr(x.dr) {}
 
    inline drawing_drawing() : dr() {}
    template<class A> inline drawing_drawing(A       &a) : dr(a) {}
    template<class A> inline drawing_drawing(const A &a) : dr(a) {}
    template<class A,class B> inline drawing_drawing(A       &a,const B &b) : dr(a,b) {}
    template<class A,class B> inline drawing_drawing(const A &a,const B &b) : dr(a,b) {}

    inline drawing_type       &drawing()       { return dr; }
    inline const drawing_type &drawing() const { return dr; }
    
    inline image_type       &image()       { return drawing().image(); }
    inline const image_type &image() const { return drawing().image(); }
    
    inline argument_type left  () const { return drawing().left  (); }
    inline argument_type right () const { return drawing().right (); }
    inline argument_type top   () const { return drawing().top   (); }
    inline argument_type bottom() const { return drawing().bottom(); }
};

template<class Draw, class V=float>
class scale_drawing : public drawing_drawing<Draw>
{
  public:
    typedef scale_drawing self;
    typedef drawing_drawing<Draw> base;
    
    typedef typename base::drawing_type  drawing_type;
    typedef typename base::image_type    image_type;
    typedef typename base::int_type      int_type;
    typedef typename base::index_type    index_type;
    
    typedef V argument_type;
      
  private:
    argument_type x0,x1;
    argument_type y0,y1;
    
  public:
    inline scale_drawing() : base() {}

    template<class A> inline scale_drawing(A &a, const argument_type &xmin, const argument_type &xmax, const argument_type &ymin, const argument_type &ymax) : base(a), x0(xmin), x1(xmax), y0(ymin), y1(ymax) {}
    template<class A,class B> inline scale_drawing(A &a, const B &b, const argument_type &xmin, const argument_type &xmax, const argument_type &ymin, const argument_type &ymax) : base(a,b), x0(xmin), x1(xmax), y0(ymin), y1(ymax) {}
    template<class A> inline scale_drawing(const A &a, const argument_type &xmin, const argument_type &xmax, const argument_type &ymin, const argument_type &ymax) : base(a), x0(xmin), x1(xmax), y0(ymin), y1(ymax) {}
    template<class A,class B> inline scale_drawing(const A &a, const B &b, const argument_type &xmin, const argument_type &xmax, const argument_type &ymin, const argument_type &ymax) : base(a,b), x0(xmin), x1(xmax), y0(ymin), y1(ymax) {}

    using base::image;
    using base::drawing;

    inline argument_type left  () const { return x0; }
    inline argument_type right () const { return x1; }
    inline argument_type top   () const { return y0; }
    inline argument_type bottom() const { return y1; }

    inline int_type x_index(argument_type x) const { return (x-x0)/(x1-x0)*image().ncols(); }
    inline int_type y_index(argument_type y) const { return (y-y0)/(y1-y0)*image().nrows(); }
    
    inline self &operator()(argument_type x, argument_type y) const { drawing()(x_index(x),y_index(y)); return *this; }
    inline self &draw      (argument_type x, argument_type y) const { drawing().draw(x_index(x),y_index(y)); return *this; }

    void line(int_type x1, int_type y1, int_type x2, int_type y2) const { line(x_index(x1),y_index(y1),x_index(x2),y_index(y2)); }
};


template<class Draw>
class line_drawing : public drawing_drawing<Draw>
{
  public:
    typedef line_drawing self;
    typedef drawing_drawing<Draw> base;
    
    typedef typename base::drawing_type  drawing_type;
    typedef typename base::image_type    image_type;
    typedef typename base::int_type      int_type;
    typedef typename base::index_type    index_type;
    typedef typename base::argument_type argument_type;        
          
  private:
    argument_type x0,y0;
    bool prev;
    
  public:
    inline line_drawing() : base(), prev(false) {}

    template<class A> inline line_drawing(A       &a) : base(a), prev(false) {}
    template<class A> inline line_drawing(const A &a) : base(a), prev(false) {}
    template<class A,class B> inline line_drawing(A       &a, const B &b) : base(a,b), prev(false) {}
    template<class A,class B> inline line_drawing(const A &a, const B &b) : base(a,b), prev(false) {}
    
    using base::image;
    using base::drawing;
    
    inline self &operator()(argument_type x, argument_type y) const { if (prev) drawing().line(x0,y0,x,y); x0=x; y0=y; prev=true; return *this; }
    inline self &draw      (argument_type x, argument_type y) const { return (*this)(x,y); }
};


template<class G>
void draw_point(Matrix<G> &X, const typename Matrix<G>::index_type &p, const typename Matrix<G>::value_type &val)
{
  point_drawing<Matrix<G> >(X,val)(p);
}

template<class G>
void draw_point(Matrix<G> &X, const typename Matrix<G>::int_type x, const typename Matrix<G>::int_type y, const typename Matrix<G>::value_type &val)
{
  point_drawing<Matrix<G> >(X,val)(x,y);
}

template<class G>
void draw_line(Matrix<G> &X, const typename Matrix<G>::index_type &p1, const typename Matrix<G>::index_type &p2, const typename Matrix<G>::value_type &val)
{
  point_drawing<Matrix<G> >(X,val).line(p1,p2);
}

template<class G>
void draw_line(Matrix<G> &X, const typename Matrix<G>::int_type &x1, const typename Matrix<G>::int_type &y1, const typename Matrix<G>::int_type &x2, const typename Matrix<G>::int_type &y2, const typename Matrix<G>::value_type &val)
{
  draw_line(X,typename Matrix<G>::index_type(y1,x1),typename Matrix<G>::index_type(y2,x2));
}

template<class G, class Brush>
void draw_brush(Matrix<G> &X, const typename Matrix<G>::index_type &p, const Brush &brush, const typename Matrix<G>::value_type &val)
{
  brush_drawing<Matrix<G>,Brush>(X,brush,val)(p);
}

template<class G, class Brush>
void draw_brush(Matrix<G> &X, const typename Matrix<G>::int_type x, const typename Matrix<G>::int_type y, const Brush &brush, const typename Matrix<G>::value_type &val)
{
  brush_drawing<Matrix<G>,Brush>(X,brush,val)(x,y);
}


template<class Draw, class G>
typename Draw::image_type &plot(const Draw &draw, const Vector<G> &F)
{
  for (typename Vector<G>::int_type i=F.lower_bound(), imax=F.upper_bound(); i<=imax; ++i)
    draw(i,F[i]);
  return draw.image();
}

template<class Draw, class Func>
typename Draw::image_type &plot(const Draw &draw, const Func &f)
{
  for (typename Draw::argument_type x=draw.left(), xmax=draw.right(), dx=(draw.right()-draw.left())/draw.image().ncols(); x<=xmax; x+=dx)
    draw(x,f(x));
  return draw.image();
}


//}

#endif


