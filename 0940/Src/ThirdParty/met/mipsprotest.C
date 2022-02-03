template <class A> class Special {};

class Container {
public:
  Container() {}
  template <class A> Container& operator=(const Special<A>& rhs) { return *this;}
  template <class Any> Container& operator=(const Any& rhs) { return *this; }
};

int main(int argc, char * argv[]) {
  Special<double> sd;
  Container c;
  
  c = sd;

  return 0;
}
