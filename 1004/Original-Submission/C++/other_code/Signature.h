#ifndef JR_SIGNATURE_H
#define JR_SIGNATURE_H

#include<array>

typedef float SignatureNumeric;

//Signature<D,M> represents the first M levels in a signature of a D dimensional path. 
//The only data stored is the actual numbers.
//It can be constructed from a displacement (representing a straight path)
//or two such signatures (Chen's formula). 
//Other classes in this file are for internal use only.
//Note that large signatures will be large objects (so don't stick them on the stack) 
//and very large signatures will cause code bloat.

template<int D, int M>
struct Signature{
	typedef Signature<D,M> Self;
	typedef Signature<D,M-1> Parent;
	static const int length = D * Parent::length;
	static const int totalLength = length+Parent::totalLength;
	Parent m_parent;
	std::array<SignatureNumeric,length> m_data;

	//construct from a single displacement
	Signature<D,M>(const std::array<SignatureNumeric,D>& a)
	  : m_parent(a){
		for(int d=0, j=0; d<D; ++d){
			for(int i=0; i!=Parent::length; ++i)
				m_data[j++]=m_parent.m_data[i]*a[d]*(1.0f/M);
		}
	}

	const Signature<D,1>& root() const{return m_parent.root();}

	//construct blank
	Signature<D,M>():m_data(){}

	//construct from two, concatenated - using Chen's identity
	Signature<D,M>(const Self& a, const Self& b);

	//a.postConcatenateWith(b)  is the same as a=Self(a,b)
	void postConcatenateWith(const Self& b);

	////a.postConcatenateWith(b)  is the same as a.postConcatenateWith(Self(b))
	//void postConcatenateWith(const std::array<SignatureNumeric,D>& a);

	void print()
	{
		m_parent.print();
		for(auto a : m_data)
			std::cout<<a<<",";
		std::cout<<std::endl;
	}
	template<class T>
	void iterateOver(T&& t) const{
		m_parent.iterateOver(t);
		for(auto a:m_data)
			t(a);
	}
	template<class T>
	void iterateOver_mutably(T&& t) {
		m_parent.iterateOver_mutably(t);
		for(auto& a:m_data)
			t(a);
	}
	template<class T, class S>
	void iterateOverWithSignal(T&& t, S&& signal) const{
		m_parent.iterateOverWithSignal(t, signal);
		signal();
		for(auto a:m_data)
			t(a);
	}

};

template<int D>
struct Signature<D,1>{
	typedef Signature<D,1> Self;
	static const int length=D;
	static const int totalLength=D;
	static const int lengthSummingOverPermutations=D;
	std::array<SignatureNumeric,D> m_data;

	//construct from a single displacement
	Signature<D,1>(const std::array<SignatureNumeric,D>& a):m_data(a){}
	const Self& root() const{return *this;}

	//construct from two, concatenated
	Signature<D,1>(const Self& a, const Self& b){
		for(int d=0;d<D;++d)
			m_data[d]=a.m_data[d]+b.m_data[d];
	}

	//construct blank
	Signature<D,1>():m_data(){}

	//a.postConcatenateWith(b)  is the same as a=Self(a,b)
	void postConcatenateWith(const Self& b){
		for(int d=0;d<D;++d)
			m_data[d]+=b.m_data[d];
	}

	////a.postConcatenateWith(b)  is the same as a.postConcatenateWith(Self(b))
	//void postConcatenateWith(const std::array<SignatureNumeric,D>& a){
	//	for(int d=0;d<D;++d)
	//		m_data[d]+=a[d];
	//}

	void print()
	{
		for(auto a : m_data)
			std::cout<<a<<",";
		std::cout<<std::endl;
	}
	template<class T>
	void iterateOver(T&& t) const{
		for(auto a:m_data)
			t(a);
	}
	template<class T>
	void iterateOver_mutably(T&& t) {
		for(auto& a:m_data)
			t(a);
	}
	template<class T, class S>
	void iterateOverWithSignal(T&& t, S&& signal) const{
		for(auto a:m_data)
			t(a);
	}
	
};

//BEGIN IMPLEMENTATION OF Self(Self,Self)
namespace{
	template<int D, int M, int MM, int Excess>
	struct Adder{
		static void addOnTo(Signature<D,MM>& target, const Signature<D,M>& me, const Signature<D,MM-M+Excess>& other) {
			Adder<D, M, MM,Excess-1>::addOnTo(target, me, other.m_parent);
		}
	};

	template<int D, int M, int MM>
	struct Adder<D, M, MM, 0>{
		static void addOnTo(Signature<D,MM>& target, const Signature<D,M>& me, const Signature<D,MM-M>& other) {
			for(int i=0, k=0; i!=me.length; ++i)
				for(int j=0; j!=other.length; ++j)
					target.m_data[k++]+=me.m_data[i]*other.m_data[j];
		}
	};

	template<int D, int MM, int Excess>
	struct Adder<D, 1, MM, Excess>{
		static void addOnTo(Signature<D,MM>& target, const Signature<D,1>& me, const Signature<D,MM-1+Excess>& other) {
			Adder<D, 1, MM,Excess-1>::addOnTo(target, me, other.m_parent);
		}
	};

	template<int D, int MM>
	struct Adder<D, 1, MM, 0>{
		static void addOnTo(Signature<D,MM>& target, const Signature<D,1>& me, const Signature<D,MM-1>& other) {
			for(int i=0, k=0; i!=D; ++i)
				for(int j=0; j!=other.length; ++j)
					target.m_data[k++]+=me.m_data[i]*other.m_data[j];
		}
	};

	template<int D, int M, int MM>
	struct ManyAdder{
		static void manyAddings(Signature<D,M>& target, const Signature<D,MM>& me, const Signature<D,M>& other){
			Adder<D, MM, M,MM>::addOnTo(target,me,other);
			ManyAdder<D, M, MM-1>::manyAddings(target, me.m_parent, other);
		}
	};
	template<int D, int M>
	struct ManyAdder<D, M, 1>{
		static void manyAddings(Signature<D,M>& target, const Signature<D,1>& me, const Signature<D,M>& other){
			Adder<D,1,M,1>::addOnTo(target,me,other);
		}
	};
}//end of anonymous namespace

template<int D, int M>
Signature<D,M>::Signature(const Self& a, const Self& b)
	: m_parent(a.m_parent,b.m_parent){
		for(int i=0; i!=length; ++i)
			m_data[i]=a.m_data[i]+b.m_data[i];
		//Parent::Adder<M,M-1>::addOnTo(*this,a.m_parent,b);
		ManyAdder<D,M,M-1>::manyAddings(*this,a.m_parent,b);
}
//END IMPLEMENTATION OF Self(Self,Self)

//BEGIN IMPLEMENTATION OF void postConcatenateWith(Self)
//use restrict?
template<int D, int M>
void Signature<D,M>::postConcatenateWith(const Self& b){
	for(int i=0; i!=length; ++i)
		m_data[i]+=b.m_data[i];
	ManyAdder<D,M,M-1>::manyAddings(*this, m_parent, b);//cant use restrict!
	m_parent.postConcatenateWith(b.m_parent);
}

//template<int D, int M>
//void Signature<D,M>::postConcatenateWith(const std::array<SignatureNumeric,D>& a){
//	m_parent.postConcatenateWith(a);
//}


//END IMPLEMENTATION OF void postConcatenateWith(Self)


#endif
