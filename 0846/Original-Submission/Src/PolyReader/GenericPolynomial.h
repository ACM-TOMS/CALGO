//GenericPolynomial.h
//Author: Xing Li
//Data:  02/20/2000


#ifndef INCLUDE_GENERICPOLYNOMIAL_H
#define INCLUDE_GENERICPOLYNOMIAL_H

#include "PolynomialException.h"
#include "ExpandedPolynomial.h"

#include <string>
#include <vector>

using namespace std;

class GenericPolynomial {
 public:
  GenericPolynomial();  //construct a zero polynomial
  GenericPolynomial(string);
  GenericPolynomial(const GenericPolynomial&);
  GenericPolynomial& operator=(const GenericPolynomial&);
  ~GenericPolynomial();
  
  ExpandedPolynomial expand();
  
  void print();
  
 private:
  struct Node;
  struct Token;
  enum NodeType{VAR, LEFT_PAREN, NUMBER, RIGHT_PAREN, COEF,
		 COMPLEX, ADD, MINUS, MULTI, COMMA, EXP, DIV};
  
  static Node* createTree(vector<Token>::iterator first, 
			  vector<Token>::iterator last, vector<string>&, vector<Node*>&)
    throw (PolynomialException);
  
  static void deleteTree(Node* aNode);
  static void printTree(const Node*, const vector<string>&);
  
  static int checkParenthesis(vector<Token>::iterator first, vector<Token>::iterator last);

  static vector<Token>::iterator partition(vector<Token>::iterator first, 
					   vector<Token>::iterator last);

  static vector<Token>::iterator findExp(vector<Token>::iterator first, 
					 vector<Token>::iterator last);

  static double toDouble(string::iterator first, string::iterator last);

  static int toInteger(string::iterator first, string::iterator last);

  static bool isVariable(string::iterator first, string::iterator last);
  static Node* copyTree(Node* const );
  static void tokenize(string poly, vector<Token>& tokens);
  static ExpandedPolynomial expand(Node* );
  static bool isop(char c)
  { return  c == '+' || c =='-' || c == '*' || c =='/' || c =='^' ; }
  
 private:
  struct Token
  {
    Token(NodeType t, string s1=string(""), string s2=string("") ): 
      m_type(t), m_str1(s1), m_str2(s2) {}
    
    NodeType m_type;
    string m_str1;
    string m_str2;
  };
  
  struct Node{
    NodeType m_type;
    union {
      Node* m_left;
      Node* m_base;
      double m_real;
      int m_var;
    };
    union {
      Node* m_right;
      int   m_exp;
      double m_imag;
    };
  };
  
  Node* m_root;
  vector<string> m_variable;
};

#endif // INCLUDE_GENERICPOLYNOMIAL_H



