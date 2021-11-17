//GenericPolynomial.cpp
//Author: Xing Li
//Data:  02/20/2000


#include "GenericPolynomial.h"
#include <iostream>
#include <iterator>
#include <algorithm>
#include <complex>

#ifdef _MSC_VER
#define for if( 0) {} else for
#endif

GenericPolynomial::GenericPolynomial():
m_root(new Node)
{
  m_root ->m_type = NUMBER;
  m_root->m_real = 0;
}

GenericPolynomial::GenericPolynomial(string poly)
{
  vector<Node*> temp;
  try{
    vector<Token> tokens;   
    tokenize(poly, tokens);
    for(vector<Token>::const_iterator it = tokens.begin(); it != tokens.end(); ++it)
      if( it->m_type == VAR )
	{
	  vector<string>::iterator it2 = m_variable.begin();
	  for( ; it2 != m_variable.end() && *it2 < it->m_str1;   ++it2)
	    ;
	  if( it2 == m_variable.end() || *it2 != it->m_str1 )
	    m_variable.insert(it2, it->m_str1);
	}
    m_root = createTree(tokens.begin(), tokens.end(), m_variable, temp);    
  } catch ( PolynomialException& e)
    {
      for(vector<Node*>::iterator it = temp.begin(); it != temp.end(); ++it)
	delete []*it;  //delocate memory;
      cout<< e.toString() << endl;
      exit(-1);
    }
  //   print();
}

GenericPolynomial::GenericPolynomial(const GenericPolynomial&  gp):
  m_variable(gp.m_variable )
{
  m_root = copyTree(gp.m_root );
}

GenericPolynomial& GenericPolynomial::operator = (const GenericPolynomial& gp)
{
  deleteTree(m_root);
  m_root = copyTree(gp.m_root);
  m_variable = gp.m_variable;
  
  return *this;
}

GenericPolynomial::~GenericPolynomial()
{
  deleteTree(m_root);
}

void GenericPolynomial::print()
{
  printTree(m_root, m_variable);
  cout<<endl;
}


GenericPolynomial::Node* 
GenericPolynomial::createTree(vector<Token>::iterator first, 
			      vector<Token>::iterator last, vector<string>& variable, vector<Node*>&  incase) 
  throw(PolynomialException)
{
  int paren = checkParenthesis(first, last);
  if( paren != 0 )
    throw PolynomialException(PolynomialException::UNPAIRED_PAREN,  ( paren > 0 )? "(" : ")" );
 
  vector<Token>::iterator it = partition(first, last);
  if(  it != last )
    {
      if(  (it == first  && ( it->m_type == MULTI  || it->m_type == DIV ) ) || ( it + 1 == last  )  )
	throw PolynomialException(PolynomialException::WRONG_PLACE_SYMBOL,  string()+ it->m_str1 );
      
      Node* temp = new Node;
      incase.push_back(temp);
      temp->m_type = it->m_type;
      temp->m_left = temp->m_right = 0;
      temp->m_left = createTree(first, it, variable, incase);
      temp->m_right = createTree(it+1, last, variable, incase);
      return temp;
    }
  
  it = findExp(first, last);
  if( it == first || it == first + (last - first -1) )
    throw PolynomialException(PolynomialException::WRONG_PLACE_SYMBOL, "^" );
  
  if( it != last )
    {
        Node*  temp = new Node;
        incase.push_back( temp);
      temp->m_type = EXP;
      temp->m_base =  0;    
      temp->m_base = createTree(first, it , variable, incase );
      temp->m_exp = toInteger( (it+1)->m_str1.begin(), (it+1)->m_str1.end());
      return temp;
    }
    
  if( first->m_type == LEFT_PAREN  && (first + (last - first -1))->m_type == RIGHT_PAREN )
   return createTree(first+1, first + (last - first -1) , variable , incase);
 
  if( first + 1 != last)
      throw PolynomialException(PolynomialException::WRONG_PLACE_SYMBOL, (first+1)->m_str1 );

  Node* temp = new Node;
  incase.push_back( temp );
  if( first->m_type == VAR )
    {
      temp->m_type = VAR;
      vector<string>::iterator p = find(variable.begin(), variable.end(), first->m_str1 );
      if( p == variable.end() )
	{
	  temp->m_var = variable.size();
	  variable.push_back ( first->m_str1 );
	}
      else 
	temp->m_var = p - variable.begin();
    }
  else if( first->m_type == NUMBER )
    {
      temp->m_type = NUMBER ;
      temp->m_real = toDouble(first->m_str1.begin(),first->m_str1.end()  );
      temp->m_imag = 0;
    }
  else if( first->m_type == COMPLEX )
    {
      temp->m_type = COMPLEX ;
      temp->m_real = toDouble(first->m_str1.begin(),first->m_str1.end()  );
      temp->m_imag = toDouble(first->m_str2.begin(),first->m_str2.end()  );
    }
  
  return temp;
}

void GenericPolynomial::printTree(const Node* aNode, const vector<string>& variable)
{
  if( aNode )
    {
      switch ( aNode->m_type )
        {
        case ADD:
        case MINUS:
	  cout<<'(';
	  printTree(aNode->m_left, variable);
	  cout << ( (aNode->m_type == ADD)? '+' : '-' );
	  printTree(aNode->m_right, variable);
	  cout<<')';
	  break;

        case MULTI:
  	  printTree(aNode->m_left, variable);
	  cout << '*';
	  printTree(aNode->m_right, variable);
	  break;
	  
        case DIV:
  	  printTree(aNode->m_left, variable);
	  cout << '/';
	  printTree(aNode->m_right, variable);
	  break;

        case EXP:
	  if( aNode->m_base->m_type != VAR && aNode->m_base->m_type != COEF )
            {
	      cout<<'(';
	      printTree(aNode->m_base, variable);
	      cout<<')';
            }
	  else
	    printTree(aNode->m_base, variable);
	  if( aNode->m_exp >= 0 )
	    cout<<'^'<<aNode->m_exp ;
	  
	  break;
	  
        case VAR:
	  cout<< variable[aNode->m_var];
	  break;
	  
        case NUMBER:
	  cout<< aNode->m_real ;
	  break;

        case COMPLEX:
	  cout<< complex<double>(aNode->m_real, aNode->m_imag) ;
	  break;

	default:
	  break;
        }
    }
}


void  GenericPolynomial::deleteTree(Node* aNode)
{
  switch ( aNode->m_type  )
    {
    case ADD:
    case MINUS:
    case MULTI:
    case DIV:
      deleteTree(aNode->m_left);
      deleteTree(aNode->m_right);
      break;

    case EXP:
        deleteTree(aNode->m_base);
        break;

    case COEF:
    case VAR:
      break;
      
    default:
      break;
    }
  
  delete aNode;
}


double GenericPolynomial::toDouble(string::iterator first, string::iterator last) 
{
  //  string::iterator ed;
  char p[100];
  copy(first, last, p);
  p[last - first]= 0;
  double  v = strtod(p, 0);
  // if( ed != last )
  // throw PolynomialException(PolynomialException::ILLEGAL_SYMBOL, string(first, last) );
  
  return v;
}

int GenericPolynomial::toInteger(string::iterator first, string::iterator last) 
{
  //  string::iterator ed;
  char p[100];
  copy(first, last, p);
  p[last - first]= 0;
  int  v = strtol(p, 0, 10);
//      if( ed != last )
//      throw PolynomialException(PolynomialException::ILLEGAL_SYMBOL, string(first, last) );
  
  return v;
}

bool GenericPolynomial::isVariable(string::iterator first, string::iterator last)
{
  if( ! isalpha( *first ) )
    return false;
  
  for(string::iterator it = first; it < last; ++it)
    if( ! isalnum( *it ) )
      return false;
  
  return true;
}

int GenericPolynomial::
checkParenthesis(vector<Token>::iterator first, vector<Token>::iterator last)
{
  int paren = 0;
  for( vector<Token>::iterator it = first ;  it != last;  ++it  )
    if( it->m_type == LEFT_PAREN )
      ++paren;
    else if(  it->m_type == RIGHT_PAREN )
      if( --paren < 0 )
	return  -1;
  
  return (paren > 0 )? 1 : 0;
}


vector<GenericPolynomial::Token>::iterator GenericPolynomial::
partition(vector<Token>::iterator first, vector<Token>::iterator last)
{
  int paren = ( first->m_type == LEFT_PAREN )? 1 : 0;
  vector<Token>::iterator it = first + 1;
  vector<Token>::iterator breakpoint = last;
  for(  ; it != last ; ++it)
  {
       switch (it->m_type)
       {
       case  LEFT_PAREN:
           ++paren;  
           break;

       case RIGHT_PAREN:
           --paren;
           break;

       case ADD:
       case MINUS:
           if( paren == 0 )
               breakpoint = it;
           break;

       default:
           break;
       }
  }
  
  if( breakpoint < last )
    return breakpoint;
  
  paren = 0; breakpoint = last;
  for( it = first;  it != last; ++it )
    if( it->m_type == LEFT_PAREN )
      ++paren;
    else if( it->m_type == RIGHT_PAREN )
      --paren;
    else if( paren == 0 && ( it->m_type == DIV || it->m_type == MULTI ) )
        breakpoint = it;
  
  return breakpoint;
}

vector<GenericPolynomial::Token>::iterator GenericPolynomial::
findExp(vector<Token>::iterator first, vector<Token>::iterator last)
{
  vector<Token>::iterator it = first;
  for(int paren = 0; (it != last) &&  ! (paren == 0 &&  it->m_type == EXP ) ; ++it )
    if( it->m_type == LEFT_PAREN )
      ++paren;
    else if( it->m_type == RIGHT_PAREN )
      --paren;
  
  return it;
}


GenericPolynomial::Node* GenericPolynomial::copyTree(Node* const T)
{
  if( ! T )
    return 0;
  
  Node* R = new Node;
  R->m_type = T->m_type;
  
  switch ( T->m_type )
    {
    case ADD:
    case MINUS:
    case MULTI:
    case DIV:
      R->m_left = copyTree(T->m_left);
      R->m_right = copyTree(T->m_right);
      break;
      
    case COMPLEX:
    case NUMBER:
      R->m_real = T->m_real;
      R->m_imag = T->m_imag;
      break;
      
    case VAR:
      R->m_var= T->m_var;
      break;
      
    case EXP:
      R->m_base = copyTree(T->m_base);
      R->m_exp = T->m_exp;
      break;

    default:
      break;
    }
  
  return R;
}

void GenericPolynomial::tokenize(string poly, vector<Token>& tokens)
{
  //delete spaces in the begining of the string
  for(string::iterator it = poly.begin(); 
      it!= poly.end() && ( *it == ' '  || *it == '\n' || *it == '\t' ); it=poly.erase(it) )
    ;
  
  //delete all other space in the string
  for(string::iterator first = poly.begin(), second=poly.begin()+1; second != poly.end(); )
    {
      if( *second == ' ' || *second == '\n' || *second == '\t' )
	{
          for( ++second; *second == ' ' || *second == '\n' || *second == '\t' ; ++second )
	    ;
          if( second != poly.end() )
	    {
              if(  (  ( isdigit(*first) || isalpha(*first) ) 
		  && (isdigit(*second) || isalpha(*second) )  )
                  || (isalnum(*first) && *second == '.') 
		  || (*first == '.' && isalnum(*second))   )
          throw PolynomialException(PolynomialException::WRONG_PLACE_SYMBOL, string(first, second) );
	    }
	  second = poly.erase(first+1, second);
	}

      if(second != poly.end() )
	{
          first = second;
          second = first + 1;
	}
    }

   vector<Token> T;
   for(string::iterator first = poly.begin(); first != poly.end(); )
    {
      if( isalpha(*first) )
        {
	  string::iterator second = first+1;
	  for( ; second != poly.end() && (isalnum(*second) ); ++second)
	    ;
	  T.push_back( Token(VAR, string(first, second) ) );
	  first= second;
        } 
      else if( isdigit(*first) || *first == '.' )
	{            
	  string::iterator second = first + 1;    
	  if( isdigit(*first) )
	    {
	      while( second != poly.end() && isdigit(*second) )
		++second;
	      if( second != poly.end() && *second == '.' )
		for( ++second; second != poly.end() && isdigit(*second); ++second)
		  ;
	      } 
	  else 
	    {
	      if( second == poly.end() || !isdigit(*second)  )
              throw PolynomialException(PolynomialException::WRONG_PLACE_SYMBOL, string(".") );

          while( second != poly.end() && isdigit(*second) )
		++second;
	    }
	  
	  if( second != poly.end() && ( *second == 'E'  || *second == 'e' ) )
	    {
	      ++second;
	      if( second == poly.end() ||  ( *second != '+' && *second != '-' ) ) 
              throw PolynomialException(PolynomialException::WRONG_PLACE_SYMBOL, string(first, second) );
	      else 
		{
		  ++second;
		  if( second == poly.end() ||  !isdigit( *second ) )
             throw PolynomialException(PolynomialException::WRONG_PLACE_SYMBOL, string(first, second) );
		}
	      
	      while(second != poly.end() &&  isdigit(*second) )
		++second;
	      if( second != poly.end() && !isop(*second) && *second != ')'   && *second != ',' )
            throw PolynomialException(PolynomialException::ILLEGAL_SYMBOL, string(first, second) );
	    }
	  else if( second != poly.end() && isalpha(*second) )
           throw PolynomialException(PolynomialException::ILLEGAL_SYMBOL, string(first, second) );

	  //Something wierd is here
	  string str = string(first, second);
	  Token tk = Token(NUMBER, str);
	  tk.m_str1.end();
	  T.push_back( tk );
	  first = second;
	}
      else 
        {
	  NodeType ty;
	  switch (*first)
	    {
	    case '(':
	      ty = LEFT_PAREN;       break;
	      
	    case ')':
	      ty = RIGHT_PAREN;     break;
	      
	    case '+':
	      ty = ADD;    break;
	      
	    case '-':
	      ty = MINUS;     break;
	      
	    case '*':
	      ty = MULTI;      break;
	      
	    case '/':
	      ty = DIV;      break;

	    case ',':
	      ty = COMMA;   break;
	      
	    case '^':
	      ty = EXP;       break;
	      
	    default:
           throw PolynomialException(PolynomialException::ILLEGAL_SYMBOL, string(first, first+1) );
	    }
	  
	  T.push_back( Token(ty,  string(first, first+1) ) );
	  ++first;
        }
    }
   
   vector<Token>::iterator it = T.begin();
   if( it->m_type == ADD )
     it++;
   else if( it->m_type == MINUS )
     {
       tokens.push_back( Token(NUMBER, string("-1") )  );
       tokens.push_back( Token(MULTI) );
       ++it;
     }     
   while(it != T.end())
     {
       if( T.end() - it > 4  && it->m_type == LEFT_PAREN   && (it+1)->m_type == NUMBER 
	     && (it+2)->m_type == COMMA &&  (it+3)->m_type == NUMBER 
	     && (it+4)->m_type == RIGHT_PAREN  )
	   {
	     tokens.push_back( Token(COMPLEX, (it+1)->m_str1, (it+3)->m_str1 ) );
	     it = it + 5;
	   }
       else if( T.end() - it > 5  && it->m_type == LEFT_PAREN   && (it+1)->m_type == MINUS
		&& (it+2)->m_type == NUMBER  && (it+3)->m_type == COMMA &&  
		(it+4)->m_type == NUMBER && (it+5)->m_type == RIGHT_PAREN  )
	   {
	     tokens.push_back( Token(COMPLEX, string("-")+(it+2)->m_str1, (it+4)->m_str1 ) );
	     it = it + 6;
	   }
       else if( T.end() - it > 6  && it->m_type == LEFT_PAREN   && (it+1)->m_type == MINUS
		&& (it+2)->m_type == NUMBER  && (it+3)->m_type == COMMA &&  
		(it+4)->m_type == MINUS && (it+5)->m_type == NUMBER 
		&& (it+6)->m_type == RIGHT_PAREN  )
	   {
	     tokens.push_back( Token(COMPLEX, string("-")+(it+2)->m_str1, 
				     string("-")+ (it+5)->m_str1 ) );
	     it = it + 7;
	   }
       else if( T.end() - it > 5  && it->m_type == LEFT_PAREN  
		&& (it+1)->m_type == NUMBER  && (it+2)->m_type == COMMA &&  
		(it+3)->m_type == MINUS && (it+4)->m_type == NUMBER 
		&& (it+5)->m_type == RIGHT_PAREN  )
	   {
	     tokens.push_back( Token(COMPLEX, (it+1)->m_str1, 
				     string("-")+ (it+4)->m_str1 ) );
	     it = it + 6;
	   }
        else if( T.end() - it > 1 && it->m_type == LEFT_PAREN && (it+1)->m_type == ADD )
	 {
	   tokens.push_back(*it);
	   it += 2;
	 }
       else if( T.end() - it > 1 && it->m_type == LEFT_PAREN && (it+1)->m_type == MINUS )
	 {
	   tokens.push_back(*it);
	   tokens.push_back( Token(NUMBER, string("-1"))  );
	   tokens.push_back(Token(MULTI, string() ) );
	   it += 2;
	 }
       else
	 tokens.push_back(*it++);
     }
//     for( vector<Token>::iterator it = tokens.begin(); it != tokens.end(); ++it)
//       {
//         string s;
//         switch (it->m_type )
//  	 {
//  	 case MULTI:
//  	   s = string("multi");
//  	   break;

//  	 case MINUS:
//  	   s = string("minus");
//  	   break;

//  	 case NUMBER:
//  	   s = string("NUMBER");
//  	   break;

//  	 case VAR:
//  	   s = string("var");
//  	   break;

//  	 default:
//  	   break;
//  	 }

//         cout<< s<<  it->m_str1;
//       }
//     cout<<endl;
}

ExpandedPolynomial GenericPolynomial::expand(Node* node)
{
  if(node->m_type == VAR)
    {
      ExpandedPolynomial::VarDegrees temp(1, ExpandedPolynomial::VarDegree(node->m_var, 1));
      return ExpandedPolynomial(ExpandedPolynomial::Pterm(ExpandedPolynomial::ValueType(1.0, 0.0), 
							temp) );
    }

  if(node->m_type == NUMBER || node->m_type == COMPLEX )
    return ExpandedPolynomial( ExpandedPolynomial::ValueType(node->m_real, node->m_imag) );

  if( node->m_type == ADD )
    return  expand(node->m_left) + expand(node->m_right);

  if( node->m_type == MINUS )
    return  expand(node->m_left) - expand(node->m_right);

  if( node->m_type == MULTI )
    return  expand(node->m_left)* expand(node->m_right);

  if( node->m_type == DIV )
    return  expand(node->m_left) / expand(node->m_right);

  if( node->m_type == EXP )
    return  expand(node->m_base)^node->m_exp;

  return ExpandedPolynomial();
}


ExpandedPolynomial GenericPolynomial::expand()
{
    ExpandedPolynomial p = expand(m_root);
    p.setVariable(m_variable);

    return p;
}
