// Author: Yong-Kang Zhu (yongkang.zhu@gmail.com)
// Code can be used only for academic purpose

#include "ExactSum.h"

ExactSum::ExactSum()
{
	set_fpu(0x27F);
	r_c = 0;
	t_s = new double[N2_EXPONENT + 1];
	t_s2 = new double[N2_EXPONENT + 1];
	Reset();
}

ExactSum::~ExactSum()
{
	delete [] t_s;
	delete [] t_s2;
}

void ExactSum::set_fpu (unsigned int mode)
{
#ifdef DOUBLE
	asm ("fldcw %0" : : "m" (*&mode));
#endif
}

void ExactSum::AddTwo(double &a, double &b)
{
	double t = a + b;
	b = ( ((str_double*)(&a))->exponent < ((str_double*)(&b))->exponent ) ?
		(b - t) + a : (a - t) + b;
	a = t;
}

int ExactSum::Round3(double s0, double s1, double s2)
{
	// To check "s1 is half-ulp", s1!=0 and all the mantissa bits=0 only work
	// for normalized numbers. But if s1 is de-normalized, then according to
	// the two places where Round3 gets called, \hat s_2 must be zero, which
	// means s0 is correctly rounded.
	if (s1 != .0 && ((str_double*)(&s1))->mantissa_high == 0 &&
		((str_double*)(&s1))->mantissa_low == 0 &&
		s1 * s2 > 0)
		return 1;
	return 0;
}	

double ExactSum::iFastSum(double *num_list, int n)
{
	if (n < 1) return .0;
	double s = 0, s_t, s1, s2, e1, e2;
	int count;          // next position in num_list to store error
	int c_n = n;          // current number of summands
	unsigned max = 0;     // the max exponent of s_t
	int i;
	double t, e_m;
	double half_ulp = .0;  // half an ulp of s

	double EPS = 1.0;
	((str_double*)(&EPS))->exponent -= 53;

	double ev_d = .0;
	((str_double*)(&ev_d))->sign = 0;

	for (i = 1; i <= n; i++)
	{
		// AddTwo, inline
		t = s + num_list[i];
		num_list[i] =
		((str_double*)(&s))->exponent
			< ((str_double*)(&(num_list[i])))->exponent ?
			(num_list[i] - t) + s : (s - t) + num_list[i];
		s = t;
	}
	
	while(1)
	{
		count = 1;
		s_t = .0;
		max = 0;
		for (i = 1; i <= c_n; i++)
		{
			// AddTwo, inline
			t = s_t + num_list[i];
			num_list[count] =
				((str_double*)(&s_t))->exponent
				< ((str_double*)(&num_list[i]))->exponent ?
				(num_list[i] - t) + s_t : (s_t - t) + num_list[i];
			s_t = t;

			if (num_list[count] != 0)			 
			{
				count++;
				if(max < ((str_double*)&s_t)->exponent)
					max = ((str_double*)&s_t)->exponent;
			}
		}
		
		// compute e_m, the estimated global error
		if (max > 0) // neither minimum exponent nor de-normalized 
		{
			((str_double*)(&ev_d))->exponent = max;
			ev_d *= EPS;
			e_m = ev_d * (count-1);
		}
		else
			e_m = .0;

		AddTwo(s, s_t);
		num_list[count] = s_t;
		c_n = count;

		// compute HalfUlp(s)
		if (((str_double*)(&s))->exponent > 0) 
		{
			((str_double*)(&half_ulp))->exponent = ((str_double*)(&s))->exponent;
			half_ulp *= EPS;
		}
		else
			half_ulp = .0;

		if (e_m < half_ulp || e_m == .0)
		{
			if (r_c > 0) return s;
			s1 = s2 = s_t;
			e1 = e_m;
			e2 = -e_m;
			AddTwo(s1, e1);
			AddTwo(s2, e2);
			if (s + s1 != s || s + s2 != s || 
				Round3(s, s1, e1) || Round3(s, s2, e2))
			{
				r_c = 1;
				double s1 = iFastSum(num_list, c_n);
				AddTwo(s, s1);
				double s2 = iFastSum(num_list, c_n);
				r_c = 0;
				if (Round3(s, s1, s2))
				{
					((str_double*)(&s1))->mantissa_low |= 0x1;   // the Magnify function
					s += s1;
				}
			}
			return s;
		}
	}
}

double ExactSum::OnlineExactSum(double *num_list, int n)
{
	if (n < 1) return .0;
	double *temp_swap;
	// create two accumulators
	double *temp_ss = new double[N2_EXPONENT + 1];
	double *temp_ss2 = new double[N2_EXPONENT + 1];
	double t;
	int i;
	unsigned exp;
	double s;
	for (i = 0; i < N2_EXPONENT; i++)
		temp_ss[i] = .0;

	int num = n > MAX_N ? MAX_N : n;
	while (1)
	{
		for (i = 1; i <= num; i++)
		{
			exp = ((str_double*)(&num_list[i]))->exponent;
			// AddTwo combined with an addition
			t = temp_ss[exp] + num_list[i];
			temp_ss[exp + N_EXPONENT] +=
					((str_double*)(&temp_ss[exp]))->exponent
					< exp ? (num_list[i] - t) + temp_ss[exp]
					: (temp_ss[exp] - t) + num_list[i];
			temp_ss[exp] = t;
		}
		if (num < n)
		{
			for (i = 0; i < N2_EXPONENT; i++)
				temp_ss2[i] = .0;
			for (i = 0; i < N2_EXPONENT; i++)
			{
				exp = ((str_double*)(&temp_ss[i]))->exponent;
				t = temp_ss2[exp] + temp_ss[i];
				temp_ss2[exp + N_EXPONENT] +=
						((str_double*)(&temp_ss2[exp]))->exponent
						< exp ? (temp_ss[i] - t) + temp_ss2[exp]
						: (temp_ss2[exp] - t) + temp_ss[i];
				temp_ss2[exp] = t;
			}
			num_list += num;
			n -= num;
			num = (n > MAX_N_AFTER_SWAP) ?  MAX_N_AFTER_SWAP : n;
			// Swap the pointers of two accumulators
			temp_swap = temp_ss;
			temp_ss = temp_ss2;
			temp_ss2 = temp_swap;
		}
		else
			break;
	}
	int j = 0;
	// extract all the non-zero accumulators
	for (i = 0; i < N2_EXPONENT; i++)
		if (temp_ss[i] != .0)
			temp_ss2[++j] = temp_ss[i];
	s = iFastSum(temp_ss2, j);
	delete [] temp_ss;
	delete [] temp_ss2;
	return s;
}

void ExactSum::AddNumber(double x)
{
	int i;
	double t;
	unsigned exp;
	if (c_num >= MAX_N)
	{
			for (i = 0; i < N2_EXPONENT; i++)
				t_s2[i] = 0;
			for (i = 0; i < N2_EXPONENT; i++)
			{
				exp = ((str_double*)(&t_s[i]))->exponent;
				// AddTwo, inline
				t = t_s2[exp] + t_s[i];
				t_s2[exp + N_EXPONENT] +=
						((str_double*)(&t_s2[exp]))->exponent
						< exp ? (t_s[i] - t) + t_s2[exp]
						: (t_s2[exp] - t) + t_s[i];
				t_s2[exp] = t;
			}
			double *t_swap = t_s;
			t_s = t_s2;
			t_s2 = t_swap;
			c_num = N2_EXPONENT;  // t_s has added N2_EXPONENT numbers
	}

	exp = ((str_double*)(&x))->exponent;
	// AddTwo, inline
	t = t_s[exp] + x;
	t_s[exp + N_EXPONENT] +=
			((str_double*)(&t_s[exp]))->exponent
			< exp ? (x - t) + t_s[exp]
			: (t_s[exp] - t) + x;
	t_s[exp] = t;
	c_num ++;
}

void ExactSum::AddArray(double *num_list, int n)
{
	if (n < 1) return ;

	// we could use AddNumber directly, but the following is more efficient
	int i;
	double t;
	unsigned exp;

	int num = (n > MAX_N - c_num) ? MAX_N - c_num : n;
	
	c_num += num;

	while (1)
	{
		for (i = 1; i <= num; i++)
		{
			exp = ((str_double*)(&num_list[i]))->exponent;
			// AddTwo combined with an addition
			t = t_s[exp] + num_list[i];
			t_s[exp + N_EXPONENT] +=
					((str_double*)(&t_s[exp]))->exponent
					< exp ? (num_list[i] - t) + t_s[exp]
					: (t_s[exp] - t) + num_list[i];
			t_s[exp] = t;
		}
		if (num < n)
		{
			for (i = 0; i < N2_EXPONENT; i++)
				t_s2[i] = .0;
			for (i = 0; i < N2_EXPONENT; i++)
			{
				exp = ((str_double*)(&t_s[i]))->exponent;
				t = t_s2[exp] + t_s[i];
				t_s2[exp + N_EXPONENT] +=
						((str_double*)(&t_s2[exp]))->exponent
						< exp ? (t_s[i] - t) + t_s2[exp]
						: (t_s2[exp] - t) + t_s[i];
				t_s2[exp] = t;
			}
			num_list += num;
			n -= num;
			num = (n > MAX_N_AFTER_SWAP) ?  MAX_N_AFTER_SWAP : n;
			// Swap the pointers of two accumulators
			double *t_swap = t_s;
			t_s = t_s2;
			t_s2 = t_swap;
			c_num = N2_EXPONENT;  // t_s has added N2_EXPONENT numbers
		}
		else
			break;
	}
}

double ExactSum::GetSum()
{
	int i, j = 0;
	// same to iHyrbidSum: use iFastSum to sum all the non-zero accumulators
	for (i = 0; i < N2_EXPONENT; i++)
		if (t_s[i] != .0)
			t_s2[++j] = t_s[i];
	return iFastSum(t_s2, j);
}

void ExactSum::Reset()
{
	c_num = 0;
	int i;
	for (i = 0; i < N2_EXPONENT; i++)
		t_s[i] = .0;
}
