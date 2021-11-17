#pragma once
#include "stdafx.h"

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

class CPointPH
{
public:
	double xPos, yPos;

	CPointPH(void)
	{
		xPos = 0.0;
		yPos = 0.0;
	}

	CPointPH(CPoint pt)
	{
		xPos = (double)pt.x;
		yPos = (double)pt.y;
	}

	CPointPH(double x, double y)
	{
		xPos = x;
		yPos = y;
	}

	~CPointPH(void)
	{
	}

	CPointPH operator +(const CPointPH &p) 
	{ 
		CPointPH temp; 
		temp.xPos = xPos+p.xPos; 
		temp.yPos = yPos+p.yPos; 
		return temp; 
	} 

	CPointPH operator +(const CPoint &p) 
	{ 
		CPointPH temp; 
		temp.xPos = xPos+(double)p.x; 
		temp.yPos = yPos+(double)p.y; 
		return temp; 
	} 

	BOOL operator ==(const CPointPH &p)
	{
		if(xPos == p.xPos && yPos == p.yPos)
			return TRUE;
		else return FALSE;
	}
};

