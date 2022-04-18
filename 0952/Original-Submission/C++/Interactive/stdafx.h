// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently,
// but are changed infrequently
#pragma once

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN            // Exclude rarely-used stuff from Windows headers
#endif
#include "targetver.h"
#define _ATL_CSTRING_EXPLICIT_CONSTRUCTORS      // some CString constructors will be explicit
// turns off MFC's hiding of some common and often safely ignored warning messages
#define _AFX_ALL_WARNINGS
#include <afxwin.h>         // MFC core and standard components
#include <afxext.h>         // MFC extensions
#include <afxdisp.h>        // MFC Automation classes
#ifndef _AFX_NO_OLE_SUPPORT
#include <afxdtctl.h>           // MFC support for Internet Explorer 4 Common Controls
#endif
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>             // MFC support for Windows Common Controls
#endif // _AFX_NO_AFXCMN_SUPPORT
#include <afxcontrolbars.h>     // MFC support for ribbons and control bars

//
#include <GL/glut.h>
#include <math.h>
const int Step = 50;
const int Dgr = 5;
const double Pi = 3.1415926535897932384;
const double Epsilon = 0.000000000001;
const int Max = 50;
const int MAXDEGREE = 100;

const int plotPoint = 100;
const int plotSpline = 101;
const int plotCtrl = 103;
const int plotAxis = 104;
const int plotOffset = 105;

const int statNew = 200;
const int statPlot = 201;
const int statEnd = 202;
const int statTrans = 203;
const int statPara = 204;
const int statArcl = 205;
const int statHermite = 206;

/*
COLOR LIST
0: Points Color
1: Activated Points Color
2: PH Curve Color
3: Cubic Curve Color
4: PH Ctrl Poly Color
5: Cubic Ctrl Poly Color
6: Analyse Color
7: Axis Color
*/
const GLfloat colors[][3] = {  
    {0.2, 0.2, 0.2},
	{0.2, 0.8, 0.2},
	{0.2, 0.2, 0.2},	
	{0.8, 0.2, 0.2},	
    {0.4, 0.4, 0.4},	
	{0.8, 0.4, 0.4},	
	{0.0, 0.0, 0.0},	
	{0.6, 0.6, 0.6}	
}; 

#ifdef _UNICODE
#if defined _M_IX86
#pragma comment(linker,"/manifestdependency:\"type='win32' name='Microsoft.Windows.Common-Controls' version='6.0.0.0' processorArchitecture='x86' publicKeyToken='6595b64144ccf1df' language='*'\"")
#elif defined _M_X64
#pragma comment(linker,"/manifestdependency:\"type='win32' name='Microsoft.Windows.Common-Controls' version='6.0.0.0' processorArchitecture='amd64' publicKeyToken='6595b64144ccf1df' language='*'\"")
#else
#pragma comment(linker,"/manifestdependency:\"type='win32' name='Microsoft.Windows.Common-Controls' version='6.0.0.0' processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")
#endif
#endif