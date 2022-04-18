#include "stdafx.h"
#include "PlanarPH_mfc.h"
#include "Analyse.h"
#include "afxdialogex.h"

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

IMPLEMENT_DYNAMIC(CAnalyse, CDialogEx)

CAnalyse::CAnalyse(PlanarPH pph, int index, int state, CWnd* pParent /*=NULL*/)
	: CDialogEx(CAnalyse::IDD, pParent)
{
	a_pph = pph;
	a_stat = state;
	a_index = index;
	margin = 30.0;
	pointing = FALSE;
}

CAnalyse::~CAnalyse()
{
}

void CAnalyse::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAnalyse, CDialogEx)
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_SIZE()
	ON_WM_PAINT()
	ON_WM_ERASEBKGND()
	ON_WM_MOUSEMOVE()
END_MESSAGE_MAP()

int CAnalyse::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CDialogEx::OnCreate(lpCreateStruct) == -1)
		return -1;

	m_oldDC = wglGetCurrentDC();
	m_oldRC = wglGetCurrentContext();

	m_pDC = new CClientDC(this);
	ASSERT(m_pDC != NULL);
	if(!SetDCPixelFormat())
		return -1;

	m_hRC= wglCreateContext(m_pDC->GetSafeHdc());  
    if(!m_hRC)  
		MessageBox("Failed to create rendering context");  
    if(!wglMakeCurrent(m_pDC->GetSafeHdc(), m_hRC))  
		MessageBox("Failed to make the rendering context the current context");

	buildFont();

	return 0;
}

BOOL CAnalyse::SetDCPixelFormat()
{
	PIXELFORMATDESCRIPTOR pfd = {    
        sizeof(PIXELFORMATDESCRIPTOR),        
        1,                                    
        PFD_DRAW_TO_WINDOW |                 
        PFD_SUPPORT_OPENGL |                
        PFD_TYPE_RGBA,                    
        24,                               
        0, 0, 0, 0, 0, 0,                  
        0,                                
        0,                               
        0,                               
        0, 0, 0, 0,                     
        32,                                 
        0,                             
        0,                             
        PFD_MAIN_PLANE,                  
        0,                              
        0, 0, 0                        
    };   

	int pixelformat; 
	if((pixelformat = ChoosePixelFormat(m_pDC->GetSafeHdc(), &pfd)) == 0)
	{
		MessageBox("Choose Pixel Format failed"); 
        return FALSE;
	}      
	if(SetPixelFormat(m_pDC->GetSafeHdc(), pixelformat, &pfd) == FALSE)
	{
		MessageBox("Set Pixel Format failed");
		return FALSE; 
	} 

	return TRUE;
}

void CAnalyse::OnDestroy()
{
	CDialogEx::OnDestroy();

	::wglMakeCurrent(NULL, NULL);   
    ::wglDeleteContext(m_hRC);  
	wglMakeCurrent(m_oldDC, m_oldRC);
}

void CAnalyse::RenderScene()
{
	int i, j;
	double t, xtemp, htemp;

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	double h_max = 0.0;
	if(a_stat == statPara)
		for(i = 1; i <= a_index; i++)
		{
			for(j = 0; j <= Step; j++)
			{
				t = 1.0/Step*j;
				if(h_max < getParaSpeed(a_pph, i, t))
					h_max = getParaSpeed(a_pph, i, t);
			}
		}
	else
		h_max = getArcLength(a_pph, a_index, 1.0);

	glClearColor(1.0, 1.0, 1.0, 1.0);   
    glClear(GL_COLOR_BUFFER_BIT);
	glLineWidth(1);
	glColor3fv(colors[6]);
	glBegin(GL_LINE_LOOP); 
		glVertex2f(margin, height-margin);
		glVertex2f(width-margin, height-margin);
		glVertex2f(width-margin, margin);
		glVertex2f(margin, margin);
	glEnd();
	glBegin(GL_LINE_STRIP);
		glColor3f(0.1, 0.1, 0.1);
		for(i = 1; i <= a_index; i++)
		{
			for(j = 0; j <= Step; j++)
			{
				t = 1.0/Step*j;
				switch (a_stat)
				{
				case statArcl:
					glVertex2f((width-2*margin)/a_index*(i-1+t)+margin, 
								margin+getArcLength(a_pph, i, t)/h_max*(height-2.0*margin));
					break;
				case statPara:
					glVertex2f((width-2*margin)/a_index*(i-1+t)+margin, 
							   margin+(height-2.0*margin)/h_max*getParaSpeed(a_pph, i, t));
					break;
				default:
					break;
				}
			}	
		}
	glEnd();

	if(pointing)
	{
		i = (int)((cursor_x-margin)*a_index/(width-2.0*margin))+1;
		t = (cursor_x-margin-(width-2.0*margin)/a_index*(i-1))/((width-2.0*margin)/a_index);
		glLineStipple(2, 0x5555);
		glEnable(GL_LINE_STIPPLE);
		glBegin(GL_LINES); 
			glVertex2f(cursor_x, margin);
			if(a_stat == statArcl)
			{
				htemp = getArcLength(a_pph, i, t);
				glVertex2f(cursor_x, margin+htemp/h_max*(height-2.0*margin));
			}
			else
			{
				htemp = getParaSpeed(a_pph, i, t);
				glVertex2f(cursor_x, margin+(height-2.0*margin)/h_max*htemp);
			}
		glEnd();
		glDisable(GL_LINE_STIPPLE);

		t += i-1;
		if(cursor_x < margin+25.0)
			xtemp = margin+25.0;
		else
			if(cursor_x > width-margin-25.0)
				xtemp = width-margin-25.0;
			else 
				xtemp = cursor_x;
		glRasterPos2f(xtemp-25.0, margin-10.0);
		printString("(%.2f)", t);
		glRasterPos2f(xtemp-25.0, height-margin+5.0);
		printString("(%.2f)", htemp);
	}

	glRasterPos2f(margin-10.0, margin-5.0);
	printString("0");
	glRasterPos2f(width-margin+5.0, margin-5.0);
	printString("%d", a_index);
	glRasterPos2f(margin-25.0, height-margin);
	printString("%d", (int)h_max);

	glFlush();
}

void CAnalyse::printString(const char* str, ...)
{
	char text[256];
	va_list ap;
	if(str == NULL)
		return;
	va_start(ap, str);
	vsprintf(text, str, ap);
	va_end(ap);
	glPushAttrib(GL_LIST_BIT);
	glListBase(base-32);
	glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
	glPopAttrib();
}

void CAnalyse::buildFont()
{
	HFONT font, oldfont;

	base = glGenLists(96);
	font = CreateFont(-12, 0, 0, 0, 
					  FW_BOLD,
					  FALSE, FALSE, FALSE,
					  ANSI_CHARSET,
					  OUT_TT_PRECIS,
					  CLIP_DEFAULT_PRECIS,
					  ANTIALIASED_QUALITY,
					  FF_DONTCARE | DEFAULT_PITCH,
					  "Courier New");
	oldfont = (HFONT)SelectObject(m_pDC->GetSafeHdc(), font);
	wglUseFontBitmaps(m_pDC->GetSafeHdc(), 32, 96, base);
	SelectObject(m_pDC->GetSafeHdc(), oldfont);
	DeleteObject(font);
}

void CAnalyse::OnSize(UINT nType, int cx, int cy)
{
	CDialogEx::OnSize(nType, cx, cy);

	width = cx;
	height = cy; 

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, width, 0.0, height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void CAnalyse::OnPaint()
{
	CPaintDC dc(this); 

	RenderScene();
}

BOOL CAnalyse::OnEraseBkgnd(CDC* pDC)
{
	return TRUE;
}

void CAnalyse::OnMouseMove(UINT nFlags, CPoint point)
{
	if(point.x >= margin && point.x <= width-margin && point.y >= margin && point.y <= height-margin)
	{
		pointing = TRUE;
		cursor_x = (float)point.x;
		RenderScene();
	}
	else 
		if(pointing)
		{
			pointing = FALSE;
			RenderScene();
		}

	CDialogEx::OnMouseMove(nFlags, point);
}
