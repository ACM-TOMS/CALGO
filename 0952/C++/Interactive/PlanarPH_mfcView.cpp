#include "stdafx.h"
#ifndef SHARED_HANDLERS
#include "Mainfrm.h"
#include "PlanarPH_mfc.h"
#endif
#include "PlanarPH_mfcDoc.h"
#include "PlanarPH_mfcView.h"
#ifdef _DEBUG
#define new DEBUG_NEW
#endif

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

IMPLEMENT_DYNCREATE(CPlanarPH_mfcView, CView)
BEGIN_MESSAGE_MAP(CPlanarPH_mfcView, CView)
	ON_WM_CREATE()   
    ON_WM_DESTROY()   
    ON_WM_SIZE()
	ON_WM_ERASEBKGND()
    ON_WM_LBUTTONDOWN() 
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONDOWN() 
	ON_WM_MOUSEMOVE()
	ON_COMMAND(ID_EDIT_NEW, OnMenuNew)
	ON_COMMAND(ID_EDIT_MULPOINTS, OnMenuEditMulpt)
	ON_COMMAND(ID_EDIT_TRANS, OnMenuTranslate)
	ON_COMMAND(ID_EDIT_OFFSET, OnMenuOffset)
	ON_COMMAND(ID_VIEW_HERMITE, OnMenuHermite)
	ON_COMMAND(ID_VIEW_PH, OnMenuPh)
	ON_COMMAND(ID_VIEW_CUBIC, OnMenuCubic)
	ON_COMMAND(ID_VIEW_PHCTRL, OnMenuPhCtrl)
	ON_COMMAND(ID_VIEW_CUBICCTRL, OnMenuCubicCtrl)
	ON_COMMAND(ID_ANALYSE_PARASPEED, OnMenuParaspeed)
	ON_COMMAND(ID_ANALYSE_ARCLENGTH, OnMenuArclength)
	ON_UPDATE_COMMAND_UI(ID_VIEW_PH, OnUpdateMenuPh)
	ON_UPDATE_COMMAND_UI(ID_VIEW_CUBIC, OnUpdateMenuCubic)
	ON_UPDATE_COMMAND_UI(ID_VIEW_CUBICCTRL, OnUpdateMenuCubicCtrl)
	ON_UPDATE_COMMAND_UI(ID_VIEW_PHCTRL, OnUpdateMenuPhCtrl)
END_MESSAGE_MAP()

CPlanarPH_mfcView::CPlanarPH_mfcView()
{
	index = 0;
	status = statNew;
	pointing = -1;
	isMenuHermite = FALSE;
	isMenuPh = TRUE;
	isMenuCubic = FALSE;
	isMenuPhCtrl = FALSE;
	isMenuCubicCtrl = FALSE;
	isLBDown = FALSE;
	offsetNum = 0;
}

CPlanarPH_mfcView::~CPlanarPH_mfcView()
{
}

BOOL CPlanarPH_mfcView::PreCreateWindow(CREATESTRUCT& cs)
{
	cs.style |= WS_CLIPCHILDREN | WS_CLIPSIBLINGS; 

	return CView::PreCreateWindow(cs);
}

void CPlanarPH_mfcView::OnDraw(CDC* pDC)
{
	CPlanarPH_mfcDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	if(!pDoc)
		return;
}

#ifdef _DEBUG
void CPlanarPH_mfcView::AssertValid() const
{
	CView::AssertValid();
}

void CPlanarPH_mfcView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CPlanarPH_mfcDoc* CPlanarPH_mfcView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CPlanarPH_mfcDoc)));
	return (CPlanarPH_mfcDoc*)m_pDocument;
}
#endif //_DEBUG

//////////////////////////////////////////////////////////////
int CPlanarPH_mfcView::OnCreate(LPCREATESTRUCT lpCreateStruct)    
{
	if(CView::OnCreate(lpCreateStruct) == -1)
		return -1;

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

BOOL CPlanarPH_mfcView::SetDCPixelFormat()
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

void CPlanarPH_mfcView::OnDestroy()    
{   
    CView::OnDestroy();   
       
    ::wglMakeCurrent(NULL, NULL);   
    ::wglDeleteContext(m_hRC);   
    if(m_pDC)   
        delete m_pDC;  
	glDeleteLists(base, 96);
}   

BOOL CPlanarPH_mfcView::OnEraseBkgnd(CDC* pDC)  
{  
	return TRUE; 
}

void CPlanarPH_mfcView::OnSize(UINT nType, int cx, int cy)    
{   
    CView::OnSize(nType, cx, cy);   
   
	width = cx;
	height = cy; 

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, width, 0.0, height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	if(status == statNew)
	{
		o_ActRF = CPointPH((double)width/2.0, (double)height/2.0);
		o_ActRFTrans = o_ActRF;
		DrawLine(0, plotAxis);
	}
	else
		DrawLine(index-1, plotSpline);
} 

void CPlanarPH_mfcView::OnMenuHermite()
{
	if(!isMenuHermite)
	{
		if(status == statTrans)
			OnMenuTranslate();
		status = statNew;
		isMenuHermite = TRUE;
		isMenuPh = FALSE;
		isMenuCubic = FALSE;
		isMenuPhCtrl = FALSE;
		isMenuCubicCtrl = FALSE;
		index = 0;
		status = statNew;
		pointing = -1;
		offsetNum = 0;
		DrawBG();
		glFlush();
	}	
}

void CPlanarPH_mfcView::OnMenuPh()
{
	if(isMenuHermite)
	{
		status = statNew;
		isMenuHermite = FALSE;
		isMenuPh = TRUE;
		DrawBG();
		glFlush();
	}
	else
	{
		isMenuPh = !isMenuPh;
		if(status == statEnd)
			DrawLine(index-1, plotSpline);
	}
}

void CPlanarPH_mfcView::OnMenuCubic()
{
	if(isMenuHermite)
	{
		status = statNew;
		isMenuHermite = FALSE;
		isMenuCubic = TRUE;
		DrawBG();
		glFlush();
	}
	else
	{
		isMenuCubic = !isMenuCubic;
		if(status == statEnd)
			DrawLine(index-1, plotSpline);
	}
}

void CPlanarPH_mfcView::OnMenuPhCtrl()
{
	if(isMenuHermite)
	{
		status = statNew;
		isMenuHermite = FALSE;
		isMenuPhCtrl = TRUE;
		DrawBG();
		glFlush();
	}
	else
	{
		isMenuPhCtrl = !isMenuPhCtrl;
		if(status == statEnd)
			DrawLine(index-1, plotSpline);
	}
}

void CPlanarPH_mfcView::OnMenuCubicCtrl()
{
	if(isMenuHermite)
	{
		status = statNew;
		isMenuHermite = FALSE;
		isMenuCubicCtrl = TRUE;
		DrawBG();
		glFlush();
	}
	else
	{
		isMenuCubicCtrl = !isMenuCubicCtrl;
		if(status == statEnd)
			DrawLine(index-1, plotSpline);
	}
}

void CPlanarPH_mfcView::OnMenuNew()
{
	if(status == statTrans)
		OnMenuTranslate();
	DrawBG();
	glFlush();
	index = 0;
	status = statNew;
	pointing = -1;
	offsetNum = 0;
}

void CPlanarPH_mfcView::OnMenuEditMulpt()
{
	if((status == statNew) | (status == statPlot))
		MessageBox("Please plot the spline first!", "Error", MB_OK);
	else
	{
		CAdjWindowA m_Awa(ptActRF, index, this);
		if(m_Awa.DoModal() == IDOK)
		{
			for(int i = 0; i < index; i++)
				ptActRF[i] = m_Awa.ReturnPoint(i);
			InitPH();
			DrawLine(index-1, plotSpline);
		}
	}	
}

void CPlanarPH_mfcView::OnMenuTranslate()
{
	if(status != statTrans)
	{
		status = statTrans;
		SetClassLong(this->GetSafeHwnd(), GCL_HCURSOR, (LONG)LoadCursor(NULL, IDC_HAND));
	}
	else
	{
		if(index == 0)
			status = statNew;
		else
			status = statEnd;
		SetClassLong(this->GetSafeHwnd(), GCL_HCURSOR, (LONG)LoadCursor(NULL, IDC_ARROW));
	}
}

void CPlanarPH_mfcView::OnMenuOffset()
{
	if(index > 2)
	{
		COffsetDlg m_OffsetDlg;
		if(m_OffsetDlg.DoModal() == IDOK)
		{
			offsetD[offsetNum] = m_OffsetDlg.d;
			if(offsetD[offsetNum] == 0.0)
				offsetNum = 0;
			else
				offsetNum++;
			DrawLine(index-1, plotSpline);
		}
	}
	else
		MessageBox("There must be more than 2 points!", "Need more points", MB_OK);
}

void CPlanarPH_mfcView::OnMenuParaspeed()
{
	if(index > 2)
	{
		if(isMenuHermite)
		{
			CAnalyse m_Analyse(pph, 1, statPara, this);
			m_Analyse.DoModal();
		}
		else
		{
			CAnalyse m_Analyse(pph, index-1, statPara, this);
			m_Analyse.DoModal();
		}
	}
	else
		MessageBox("There must be more than 2 points!", "Need more points", MB_OK);
}

void CPlanarPH_mfcView::OnMenuArclength()
{
	if(index > 2)
	{
		if(isMenuHermite)
		{
			CAnalyse m_Analyse(pph, 1, statArcl, this);
			m_Analyse.DoModal();
		}
		else
		{
			CAnalyse m_Analyse(pph, index-1, statArcl, this);
			m_Analyse.DoModal();
		}
	}
	else
		MessageBox("There must be more than 2 points!", "Need more points", MB_OK);
}

void CPlanarPH_mfcView::OnUpdateMenuPh(CCmdUI* pCmdUI)
{
	pCmdUI->SetCheck(isMenuPh);
}

void CPlanarPH_mfcView::OnUpdateMenuCubic(CCmdUI* pCmdUI)
{
	pCmdUI->SetCheck(isMenuCubic);
}

void CPlanarPH_mfcView::OnUpdateMenuPhCtrl(CCmdUI* pCmdUI)
{
	pCmdUI->SetCheck(isMenuPhCtrl);
}

void CPlanarPH_mfcView::OnUpdateMenuCubicCtrl(CCmdUI* pCmdUI)
{
	pCmdUI->SetCheck(isMenuCubicCtrl);
}

void CPlanarPH_mfcView::DrawBG()   
{
	glClearColor(1.0, 1.0, 1.0, 1.0);   
    glClear(GL_COLOR_BUFFER_BIT);

	glBegin(GL_LINES);
		glColor3fv(colors[7]);
			glVertex2f(0.0, (float)height-(float)o_ActRF.yPos);
			glVertex2f(width, (float)height-(float)o_ActRF.yPos);
			glVertex2f((float)o_ActRF.xPos, 0.0);
			glVertex2f((float)o_ActRF.xPos, (float)height);
		glEnd();

	glRasterPos2f(0.0, 0.0);
	if(isMenuHermite)
		printString("First-order Hermite");
	else
		printString("Planar PH");
}

/*
@FN:		PRETRANSLATEMESSAGE
@BRIEF:		INTERCEPT THE MESSAGE FROM KEYBOARD
			FUNCTION CALLED WHEN SYSTEM GET A MESSAGE FROM KEYBOARD
			N:	SAME AS NEW MENU OPTION;
			CTRL+E:	SAME AS EDIT MULTIPLE POINTS MENU OPTION;
			H:	SAME AS OFFSET MENU OPTION;
			C:	CLOSE THE SPLINE, END PLOT STATUS, PLOT THE CLOSE CURVE
@PARAM:		PMSG<MSG *>
@RETURN:	<BOOL>
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		DEC-22-2013
*/
BOOL CPlanarPH_mfcView::PreTranslateMessage(MSG*pMsg)
{
	if(pMsg->message == WM_KEYDOWN)
	{
		switch(pMsg->wParam)
		{
		case 'N':
			OnMenuNew();
			return TRUE;
		case 'E':
			if(GetKeyState(VK_CONTROL) < 0)
				OnMenuEditMulpt();
			return TRUE;
		case 'H':
			OnMenuHermite();
			return TRUE;
		case 'T':
			OnMenuTranslate();
			return TRUE;
		case 'O':
			OnMenuOffset();
			return TRUE;
		case 'A':
			if(GetKeyState(VK_CONTROL) < 0)
				OnMenuParaspeed();
			else
				OnMenuArclength();
			return TRUE;
		case 'P':
			if(GetKeyState(VK_CONTROL) < 0)
				OnMenuPhCtrl();
			else
				OnMenuPh();
			return TRUE;
		case 'B':
			if(GetKeyState(VK_CONTROL) < 0)
				OnMenuCubicCtrl();
			else
				OnMenuCubic();
			return TRUE;
		case 'C':
			if(status == statPlot)
			{
				status = statEnd;
				ptActRF[index] = ptActRF[0];
				index++;
				InitPH();
				DrawLine(index-1, plotSpline);
			}
			return TRUE;
		}
	}
	return CView::PreTranslateMessage(pMsg);
}

/*
@FN:		ONMOUSEMOVE
@BRIEF:		FUNCTION CALLED WHEN MOVE MOUSE
			STATUS BAR ALWAYS SHOWS THE COORDINATE OF CURSOR
			WHEN STATUS == PLOT		PREVIEW THE SPLINE WITH THE POINT WHICH CURSOR IS POINTING
			WHEN STATUS == OFFSET	MOVE THE WHOLE SPLINE
			WHEN STATUS == END, INDEX != 0
									DETECT ALL OF CONTROL POINTS, ACTIVE THE POINT WHICH IS POINTED BY CURSOR
@PARAM:		NFLAGS<UINT>
			POINT<CPOINT>
@RETURN:
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		DEC-23-2013
*/
void CPlanarPH_mfcView::OnMouseMove(UINT nFlags, CPoint point)
{
	int i;
	CMainFrame* pFrame = (CMainFrame*)AfxGetApp()->m_pMainWnd;
	CStatusBar* pStatusBar = (CStatusBar*)&pFrame->m_wndStatusBar;
	CString str;   
	str.Format("(%.2f, %.2f)", getActCoord(point).xPos, getActCoord(point).yPos);

	switch(status)
	{
	case statPlot:
		ptActRF[index] = getActCoord(point);
		index++;
		if(isMenuHermite)
		{
			if(index == 4)
				InitPH();
		}
		else
			InitPH();
		DrawLine(index-1, plotSpline);
		index--;
		if(isMenuHermite)
			str.Format("(%.2f, %.2f)", getActCoord(point).xPos, getActCoord(point).yPos);
		else
			str.Format("(%.2f, %.2f), ITER:%d", getActCoord(point).xPos, getActCoord(point).yPos, getIter(pph));
		break;
	case statTrans:
		if(isLBDown)
		{
			o_ActRF = CPointPH(o_ActRFTrans.xPos+(double)point.x-LBDownPt.xPos, 
							   o_ActRFTrans.yPos+(double)point.y-LBDownPt.yPos);
			if(index == 0)
			{
				DrawBG();
				glFlush();
			}
			else
				DrawLine(index-1, plotSpline);
		}
		break;
	case statEnd:
		for(i = 0; i < index; i++)
			if(Distance(point, getAbsCoord(ptActRF[i])) < 15.0)
			{
				DrawLine(i, plotPoint);
				pointing = i;
				str.Format("P%d(%.2f, %.2f)", pointing, ptActRF[pointing].xPos, ptActRF[pointing].yPos);
				break;
			}
		if(i == index && pointing != -1)
		{
			DrawLine(index-1, plotSpline);
			pointing = -1;
		}
		if(pointing == -1)
		{
			if(isMenuHermite)
				str.Format("(%.2f, %.2f), Bending Energy:%lf", 
						   getActCoord(point).xPos, getActCoord(point).yPos, getEnergy(pph, 1));
			else
				str.Format("(%.2f, %.2f), ITER:%d, Bending Energy:%lf", 
						   getActCoord(point).xPos, getActCoord(point).yPos, getIter(pph), pph.getEnergy());
		}
		break;
	default:
		break;
	}
	pStatusBar->SetWindowText(str); 

	CView::OnMouseMove(nFlags, point);
}

/*
@FN:		ONLBUTTONDOWN
@BRIEF:		FUNCTION CALLED WHEN PRESS THE LEFT BUTTON OF MOUSE
			WHEN STATUS == OFFSET	BEGIN TO MOVE THE WHOLE SPLINE
			WHEN INDEX == 0 OR STATUS == PLOT
									KEEP THE PLOT STATUS TRUE
									SET UP THE POINT AS A NEW ONE IN CONTROL POINTS ARRAY
@PARAM:		NFLAGS<UINT>
			POINT<CPOINT>
@RETURN:
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		DEC-22-2013
*/
void CPlanarPH_mfcView::OnLButtonDown(UINT nFlags, CPoint point)    
{   
	switch(status)
	{
	case statTrans:
		isLBDown = TRUE;
		LBDownPt = point;
		break;
	case statNew:
		status = statPlot;
	case statPlot:
		ptActRF[index] = getActCoord(point);
		index++;
		if(isMenuHermite)
		{
			if(index == 4)
			{
				status = statEnd;
				InitPH();
			}
		}
		else
			if(index > 2)
				InitPH();
		DrawLine(index-1, plotPoint);
		break;
	default:
		break;
	}

    CView::OnLButtonDown(nFlags, point);   
}

/*
@FN:		ONLBUTTONUP
@BRIEF:		FUNCTION CALLED WHEN OFFSET STOP, WORKS ONLY IF STATUS == OFFSET
			MAKE THE LBUTTONDOWN STATUS FALSE
			SAVE THE ORIGIN DATA IN O_ACTRFOFFSET FOR NEXT OFFSET
@PARAM:		NFLAGS<UINT>
			POINT<CPOINT>
@RETURN:
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		DEC-22-2013
*/
void CPlanarPH_mfcView::OnLButtonUp(UINT nFlags, CPoint point)   
{
	if(status == statTrans)
	{
		isLBDown = FALSE;
		o_ActRFTrans = o_ActRF;
	}
	
	CView::OnLButtonUp(nFlags, point);   
}

/*
@FN:		ONRBUTTONDOWN
@BRIEF:		FUNCTION CALLED WHEN CLICK THE RIGHT BUTTON OF MOUSE
			WHEN STATUS == PLOT		FINISH PLOTTING AND PLOT THE FINAL SPLINE;
			WHEN STATUS == END		IF THE CURSOR IS POINTING ON ONE POINT, POPUP THE WINDOW 
									FOR ADJUSTING THIS SINGLE POINT, GET THE NEW COORDINATE AND RENEW THE SPLINE
@PARAM:		NFLAGS<UINT>
			POINT<CPOINT>
@RETURN:
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		DEC-22-2013
*/
void CPlanarPH_mfcView::OnRButtonDown(UINT nFlags, CPoint point)    
{
	switch(status)
	{
	case statPlot:
		if(!isMenuHermite)
		{
			status = statEnd;
			InitPH();
			DrawLine(index-1, plotSpline);
		}
		break;
	case statEnd:
		if(pointing > -1)
		{
			CAdjWindowS m_Aws(ptActRF[pointing], this);				//Adjust Window For Single Point
			if(m_Aws.DoModal() == IDOK)
			{
				ptActRF[pointing] = m_Aws.ReturnPoint();
				InitPH();
				DrawLine(index-1, plotSpline);
			}
		}
		break;
	default:
		break;
	}
	
    CView::OnRButtonDown(nFlags, point);   
}

void CPlanarPH_mfcView::printString(const char* str)
{
//    static int firstCall = 1;

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

void CPlanarPH_mfcView::buildFont()
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

/*
@FN:		DRAWLINE
@BRIEF:		PLOT CURVES AND POINTS USING OPENGL
@PARAM:		PT<CPOINT []>:
				CONTROL POINTS ARRAY
			NINDEX<INT>:
				THE POINT NUMBER WHEN FLAG == 1; HOW MANY POINTS WHEN FLAG == 2
			FLAG<INT>:
				JUST PLOT THE POINT WHEN FLAG == 1; PLOT THE WHOLE SPLINE AND CONTROL POINTS WHEN FLAG == 2
@RETURN:
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		NOV-26-2013
*/
void CPlanarPH_mfcView::DrawLine(int nIndex, int flag)
{	
	int i, j;
	double t;
	CPointPH ptTemp;

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glLineWidth(1.5);
	glPointSize(8);
	
	switch(flag)
	{
	case plotPoint:
		if(index == 0)
		{
			DrawBG();
		}
		glBegin(GL_POINTS); 
		if(status == statPlot)			
				glColor3fv(colors[0]);
			else
				glColor3fv(colors[1]);
			glVertex2f(getAbsCoord(ptActRF[nIndex]).xPos, (float)height-getAbsCoord(ptActRF[nIndex]).yPos);
		glEnd();
		break;
	case plotSpline:
		DrawBG();
		if(isMenuHermite)
		{
			switch(nIndex)
			{
			case 0:
				glColor3fv(colors[0]);
				glBegin(GL_POINTS); 
					glVertex2f(getAbsCoord(ptActRF[0]).xPos, (float)height-getAbsCoord(ptActRF[0]).yPos);
				glEnd();
				break;
			case 1:
			case 2:
				glColor3fv(colors[4]);
				glBegin(GL_LINES); 
					glVertex2f(getAbsCoord(ptActRF[0]).xPos, (float)height-getAbsCoord(ptActRF[0]).yPos);
					glVertex2f(getAbsCoord(ptActRF[1]).xPos, (float)height-getAbsCoord(ptActRF[1]).yPos);
				glEnd();
				glColor3fv(colors[0]);
				glBegin(GL_POINTS); 
					glVertex2f(getAbsCoord(ptActRF[0]).xPos, (float)height-getAbsCoord(ptActRF[0]).yPos);
					glVertex2f(getAbsCoord(ptActRF[1]).xPos, (float)height-getAbsCoord(ptActRF[1]).yPos);
					if(nIndex == 2)
						glVertex2f(getAbsCoord(ptActRF[2]).xPos, (float)height-getAbsCoord(ptActRF[2]).yPos);
				glEnd();
				break;
			case 3:
				glColor3fv(colors[4]);
				glBegin(GL_LINES); 
					for(i = 0; i <= nIndex; i++)
						glVertex2f(getAbsCoord(ptActRF[i]).xPos, (float)height-getAbsCoord(ptActRF[i]).yPos);
				glEnd();
				glColor3fv(colors[2]);		
				glBegin(GL_LINE_STRIP);							
					for(i = 0; i <= Step; i++)
					{
						t = 1.0/Step*i;
						ptTemp = CPointPH(real(beval(Dgr, pph.spline[1].p, t)), imag(beval(Dgr, pph.spline[1].p, t)));
						glVertex2f(getAbsCoord(ptTemp).xPos, (float)height-getAbsCoord(ptTemp).yPos);
					}
				glEnd();
				glColor3fv(colors[0]);
				glBegin(GL_POINTS); 
					for(i = 0; i <= nIndex; i++)
						glVertex2f(getAbsCoord(ptActRF[i]).xPos, (float)height-getAbsCoord(ptActRF[i]).yPos);
				glEnd();
				break;
			default:
				break;
			}
		}
		else
		{
		if(nIndex > 1)
		{
			if(isMenuPh)
			{
				if(isMenuPhCtrl)
				{
					glLineWidth(1);
					glPointSize(5);
					glBegin(GL_LINE_STRIP);
						glColor3fv(colors[4]);
						for(i = 1; i <= nIndex; i++)
							for(j = 0; j <= Dgr; j++)
							{
								ptTemp = CPointPH(real(getCtrlPt(pph, i, j)), imag(getCtrlPt(pph, i, j)));
								glVertex2f(getAbsCoord(ptTemp).xPos, height-getAbsCoord(ptTemp).yPos);
							}
					glEnd();
					glBegin(GL_POINTS);   
						glColor3fv(colors[4]);
						for(i = 1; i <= nIndex; i++)
							for(j = 0; j <= Dgr; j++)
							{
								ptTemp = CPointPH(real(getCtrlPt(pph, i, j)), imag(getCtrlPt(pph, i, j)));
								glVertex2f(getAbsCoord(ptTemp).xPos, height-getAbsCoord(ptTemp).yPos);
							}
					glEnd(); 
					glLineWidth(1.5);
					glPointSize(8);
				}
				glBegin(GL_LINE_STRIP);
					glColor3fv(colors[2]);									
					for(i = 1; i <= nIndex; i++)
					{
						for(j = 0; j <= Step; j++)
						{
							t = 1.0/Step*j;
							ptTemp = CPointPH(real(beval(Dgr, pph.spline[i].p, t)), imag(beval(Dgr, pph.spline[i].p, t)));
							glVertex2f(getAbsCoord(ptTemp).xPos, (float)height-getAbsCoord(ptTemp).yPos);
						}
					} 
				glEnd();
			}
			if(isMenuCubic)
			{
				CCubic m_cubic(ptActRF, nIndex);

				if(isMenuCubicCtrl)
				{
					glLineWidth(1);
					glPointSize(5);
					glBegin(GL_LINE_STRIP);
						glColor3fv(colors[5]);
						for(i = 0; i <= nIndex; i++)
						{
							if(ptActRF[0] == ptActRF[nIndex] && i == nIndex)
								ptTemp = CPointPH(real(m_cubic.dp[0]), imag(m_cubic.dp[0]));
							else
								ptTemp = CPointPH(real(m_cubic.dp[i]), imag(m_cubic.dp[i]));
							glVertex2f(getAbsCoord(ptTemp).xPos, height-getAbsCoord(ptTemp).yPos);
						}
					glEnd();
					glBegin(GL_POINTS);   
						glColor3fv(colors[5]);
						for(i = 0; i <= nIndex; i++)
						{
							if(ptActRF[0] == ptActRF[nIndex] && i == nIndex)
								ptTemp = CPointPH(real(m_cubic.dp[0]), imag(m_cubic.dp[0]));
							else
								ptTemp = CPointPH(real(m_cubic.dp[i]), imag(m_cubic.dp[i]));
							glVertex2f(getAbsCoord(ptTemp).xPos, height-getAbsCoord(ptTemp).yPos);
						}
					glEnd(); 
					glLineWidth(1.5);
					glPointSize(8);
				}
				glBegin(GL_LINE_STRIP);
					glColor3fv(colors[3]);
					if(ptActRF[0] == ptActRF[nIndex])
						for(i = 3; i <= nIndex+2; i++)
							for(j = 0; j <= Step; j++)
							{
								t= (i-3.0)+1.0/Step*j;
								ptTemp = CPointPH(real(getSpline(m_cubic, i, t)), imag(getSpline(m_cubic, i, t)));
								glVertex2f(getAbsCoord(ptTemp).xPos, (float)height-getAbsCoord(ptTemp).yPos);
							}
					else
						for(i = 1; i <= nIndex; i++)
							for(j = 0; j <= Step; j++)
							{
								t = 1.0/Step*j;
								ptTemp = CPointPH(real(getSpline(m_cubic, i, t)), imag(getSpline(m_cubic, i, t)));
								glVertex2f(getAbsCoord(ptTemp).xPos, (float)height-getAbsCoord(ptTemp).yPos);
							}
				glEnd();
			}
		}
		if(nIndex == 1)
		{
			glBegin(GL_LINE_STRIP);
				glColor3fv(colors[2]);
				glVertex2f(getAbsCoord(ptActRF[0]).xPos, (float)height-getAbsCoord(ptActRF[0]).yPos);
				glVertex2f(getAbsCoord(ptActRF[1]).xPos, (float)height-getAbsCoord(ptActRF[1]).yPos);
			glEnd();
		}
		
		glBegin(GL_POINTS);   
			glColor3fv(colors[2]);
			for(i = 0; i < nIndex; i++)
				glVertex2f(getAbsCoord(ptActRF[i]).xPos, (float)height-getAbsCoord(ptActRF[i]).yPos);
			if(status == statPlot)
				glColor3fv(colors[1]);
			glVertex2f(getAbsCoord(ptActRF[nIndex]).xPos, (float)height-getAbsCoord(ptActRF[nIndex]).yPos);
		glEnd();
		}
		if(offsetNum != 0)
			for(i = 0; i < offsetNum; i++)
				DrawLine(i, plotOffset);
		break;
	case plotOffset:
		if(isMenuPh)
		{
			if(isMenuPhCtrl)
			{
				glLineWidth(1);
				glPointSize(5);
				glBegin(GL_LINE_STRIP);
					glColor3fv(colors[4]);
					for(i = 1; i <= index-1; i++)
						for(j = 0; j <= 2*Dgr-1; j++)
						{
							ptTemp = CPointPH(real(getOffset(pph, i, offsetD[nIndex])[j]), imag(getOffset(pph, i, offsetD[nIndex])[j]));
							glVertex2f(getAbsCoord(ptTemp).xPos, height-getAbsCoord(ptTemp).yPos);
						}
				glEnd();
				glBegin(GL_POINTS);   
					glColor3fv(colors[4]);
					for(i = 1; i <= index-1; i++)
						for(j = 0; j <= 2*Dgr-1; j++)
						{
							ptTemp = CPointPH(real(getOffset(pph, i, offsetD[nIndex])[j]), imag(getOffset(pph, i, offsetD[nIndex])[j]));
							glVertex2f(getAbsCoord(ptTemp).xPos, height-getAbsCoord(ptTemp).yPos);
						}
				glEnd(); 
				glLineWidth(1.5);
				glPointSize(8);
			}
			glBegin(GL_LINE_STRIP);
				glColor3fv(colors[2]);									
				for(i = 1; i <= index-1; i++)
				{
					for(j = 0; j <= Step; j++)
					{
						t = 1.0/Step*j;
						ptTemp = CPointPH(real(beval(2*Dgr-1, getOffset(pph, i, offsetD[nIndex]), t)), imag(beval(2*Dgr-1, getOffset(pph, i, offsetD[nIndex]), t)));
						glVertex2f(getAbsCoord(ptTemp).xPos, (float)height-getAbsCoord(ptTemp).yPos);
					}
				} 
			glEnd();
		}
		break;
	default:
		DrawBG();
		break;
	}

	glFlush();
}

/*
@FN:		DISTANCE
@BRIEF:		CULCULATE THE DISTANCE BETWEEN TWO POINTS
@PARAM:		P1<CPOINT>:
			P2<CPOINT>:
				POINTS NEED TO BE CULCULATED
@RETURN:	<DOUBLE>
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		NOV-15-2013
*/
double CPlanarPH_mfcView::Distance(CPointPH p1, CPointPH p2)
{
	return sqrt(pow(p1.xPos-p2.xPos, 2)+pow(p1.yPos-p2.yPos, 2));
}

/*
@FN:		GETABSCOORD
@BRIEF:		GET ABSOLUTE COORDINATE OF CERTAIN POINT
@PARAM:		PT<CPOINT>:
@RETURN:	absCoord<CPOINT>
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		DEC-23-2013
*/
CPointPH CPlanarPH_mfcView::getAbsCoord(CPointPH pt)
{
	CPointPH absCoord(pt.xPos+o_ActRF.xPos, o_ActRF.yPos-pt.yPos);
	return absCoord;
}

/*
@FN:		GETACTCOORD
@BRIEF:		GET ACTUAL COORDINATE OF CERTAIN POINT
@PARAM:		PT<CPOINT>:
@RETURN:	actCoord<CPOINT>
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		DEC-23-2013
*/
CPointPH CPlanarPH_mfcView::getActCoord(CPointPH pt)
{
	CPointPH actCoord(pt.xPos-o_ActRF.xPos, o_ActRF.yPos-pt.yPos);
	return actCoord;
}

/*
@FN:		CALPHCTRL
@BRIEF:		CALCULATE PH CONTROL POINTS
@PARAM:		
@RETURN:	
@WARNING:
@AUTHOR:	BOHAN DONG
@DATE:		DEC-23-2013
*/
void CPlanarPH_mfcView::InitPH()
{
	if(isMenuHermite)
		pph = PlanarPH(ptActRF);
	else
		pph = PlanarPH(ptActRF, index-1);
}