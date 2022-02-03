//AnalyseView.cpp:
//

#include "stdafx.h"
#include "PlanarPH_mfc.h"
#include "AnalyseView.h"
#include "afxdialogex.h"

//CAnalyseView Dialog

IMPLEMENT_DYNAMIC(CAnalyseView, CDialogEx)

CAnalyseView::CAnalyseView(CWnd* pParent /*=NULL*/)
	: CDialogEx(CAnalyseView::IDD, pParent)
{
}

CAnalyseView::CAnalyseView(CPoint pt[], int nIndex, CWnd* pParent /*=NULL*/)
	: CDialogEx(CAnalyseView::IDD, pParent)
{
	index_Analyse = nIndex;
	for(int i = 0; i <= index_Analyse; i++)
		pt_Analyse[i] = pt[i];
}

CAnalyseView::~CAnalyseView()
{
}

void CAnalyseView::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAnalyseView, CDialogEx)
	ON_WM_PAINT()
    ON_WM_DESTROY()   
	ON_WM_ERASEBKGND()
END_MESSAGE_MAP()

//CAnalyseView
BOOL CAnalyseView::OnInitDialog()
{
	InitOpengl();

	return CDialog::OnInitDialog();
}

BOOL CAnalyseView::InitOpengl()
{
	m_pDC = new CClientDC(GetDlgItem(IDC_PIC));  
    if(m_pDC == NULL)  
        return FALSE;  
    if(!SetupPixelFormat())  
        return FALSE;   
    m_hRC = wglCreateContext(m_pDC->GetSafeHdc());  
    if(m_hRC == 0)  
		return FALSE;   
    if(wglMakeCurrent(m_pDC->GetSafeHdc(), m_hRC) == FALSE)  
		return FALSE;  
  
    return TRUE;  
}

BOOL CAnalyseView::SetupPixelFormat() 
{
	static PIXELFORMATDESCRIPTOR pfd = {  
        sizeof(PIXELFORMATDESCRIPTOR),  
        1,  
        PFD_DRAW_TO_WINDOW|  
        PFD_SUPPORT_OPENGL|   
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
        0,0,0  
    };  

    int m_nPixelFormat = ChoosePixelFormat(m_pDC->GetSafeHdc(), &pfd);
    if(m_nPixelFormat==0)  
	{
		MessageBox("Choose Pixel Format failed"); 
		return FALSE;  
	}
    if(SetPixelFormat(m_pDC->GetSafeHdc(), m_nPixelFormat, &pfd) == FALSE)
	{
		MessageBox("Set Pixel Format failed");
		return FALSE;  
	}
 
    return TRUE;  
}

void CAnalyseView::OnDestroy()    
{   
    CDialogEx::OnDestroy();  
  
    if(wglMakeCurrent(NULL, NULL) == FALSE) 
		MessageBox("Cannot release RC");  
    if(wglDeleteContext(m_hRC) == FALSE)
		MessageBox("Cannot delete RC");  
    if(m_pDC)  
		delete m_pDC;  
    m_pDC=NULL;
}   

BOOL CAnalyseView::OnEraseBkgnd(CDC* pDC)  
{  
	return TRUE; 
}

void CAnalyseView::OnPaint()
{
	DrawScene();

	CDialogEx::OnPaint();
}

void CAnalyseView::DrawScene()
{
	int i, j;
	double t;
	PlanarPH m_pph(pt_Analyse, index_Analyse);
	if(pt_Analyse[index_Analyse] == pt_Analyse[0])
		m_pph.pphControl(TRUE);
	else m_pph.pphControl(FALSE);

	glClearColor(1.0, 1.0, 1.0, 1.0); 
	glClear(GL_COLOR_BUFFER_BIT);  
    glLoadIdentity();  

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glLineWidth(1.5);
      
    glBegin(GL_LINE_STRIP);  
        glColor3f(0.0, 0.0, 0.0);  
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(1.0, 0.0, 0.0);
		glVertex3f(0.0, 1.0, 0.0);
        for(i = 1; i <= index_Analyse; i++)
			for(j = 1; j < 200; j++)
			{
				t = 1.0/200*j;
				glVertex2f(200*(i-1)+j, m_pph.pphParaspeed(i, t)+100);
			}
    glEnd();  

	glFlush();
}

/*
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

	if(index != 0)
		DrawLine(m_ptDisplay, index-1, 2);
} 
*/