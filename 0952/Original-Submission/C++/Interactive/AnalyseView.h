#pragma once

//CAnalyseView Dialog
#include "PlanarPH.h"
#include <GL/glut.h>

class CAnalyseView : public CDialogEx
{
	DECLARE_DYNAMIC(CAnalyseView)

public:
	CAnalyseView(CWnd* pParent = NULL); 
	CAnalyseView(CPoint pt[], int nIndex, CWnd* pParent = NULL);  
	virtual ~CAnalyseView();

	enum { IDD = IDD_DIALOG4 };

protected:
	HGLRC m_hRC;					//OpenGL Rendering Context 
	CClientDC* m_pDC;				//OpenGL Device Context
	GLsizei width, height;  
    GLfloat aspect, m_xPos, m_yPos;
	CPoint pt_Analyse[N+1];
	int index_Analyse;

	afx_msg void OnDestroy(); 
	afx_msg void OnPaint();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
//	afx_msg void OnSize(UINT nType, int cx, int cy);
	BOOL InitOpengl();
	
	virtual void DoDataExchange(CDataExchange* pDX);   
	virtual BOOL OnInitDialog();
	BOOL SetupPixelFormat();
	void DrawScene();

	DECLARE_MESSAGE_MAP()
public:
	
};
