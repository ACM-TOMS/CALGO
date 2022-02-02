#pragma once
#include "PlanarPH.h"
#include <GL/glut.h>

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

class CAnalyse : public CDialogEx
{
	DECLARE_DYNAMIC(CAnalyse)
protected:
	HGLRC m_oldRC, m_hRC;			//OpenGL Rendering Context 
	CClientDC* m_pDC;				//OpenGL Device Context
	HDC m_oldDC, m_hDC;
	PlanarPH a_pph;
	int a_stat, a_index;
	GLfloat margin;
	GLsizei width, height;
	BOOL pointing;
	GLfloat cursor_x;
	GLuint base;

	BOOL SetDCPixelFormat();
	virtual void DoDataExchange(CDataExchange* pDX); 
	void RenderScene();
	void printString(const char* str, ...);
	void buildFont();

	DECLARE_MESSAGE_MAP()

public:
	CAnalyse(PlanarPH pph, int index, int state, CWnd* pParent = NULL);   
	virtual ~CAnalyse();
	enum {IDD = IDD_DIALOG3};

	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnPaint();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
};
