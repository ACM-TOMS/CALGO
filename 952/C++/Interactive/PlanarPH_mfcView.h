#pragma once
#include "PlanarPH.h"
#include "PointPH.h"
#include "Cubic.h"
#include "AdjWindowS.h"
#include "AdjWindowA.h"
#include "Analyse.h"
#include "OffsetDlg.h"
#include <stdarg.h>

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

class CPlanarPH_mfcView : public CView
{
protected:
	CPlanarPH_mfcView();
	DECLARE_DYNCREATE(CPlanarPH_mfcView)

	virtual ~CPlanarPH_mfcView();

	HGLRC m_hRC;					//OpenGL Rendering Context 
	CClientDC* m_pDC;				//OpenGL Device Context
	GLsizei width, height; 
	PlanarPH pph;
	int index, pointing, status;		//Pointing is pointing on which point. = -1 if pointing nothing
	CPointPH o_ActRF, o_ActRFTrans;		//Origin of Actual Reference Frame
	CPointPH ptActRF[Max+1];			//Points in Actual Reference Frame
	BOOL isMenuHermite, isMenuPh, isMenuCubic, isMenuPhCtrl, isMenuCubicCtrl;		//Status of Menu
	CPointPH LBDownPt;
	BOOL isLBDown;
	double offsetD[10];
	int offsetNum;
	GLuint base;
	
	void DrawBG();
	void DrawLine(int nIndex, int flag);
	void InitPH();
	BOOL SetDCPixelFormat();
	void printString(const char* str);
	void buildFont();
	double Distance(CPointPH p1, CPointPH p2);
	CPointPH getAbsCoord(CPointPH pt);
	CPointPH getActCoord(CPointPH pt);

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	DECLARE_MESSAGE_MAP()

	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct); 
	afx_msg void OnDestroy(); 
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);

	void OnMenuHermite();
	void OnMenuPh();
	void OnMenuCubic();
	void OnMenuPhCtrl();
	void OnMenuCubicCtrl();
	void OnMenuNew();
	void OnMenuEditMulpt();
	void OnMenuTranslate();
	void OnMenuOffset();
	void OnMenuParaspeed();
	void OnMenuArclength();
	void OnUpdateMenuPh(CCmdUI* pCmdUI);	//Function for Ph Curve Menu
	void OnUpdateMenuCubic(CCmdUI* pCmdUI);	//Function for Cubic Curve Menu
	void OnUpdateMenuPhCtrl(CCmdUI* pCmdUI);
	void OnUpdateMenuCubicCtrl(CCmdUI* pCmdUI);
	BOOL PreTranslateMessage(MSG* pMsg);

public:
	CPlanarPH_mfcDoc* GetDocument() const;

	virtual void OnDraw(CDC* pDC);			// overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
};

#ifndef _DEBUG  // debug version in PlanarPH_mfcView.cpp
inline CPlanarPH_mfcDoc* CPlanarPH_mfcView::GetDocument() const
   { return reinterpret_cast<CPlanarPH_mfcDoc*>(m_pDocument); }
#endif