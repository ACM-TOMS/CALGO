#pragma once
#include "afxext.h"

class CMainFrame : public CFrameWnd
{
protected:
	CMainFrame();
	DECLARE_DYNCREATE(CMainFrame)

public:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	virtual ~CMainFrame();
	CStatusBar m_wndStatusBar;

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	CToolBar m_wndToolBar;

	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

	DECLARE_MESSAGE_MAP()
};


