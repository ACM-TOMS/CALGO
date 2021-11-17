#pragma once
#include "AnalyseView.h"

//CAnalyseDlg Dialog

class CAnalyseDlg : public CDialogEx
{
	DECLARE_DYNAMIC(CAnalyseDlg)

public:
	CAnalyseDlg(CWnd* pParent = NULL);   
	virtual ~CAnalyseDlg();

	enum { IDD = IDD_DIALOG3 };

protected:
	CTabCtrl m_tab;
	
	int m_CurSelTab; 
	CDialog* pDialog[2];

	virtual void DoDataExchange(CDataExchange* pDX);  
	virtual BOOL OnInitDialog();
	void OnTcnSelchangeTab(NMHDR *pNMHDR, LRESULT *pResult);

	DECLARE_MESSAGE_MAP()
};
