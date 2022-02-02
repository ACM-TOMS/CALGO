//AnalyseDlg.cpp:
//

#include "stdafx.h"
#include "PlanarPH_mfc.h"
#include "AnalyseDlg.h"
#include "afxdialogex.h"

// CAnalyseDlg Dialog

IMPLEMENT_DYNAMIC(CAnalyseDlg, CDialogEx)

CAnalyseDlg::CAnalyseDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CAnalyseDlg::IDD, pParent)
{
}

CAnalyseDlg::~CAnalyseDlg()
{
}

void CAnalyseDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);

	DDX_Control(pDX, IDC_TAB1, m_tab);
}

BEGIN_MESSAGE_MAP(CAnalyseDlg, CDialogEx)
	ON_NOTIFY(TCN_SELCHANGE, IDC_TAB1, OnTcnSelchangeTab)
END_MESSAGE_MAP()

//CAnalyseDlg
BOOL CAnalyseDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	CAnalyseView m_page1;
	CAnalyseView m_page2;

	m_tab.InsertItem(0, "Para Speed");
	m_tab.InsertItem(1, "Arc Length");

	CRect rc;
	m_tab.GetClientRect(rc);
	rc.top += 20;
	rc.bottom -= 0;
	rc.left += 0;
	rc.right -= 0;
	m_page1.MoveWindow(&rc);
	m_page2.MoveWindow(&rc);

	m_CurSelTab = 0;

	pDialog[0] = &m_page1;
	pDialog[1] = &m_page2;
	pDialog[0]->ShowWindow(SW_SHOW);
	pDialog[1]->ShowWindow(SW_HIDE);

	return TRUE;
}

void CAnalyseDlg::OnTcnSelchangeTab(NMHDR *pNMHDR, LRESULT *pResult)
{
    pDialog[m_CurSelTab]->ShowWindow(SW_HIDE);
    m_CurSelTab = m_tab.GetCurSel();
    pDialog[m_CurSelTab]->ShowWindow(SW_SHOW);
    *pResult = 0;
}