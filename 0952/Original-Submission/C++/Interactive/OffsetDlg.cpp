#include "stdafx.h"
#include "PlanarPH_mfc.h"
#include "OffsetDlg.h"
#include "afxdialogex.h"

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

IMPLEMENT_DYNAMIC(COffsetDlg, CDialogEx)

COffsetDlg::COffsetDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(COffsetDlg::IDD, pParent)
{
}

COffsetDlg::~COffsetDlg()
{
}

void COffsetDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(COffsetDlg, CDialogEx)
	ON_BN_CLICKED(IDOK, &COffsetDlg::OnBnClickedOk)
END_MESSAGE_MAP()

void COffsetDlg::OnBnClickedOk()
{
	o_ScaleNum = (CEdit *)GetDlgItem(IDC_EDIT1);
	o_ScaleNum->GetWindowText(o_str);
	d = atof(o_str);
	
	CDialogEx::OnOK();
}