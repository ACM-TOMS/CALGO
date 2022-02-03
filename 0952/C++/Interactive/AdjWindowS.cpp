#include "stdafx.h"
#include "PlanarPH_mfc.h"
#include "PlanarPH_mfcDoc.h"
#include "PlanarPH_mfcView.h"
#include "afxdialogex.h"
#include "AdjWindowS.h"

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

IMPLEMENT_DYNAMIC(CAdjWindowS, CDialogEx)

CAdjWindowS::CAdjWindowS(CPointPH m_pt, CWnd* pParent /*=NULL*/)
	: CDialogEx(CAdjWindowS::IDD, pParent)
{
	s_pt = m_pt;
}

CAdjWindowS::~CAdjWindowS()
{
}

void CAdjWindowS::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAdjWindowS, CDialogEx)
	ON_BN_CLICKED(IDOK, &CAdjWindowS::OnBnClickedOk)
END_MESSAGE_MAP()


BOOL CAdjWindowS::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	s_str.Format("%.2f", s_pt.xPos);
	s_ScaleName = (CEdit *)GetDlgItem(IDC_EDIT1);
	s_ScaleName->SetWindowTextA(s_str);
	s_str.Format("%.2f", s_pt.yPos);
	s_ScaleName = (CEdit *)GetDlgItem(IDC_EDIT2);
	s_ScaleName->SetWindowTextA(s_str);
	s_ScaleName = (CEdit *)GetDlgItem(IDC_EDIT3);
	s_ScaleName->SetWindowTextA("X:");
	s_ScaleName = (CEdit *)GetDlgItem(IDC_EDIT4);
	s_ScaleName->SetWindowTextA("Y:");

	return TRUE; 
}

void CAdjWindowS::OnBnClickedOk()
{
	double xNew, yNew;

	s_ScaleNum = (CEdit *)GetDlgItem(IDC_EDIT1);
	s_ScaleNum->GetWindowText(s_str);
	xNew = atof(s_str);
	s_ScaleNum = (CEdit *)GetDlgItem(IDC_EDIT2);
	s_ScaleNum->GetWindowText(s_str);
	yNew = atof(s_str);

	s_pt = CPointPH(xNew, yNew);

	CDialogEx::OnOK();
}

CPointPH CAdjWindowS::ReturnPoint()
{
	return s_pt;
}