#pragma once

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

class COffsetDlg : public CDialogEx
{
	DECLARE_DYNAMIC(COffsetDlg)
protected:
	virtual void DoDataExchange(CDataExchange* pDX);

	DECLARE_MESSAGE_MAP()

public:
	COffsetDlg(CWnd* pParent = NULL);
	virtual ~COffsetDlg();

	CEdit *o_ScaleNum;
	CString o_str;
	double d;

	afx_msg void OnBnClickedOk();

	enum { IDD = IDD_DIALOG4 };
};
