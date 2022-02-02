#pragma once

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

class CAdjWindowS : public CDialogEx
{
	DECLARE_DYNAMIC(CAdjWindowS)
protected:
	virtual void DoDataExchange(CDataExchange* pDX); 
	DECLARE_MESSAGE_MAP()

public:
	CAdjWindowS(CPointPH m_pt, CWnd* pParent = NULL); 
	virtual ~CAdjWindowS();
	enum { IDD = IDD_DIALOG1 };

	CEdit *s_ScaleName, *s_ScaleNum;
	CPointPH s_pt;
	CString s_str;
	int s_index;

	CPointPH ReturnPoint();
	virtual BOOL OnInitDialog();
	afx_msg void OnBnClickedOk();
};
