#pragma once
#include "afxcmn.h"
#include "PointPH.h"
#include "PlanarPH.h"

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

class CAdjWindowA : public CDialogEx
{
	DECLARE_DYNAMIC(CAdjWindowA)
protected:
	virtual void DoDataExchange(CDataExchange* pDX); 
	afx_msg void OnOK();
	DECLARE_MESSAGE_MAP()

public:
	CAdjWindowA(CPointPH m_pt[], int index, CWnd* pParent = NULL);  
	virtual ~CAdjWindowA();
	enum { IDD = IDD_DIALOG2 };

	CListCtrl a_ctllist;
	CPointPH a_pt[Max+1];
	int a_index, nItem, nSubItem;
	CString a_str, a_strx, a_stry;

	virtual BOOL OnInitDialog();
	void InsertItems();
	CPointPH ReturnPoint(int i);
	void SetCell(HWND hWnd, CString value, int nRow, int nCol);
	CString GetItemText(HWND hWnd, int nItem, int nSubItem)const;
	afx_msg void OnNMClickList1(NMHDR *pNMHDR, LRESULT *pResult);
};
