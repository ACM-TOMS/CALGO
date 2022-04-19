#include "stdafx.h"
#include "PlanarPH_mfc.h"
#include "AdjWindowA.h"
#include "afxdialogex.h"

/*
Author: Bohan Dong, University of California, Davis(2014)
This code is intended only for experimental and research purposes. 
No commercial use may be made without the author's express permision. 
*/

IMPLEMENT_DYNAMIC(CAdjWindowA, CDialogEx)

CAdjWindowA::CAdjWindowA(CPointPH m_pt[], int index, CWnd* pParent /*=NULL*/)
	: CDialogEx(CAdjWindowA::IDD, pParent)
{
	a_index = index;
	for(int i = 0; i < a_index; i++)
		a_pt[i] = m_pt[i];
}

CAdjWindowA::~CAdjWindowA()
{
}

void CAdjWindowA::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_LIST1, a_ctllist);
}

BEGIN_MESSAGE_MAP(CAdjWindowA, CDialogEx)
	ON_NOTIFY(NM_CLICK, IDC_LIST1, OnNMClickList1)
END_MESSAGE_MAP()

BOOL CAdjWindowA::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	ListView_SetExtendedListViewStyle(::GetDlgItem(m_hWnd, IDC_LIST1), LVS_EX_FULLROWSELECT | LVS_EX_GRIDLINES); 
	InsertItems();  
	::ShowWindow(::GetDlgItem(m_hWnd, IDC_EDIT1), SW_HIDE); 

	return TRUE;
}

void CAdjWindowA::InsertItems()
{
	HWND hWnd = ::GetDlgItem(m_hWnd, IDC_LIST1);
	LVCOLUMN list;
	list.mask =	LVCF_TEXT | LVCF_WIDTH | LVCF_FMT | LVCF_SUBITEM;
	list.fmt = LVCFMT_LEFT;
	list.cx = 40; 
	list.pszText = "Point";
	list.iSubItem = 0;
	::SendMessage(hWnd, LVM_INSERTCOLUMN, (WPARAM)0, (WPARAM)&list);
	list.cx = 80;  
	list.pszText = "X"; 
	list.iSubItem = 1;      
	::SendMessage(hWnd, LVM_INSERTCOLUMN, (WPARAM)1, (WPARAM)&list);
	list.cx = 80;  
	list.pszText = "Y";
	list.iSubItem = 2;
	::SendMessage(hWnd, LVM_INSERTCOLUMN, (WPARAM)2, (WPARAM)&list);

	for(int i = 0; i < a_index; i++)
	{
		a_str.Format("P%d", i);
		a_strx.Format("%.2f", a_pt[i].xPos);
		a_stry.Format("%.2f", a_pt[i].yPos);
		SetCell(hWnd, a_str, i, 0);  
		SetCell(hWnd, a_strx, i, 1); 
		SetCell(hWnd, a_stry, i, 2); 
	}
}

void CAdjWindowA::SetCell(HWND hWnd, CString value, int nRow, int nCol)
{
	TCHAR szString[256]; 
	wsprintf(szString, value, 0);
	LVITEM lvItem;  
	lvItem.mask = LVIF_TEXT; 
	lvItem.iItem = nRow; 
	lvItem.pszText = szString; 
	lvItem.iSubItem = nCol; 
	if(nCol > 0) 
		::SendMessage(hWnd, LVM_SETITEM, (WPARAM)0, (WPARAM)&lvItem); 
	else 
		ListView_InsertItem(hWnd, &lvItem);
}

CString CAdjWindowA::GetItemText(HWND hWnd, int nItem, int nSubItem)const 
{
	LVITEM lvi;  
	memset(&lvi, 0, sizeof(LVITEM)); 
	lvi.iSubItem = nSubItem; 
	CString str; 
	int nLen = 128; 
	int nRes; 
	
	nLen *= 2;  
	lvi.cchTextMax = nLen;  
	lvi.pszText = str.GetBufferSetLength(nLen); 
	nRes = (int)::SendMessage(hWnd, LVM_GETITEMTEXT, (WPARAM)nItem, (LPARAM)&lvi); 
	str.ReleaseBuffer(); 
	return str; 
}

void CAdjWindowA::OnNMClickList1(NMHDR *pNMHDR, LRESULT *pResult)
{
	Invalidate();
	HWND hWnd1 = ::GetDlgItem(m_hWnd, IDC_LIST1);
	LPNMITEMACTIVATE temp = (LPNMITEMACTIVATE)pNMHDR;
	RECT rect; 
	nItem = temp->iItem; 
	nSubItem = temp->iSubItem;  
	if(nSubItem == 0 || nSubItem == -1 || nItem == -1) 
		return;
	CString str = GetItemText(hWnd1, nItem , nSubItem);
	RECT rect1, rect2; 
	ListView_GetSubItemRect(hWnd1, temp->iItem, temp->iSubItem, LVIR_BOUNDS, &rect); 
	::GetWindowRect(temp->hdr.hwndFrom, &rect1);         
	::GetWindowRect(m_hWnd, &rect2);   
	int x = rect1.left-rect2.left; 
	int y = rect1.top-rect2.top;   
	if(nItem != -1)
		::SetWindowPos(::GetDlgItem(m_hWnd, IDC_EDIT1), HWND_TOP, rect.left+x-6, rect.top+13, rect.right-rect.left, rect.bottom-rect.top, NULL);      
	::ShowWindow(::GetDlgItem(m_hWnd, IDC_EDIT1), SW_SHOW);     
	::SetFocus(::GetDlgItem(m_hWnd, IDC_EDIT1));
	::Rectangle(::GetDC(temp->hdr.hwndFrom), rect.left, rect.top, rect.right, rect.bottom);
	::SetWindowText(::GetDlgItem(m_hWnd, IDC_EDIT1), str); 
	
	*pResult = 0;
}

void CAdjWindowA::OnOK()
{
	CWnd* pwndCtrl = GetFocus(); 
	int ctrl_ID = pwndCtrl->GetDlgCtrlID(); 
	CString str;  
	switch(ctrl_ID)     
	{ 
	case IDC_EDIT1: 
		GetDlgItemText(IDC_EDIT1, str);      
		SetCell(::GetDlgItem(m_hWnd, IDC_LIST1), str, nItem, nSubItem);
		if(nSubItem == 1)
			a_pt[nItem].xPos = atof(str);
		if(nSubItem == 2)
			a_pt[nItem].yPos = atof(str);
		::SendDlgItemMessage(m_hWnd, IDC_EDIT1, WM_KILLFOCUS, 0, 0);
		::ShowWindow(::GetDlgItem(m_hWnd, IDC_EDIT1), SW_HIDE); 
		break;      
	default: 
		CDialogEx::OnOK();
		break;     
	} 
}

CPointPH CAdjWindowA::ReturnPoint(int i)
{
	if(i < a_index)
		return a_pt[i];
	else
		return (CPointPH)(0.0, 0.0);
}