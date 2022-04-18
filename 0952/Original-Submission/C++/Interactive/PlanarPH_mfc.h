#pragma once
#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif
#include "resource.h"

class CPlanarPH_mfcApp : public CWinApp
{
public:
	CPlanarPH_mfcApp();

	virtual BOOL InitInstance();
	virtual int ExitInstance();

	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};
extern CPlanarPH_mfcApp theApp;
