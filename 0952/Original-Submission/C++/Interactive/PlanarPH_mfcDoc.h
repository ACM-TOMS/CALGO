#pragma once

class CPlanarPH_mfcDoc : public CDocument
{
protected:
	CPlanarPH_mfcDoc();
	DECLARE_DYNCREATE(CPlanarPH_mfcDoc)
	DECLARE_MESSAGE_MAP()

public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif 
	virtual ~CPlanarPH_mfcDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif
#ifdef SHARED_HANDLERS
	void SetSearchContent(const CString& value);
#endif 
};
