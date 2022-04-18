#ifndef GDIPLUS_H
#define GDIPLUS_H

#ifdef WIN32

#include <windows.h>
#undef min
#undef max
namespace gdi
{
  #include <GdiPlus.h>
  using namespace Gdiplus;
}

static class GdiDeclare
{
  public:
    GdiDeclare()  { gdi::GdiplusStartup(&m_gdiplusToken, &gdiplusStartupInput, NULL);}
    ~GdiDeclare() { gdi::GdiplusShutdown(m_gdiplusToken);                            }
    ULONG_PTR m_gdiplusToken;
    gdi::GdiplusStartupInput gdiplusStartupInput;
    
} Gdideclare;


inline void assert_gdi(gdi::Status status)  { assert(status==gdi::Ok); }
inline void check_gdi (gdi::Status status,string errormessage)  { if(status!=gdi::Ok) throw error(errormessage); }

//Retrieves the encoder information for saving files of the formats: jpg,bmp,gif,tiff,png
inline void GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
   unsigned int num  = 0; // number of image encoders
   unsigned int size = 0; // size of the image encoder array in bytes

   gdi::ImageCodecInfo* pImageCodecInfo = NULL;

   gdi::GetImageEncodersSize(&num, &size);
   if(size == 0) throw error("ImageWrite: No Encoder found!");

   pImageCodecInfo = (gdi::ImageCodecInfo*)(malloc(size));
   if(pImageCodecInfo == NULL)  throw error("ImageWrite: Could not allocate Codec!");  // Failure

   gdi::GetImageEncoders(num, size, pImageCodecInfo);

   for(unsigned int j = 0; j < num; ++j)
      if( wcscmp(pImageCodecInfo[j].MimeType, format) == 0 )
      {
         *pClsid = pImageCodecInfo[j].Clsid;
         free(pImageCodecInfo);
         return;  // Success
      }    

   free(pImageCodecInfo);
   throw error("ImageWrite: Encoder not found!");
}



#endif // WIN32

#endif //GDIPLUS_H

