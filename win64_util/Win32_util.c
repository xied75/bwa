#include "Win32_util.h"

//Created 2012-04-14, not tested
LPCWSTR str2Wide(char * const str)
{
    const int InitSize = 200;
    WCHAR * out; //I reckon the file names could be very long.
    int len = 0;

    out = (WCHAR*)calloc(200, sizeof(WCHAR));
    len = MultiByteToWideChar(CP_ACP, 0, str, -1, out, 0); //Please read MSDN on this function

    if (len > 0)
    {
        if (len > InitSize) out = (WCHAR*)realloc(out, len);
        MultiByteToWideChar(CP_ACP, 0, str, -1, out, len);
    }
    else
    {
        //error handling, not implemented yet.
    }

    return out;
}


LPVOID mmf_open(char * const fname)
{
    LPVOID lpView;
    HANDLE hFile, hMapFile;
    DWORD err;
    
    //LPCWSTR lpFileName = str2Wide(fname);

    //dwShareMode: 0, exclusive; 
    hFile = CreateFile(fname, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);

    if (hFile == INVALID_HANDLE_VALUE)
    {
        //error handling, not implemented yet.
        return NULL;
    }

    hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);

    if (hMapFile == NULL)
    {
        //error handling, not implemented yet.
        return NULL;
    }

    lpView = MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0);

    if (lpView == NULL)
    {
        //error handling, not implemented yet.
        return NULL;
    }

    return lpView;
}

BOOL mmf_close(LPVOID lpView)
{
    if (lpView != NULL) UnmapViewOfFile(lpView);

    return TRUE;
}

LPVOID mmf_open_ex(char * const fname, LARGE_INTEGER *lpFileSize)
{
    LPVOID lpView;
    HANDLE hFile, hMapFile;
    DWORD err;
    BOOL b;
    
    //LPCWSTR lpFileName = str2Wide(fname);

    //dwShareMode: 0, exclusive; 
    hFile = CreateFile(fname, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);

    if (hFile == INVALID_HANDLE_VALUE)
    {
        //error handling, not implemented yet.
        return NULL;
    }

    b = GetFileSizeEx(hFile, lpFileSize);

    if (!b)
    {
        //error handling, not implemented yet.
        return NULL;
    }

    hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);

    if (hMapFile == NULL)
    {
        //error handling, not implemented yet.
        return NULL;
    }

    //dwNumberOfBytesToMap=0: 
    lpView = MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0);

    if (lpView == NULL)
    {
        //error handling, not implemented yet.
        return NULL;
    }

    return lpView;
}