#include <Windows.h>
#include <stdio.h>
#include <stdint.h>

//Open a memory mapped file from filename and return a pointer
LPVOID mmf_open(char * const fname);
//Open a memory mapped file and also out FileSize
LPVOID mmf_open_ex(char * const fname, LARGE_INTEGER *lpFileSize);

//Close a memory mapped file from pointer
BOOL mmf_close(LPVOID lpView);