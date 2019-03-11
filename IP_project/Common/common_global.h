#pragma once

/*
DLL configuration.
*/
#ifdef _COMMON_EXPORTS_
#define MC_COMMON_EXPORTS	__declspec(dllexport)
#else
#define MC_COMMON_EXPORTS	__declspec(dllimport)
#endif

