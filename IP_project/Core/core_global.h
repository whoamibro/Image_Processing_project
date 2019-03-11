#pragma once

/*
DLL configuration.
*/
#ifdef _CORE_EXPORTS_
#define MC_CORE_EXPORTS		__declspec(dllexport)
#else
#define	MC_CORE_EXPORTS		__declspec(dllimport)
#endif

