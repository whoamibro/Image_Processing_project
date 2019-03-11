#pragma once
/*=========================================================================

File:		Memory macros definitions.
Language:	Modern C++ (C++11/C++14)
Library:	Standard C++ Library
Author:		Minyoung Chung
Date:		2016-12-14
Version:	1.0
Mail:		chungmy.freddie@gmail.com

Copyright (c) 2017 All rights reserved by Bokkmani Studio corp.

=========================================================================*/
#include <new>
#include <cassert>
#include <cstring>

namespace mc
{
	namespace memory
	{
		/*
		* inline function for allocation of memory.
		* noexcept.
		*/
		extern "C++"
			template <typename TYPE>
		__forceinline TYPE* allocate_1D(size_t size) noexcept
		{
			TYPE *p = new TYPE[size];
			return p;
		}
		/*
		* inline function for allocation of memory.
		* noexcept.
		*/
		extern "C++"
			template <typename TYPE>
		__forceinline TYPE** allocate_2D(size_t depth, size_t szSlice) noexcept
		{
			TYPE **p = new TYPE*[depth];
			for (int i = 0; i < static_cast<int>(depth); ++i) {
				p[i] = new TYPE[szSlice];
			}
			return p;
		}
	}
}

/*
* MACRO : safe allocation of memory.
* 1D
*/
#define SAFE_ALLOC_1D( TYPE, size )		mc::memory::allocate_1D<TYPE>((size))

/*
* MACRO : safe allocation of memory.
* "Volume allocation"
* OUTPUT : double pointer
*/
#define SAFE_ALLOC_VOLUME( TYPE, depth, szSlice )	mc::memory::allocate_2D<TYPE>((depth), (szSlice))

/*
* MACRO : safe deallocation of memory.
*/
#define SAFE_DELETE_OBJECT(p)	{ if(p) delete(p);		p = nullptr; }
#define SAFE_DELETE_ARRAY(p)	{ if(p) delete[](p);	p = nullptr; }
#define SAFE_DELETE_VOLUME(p, depth) { if(p) {for (int i=0; i<depth; ++i)	SAFE_DELETE_ARRAY(p[i]); SAFE_DELETE_ARRAY(p);} }

/*
* MACRO : copy & set memory.
*/
#define ZERO_VOLUME(TYPE, p, depth, szSlice)			{ if(p)		{for (int i=0; i<depth; ++i)	std::memset(p[i], 0, sizeof(TYPE)*szSlice);} }
#define MEMCPY_VOLUME(TYPE, dst, src, depth, szSlice)	{ if(dst)	{for (int i=0; i<depth; ++i)	std::memcpy(dst[i], src[i], sizeof(TYPE)*szSlice);} }
