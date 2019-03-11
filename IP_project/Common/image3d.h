#pragma once
/*=========================================================================

File:		class image3d
Language:	Modern C++ (C++11/C++14)
Library:	Standard C++ Library
Author:		Minyoung Chung
Date:		2016-12-14
Version:	1.0
Mail:		chungmy.freddie@gmail.com

Copyright (c) 2017 All rights reserved by Bokkmani Studio corp.

=========================================================================*/
#include "common_global.h"

namespace mc
{
	/*class MC_COMMON_EXPORTS image3d_base
	{
	public:
	image3d_base() {}
	virtual ~image3d_base() {}

	};*/
	/*
	* image3d class description.
	- handle 3D volume data as 2D array.

	* Instantiated template types (available for use) :
	- char
	- unsigned char
	- short
	- unsigned short
	- float

	* Examples.
	auto Image1 = new image3d<short>(w, h, d);
	auto Image2 = new image3d<float>(w, h, d);
	*/
	template <typename TYPE>
	class MC_COMMON_EXPORTS image3d// : public image3d_base
	{
	public:
		explicit image3d(const int w, const int h, const int d) noexcept;
		~image3d() noexcept;
		image3d(const image3d& rhs) noexcept;				// copy-constructor.
		image3d& operator=(const image3d& rhs) noexcept;	// copy assignment.
		image3d(image3d&& rhs) noexcept;					// move-constructor.
		image3d& operator=(image3d&& rhs) noexcept;			// move assignment.

															/*
															* direct access to volume data.
															*/
		__forceinline TYPE** data() const noexcept { return m_ppData; }
		__forceinline TYPE* data(int sliceNum) const noexcept { return m_ppData[sliceNum]; }
		__forceinline TYPE get(int x, int y, int z) const noexcept {
			return m_ppData[z][x + y*m_iWidth];
		}

		/*
		* reset buffer. (zero_memory)
		*/
		void zeroImage(void) noexcept;

		__forceinline int width() const noexcept { return m_iWidth; }
		__forceinline int height() const noexcept { return m_iHeight; }
		__forceinline int depth() const noexcept { return m_iDepth; }

		// size of 2D slice.
		__forceinline int sizeSlice() const noexcept { return m_iWidth*m_iHeight; }
		// size of 3D volume.
		__forceinline int sizeVol() const noexcept { return m_iWidth*m_iHeight*m_iDepth; }

		// pixel spacing.
		__forceinline float pixelSpacing() const noexcept { return m_fPixelSpacing; }
		// slice spacing.
		__forceinline float sliceSpacing() const noexcept { return m_fSliceSpacing; }

		// set pixel spacing.
		__forceinline void setPixelSpacing(const float ps) noexcept { m_fPixelSpacing = ps; }
		// set slice spacing.
		__forceinline void setSliceSpacing(const float ss) noexcept { m_fSliceSpacing = ss; }

	private:
		// private member fields.
		TYPE**		m_ppData;			// volume data.
		int			m_iWidth;			// volume width.
		int			m_iHeight;			// volume height.
		int			m_iDepth;			// volume depth.
		float		m_fPixelSpacing;	// pixel spacing. (on-plane spacing)
		float		m_fSliceSpacing;	// slice spacing. (Z-axis spacing)
	};
}