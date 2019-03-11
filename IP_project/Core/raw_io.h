#pragma once
/*=========================================================================

File:		class raw_io
Language:	Modern C++ (C++11/C++14)
Library:	Standard C++ Library
Author:		Minyoung Chung
Date:		2016-12-14
Version:	1.0
Mail:		chungmy.freddie@gmail.com

Copyright (c) 2017 All rights reserved by Bokkmani Studio corp.

=========================================================================*/
#include "core_global.h"
#include <string>

namespace mc
{
	template <typename T> class image3d;
	/*
	* class raw_io.
	- Raw data i/o class.
	- configuration settings are optional.

	Default settings : Little-endian, Short type.

	* Instantiated template types (available image types) :
	- char
	- unsigned char
	- short
	- unsigned short
	- float

	* Examples.
	using namespace mc;
	typedef [char | unsigned char | short | unsigned short | float] IMAGETYPE;

	// 1. w/ default configurations.
	try {
	auto Image = raw_io<IMAGETYPE>("path").read(width, height, depth);
	}catch(raw_io_exception& e) {
	//  error handlings.
	// ...
	}
	// 2. w/ configuration settings.
	try {
	raw_io<IMAGETYPE> io("path");
	io.setEndianType(raw_io<IMAGETYPE>::EENDIAN_TYPE::BIG);
	io.setDataType(raw_io<IMAGETYPE>::EDATA_TYPE::FLOAT);
	auto Image = io.read(width, height, depth);
	}catch(raw_io_exception& e) {
	//  error handlings.
	// ...
	}
	*/
	template <typename IMAGETYPE>
	class MC_CORE_EXPORTS raw_io
	{
	public:
		enum class EENDIAN_TYPE {
			LITTLE,
			BIG
		};
		enum class EDATA_TYPE {
			SHORT,
			USHORT,
			FLOAT,
			CHAR,
			UCHAR
		};
		raw_io(const char* const path, EENDIAN_TYPE eEndianT = EENDIAN_TYPE::LITTLE, EDATA_TYPE eDataT = EDATA_TYPE::SHORT) noexcept;
		~raw_io() noexcept;

		// PUBLIC INTERFACES.
		/*
		* set endian type of target file.
		*/
		void setEndianType(EENDIAN_TYPE eType) noexcept { m_eEndianType = eType; }
		/*
		* set data type of target file.
		*/
		void setDataType(EDATA_TYPE eType) noexcept { m_eDataType = eType; }
		/*
		* read data and allocate image3d.
		*/
		mc::image3d<IMAGETYPE>* read(unsigned int w, unsigned int h, unsigned int d);

	private:
		template <typename DATATYPE>
		void convertToBidEndian(DATATYPE *data, int iSize) noexcept;
		template <typename DATATYPE>
		void copyData(mc::image3d<IMAGETYPE>* vol, unsigned char *buf) noexcept;

	private:
		const char* const	m_rStrPath;
		EENDIAN_TYPE		m_eEndianType;
		EDATA_TYPE			m_eDataType;
	};
}