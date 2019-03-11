#pragma once
/*=========================================================================

File:		class raw_io_exception
Language:	Modern C++ (C++11/C++14)
Library:	Standard C++ Library
Author:		Minyoung Chung
Date:		2016-12-14
Version:	1.0
Mail:		chungmy.freddie@gmail.com

Copyright (c) 2017 All rights reserved by Bokkmani Studio corp.

=========================================================================*/
#include "core_global.h"

/*
* raw_io_exception class's exception class.
[ERROR_ID]
FILE_OPEN
FILE_READ
FILE_WRITE
UNKOWN
*/
namespace mc
{
	/*
	* exception class for raw_io.
	* similar interface with std::exception.
	*/
	class MC_CORE_EXPORTS raw_io_exception
	{
	public:
		// define exception IDs.
		enum class EID
		{
			FILE_OPEN,
			FILE_READ,
			FILE_WRITE,
			UNKNOWN
		};
	public:
		raw_io_exception(EID eid = EID::UNKNOWN, const char* str = nullptr) :
			m_eID(eid), m_rStrErr(str)
		{}
		inline const char* what() const noexcept { return m_rStrErr; }
		inline const EID ErrorType() const noexcept { return m_eID; }

	private:
		EID				m_eID;
		const char*		m_rStrErr;
	};
}

