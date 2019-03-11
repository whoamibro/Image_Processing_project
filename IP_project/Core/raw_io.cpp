#include "raw_io.h"
#include "raw_io_exception.h"
#include "../Common/image3d.h"
#include "../Common/memory.h"
#include <fstream>
#include <memory>

using mc::raw_io;
using mc::image3d;

/*
* Explicit Template Instantiation.
*/
template class __declspec(dllexport) raw_io<char>;
template class __declspec(dllexport) raw_io<unsigned char>;
template class __declspec(dllexport) raw_io<short>;
template class __declspec(dllexport) raw_io<unsigned short>;
template class __declspec(dllexport) raw_io<float>;

template <typename IMAGETYPE>
raw_io<IMAGETYPE>::raw_io(const char* const path, EENDIAN_TYPE eEndianT, EDATA_TYPE eDataT) noexcept :
m_rStrPath(path),
m_eEndianType(eEndianT),
m_eDataType(eDataT)
{}

template <typename IMAGETYPE>
raw_io<IMAGETYPE>::~raw_io() noexcept
{
}

template <typename IMAGETYPE>
image3d<IMAGETYPE>* raw_io<IMAGETYPE>::read(unsigned int width, unsigned int height, unsigned int depth)
{
	//typedef MyRawIOException EXCEPTION;
	using EXCEPTION = mc::raw_io_exception;

	const int iBufSize = width * height * depth;
	int szType;
	switch (m_eDataType)
	{
	case EDATA_TYPE::SHORT:
		szType = sizeof(short);
		break;
	case EDATA_TYPE::USHORT:
		szType = sizeof(unsigned short);
		break;
	case EDATA_TYPE::FLOAT:
		szType = sizeof(float);
		break;
	case EDATA_TYPE::UCHAR:
		szType = sizeof(unsigned char);
		break;
	}

	// input file stream creation.
	std::ifstream inFile(m_rStrPath, std::ios::binary | std::ios::in);
	if (!inFile.is_open())
		throw EXCEPTION(EXCEPTION::EID::FILE_OPEN, "file open error.");

	// create buffer for reading.
	std::unique_ptr<unsigned char> readBuffer(new unsigned char[iBufSize * szType]);
	std::memset(readBuffer.get(), 0, sizeof(char)*iBufSize * szType);

	// read.
	inFile.read(reinterpret_cast<char*>(readBuffer.get()), iBufSize * szType);
	if (m_eEndianType == EENDIAN_TYPE::BIG)
	{
		switch (m_eDataType)
		{
		case EDATA_TYPE::SHORT:
			convertToBidEndian<short>(reinterpret_cast<short*>(readBuffer.get()), iBufSize);
			break;
		case EDATA_TYPE::USHORT:
			convertToBidEndian<unsigned short>(reinterpret_cast<unsigned short*>(readBuffer.get()), iBufSize);
			break;
		case EDATA_TYPE::FLOAT:
			convertToBidEndian<float>(reinterpret_cast<float*>(readBuffer.get()), iBufSize);
			break;
		}
	}

	// create volume.
	auto pVol = new image3d<IMAGETYPE>(width, height, depth);

	// copy.
	switch (m_eDataType)
	{
	case EDATA_TYPE::SHORT:
		this->copyData<short>(pVol, readBuffer.get());
		break;
	case EDATA_TYPE::USHORT:
		this->copyData<unsigned short>(pVol, readBuffer.get());
		break;
	case EDATA_TYPE::FLOAT:
		this->copyData<float>(pVol, readBuffer.get());
		break;
	case EDATA_TYPE::UCHAR:
		this->copyData<unsigned char>(pVol, readBuffer.get());
		break;
	}

	// default settings.
	pVol->setPixelSpacing(1.0f);
	pVol->setSliceSpacing(1.0f);

	inFile.close();
	return pVol;
}

template <class IMAGETYPE> template <class DATATYPE>
void raw_io<IMAGETYPE>::copyData(image3d<IMAGETYPE> *vol, unsigned char *buf) noexcept
{
	auto bufPtr = reinterpret_cast<DATATYPE*>(buf);
	auto slicePtr = vol->data();
	for (int slice = 0; slice < vol->depth(); ++slice)
	{
		auto ptr = *(slicePtr++);
		for (int itr = 0; itr < vol->sizeSlice(); ++itr)
			//*(ptr++) = static_cast<short>(reinterpret_cast<TYPE*>(buf)[idx++]);
			*(ptr++) = static_cast<IMAGETYPE>(*(bufPtr++));
	}
}

template <class IMAGETYPE> template <typename DATATYPE>
void raw_io<IMAGETYPE>::convertToBidEndian(DATATYPE *data, int iSize) noexcept
{
	union {
		DATATYPE f;
		unsigned char c8[sizeof(DATATYPE)];
	} source, dest;
	for (int i = 0; i<iSize; ++i)
	{
		// convert
		source.f = data[i];
		for (size_t k = 0; k < sizeof(DATATYPE); ++k)
			dest.c8[k] = source.c8[sizeof(DATATYPE) - k - 1];
		data[i] = dest.f;
	}
}
